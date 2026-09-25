//
// Log-likelihood of every fragment under all correlation patterns, in one pass.
//
// The original ran one job per (sample, phi, s) and wrote one file each, so the
// 21 patterns cost 21 full passes over the fragments and 21 files per sample.
// Everything those passes recomputed -- parsing the fragment, deciding which
// CpGs are usable, the probit thresholds, the sign flips, the cumulative
// positions -- depends only on the fragment, not on the pattern.  Here it is
// computed once and only the correlation matrix is rebuilt per pattern.
//
// Two further shortcuts, both exact rather than approximate:
//   K == 1  the answer is a Bernoulli log-probability and is the SAME for every
//           pattern; it is computed once and copied across the row.
//   K == 2  a bivariate normal orthant probability, available in closed form.
//           The EP solver is only needed from K >= 3.  This one is OFF by
//           default: the EP solver is approximate, and measured against
//           pbv::pbvnorm it carries up to 3.1e-02 absolute error in probability
//           at K = 2, so switching the exact form on changes the numbers the
//           original produced.  It is more accurate, not less -- but it is a
//           change, so it is opt-in.

#include <RcppArmadillo.h>
#include "ep_core.h"
// A_aux's helpers are pulled in as source rather than linked, so this file is a
// single self-contained translation unit.  Linking a separate object left a
// stale .o from another platform in the directory, which sourceCpp picked up.
#include "A_aux_impl.h"
#include <vector>
#include <string>
#include <cmath>

using namespace Rcpp;

// Split a comma-separated string of integers.  Same parsing the original used.
static inline std::vector<int> split_ints_local(const std::string& s) {
  std::vector<int> out;
  if (s.empty()) return out;
  const char* p = s.c_str();
  const char* end = p + s.size();
  while (p < end) {
    char* stop;
    long v = std::strtol(p, &stop, 10);
    if (stop == p) break;
    out.push_back((int)v);
    p = stop;
    while (p < end && (*p == ',' || *p == ' ')) ++p;
  }
  return out;
}

// Bivariate standard normal CDF, P(Z1 <= h, Z2 <= k) with correlation r.
// Drezner-Wesolowsky with 20-point Gauss-Legendre: accurate to ~1e-14, far
// beyond the 3 decimal places the pipeline keeps.
static double pbivnorm(double h, double k, double r) {
  if (r > 1.0) r = 1.0;
  if (r < -1.0) r = -1.0;
  const double Ph = R::pnorm5(h, 0.0, 1.0, 1, 0);
  const double Pk = R::pnorm5(k, 0.0, 1.0, 1, 0);
  if (std::fabs(r) < 1e-14) return Ph * Pk;
  if (r > 1.0 - 1e-14) return std::min(Ph, Pk);
  if (r < -1.0 + 1e-14) return std::max(0.0, Ph + Pk - 1.0);

  static const double x[20] = {
    -0.9931285991850949, -0.9639719272779138, -0.9122344282513259,
    -0.8391169718222188, -0.7463319064601508, -0.6360536807265150,
    -0.5108670019508271, -0.3737060887154195, -0.2277858511416451,
    -0.0765265211334973,  0.0765265211334973,  0.2277858511416451,
     0.3737060887154195,  0.5108670019508271,  0.6360536807265150,
     0.7463319064601508,  0.8391169718222188,  0.9122344282513259,
     0.9639719272779138,  0.9931285991850949 };
  static const double w[20] = {
    0.0176140071391521, 0.0406014298003869, 0.0626720483341091,
    0.0832767415767048, 0.1019301198172404, 0.1181945319615184,
    0.1316886384491766, 0.1420961093183820, 0.1491729864726037,
    0.1527533871307258, 0.1527533871307258, 0.1491729864726037,
    0.1420961093183820, 0.1316886384491766, 0.1181945319615184,
    0.1019301198172404, 0.0832767415767048, 0.0626720483341091,
    0.0406014298003869, 0.0176140071391521 };

  // integrate the density of the correlation from 0 to r
  const double a = 0.0, b = std::asin(r);
  const double c1 = 0.5 * (b - a), c2 = 0.5 * (b + a);
  double acc = 0.0;
  for (int i = 0; i < 20; ++i) {
    const double th = c1 * x[i] + c2;
    const double st = std::sin(th);
    const double denom = 2.0 * (1.0 - st * st);
    acc += w[i] * std::exp(-(h * h + k * k - 2.0 * h * k * st) / denom);
  }
  acc *= c1 / (2.0 * M_PI);
  double v = Ph * Pk + acc;
  if (v < 0.0) v = 0.0;
  if (v > 1.0) v = 1.0;
  return v;
}

//' Log-likelihood of each fragment under each correlation pattern
//'
//' @param fragment  character, the C/T/. string of each fragment
//' @param dist_sub  character, comma-separated bp gaps between consecutive CpGs
//' @param M_str,N_str character, comma-separated methylated / total counts
//' @param phi,s     the patterns; rho(d) = phi * exp(-d / s)
//' @param eps       clamp on the per-site methylation rate
//' @param min_n     minimum read depth for a site to be used
//' @param use_bvn   closed-form bivariate normal when exactly two sites are used
//' @return a matrix, one row per fragment and one column per pattern
// [[Rcpp::export]]
NumericMatrix loglik_all_patterns(CharacterVector fragment,
                                  CharacterVector dist_sub,
                                  CharacterVector M_str,
                                  CharacterVector N_str,
                                  NumericVector phi,
                                  NumericVector s,
                                  double eps = 1e-6,
                                  int min_n = 5,
                                  bool use_bvn = false,
                                  int print_every = 200000) {
  const int n = fragment.size();
  const int P = phi.size();
  if (s.size() != P) stop("phi and s must have the same length");

  NumericMatrix out(n, P);
  std::fill(out.begin(), out.end(), NA_REAL);

  // rho(d) = phi * exp(-d/s), tabulated once per pattern over the distances
  // that can occur.  The original built a 7001-long table per pattern; the same
  // cap is kept so distances beyond it saturate identically.
  const int MAXD = 7000;
  std::vector< std::vector<double> > rho_tab(P, std::vector<double>(MAXD + 1));
  for (int p = 0; p < P; ++p)
    for (int d = 0; d <= MAXD; ++d)
      rho_tab[p][d] = phi[p] * std::exp(-(double)d / s[p]);

  std::vector<int> keep_idx, x_full, cum;
  std::vector<double> mu, sgn, upper;

  for (int r = 0; r < n; ++r) {
    if (print_every > 0 && r > 0 && (r % print_every == 0))
      Rcpp::Rcout << "  " << r << " / " << n << "\n";

    if (fragment[r] == NA_STRING || M_str[r] == NA_STRING || N_str[r] == NA_STRING)
      continue;

    const std::string frag = as<std::string>(fragment[r]);
    const int L = (int)frag.size();
    std::vector<int> Mv = split_ints_local(as<std::string>(M_str[r]));
    std::vector<int> Nv = split_ints_local(as<std::string>(N_str[r]));
    if ((int)Mv.size() != L || (int)Nv.size() != L) continue;

    // ---- which sites are usable, and their probit thresholds ----------------
    keep_idx.clear(); x_full.assign(L, 0);
    for (int i = 0; i < L; ++i) {
      const char c = frag[i];
      if (c == 'C' || c == 'c') x_full[i] = +1;
      else if (c == 'T' || c == 't') x_full[i] = -1;
      if (x_full[i] != 0 && Nv[i] >= min_n) keep_idx.push_back(i);
    }
    const int K = (int)keep_idx.size();
    if (K == 0) continue;

    bool bad = false;
    mu.assign(K, 0.0);
    for (int k = 0; k < K; ++k) {
      const int j = keep_idx[k];
      if (Nv[j] == 0) { bad = true; break; }
      double p = (double)Mv[j] / (double)Nv[j];
      if (p <= eps) p = eps;
      if (p >= 1.0 - eps) p = 1.0 - eps;
      mu[k] = R::qnorm5(p, 0.0, 1.0, 1, 0);
    }
    if (bad) continue;

    // ---- K == 1: independent of the pattern, so fill the row and move on ----
    if (K == 1) {
      const int j = keep_idx[0];
      double p = (double)Mv[j] / (double)Nv[j];
      if (p <= eps) p = eps;
      if (p >= 1.0 - eps) p = 1.0 - eps;
      const double lp = (x_full[j] == +1) ? std::log(p) : std::log(1.0 - p);
      const double v = std::round(lp * 1000.0) / 1000.0;
      for (int p_i = 0; p_i < P; ++p_i) out(r, p_i) = v;
      continue;
    }

    // ---- cumulative bp positions of the fragment's CpGs ----------------------
    if (dist_sub[r] == NA_STRING) continue;
    std::vector<int> dvec = split_ints_local(as<std::string>(dist_sub[r]));
    if ((int)dvec.size() != L - 1) continue;
    cum.assign(L, 0);
    for (int i = 1; i < L; ++i) cum[i] = cum[i - 1] + dvec[i - 1];

    // ---- signs and integration limits (pattern-independent) ------------------
    sgn.assign(K, 0.0); upper.assign(K, 0.0);
    for (int k = 0; k < K; ++k) {
      const int j = keep_idx[k];
      sgn[k]   = (x_full[j] == +1) ? -1.0 : 1.0;
      upper[k] = (x_full[j] == +1) ? mu[k] : -mu[k];
    }

    // ---- K == 2: bivariate normal in closed form ----------------------------
    if (K == 2 && use_bvn) {
      const int d = std::abs(cum[keep_idx[1]] - cum[keep_idx[0]]);
      const int dd = d > MAXD ? MAXD : d;
      for (int p_i = 0; p_i < P; ++p_i) {
        const double rho = rho_tab[p_i][dd] * sgn[0] * sgn[1];
        const double pr = pbivnorm(upper[0], upper[1], rho);
        const double lp = (pr <= 0.0) ? R_NegInf : std::log(pr);
        out(r, p_i) = R_finite(lp) ? std::round(lp * 1000.0) / 1000.0 : NA_REAL;
      }
      continue;
    }

    // ---- K >= 3: the EP solver, once per pattern ----------------------------
    NumericVector upperR(K);
    for (int k = 0; k < K; ++k) upperR[k] = upper[k];

    for (int p_i = 0; p_i < P; ++p_i) {
      NumericMatrix Sigma(K, K);
      for (int i = 0; i < K; ++i) {
        Sigma(i, i) = 1.0;
        for (int j = i + 1; j < K; ++j) {
          int d = std::abs(cum[keep_idx[j]] - cum[keep_idx[i]]);
          if (d > MAXD) d = MAXD;
          const double v = rho_tab[p_i][d] * sgn[i] * sgn[j];
          Sigma(i, j) = v;
          Sigma(j, i) = v;
        }
      }
      try {
        const double logcdf = FAStCDF_algo1_cpp(upperR, Sigma, true, 2.0, 1e-4, 20000);
        out(r, p_i) = R_finite(logcdf) ? std::round(logcdf * 1000.0) / 1000.0 : NA_REAL;
      } catch (...) {
        out(r, p_i) = NA_REAL;
      }
    }
  }
  return out;
}
