// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <string>
#include <sstream>
#include <cstring>
using namespace Rcpp;

#include "A_aux.h"

// [[Rcpp::export]]
Rcpp::List getParamsEP_priorMean_algo1_cpp(
    arma::mat X,
    arma::vec y,
    arma::vec priorMean,
    double nu2,
    double tolerance,
    int maxIter=1e4){
  
  int n = X.n_rows;
  int p = X.n_cols;
  
  double log_tolerance = log(tolerance);
  
  arma::mat Xt = X.t(); 
  double diff = 1;
  
  int count =0;
  arma::vec r = priorMean/nu2;
  double logDetinvQ = -p*log(nu2);
  
  // Containers
  arma::colvec logZ(n);
  arma::colvec k(n);
  arma::mat repl_k(n,p);
  arma::colvec m(n);
  
  arma::vec meanBeta(p);
  
  
  arma::mat invQ = arma::eye(p,p) * nu2 ; 
  
  int nnIter = 0;
  
  for(int nIter = 0; nIter<maxIter; nIter++){
    Rcpp::checkUserInterrupt();
    if(diff < log_tolerance){nnIter = nIter; 
      break; }
    
    diff = - arma::datum::inf;
    count = 0;
    
    for(int i=0; i<n; i++){
      
      arma::vec xi = (Xt.col(i));
      arma::vec r_i = r - m[i]*xi;
      
      arma::mat Oxi = invQ * xi;
      double xitOxi = arma::dot(xi,Oxi);
      /////////////////////////////////////////
      arma::mat Oi = invQ + (Oxi*Oxi.t())*k[i] / (1.0 - k[i]*xitOxi);
      /////////////////////////////////////////
      arma::mat Oixi = Oi*xi;
      double xiOixi = arma::dot(xi, Oixi);
      
      
      if(xiOixi>0){
        
        double r_iOixi = arma::dot(r_i,Oixi);
        
        double s = (2.0 * y[i] - 1.0) / sqrt(1.0 + xiOixi);
        double tau = s * r_iOixi;
        
        double z1 = zeta1(tau);
        double z2 = zeta2(tau,z1);
        
        double kNew = -z2/(1.0 + xiOixi + z2 * xiOixi);
        double delta_k = kNew - k[i];
        double mNew = (z1*s + kNew * r_iOixi + kNew * z1 * s * xiOixi);
        double delta_m = mNew - m[i];
        
        double rel_logdelta_k = log(std::fabs(delta_k)) - log(std::fabs(k[i]));
        double rel_logdelta_m = log(std::fabs(delta_m)) - log(std::fabs(m[i]));
        
        
        k[i] = kNew;
        m[i] = mNew;
        
        double prev_logZ = logZ[i];
        logZ[i] = (( 2 * m[i] * r_iOixi + 
          pow(m[i],2) * xiOixi - k[i]*pow(r_iOixi,2))/
            ( 1 + k[i] * xiOixi ) - 
              log(1.0 + k[i] * xiOixi ) )*.5 - 
              log(arma::normcdf(tau));
        double delta_Z = exp(logZ[i]) - exp(prev_logZ);
        double rel_logdelta_Z = log(std::fabs(delta_Z)) - logZ[i];
        
        
        r = r_i + m[i]*xi;
        
        double denominator = (1.+(delta_k) * quadform(xi,invQ) );
        logDetinvQ += log(denominator);
        invQ = Oi + z2 * pow(s,2)* (Oixi*Oixi.t());
        
        diff = check_convergence_diffs(rel_logdelta_k,
                                       rel_logdelta_m,
                                       rel_logdelta_Z,
                                       diff);
        
      }else{
        count = count+1;
        Rcpp::Rcout << count << " units skipped\n";;
      }
      
      
      
    }
    
    
  }
  
  meanBeta  = invQ * r;
  
  
  double logML = (arma::dot(r,meanBeta) - logDetinvQ - p*log(nu2) - 
                  arma::dot(priorMean,priorMean)/nu2) *.5 - arma::accu(logZ);
  
  arma::colvec id = invQ.diag();
  List results = List::create(_["meanBeta"] = meanBeta,
                              _["diagOmega"] = id,
                              _["logML"] = logML,
                              _["nIter"] = nnIter,
                              _["kEP"] = k,
                              _["mEP"] = m);
  return(results);
  
}


// [[Rcpp::export]]
double FAStCDF_algo1_cpp(const NumericVector& x,
                         const NumericMatrix& Sigma,
                         bool logp = false,
                         double eps = 2.0,
                         double tolerance = 1e-5,
                         int maxIter = 1000000)
{
  // 1) handle Inf: drop the corresponding dimensions
  std::vector<int> keep;
  keep.reserve(x.size());
  for (int i = 0; i < x.size(); ++i) {
    if (R_finite(x[i])) keep.push_back(i);
  }
  if (keep.empty()) return logp ? 0.0 : 1.0;  // log(1)=0
  
  const int n = (int)keep.size();
  
  arma::vec xk(n);
  arma::mat Sk(n, n);
  for (int i = 0; i < n; ++i) {
    xk(i) = x[ keep[i] ];
    for (int j = 0; j < n; ++j) {
      Sk(i, j) = Sigma( keep[i], keep[j] );
    }
  }
  
  // 2) x <- x - mu (mu is the zero vector, so x is unchanged) -- kept for parity with the R version:
  //    xk is used directly here
  
  // 3) eigen & rescale
  arma::vec eigvals;
  arma::eig_sym(eigvals, Sk);              // ascending
  double smallestEig = eigvals(0);
  
  // optional safeguard: if the smallest eigenvalue <= 0, lift it slightly (left off by default, for parity with the R version)
  // if (smallestEig <= 0) smallestEig = std::max(smallestEig, 1e-12);
  
  const double scalingOfEig = eps;         // as in the R version
  arma::mat tildeSigma = Sk - (smallestEig / scalingOfEig) * arma::eye(n, n);
  
  // 4) Cholesky (lower)
  arma::mat L;
  bool ok = arma::chol(L, tildeSigma, "lower");  // L * L^T = tildeSigma
  if (!ok) {
    // on extreme numerical trouble, add a small lift and retry
    tildeSigma.diag() += 1e-12;
    if (!arma::chol(L, tildeSigma, "lower"))
      Rcpp::stop("Cholesky failed in FAStCDF_algo1_cpp.");
  }
  
  // 5) priorMean = sqrt(eps / smallestEig) * solve(L, xk)
  const double nu2 = scalingOfEig / smallestEig;
  const double scale = std::sqrt(nu2);
  // solve L * z = xk  -> z = solve(L, xk) (lower triangular)
  arma::vec z = arma::solve(arma::trimatl(L), xk);
  arma::vec priorMean = scale * z;
  
  // 6) call the main EP routine
  arma::vec yvec(n, arma::fill::ones);
  Rcpp::List paramsEP = getParamsEP_priorMean_algo1_cpp(
    L,              // X
    yvec,           // y = 1
    priorMean,      // priorMean
    nu2,            // nu2
    tolerance,
    maxIter
  );
  
  double logML = Rcpp::as<double>(paramsEP["logML"]);
  return logp ? logML : std::exp(logML);
}


static inline std::vector<int> split_ints(const std::string& s) {
  std::vector<int> out;
  int n = s.size();
  int i = 0;
  while (i < n) {
    while (i < n && (s[i] == ' ' || s[i] == ',')) i++;
    if (i >= n) break;
    int sign = 1;
    if (s[i] == '+') { sign = 1; i++; }
    else if (s[i] == '-') { sign = -1; i++; }
    long long v = 0; bool any = false;
    while (i < n && s[i] >= '0' && s[i] <= '9') {
      v = v*10 + (s[i]-'0');
      i++; any = true;
    }
    if (any) out.push_back((int)(sign*v));
    while (i < n && s[i] != ',') i++;
  }
  return out;
}




static inline void fast_split_ints(const char* s, std::vector<int>& out) {
  out.clear();
  if (!s) return;
  
  int sign = 1, val = 0;
  bool in_num = false;
  
  for (const char* p = s; *p; ++p) {
    unsigned char c = (unsigned char)(*p);
    if (c == '-') {
      sign = -1;
      val = 0;
      in_num = true;
    } else if (c >= '0' && c <= '9') {
      if (!in_num) {
        sign = 1;
        val = 0;
        in_num = true;
      }
      val = val * 10 + (c - '0');
    } else {
      if (in_num) {
        out.push_back(sign * val);
        in_num = false;
        sign = 1;
        val = 0;
      }
    }
  }
  if (in_num) out.push_back(sign * val);
}

// [[Rcpp::export]]
NumericVector compute_loglik_cpp_coord(CharacterVector fragment,
                                       CharacterVector dist_sub,
                                       List M_list,
                                       List N_list,
                                       NumericVector exp_table,
                                       double eps = 1e-6,
                                       int print_every = 100000) {
  
  const int n = fragment.size();
  NumericVector out(n, NA_REAL);
  const int max_d = exp_table.size() - 1;
  
  std::vector<int> dvec;
  std::vector<int> keep_idx;
  std::vector<int> x;
  std::vector<double> sgn;
  std::vector<int> cum_keep;
  std::vector<double> upper;
  
  for (int r = 0; r < n; ++r) {
    if (print_every > 0 && (r % print_every == 0)) {
      Rcpp::Rcout << "Processing row " << (r + 1) << " / " << n << std::endl;
    }
    
    if (fragment[r] == NA_STRING) continue;
    
    Rcpp::String fragS(fragment[r]);
    const char* frag = fragS.get_cstring();
    const int L = (int)std::strlen(frag);
    
    const IntegerVector Mv = M_list[r];
    const IntegerVector Nv = N_list[r];
    if (Mv.size() != L || Nv.size() != L) continue;
    
    if (L == 1) {
      const char c = frag[0];
      if (!((c == 'C' || c == 'c' || c == 'T' || c == 't') && Nv[0] >= 5)) continue;
      if (Nv[0] == 0) continue;
      double p = (double)Mv[0] / (double)Nv[0];
      if (p <= eps) p = eps;
      if (p >= 1.0 - eps) p = 1.0 - eps;
      double logprob = (c == 'C' || c == 'c') ? std::log(p) : std::log(1.0 - p);
      out[r] = std::round(logprob * 1000.0) / 1000.0;
      continue;
    }
    
    if (dist_sub[r] == NA_STRING) continue;
    
    Rcpp::String dS(dist_sub[r]);
    const char* dstr = dS.get_cstring();
    fast_split_ints(dstr, dvec);
    if ((int)dvec.size() != L - 1) continue;
    
    keep_idx.clear();
    keep_idx.reserve(L);
    
    for (int i = 0; i < L; ++i) {
      const char c = frag[i];
      if ((c == 'C' || c == 'c' || c == 'T' || c == 't') && Nv[i] >= 5) {
        keep_idx.push_back(i);
      }
    }
    
    const int K = (int)keep_idx.size();
    if (K == 0) { out[r] = 0.0; continue; }
    
    if (K == 1) {
      const int j = keep_idx[0];
      if (Nv[j] == 0) continue;
      double p = (double)Mv[j] / (double)Nv[j];
      if (p <= eps) p = eps;
      if (p >= 1.0 - eps) p = 1.0 - eps;
      const char c = frag[j];
      double logprob = (c == 'C' || c == 'c') ? std::log(p) : std::log(1.0 - p);
      out[r] = std::round(logprob * 1000.0) / 1000.0;
      continue;
    }
    
    x.resize(K);
    sgn.resize(K);
    cum_keep.resize(K);
    upper.resize(K);
    
    int cum = 0;
    cum_keep[0] = 0;
    int kk = 1;
    for (int i = 1; i < L; ++i) {
      cum += dvec[i - 1];
      if (kk < K && i == keep_idx[kk]) {
        cum_keep[kk] = cum;
        ++kk;
      }
    }
    
    bool valid = true;
    NumericVector mu(K);
    for (int k = 0; k < K; ++k) {
      const int j = keep_idx[k];
      if (Nv[j] == 0) { valid = false; break; }
      double p = (double)Mv[j] / (double)Nv[j];
      if (p <= eps) p = eps;
      if (p >= 1.0 - eps) p = 1.0 - eps;
      
      mu[k] = R::qnorm5(p, 0.0, 1.0, 1, 0);
      
      const char c = frag[j];
      const int xi = (c == 'C' || c == 'c') ? +1 : -1;
      x[k] = xi;
      sgn[k] = (xi == +1) ? -1.0 : 1.0;
      upper[k] = (xi == +1) ? mu[k] : -mu[k];
    }
    if (!valid) continue;
    
    NumericMatrix Sigma(K, K);
    for (int i = 0; i < K; ++i) {
      Sigma(i, i) = 1.0;
      const double si = sgn[i];
      const int ci = cum_keep[i];
      for (int j = i + 1; j < K; ++j) {
        int d = std::abs(cum_keep[j] - ci);
        if (d > max_d) d = max_d;
        const double rho = exp_table[d] * si * sgn[j];
        Sigma(i, j) = rho;
        Sigma(j, i) = rho;
      }
    }
    
    try {
      NumericVector upperR(K);
      for (int i = 0; i < K; ++i) upperR[i] = upper[i];
      
      double logcdf = FAStCDF_algo1_cpp(
        upperR, Sigma, true, 2.0, 1e-4, 20000
      );
      if (!R_finite(logcdf)) out[r] = NA_REAL;
      else out[r] = std::round(logcdf * 1000.0) / 1000.0;
    } catch (...) {
      out[r] = NA_REAL;
    }
  }
  
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector compute_loglik_cpp_coord_include(
    Rcpp::CharacterVector fragment,
    Rcpp::CharacterVector dist_sub,
    Rcpp::List M_list,
    Rcpp::List N_list,
    Rcpp::NumericVector exp_table,
    Rcpp::IntegerVector start_idx,              // NEW: length n, 1-based global new_index start
    Rcpp::IntegerVector include_by_new_index,   // NEW: length max_new, 1=include 0=exclude
    double eps = 1e-6,
    int print_every = 100000) {
  
  const int n = fragment.size();
  Rcpp::NumericVector out(n, NA_REAL);
  const int max_d = exp_table.size() - 1;
  
  std::vector<int> dvec;
  std::vector<int> keep_idx;
  std::vector<int> x;
  std::vector<double> sgn;
  std::vector<int> cum_keep;
  std::vector<double> upper;
  
  for (int r = 0; r < n; ++r) {
    
    if (print_every > 0 && (r % print_every == 0)) {
      Rcpp::Rcout << "Processing row " << (r + 1) << " / " << n << std::endl;
    }
    
    if (fragment[r] == NA_STRING) continue;
    
    Rcpp::String fragS(fragment[r]);
    const char* frag = fragS.get_cstring();
    const int L = (int)std::strlen(frag);
    
    // start_idx sanity
    if (start_idx.size() != n) continue;
    const int st = start_idx[r]; // 1-based
    if (st == NA_INTEGER || st <= 0) continue;
    
    const Rcpp::IntegerVector Mv = M_list[r];
    const Rcpp::IntegerVector Nv = N_list[r];
    if (Mv.size() != L || Nv.size() != L) continue;
    
    // Special case L==1
    if (L == 1) {
      // include check: global idx = st
      if (st < 1 || st > include_by_new_index.size()) continue;
      if (include_by_new_index[st - 1] == 0) continue; // excluded => treat as '.'
      
      const char c = frag[0];
      if (!((c == 'C' || c == 'c' || c == 'T' || c == 't') && Nv[0] >= 5)) continue;
      if (Nv[0] == 0) continue;
      
      double p = (double)Mv[0] / (double)Nv[0];
      if (p <= eps) p = eps;
      if (p >= 1.0 - eps) p = 1.0 - eps;
      
      double logprob = (c == 'C' || c == 'c') ? std::log(p) : std::log(1.0 - p);
      out[r] = std::round(logprob * 1000.0) / 1000.0;
      continue;
    }
    
    if (dist_sub[r] == NA_STRING) continue;
    
    // parse distances
    Rcpp::String dS(dist_sub[r]);
    const char* dstr = dS.get_cstring();
    fast_split_ints(dstr, dvec);
    if ((int)dvec.size() != L - 1) continue;
    
    // choose kept CpGs: must be (C/T), depth>=5, AND include_by_new_index==1
    keep_idx.clear();
    keep_idx.reserve(L);
    
    for (int i = 0; i < L; ++i) {
      const int g = st + i; // 1-based global new_index for position i (0-based)
      if (g < 1 || g > include_by_new_index.size()) continue;
      if (include_by_new_index[g - 1] == 0) continue; // excluded => treat as '.'
      
      const char c = frag[i];
      if ((c == 'C' || c == 'c' || c == 'T' || c == 't') && Nv[i] >= 5) {
        keep_idx.push_back(i);
      }
    }
    
    const int K = (int)keep_idx.size();
    
    // If nothing kept, match your previous behavior
    if (K == 0) { out[r] = 0.0; continue; }
    
    if (K == 1) {
      const int j = keep_idx[0];
      if (Nv[j] == 0) continue;
      
      double p = (double)Mv[j] / (double)Nv[j];
      if (p <= eps) p = eps;
      if (p >= 1.0 - eps) p = 1.0 - eps;
      
      const char c = frag[j];
      double logprob = (c == 'C' || c == 'c') ? std::log(p) : std::log(1.0 - p);
      out[r] = std::round(logprob * 1000.0) / 1000.0;
      continue;
    }
    
    // allocate per-row
    x.resize(K);
    sgn.resize(K);
    cum_keep.resize(K);
    upper.resize(K);
    
    // cum distances for kept positions
    int cum = 0;
    cum_keep[0] = 0;
    int kk = 1;
    for (int i = 1; i < L; ++i) {
      cum += dvec[i - 1];
      if (kk < K && i == keep_idx[kk]) {
        cum_keep[kk] = cum;
        ++kk;
      }
    }
    
    bool valid = true;
    Rcpp::NumericVector mu(K);
    
    for (int k = 0; k < K; ++k) {
      const int j = keep_idx[k];
      if (Nv[j] == 0) { valid = false; break; }
      
      double p = (double)Mv[j] / (double)Nv[j];
      if (p <= eps) p = eps;
      if (p >= 1.0 - eps) p = 1.0 - eps;
      
      mu[k] = R::qnorm5(p, 0.0, 1.0, 1, 0);
      
      const char c = frag[j];
      const int xi = (c == 'C' || c == 'c') ? +1 : -1;
      x[k] = xi;
      sgn[k] = (xi == +1) ? -1.0 : 1.0;
      upper[k] = (xi == +1) ? mu[k] : -mu[k];
    }
    if (!valid) continue;
    
    Rcpp::NumericMatrix Sigma(K, K);
    for (int i = 0; i < K; ++i) {
      Sigma(i, i) = 1.0;
      const double si = sgn[i];
      const int ci = cum_keep[i];
      for (int j = i + 1; j < K; ++j) {
        int d = std::abs(cum_keep[j] - ci);
        if (d > max_d) d = max_d;
        const double rho = exp_table[d] * si * sgn[j];
        Sigma(i, j) = rho;
        Sigma(j, i) = rho;
      }
    }
    
    try {
      Rcpp::NumericVector upperR(K);
      for (int i = 0; i < K; ++i) upperR[i] = upper[i];
      
      double logcdf = FAStCDF_algo1_cpp(upperR, Sigma, true, 2.0, 1e-4, 20000);
      if (!R_finite(logcdf)) out[r] = NA_REAL;
      else out[r] = std::round(logcdf * 1000.0) / 1000.0;
      
    } catch (...) {
      out[r] = NA_REAL;
    }
  }
  
  return out;
}


// [[Rcpp::export]]
NumericVector compute_loglik_cpp_coord_0212(CharacterVector fragment,
                                            CharacterVector dist_sub,
                                            List M_list,              // list of IntegerVector
                                            List N_list,              // list of IntegerVector
                                            NumericVector exp_table,
                                            double eps = 1e-6,
                                            int print_every = 100000) {
  
  const int n = fragment.size();
  NumericVector out(n, NA_REAL);
  const int max_d = exp_table.size() - 1;
  
  for (int r = 0; r < n; ++r) {
    if (print_every > 0 && (r % print_every == 0)) {
      Rcpp::Rcout << "Processing row " << (r + 1) << " / " << n << std::endl;
    }
    
    if (fragment[r] == NA_STRING) continue;
    
    const std::string frag = as<std::string>(fragment[r]);
    const int L = (int)frag.size();
    
    // take the vectors straight from the list
    IntegerVector Mv_r = M_list[r];
    IntegerVector Nv_r = N_list[r];
    if (Mv_r.size() != L || Nv_r.size() != L) continue;
    
    std::vector<int> Mv(Mv_r.begin(), Mv_r.end());
    std::vector<int> Nv(Nv_r.begin(), Nv_r.end());
    
    // ---------- L == 1 special case ----------
    if (L == 1) {
      const char c = frag[0];
      if (!(c == 'C' || c == 'c' || c == 'T' || c == 't') || Nv[0] < 5) continue;
      if (Nv[0] == 0) continue;
      double p = (double)Mv[0] / (double)Nv[0];
      if (p <= eps) p = eps;
      if (p >= 1.0 - eps) p = 1.0 - eps;
      
      double logprob = (c == 'C' || c == 'c') ? std::log(p) : std::log(1.0 - p);
      out[r] = std::round(logprob * 1000.0) / 1000.0;
      continue;
    }
    
    // ---------- L >= 2 needs dist_sub ----------
    if (dist_sub[r] == NA_STRING) continue;
    const std::string dstr = as<std::string>(dist_sub[r]);
    std::vector<int> dvec = split_ints(dstr);
    if ((int)dvec.size() != L - 1) continue;
    
    // select the usable sites
    std::vector<int> keep_idx; keep_idx.reserve(L);
    std::vector<int> x_full(L, 0);
    for (int i = 0; i < L; ++i) {
      const char c = frag[i];
      if (c == 'C' || c == 'c') x_full[i] = +1;
      else if (c == 'T' || c == 't') x_full[i] = -1;
      else x_full[i] = 0;
      if (x_full[i] != 0 && Nv[i] >= 5) keep_idx.push_back(i);
    }
    const int K = (int)keep_idx.size();
    if (K == 0) { out[r] = 0.0; continue; }
    
    // K == 1 special case
    if (K == 1) {
      const int j = keep_idx[0];
      if (Nv[j] == 0) continue;
      double p = (double)Mv[j] / (double)Nv[j];
      if (p <= eps) p = eps;
      if (p >= 1.0 - eps) p = 1.0 - eps;
      double logprob = (x_full[j] == +1) ? std::log(p) : std::log(1.0 - p);
      out[r] = std::round(logprob * 1000.0) / 1000.0;
      continue;
    }
    
    // K >= 2
    std::vector<int> cum(L, 0);
    for (int i = 1; i < L; ++i) cum[i] = cum[i-1] + dvec[i-1];
    
    NumericVector mu(K);
    std::vector<int> x(K);
    bool valid = true;
    for (int k = 0; k < K; ++k) {
      const int j = keep_idx[k];
      if (Nv[j] == 0) { valid = false; break; }
      double p = (double)Mv[j] / (double)Nv[j];
      if (p <= eps) p = eps;
      if (p >= 1.0 - eps) p = 1.0 - eps;
      mu[k] = R::qnorm5(p, 0.0, 1.0, 1, 0);
      x[k]  = x_full[j];
    }
    if (!valid) continue;
    
    std::vector<double> sgn(K);
    for (int i = 0; i < K; ++i) sgn[i] = (x[i] == 1) ? -1.0 : 1.0;
    
    NumericMatrix Sigma(K, K);
    for (int i = 0; i < K; ++i) {
      Sigma(i,i) = 1.0;
      const int ii = keep_idx[i];
      for (int j = i+1; j < K; ++j) {
        const int jj = keep_idx[j];
        int d = std::abs(cum[jj] - cum[ii]);
        if (d > max_d) d = max_d;
        const double rho = exp_table[d];
        Sigma(i,j) = rho;
        Sigma(j,i) = rho;
      }
    }
    
    for (int i = 0; i < K; ++i) {
      const double si = sgn[i];
      for (int j = 0; j < K; ++j) Sigma(i,j) *= si;
    }
    for (int j = 0; j < K; ++j) {
      const double sj = sgn[j];
      for (int i = 0; i < K; ++i) Sigma(i,j) *= sj;
    }
    
    NumericVector upper(K);
    for (int i = 0; i < K; ++i) upper[i] = (x[i] == 1) ? mu[i] : -mu[i];
    
    try {
      double logcdf = FAStCDF_algo1_cpp(
        upper, Sigma, true, 2.0, 1e-5, 100000
      );
      if (!R_finite(logcdf)) out[r] = NA_REAL;
      else out[r] = std::round(logcdf * 1000.0) / 1000.0;
    } catch (...) {
      out[r] = NA_REAL;
    }
  }
  
  return out;
}

// [[Rcpp::export]]
NumericVector compute_loglik_cpp_coord_string(CharacterVector fragment,
                                       CharacterVector dist_sub,
                                       CharacterVector M_str,
                                       CharacterVector N_str,
                                       NumericVector exp_table,
                                       double eps = 1e-6,
                                       int print_every = 100000) {
  const int n = fragment.size();
  NumericVector out(n, NA_REAL);
  const int max_d = exp_table.size() - 1;
  
  for (int r = 0; r < n; ++r) {
    if (print_every > 0 && (r % print_every == 0)) {
      Rcpp::Rcout << "Processing row " << (r + 1) << " / " << n << std::endl;
    }
    if (fragment[r] == NA_STRING ||
        M_str[r] == NA_STRING || N_str[r] == NA_STRING) {
      continue;
    }
    const std::string frag = as<std::string>(fragment[r]);
    const std::string mstr = as<std::string>(M_str[r]);
    const std::string nstr = as<std::string>(N_str[r]);
    const int L = (int)frag.size();
    
    std::vector<int> Mv = split_ints(mstr);
    std::vector<int> Nv = split_ints(nstr);
    if ((int)Mv.size() != L || (int)Nv.size() != L) {
      continue;
    }
    
    if (L == 1) {
      const char c = frag[0];
      if (!(c == 'C' || c == 'c' || c == 'T' || c == 't') || Nv[0] < 5) {
        continue;
      }
      if (Nv[0] == 0) {
        continue;
      }
      double p = (double)Mv[0] / (double)Nv[0];
      if (p <= eps) p = eps;
      if (p >= 1.0 - eps) p = 1.0 - eps;
      double logprob;
      if (c == 'C' || c == 'c') {
        logprob = std::log(p);
      } else {
        logprob = std::log(1.0 - p);
      }
      out[r] = std::round(logprob * 1000.0) / 1000.0;
      continue;
    }
    
    if (dist_sub[r] == NA_STRING) {
      continue;
    }
    const std::string dstr = as<std::string>(dist_sub[r]);
    std::vector<int> dvec = split_ints(dstr);
    if ((int)dvec.size() != L - 1) {
      continue;
    }
    
    std::vector<int> keep_idx; keep_idx.reserve(L);
    std::vector<int> x_full(L, 0);
    for (int i = 0; i < L; ++i) {
      const char c = frag[i];
      if (c == 'C' || c == 'c') x_full[i] = +1;
      else if (c == 'T' || c == 't') x_full[i] = -1;
      else x_full[i] = 0;
      if (x_full[i] != 0 && Nv[i] >= 5) keep_idx.push_back(i);
    }
    const int K = (int)keep_idx.size();
    if (K == 0) {
      continue;
    }
    
    if (K == 1) {
      const int j = keep_idx[0];
      if (Nv[j] == 0) {
        continue;
      }
      double p = (double)Mv[j] / (double)Nv[j];
      if (p <= eps) p = eps;
      if (p >= 1.0 - eps) p = 1.0 - eps;
      double logprob;
      if (x_full[j] == +1) {
        logprob = std::log(p);
      } else if (x_full[j] == -1) {
        logprob = std::log(1.0 - p);
      } else {
        continue;
      }
      out[r] = std::round(logprob * 1000.0) / 1000.0;
      continue;
    }
    
    std::vector<int> cum(L, 0);
    for (int i = 1; i < L; ++i) cum[i] = cum[i-1] + dvec[i-1];
    
    NumericVector mu(K);
    std::vector<int> x(K);
    bool valid = true;
    for (int k = 0; k < K; ++k) {
      const int j = keep_idx[k];
      if (Nv[j] == 0) { valid = false; break; }
      double p = (double)Mv[j] / (double)Nv[j];
      if (p <= eps) p = eps;
      if (p >= 1.0 - eps) p = 1.0 - eps;
      mu[k] = R::qnorm5(p, 0.0, 1.0, 1, 0);
      x[k] = x_full[j];
    }
    if (!valid) continue;
    
    std::vector<double> sgn(K);
    for (int i = 0; i < K; ++i) sgn[i] = (x[i] == 1) ? -1.0 : 1.0;
    
    NumericMatrix Sigma(K, K);
    for (int i = 0; i < K; ++i) {
      Sigma(i,i) = 1.0;
      const int ii = keep_idx[i];
      for (int j = i+1; j < K; ++j) {
        const int jj = keep_idx[j];
        int d = std::abs(cum[jj] - cum[ii]);
        if (d > max_d) d = max_d;
        const double rho = exp_table[d];
        Sigma(i,j) = rho;
        Sigma(j,i) = rho;
      }
    }
    
    for (int i = 0; i < K; ++i) {
      const double si = sgn[i];
      for (int j = 0; j < K; ++j) Sigma(i,j) *= si;
    }
    for (int j = 0; j < K; ++j) {
      const double sj = sgn[j];
      for (int i = 0; i < K; ++i) Sigma(i,j) *= sj;
    }
    
    NumericVector upper(K);
    for (int i = 0; i < K; ++i) upper[i] = (x[i] == 1) ? mu[i] : -mu[i];
    
    try {
      double logcdf = FAStCDF_algo1_cpp(
        upper, Sigma, true, 2.0, 1e-4, 20000
      );
      if (!R_finite(logcdf)) out[r] = NA_REAL;
      else out[r] = std::round(logcdf * 1000.0) / 1000.0;
    } catch (...) {
      out[r] = NA_REAL;
    }
  }
  
  return out;
}


// [[Rcpp::export]]
NumericVector compute_loglik_cpp_coord_list(CharacterVector fragment,
                                       List dist_sub_list,
                                       List M_list,
                                       List N_list,
                                       NumericVector exp_table,
                                       double eps = 1e-6,
                                       int print_every = 100000) {
  
  const int n = fragment.size();
  NumericVector out(n, NA_REAL);
  const int max_d = exp_table.size() - 1;
  
  for (int r = 0; r < n; ++r) {
    if (print_every > 0 && (r % print_every == 0)) {
      Rcpp::Rcout << "Processing row " << (r + 1) << " / " << n << std::endl;
    }
    if (fragment[r] == NA_STRING) continue;
    
    const std::string frag = as<std::string>(fragment[r]);
    const int L = (int)frag.size();
    
    IntegerVector Mv_r = M_list[r];
    IntegerVector Nv_r = N_list[r];
    if (Mv_r.size() != L || Nv_r.size() != L) continue;
    
    if (L == 1) {
      const char c = frag[0];
      if (!(c == 'C' || c == 'c' || c == 'T' || c == 't') || Nv_r[0] < 5) continue;
      if (Nv_r[0] == 0) continue;
      double p = (double)Mv_r[0] / (double)Nv_r[0];
      if (p <= eps) p = eps;
      if (p >= 1.0 - eps) p = 1.0 - eps;
      
      double logprob = (c == 'C' || c == 'c') ? std::log(p) : std::log(1.0 - p);
      out[r] = std::round(logprob * 1000.0) / 1000.0;
      continue;
    }
    
    SEXP dvec_sexp = dist_sub_list[r];
    if (dvec_sexp == R_NilValue) continue;
    
    IntegerVector dvec_r(dvec_sexp);
    if (dvec_r.size() != L - 1) continue;
    
    std::vector<int> keep_idx; keep_idx.reserve(L);
    std::vector<int> x_full(L, 0);
    for (int i = 0; i < L; ++i) {
      const char c = frag[i];
      if (c == 'C' || c == 'c') x_full[i] = +1;
      else if (c == 'T' || c == 't') x_full[i] = -1;
      else x_full[i] = 0;
      if (x_full[i] != 0 && Nv_r[i] >= 5) keep_idx.push_back(i);
    }
    const int K = (int)keep_idx.size();
    if (K == 0) continue;
    
    if (K == 1) {
      const int j = keep_idx[0];
      if (Nv_r[j] == 0) continue;
      double p = (double)Mv_r[j] / (double)Nv_r[j];
      if (p <= eps) p = eps;
      if (p >= 1.0 - eps) p = 1.0 - eps;
      double logprob = (x_full[j] == +1) ? std::log(p) : std::log(1.0 - p);
      out[r] = std::round(logprob * 1000.0) / 1000.0;
      continue;
    }
    
    std::vector<int> cum(L, 0);
    for (int i = 1; i < L; ++i) cum[i] = cum[i-1] + dvec_r[i-1];
    
    NumericVector mu(K);
    std::vector<int> x(K);
    bool valid = true;
    for (int k = 0; k < K; ++k) {
      const int j = keep_idx[k];
      if (Nv_r[j] == 0) { valid = false; break; }
      double p = (double)Mv_r[j] / (double)Nv_r[j];
      if (p <= eps) p = eps;
      if (p >= 1.0 - eps) p = 1.0 - eps;
      mu[k] = R::qnorm5(p, 0.0, 1.0, 1, 0);
      x[k]  = x_full[j];
    }
    if (!valid) continue;
    
    std::vector<double> sgn(K);
    for (int i = 0; i < K; ++i) sgn[i] = (x[i] == 1) ? -1.0 : 1.0;
    
    NumericMatrix Sigma(K, K);
    for (int i = 0; i < K; ++i) {
      Sigma(i,i) = 1.0;
      const int ii = keep_idx[i];
      for (int j = i+1; j < K; ++j) {
        const int jj = keep_idx[j];
        int d = std::abs(cum[jj] - cum[ii]);
        if (d > max_d) d = max_d;
        const double rho = exp_table[d];
        Sigma(i,j) = rho;
        Sigma(j,i) = rho;
      }
    }
    
    for (int i = 0; i < K; ++i) {
      const double si = sgn[i];
      for (int j = 0; j < K; ++j) Sigma(i,j) *= si;
    }
    for (int j = 0; j < K; ++j) {
      const double sj = sgn[j];
      for (int i = 0; i < K; ++i) Sigma(i,j) *= sj;
    }
    
    NumericVector upper(K);
    for (int i = 0; i < K; ++i) upper[i] = (x[i] == 1) ? mu[i] : -mu[i];
    
    try {
      double logcdf = FAStCDF_algo1_cpp(upper, Sigma, true, 2.0, 1e-5, 100000);
      if (!R_finite(logcdf)) out[r] = NA_REAL;
      else out[r] = std::round(logcdf * 1000.0) / 1000.0;
    } catch (...) {
      out[r] = NA_REAL;
    }
  }
  
  return out;
}


// [[Rcpp::export]]
CharacterVector extract_dist_sub(CharacterVector distance,
                                 IntegerVector startCpG,
                                 IntegerVector endCpG,
                                 IntegerVector markerStartCpG) {
  int n = distance.size();
  CharacterVector out(n);
  
  for (int i = 0; i < n; i++) {
    if (distance[i] == NA_STRING) {
      out[i] = NA_STRING;
      continue;
    }
    
    std::string d_str = as<std::string>(distance[i]);
    
    // split by comma
    std::vector<std::string> d_vec;
    std::stringstream ss(d_str);
    std::string token;
    while (std::getline(ss, token, ',')) {
      d_vec.push_back(token);
    }
    
    int s = startCpG[i] - markerStartCpG[i];  // 0-based
    int e = endCpG[i] - 1 - markerStartCpG[i]; // 0-based
    
    if (s < 0) s = 0;
    if (e >= (int)d_vec.size()) e = d_vec.size() - 1;
    
    if (s > e || d_vec.empty()) {
      out[i] = NA_STRING;
    } else {
      std::string result = d_vec[s];
      for (int j = s + 1; j <= e; j++) {
        result += "," + d_vec[j];
      }
      out[i] = result;
    }
  }
  
  return out;
}


// [[Rcpp::export]]
CharacterVector list_to_string(List x) {
  int n = x.size();
  CharacterVector out(n);
  
  for (int i = 0; i < n; ++i) {
    if (Rf_isNull(x[i])) {
      out[i] = NA_STRING;
      continue;
    }
    
    IntegerVector vec = x[i];
    int len = vec.size();
    if (len == 0) {
      out[i] = "";
      continue;
    }
    
    std::ostringstream oss;
    oss << vec[0];
    for (int j = 1; j < len; ++j) {
      oss << "," << vec[j];
    }
    out[i] = oss.str();
  }
  
  return out;
}


// [[Rcpp::export]]
CharacterVector mask_fragment_by_idx(CharacterVector fragment,
                                     IntegerVector startCpG,
                                     IntegerVector endCpG,
                                     IntegerVector idx_keep,
                                     int print_every = 100000,
                                     bool strict_length = true) {
  const int n = fragment.size();
  if (startCpG.size() != n || endCpG.size() != n) {
    stop("startCpG/endCpG length must match fragment length.");
  }
  
  // hash set for CpG to keep
  std::unordered_set<int> keep;
  keep.reserve(idx_keep.size() * 2 + 1);
  for (int v : idx_keep) keep.insert(v);
  
  CharacterVector out(n);
  
  for (int i = 0; i < n; ++i) {
    
    // progress printing
    if (print_every > 0 && (i + 1) % print_every == 0) {
      Rcpp::Rcout << "[mask_fragment] processed "
                  << (i + 1) << " / " << n << " rows\n";
    }
    
    if (fragment[i] == NA_STRING) {
      out[i] = NA_STRING;
      continue;
    }
    
    std::string s = as<std::string>(fragment[i]);
    const int L = (int)s.size();
    const int expected = endCpG[i] - startCpG[i] + 1;
    
    if (strict_length && expected != L) {
      out[i] = NA_STRING;
      continue;
    }
    
    const int useL = std::min(L, expected);
    int cpg = startCpG[i];
    
    for (int k = 0; k < useL; ++k, ++cpg) {
      if (keep.find(cpg) == keep.end()) s[k] = '.';
    }
    
    if (!strict_length && L > expected) {
      for (int k = expected; k < L; ++k) s[k] = '.';
    }
    
    out[i] = s;
  }
  
  // final newline for cleanliness
  Rcpp::Rcout << "[mask_fragment] done. total rows = " << n << "\n";
  
  return out;
}





// [[Rcpp::export]]
List process_fragments_by_region(CharacterVector chr, IntegerVector start_idx, 
                                 CharacterVector sequence, IntegerVector count,
                                 IntegerVector ref_old_idx, IntegerVector ref_new_idx,
                                 IntegerVector ref_region_id) {
  
  int n_fragments = chr.size();
  int n_ref = ref_old_idx.size();
  
  // build the map: old_idx -> (new_idx, region_id)
  std::unordered_map<int, int> old_to_new;
  std::unordered_map<int, int> old_to_region;
  std::unordered_set<int> ref_set;
  
  for(int i = 0; i < n_ref; i++) {
    int old_idx = ref_old_idx[i];
    ref_set.insert(old_idx);
    old_to_new[old_idx] = ref_new_idx[i];
    old_to_region[old_idx] = ref_region_id[i];
  }
  
  // result vectors
  std::vector<std::string> out_chr;
  std::vector<int> out_start;
  std::vector<std::string> out_seq;
  std::vector<int> out_count;
  std::vector<int> out_region;
  
  // process each fragment
  for(int i = 0; i < n_fragments; i++) {
    if(i % 500000 == 0) {
      Rcpp::checkUserInterrupt();
      if(i % 5000000 == 0 && i > 0) {
        Rcpp::Rcout << "Processed " << i << " fragments..." << std::endl;
      }
    }
    
    int frag_start = start_idx[i];
    std::string seq = Rcpp::as<std::string>(sequence[i]);
    int seq_len = seq.length();
    
    // find the first position present in ref, which fixes the region
    int target_region = -1;
    for(int j = 0; j < seq_len; j++) {
      int current_idx = frag_start + j;
      if(ref_set.count(current_idx) > 0) {
        target_region = old_to_region[current_idx];
        break;
      }
    }
    
    // if no position is present in ref, skip this fragment
    if(target_region == -1) continue;
    
    // keep only the CpGs that belong to target_region
    std::string new_seq = "";
    int new_start = -1;
    
    for(int j = 0; j < seq_len; j++) {
      int current_idx = frag_start + j;
      
      if(ref_set.count(current_idx) > 0 && old_to_region[current_idx] == target_region) {
        if(new_start == -1) {
          new_start = old_to_new[current_idx];
        }
        new_seq += seq[j];
      }
    }
    
    // if any CpG was kept, add it to the result
    if(new_seq.length() > 0) {
      out_chr.push_back(Rcpp::as<std::string>(chr[i]));
      out_start.push_back(new_start);
      out_seq.push_back(new_seq);
      out_count.push_back(count[i]);
      out_region.push_back(target_region);
    }
  }
  
  return List::create(
    Named("chr") = out_chr,
    Named("start_idx") = out_start,
    Named("sequence") = out_seq,
    Named("count") = out_count,
    Named("region_id") = out_region
  );
}

// [[Rcpp::export]]
List process_fragments_by_region_with_celltype(CharacterVector chr, IntegerVector start_idx, 
                                               CharacterVector sequence, IntegerVector count,
                                               CharacterVector celltype,
                                               IntegerVector ref_old_idx, IntegerVector ref_new_idx,
                                               IntegerVector ref_region_id) {
  
  int n_fragments = chr.size();
  int n_ref = ref_old_idx.size();
  
  std::unordered_map<int, int> old_to_new;
  std::unordered_map<int, int> old_to_region;
  std::unordered_set<int> ref_set;
  
  for (int i = 0; i < n_ref; i++) {
    int old_idx = ref_old_idx[i];
    ref_set.insert(old_idx);
    old_to_new[old_idx] = ref_new_idx[i];
    old_to_region[old_idx] = ref_region_id[i];
  }
  
  std::vector<std::string> out_chr;
  std::vector<int> out_start;
  std::vector<std::string> out_seq;
  std::vector<int> out_count;
  std::vector<int> out_region;
  std::vector<std::string> out_celltype;
  
  for (int i = 0; i < n_fragments; i++) {
    if (i % 500000 == 0) {
      Rcpp::checkUserInterrupt();
      if (i % 5000000 == 0 && i > 0) {
        Rcpp::Rcout << "Processed " << i << " fragments..." << std::endl;
      }
    }
    
    int frag_start = start_idx[i];
    std::string seq = Rcpp::as<std::string>(sequence[i]);
    int seq_len = seq.length();
    
    int target_region = -1;
    for (int j = 0; j < seq_len; j++) {
      int current_idx = frag_start + j;
      if (ref_set.count(current_idx) > 0) {
        target_region = old_to_region[current_idx];
        break;
      }
    }
    
    if (target_region == -1) continue;
    
    std::string new_seq = "";
    int new_start = -1;
    
    for (int j = 0; j < seq_len; j++) {
      int current_idx = frag_start + j;
      
      if (ref_set.count(current_idx) > 0 && old_to_region[current_idx] == target_region) {
        if (new_start == -1) {
          new_start = old_to_new[current_idx];
        }
        new_seq += seq[j];
      }
    }
    
    if (new_seq.length() > 0) {
      out_chr.push_back(Rcpp::as<std::string>(chr[i]));
      out_start.push_back(new_start);
      out_seq.push_back(new_seq);
      out_count.push_back(count[i]);
      out_region.push_back(target_region);
      out_celltype.push_back(Rcpp::as<std::string>(celltype[i]));
    }
  }
  
  return List::create(
    Named("chr") = out_chr,
    Named("start_idx") = out_start,
    Named("sequence") = out_seq,
    Named("count") = out_count,
    Named("region_id") = out_region,
    Named("celltype") = out_celltype
  );
}


// [[Rcpp::export]]
List prepare_fragment_data_cpp(IntegerVector start_idx, 
                               CharacterVector sequence,
                               IntegerMatrix beta,
                               IntegerVector ref_new_index,
                               IntegerVector ref_pos) {
  
  int n_frags = start_idx.size();
  int n_ref = ref_new_index.size();
  
  // map ref's new_index to pos
  std::unordered_map<int, int> index_to_pos;
  for(int i = 0; i < n_ref; i++) {
    index_to_pos[ref_new_index[i]] = ref_pos[i];
  }
  
  // output
  List M_list(n_frags);
  List N_list(n_frags);
  CharacterVector dist_sub(n_frags);
  
  for(int i = 0; i < n_frags; i++) {
    if(i % 500000 == 0) {
      Rcpp::checkUserInterrupt();
      Rcpp::Rcout << "Prepared " << i << " / " << n_frags << " fragments" << std::endl;
    }
    
    int start = start_idx[i];
    std::string seq = Rcpp::as<std::string>(sequence[i]);
    int frag_len = seq.length();
    
    // extract M and N
    IntegerVector M_vec(frag_len);
    IntegerVector N_vec(frag_len);
    
    for(int j = 0; j < frag_len; j++) {
      int idx = start + j - 1;  // R is 1-based, C++ is 0-based
      if(idx >= 0 && idx < beta.nrow()) {
        M_vec[j] = beta(idx, 0);
        N_vec[j] = beta(idx, 1);
      }
    }
    
    M_list[i] = M_vec;
    N_list[i] = N_vec;
    
    // compute distances
    if(frag_len > 1) {
      std::vector<int> positions;
      for(int j = 0; j < frag_len; j++) {
        int cpg_idx = start + j;
        if(index_to_pos.count(cpg_idx) > 0) {
          positions.push_back(index_to_pos[cpg_idx]);
        }
      }
      
      if(positions.size() == frag_len) {
        std::ostringstream oss;
        for(size_t j = 1; j < positions.size(); j++) {
          if(j > 1) oss << ",";
          oss << (positions[j] - positions[j-1]);
        }
        dist_sub[i] = oss.str();
      } else {
        dist_sub[i] = NA_STRING;
      }
    } else {
      dist_sub[i] = "";
    }
  }
  
  return List::create(
    Named("M_list") = M_list,
    Named("N_list") = N_list,
    Named("dist_sub") = dist_sub
  );
}



// [[Rcpp::export]]
CharacterVector mask_seq_batch(CharacterVector seq,
                               IntegerVector start_idx,
                               IntegerVector region_id,
                               List include_list) {
  int n = seq.size();
  CharacterVector out(n);
  
  for (int i = 0; i < n; i++) {
    if (seq[i] == NA_STRING || start_idx[i] == NA_INTEGER || region_id[i] == NA_INTEGER) {
      out[i] = seq[i];
      continue;
    }
    
    std::string s = Rcpp::as<std::string>(seq[i]);
    int L = (int)s.size();
    if (L == 0) { out[i] = seq[i]; continue; }
    
    std::string rid = std::to_string(region_id[i]);
    if (!include_list.containsElementNamed(rid.c_str())) {
      // missing region -> mask all
      for (int t = 0; t < L; t++) s[t] = '.';
      out[i] = s;
      continue;
    }
    
    LogicalVector inc = include_list[rid]; // length = max_new_index for this region
    int inc_len = inc.size();
    
    int st = start_idx[i];          // 1-based new_index start
    int s0 = st - 1;                // convert to 0-based offset
    for (int t = 0; t < L; t++) {
      int idx = s0 + t;             // 0-based index into inc
      bool keep = false;
      if (idx >= 0 && idx < inc_len) {
        keep = inc[idx] == TRUE;
      }
      if (!keep) s[t] = '.';
    }
    
    out[i] = s;
  }
  
  return out;
}


static inline bool is_finite_cpp(double x) {
  return std::isfinite(x);
}

static NumericVector normalize_simplex_cpp(const NumericVector& x) {
  int p = x.size();
  NumericVector out(p);
  double s = 0.0;
  
  for (int j = 0; j < p; ++j) {
    double v = x[j];
    if (!is_finite_cpp(v) || v < 0.0) v = 0.0;
    out[j] = v;
    s += v;
  }
  
  if (s <= 0.0) {
    for (int j = 0; j < p; ++j) out[j] = 1.0 / (double)p;
    return out;
  }
  
  for (int j = 0; j < p; ++j) out[j] /= s;
  return out;
}

static NumericVector safe_init_cpp(const NumericVector& x_init, int p) {
  if (x_init.size() == p) {
    return normalize_simplex_cpp(x_init);
  } else {
    NumericVector x(p);
    for (int j = 0; j < p; ++j) x[j] = 1.0 / (double)p;
    return x;
  }
}

static double obj_logmix_avg_cpp(const NumericMatrix& logL,
                                 const NumericVector& w,
                                 const NumericVector& x) {
  int n = logL.nrow();
  int p = logL.ncol();
  
  double out = 0.0;
  double wsum = 0.0;
  
  for (int i = 0; i < n; ++i) {
    double wi = w[i];
    if (!is_finite_cpp(wi) || wi <= 0.0) continue;
    wsum += wi;
    
    double mx = R_NegInf;
    for (int j = 0; j < p; ++j) {
      double lij = logL(i, j);
      if (!is_finite_cpp(lij) || x[j] <= 0.0) continue;
      double v = std::log(x[j]) + lij;
      if (v > mx) mx = v;
    }
    
    if (!is_finite_cpp(mx)) continue;
    
    double s = 0.0;
    for (int j = 0; j < p; ++j) {
      double lij = logL(i, j);
      if (!is_finite_cpp(lij) || x[j] <= 0.0) continue;
      s += std::exp(std::log(x[j]) + lij - mx);
    }
    
    if (s > 0.0) out += wi * (mx + std::log(s));
  }
  
  if (wsum <= 0.0) return R_NegInf;
  return out / wsum;
}

static NumericVector grad_logmix_avg_cpp(const NumericMatrix& logL,
                                         const NumericVector& w,
                                         const NumericVector& x) {
  int n = logL.nrow();
  int p = logL.ncol();
  
  NumericVector g(p);
  std::fill(g.begin(), g.end(), 0.0);
  
  double wsum = 0.0;
  
  for (int i = 0; i < n; ++i) {
    double wi = w[i];
    if (!is_finite_cpp(wi) || wi <= 0.0) continue;
    wsum += wi;
    
    double mx = R_NegInf;
    for (int j = 0; j < p; ++j) {
      double lij = logL(i, j);
      if (!is_finite_cpp(lij) || x[j] <= 0.0) continue;
      double v = std::log(x[j]) + lij;
      if (v > mx) mx = v;
    }
    
    if (!is_finite_cpp(mx)) continue;
    
    double denom = 0.0;
    for (int j = 0; j < p; ++j) {
      double lij = logL(i, j);
      if (!is_finite_cpp(lij) || x[j] <= 0.0) continue;
      denom += std::exp(std::log(x[j]) + lij - mx);
    }
    
    if (denom <= 0.0) continue;
    
    for (int j = 0; j < p; ++j) {
      double lij = logL(i, j);
      if (!is_finite_cpp(lij) || x[j] <= 0.0) continue;
      
      double num = std::exp(std::log(x[j]) + lij - mx);
      double rij = num / denom;
      g[j] += wi * (rij / x[j]);
    }
  }
  
  if (wsum > 0.0) {
    for (int j = 0; j < p; ++j) g[j] /= wsum;
  }
  
  return g;
}

// group-wise Dirichlet / log-barrier on q_g = (x_j + eps)/(s_g + m eps)
// P_g = - sum_j log q_gj
//     = - sum_j log(x_j + eps) + m log(s_g + m eps)
static List penalty_groupDirichlet_q_cpp(const NumericVector& x,
                                         const List& groups,
                                         const LogicalVector& penalize_group,
                                         double eps = 1e-12) {
  int p = x.size();
  NumericVector grad(p);
  std::fill(grad.begin(), grad.end(), 0.0);
  
  double P = 0.0;
  
  for (int g = 0; g < groups.size(); ++g) {
    if (g < penalize_group.size() && !penalize_group[g]) continue;
    
    IntegerVector idx = groups[g];
    int m = idx.size();
    if (m <= 1) continue;
    
    double s = 0.0;
    for (int k = 0; k < m; ++k) {
      int j = idx[k] - 1;
      s += x[j];
    }
    
    double denom = s + m * eps;
    
    for (int k = 0; k < m; ++k) {
      int j = idx[k] - 1;
      P += -std::log(x[j] + eps) + std::log(denom);
      grad[j] += -1.0 / (x[j] + eps) + ((double)m / denom);
    }
  }
  
  return List::create(
    _["penalty"] = P,
    _["grad"] = grad
  );
}

static NumericVector mirror_update_cpp(const NumericVector& x,
                                       const NumericVector& grad,
                                       double step) {
  int p = x.size();
  NumericVector out(p);
  
  double mx = R_NegInf;
  for (int j = 0; j < p; ++j) {
    double v = std::log(std::max(x[j], 1e-300)) + step * grad[j];
    out[j] = v;
    if (v > mx) mx = v;
  }
  
  double s = 0.0;
  for (int j = 0; j < p; ++j) {
    out[j] = std::exp(out[j] - mx);
    s += out[j];
  }
  
  if (s <= 0.0) {
    for (int j = 0; j < p; ++j) out[j] = 1.0 / (double)p;
    return out;
  }
  
  for (int j = 0; j < p; ++j) out[j] /= s;
  return out;
}

// [[Rcpp::export]]
List mix_groupDirichlet_penalized_cpp(const NumericMatrix& logL,
                                      const NumericVector& w,
                                      const List& groups,
                                      const LogicalVector& penalize_group,
                                      const NumericVector& x_init,
                                      double lambda = 0.0,
                                      double eps = 1e-12,
                                      int maxit = 30,
                                      double tol = 1e-6,
                                      double step0 = 0.5,
                                      double backtrack = 0.5,
                                      int bt_max = 20,
                                      bool verbose = false) {
  NumericVector x = safe_init_cpp(x_init, logL.ncol());
  
  List pen0 = penalty_groupDirichlet_q_cpp(x, groups, penalize_group, eps);
  double f0 = obj_logmix_avg_cpp(logL, w, x);
  double P0 = as<double>(pen0["penalty"]);
  double obj_prev = f0 - lambda * P0;
  
  bool converged = false;
  int iter = 0;
  
  for (iter = 1; iter <= maxit; ++iter) {
    NumericVector g_ll = grad_logmix_avg_cpp(logL, w, x);
    List pen = penalty_groupDirichlet_q_cpp(x, groups, penalize_group, eps);
    NumericVector g_pen = pen["grad"];
    
    NumericVector g_tot(x.size());
    for (int j = 0; j < x.size(); ++j) {
      g_tot[j] = g_ll[j] - lambda * g_pen[j];
    }
    
    double step = step0;
    NumericVector x_new(x.size());
    double obj_new = R_NegInf;
    double f_new = R_NegInf;
    double P_new = R_NegInf;
    bool accepted = false;
    
    for (int bt = 0; bt < bt_max; ++bt) {
      x_new = mirror_update_cpp(x, g_tot, step);
      
      List pen_new = penalty_groupDirichlet_q_cpp(x_new, groups, penalize_group, eps);
      f_new = obj_logmix_avg_cpp(logL, w, x_new);
      P_new = as<double>(pen_new["penalty"]);
      obj_new = f_new - lambda * P_new;
      
      if (is_finite_cpp(obj_new) && obj_new >= obj_prev) {
        accepted = true;
        break;
      }
      
      step *= backtrack;
    }
    
    if (!accepted) {
      if (verbose) {
        Rcout << "iter " << iter << "  no accepted step, stopping." << std::endl;
      }
      break;
    }
    
    double rel = std::fabs(obj_new - obj_prev) / (std::fabs(obj_prev) + 1e-12);
    
    if (verbose) {
      Rcout << "iter " << iter
            << "  obj=" << obj_new
            << "  avg_loglik=" << f_new
            << "  penalty=" << P_new
            << "  step=" << step
            << "  rel=" << rel
            << std::endl;
    }
    
    x = x_new;
    obj_prev = obj_new;
    
    if (rel < tol) {
      converged = true;
      break;
    }
  }
  
  List penF = penalty_groupDirichlet_q_cpp(x, groups, penalize_group, eps);
  double fF = obj_logmix_avg_cpp(logL, w, x);
  double PF = as<double>(penF["penalty"]);
  
  return List::create(
    _["x"] = x,
    _["objective"] = fF - lambda * PF,
    _["avg_loglik_part"] = fF,
    _["penalty"] = PF,
    _["lambda"] = lambda,
    _["converged"] = converged,
    _["iter"] = iter
  );
}