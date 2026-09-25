// ep_core.h -- the EP multivariate-normal orthant probability.
//
// Expectation propagation for the multivariate normal CDF, from the EPmvnCDF
// package (Francesco Denti and Augusto Fasano; MIT licence, see
// inst/COPYRIGHTS), in the C++ form used by the original implementation of this
// method (validation/reference/fragment_lik.cpp).  Kept identical to that code so
// the likelihood is numerically identical to the original: do not reformat, any
// change here changes the numbers.
#ifndef SPATCORR_EP_CORE_H
#define SPATCORR_EP_CORE_H

#include <RcppArmadillo.h>
#include "A_aux.h"
using namespace Rcpp;

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

#endif
