#include <Rcpp.h>
#include <cmath>
#include <Rmath.h>
#include <limits>
using namespace Rcpp;

// from https://github.com/myajima/bppkgx/blob/master/bppkgx/src/utilarmadillo.cpp
double dbetabinom_cpp(double x, int n, double a, double b, int give_log) {
  //assert(a > 0 && b > 0);
  //assert(n > 0 && n > x && x >=0);
  if(a <= 0) Rf_error("a must be > 0");
  if(b <= 0) Rf_error("b must be > 0");
  if(n <= 0) Rf_error("n must be > 0");
  if(x <  0) Rf_error("x must be >= 0");
  if(n < x) return 0;
  return( give_log?
            Rf_lchoose(n, x) + Rf_lbeta(x + a, n - x + b) - Rf_lbeta(a, b) :
            Rf_choose(n, x) * Rf_beta(x + a, n - x + b) / Rf_beta(a, b ) );
}

// [[Rcpp::export]]
double test(NumericVector p, IntegerVector ct, NumericVector sf, IntegerVector y, NumericVector z, NumericMatrix DO_par){
  int i = 0;
//  int j = 0;
  double pi0 = exp(p[0])/(1 + exp(p[0]));
  double size = exp(-p[p.size()-1]);
  return (pi0 + (1-pi0)*dnbinom_mu(0, size, exp(p[1] + p[ct[i]] * (ct[i] > 1)) * sf[i], 0));
}

// [[Rcpp::export]]
double loglIBBCpp(NumericVector p, IntegerVector ct, NumericVector sf, IntegerVector y, NumericVector z, NumericMatrix DO_par) {
  int i, j;
  int n = ct.size();
  int m = y.size();
  double mu;
  double pi0 = exp(p[0])/(1 + exp(p[0]));
  double size = exp( -p[p.size()-1] );
  
  NumericVector PZ(n);
  double prob;
  for (i = 0; i < n; i++) {
    PZ[i] = 0;
    mu = exp(p[1] + p[ct[i]] * (ct[i] > 1)) * sf[i];
    for (j = 0; j < m; j++) {
      prob = exp(DO_par(i, 0) + DO_par(i, 1)*log(y[j]+1)) / ( 1 + exp(DO_par(i, 0) + DO_par(i, 1)*log(y[j]+1)) );
      PZ[i] += dbetabinom_cpp(z[i], y[j], 4*prob, 4*(1-prob), 0) * 
               dnbinom_mu(y[j], size, mu, 0) * (1-pi0);
    }
    if (z[i] == 0){
      PZ[i] += ( pi0 + (1-pi0)*dnbinom_mu(0, size, mu, 0) );
    }
    if (PZ[i] == 0){
      PZ[i] = std::numeric_limits<double>::min();
    }
  }
  return( -sum(log(PZ)) );
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
*/
