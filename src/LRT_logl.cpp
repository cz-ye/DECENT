// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <Rcpp.h>
#include <cmath>
#include <Rmath.h>
#include <limits>
#include <RcppNumerical.h>
using namespace Numer;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::IntegerVector;
typedef Eigen::Map<Eigen::MatrixXd> MapMat;
typedef Eigen::Map<Eigen::VectorXd> MapVec;

// Adapted from https://github.com/myajima/bppkgx/blob/master/bppkgx/src/utilarmadillo.cpp
double dbetabinom_cpp(double x, int n, double a, double b) {
  //assert(a > 0 && b > 0);
  //assert(n > 0 && n > x && x >=0);
  if(a <= 0) Rf_error("a must be > 0");
  if(b <= 0) Rf_error("b must be > 0");
  if(n <= 0) Rf_error("n must be > 0");
  if(x <  0) Rf_error("x must be >= 0");
  if(n < x) return 0;
  return( Rf_choose(n, x) * Rf_beta(x + a, n - x + b) / Rf_beta(a, b ) );
}

class OptimLRTCpp: public MFuncGrad
{
private:
  const IntegerVector ct;
  const NumericVector sf;
  const IntegerVector y;
  const NumericVector z;
  const NumericMatrix DO_par;
  double loglIBBCpp(NumericVector p);
  NumericVector loglIBBGradCpp(NumericVector p);
  NumericVector zinbGradCpp(double y, NumericVector p,  double sf, int ct);
  NumericVector zinbGrad0Cpp(NumericVector p,  double sf, int ct);

public:
  OptimLRTCpp(const IntegerVector ct_, const NumericVector sf_, const IntegerVector y_, const NumericVector z_,
              const NumericMatrix DO_par_) : ct(ct_), sf(sf_), y(y_), z(z_), DO_par(DO_par_) {}

  double f_grad(Constvec& par, Refvec grad) {
    NumericVector p(par.size());
    for (int i = 0; i < par.size(); i++) {
      p[i] = par[i];
    }
    grad = Rcpp::as<MapVec>(loglIBBGradCpp(p)); //
    return(loglIBBCpp(p));
  }
};

// [[Rcpp::export]]
Rcpp::List optimLRTCpp(const NumericVector p, const IntegerVector ct, const NumericVector sf, const IntegerVector y,
                       const NumericVector z, const NumericMatrix DO_par) {
  OptimLRTCpp f(ct, sf, y, z, DO_par);
  Eigen::VectorXd par(p.size());
  for (int i = 0; i < p.size(); i++) {
    par[i] = p[i];
  }
  double fopt;
  int status = optim_lbfgs(f, par, fopt);
  return Rcpp::List::create(
    Rcpp::Named("par") = par,
    Rcpp::Named("fopt") = fopt,
    Rcpp::Named("status") = status
  );
}

// Cost function (Negative log-likelihood)
double OptimLRTCpp::loglIBBCpp(NumericVector p) {
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
      PZ[i] += dbetabinom_cpp(z[i], y[j], 4*prob, 4*(1-prob)) *
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

NumericVector OptimLRTCpp::zinbGradCpp(double y, NumericVector p,  double sf, int ct) {
  int nct = p.size() - 2;
  NumericVector dmu(nct); NumericVector mu(nct);
  if (nct == 1) {
    dmu[0] = 1;
    mu[0] = exp(p[1]);
  } else {
    dmu[0] = 1; dmu[1] = exp(p[2]);
    mu = exp(p[1]) * dmu;
  }
  double theta = exp(-p[nct+1]);

  ct = ct - 1; // change of indexing in cpp

  NumericVector grad_mat(nct+2);
  grad_mat[0] = 1 / (1 + exp(-p[0]));
  grad_mat[1] = -mu[0]*(y/mu[0] - (y + theta)*dmu[ct]*sf/(theta + sf*mu[ct]));
  if (nct != 1) {
    grad_mat[2] = -dmu[ct] * (ct > 0) * (y/dmu[ct] - (y + theta)*sf*mu[0]/(theta + sf*mu[ct]));
  }
  grad_mat[nct+1] = theta * (Rf_digamma(y + theta) - Rf_digamma(theta) + log(theta) + 1 -
    log(sf*mu[ct] + theta) - (y + theta)/(sf*mu[ct] + theta) );
  return(grad_mat);
}

NumericVector OptimLRTCpp::zinbGrad0Cpp(NumericVector p,  double sf, int ct) {
  int nct = p.size() - 2;
  NumericVector dmu(nct); NumericVector mu(nct);
  if (nct == 1) {
    dmu[0] = 1;
    mu[0] = exp(p[1]); //* dmu
  } else {
    dmu[0] = 1; dmu[1] = exp(p[2]);
    mu = exp(p[1]) * dmu;
  }
  double pi0 = 1 / (1 + exp(-p[0]));
  double theta = exp(-p[nct+1]);
  ct = ct - 1; // change of indexing in cpp

  NumericVector grad_mat(nct+2);
  grad_mat[0]  = -(1 - pow(theta/(sf*mu[ct] + theta), theta)) /
    (pi0 + (1 - pi0) * pow(theta/(sf*mu[ct] + theta), theta)) * pi0 * (1 - pi0);
  grad_mat[1] = -mu[0] * (-dmu[ct]*sf*(1 - pi0) * pow(theta/(sf*mu[ct] + theta), (theta + 1)) /
    (pi0 + (1-pi0) * pow(theta/(sf*mu[ct] + theta), theta)));
  if (nct != 1) {
    grad_mat[2] = -dmu[ct]*(ct > 0) *
      (sf*mu[0] *(-(1 - pi0) * pow(theta/(sf*mu[ct] + theta), (theta + 1)) /
        (pi0 + (1 - pi0)*pow(theta/(sf*mu[ct] + theta), theta))));
  }
  grad_mat[nct+1] = theta * ((log(theta) + 1 - log(sf*mu[ct]+theta) - theta/(sf*mu[ct]+theta)) /
    (1 + pi0/(1 - pi0) * pow((sf*mu[ct] + theta)/theta, theta)));
  return(grad_mat);
}

// Gradient
NumericVector OptimLRTCpp::loglIBBGradCpp(NumericVector p) {
  int i, j, k;
  int n = ct.size();
  int m = y.size();
  int nct = max(ct);
  double pi0 = exp(p[0])/(1 + exp(p[0]));
  double size = exp(-p[p.size()-1]);
  double mu;

  NumericVector grad_mat(nct+2);
  NumericVector grad_vec(nct+2);
  double prob, DO_prob, NB_prob;
  double PZ, f0;
  for (i = 0; i < n; i++) {
    PZ = 0;
    for (k = 0; k < nct+2; k++) { // initialize every iteration
      grad_mat[k] = 0;
    }
    mu = exp(p[1] + p[ct[i]] * (ct[i] > 1)) * sf[i];
    f0 = pi0 + (1-pi0)*dnbinom_mu(0, size, mu, 0);
    for (j = 0; j < m; j++) {
      prob = exp(DO_par(i, 0) + DO_par(i, 1)*log(y[j]+1)) / ( 1 + exp(DO_par(i, 0) + DO_par(i, 1)*log(y[j]+1)) );
      DO_prob = dbetabinom_cpp(z[i], y[j], 4*prob, 4*(1-prob));
      NB_prob = dnbinom_mu(y[j], size, mu, 0) * (1-pi0);
      PZ +=  DO_prob * NB_prob ;
      for (k = 0; k < nct+2; k++) {
        grad_mat[k] += DO_prob * zinbGradCpp(y[j], p, sf[i], ct[i])[k] * NB_prob;
      }
    }
    if (z[i] == 0){
      PZ += f0;
      grad_mat = grad_mat + f0 * zinbGrad0Cpp(p, sf[i], ct[i]);
    }
    if (PZ == 0){
      PZ = std::numeric_limits<double>::min();
    }
    for (k = 0; k < nct+2; k++) {
      grad_vec[k] += grad_mat[k]/PZ;
    }
  }
  return(grad_vec);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
*/
