#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
#include "dmvt.h"

using namespace Rcpp;

//' Calculate the multivariate logistic density of a vector
//'
//' @param x numeric vector containing the random numbers
//' @param mean numeric vector
//' @param R numeric scale matrix with 1's on the diagonal
//' @param df degrees of freedom
//' @param logd boolean indicating whether the log-transformed density
//'   should be returned
// [[Rcpp::export]]
double dmvlogis(const arma::vec& x,
                const arma::vec& mean,
                const arma::mat& R,
                double df = 7.3,
                bool logd = false) {

  if((x.n_elem != mean.n_elem) || (x.n_elem != R.n_rows))
    throw Rcpp::exception("Dimensions of x, mean and R must match");

  if(R.n_rows != R.n_cols)
    throw Rcpp::exception("R is not a square matrix");

  if(!R.is_sympd())
    throw Rcpp::exception("R is not symmetric positive definite");

  if(!all(R.diag() == 1))
    throw Rcpp::exception("R should only have 1's on the diagonal");

  int q = x.n_elem;
  double dens = 0.0;

  NumericVector qe =
    as<NumericVector>(wrap(arma::exp(x - mean)/(1 + arma::exp(x - mean))));
  NumericVector e = qt(qe, df);
  double mvtdens = dmvt(e, arma::zeros(q), R, df, true);
  double logisdenssum = 0.0;
  for (int i = 0; i < q; i++) {
    logisdenssum += R::dlogis(x(i), mean(i), 1.0, true);
  }

  dens = mvtdens + logisdenssum - sum(Rcpp::dt(e, df, true));

  if (!logd) {
    dens = exp(dens);
  }

  return(dens);
}