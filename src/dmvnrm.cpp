#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

// Adapted from:
// https://gallery.rcpp.org/articles/dmvnorm_arma/

const double log2pi = std::log(2.0 * M_PI);

//' Calculate the multivariate normal density of a vector
//'
//' @param x numeric vector containing the random numbers
//' @param mean numeric vector
//' @param sigma numeric covariance matrix
//' @param logd boolean indicating whether the log-transformed density
//'   should be returned
// [[Rcpp::export]]
double dmvnrm(const arma::vec& x,
              const arma::vec& mean,
              const arma::mat& sigma,
              bool logd = false) {

  if((x.n_elem != mean.n_elem) || (x.n_elem != sigma.n_rows))
    throw Rcpp::exception("Dimensions of x, mean and sigma must match");

  if(sigma.n_rows != sigma.n_cols)
    throw Rcpp::exception("sigma is not a square matrix");

  if(!sigma.is_sympd())
    throw Rcpp::exception("sigma is not symmetric positive definite");

  int q = x.n_elem;
  double dens;

  arma::mat Qinv = (inv(trimatu(chol(sigma)))).t();
  double sumlogQinv = arma::sum(log(Qinv.diag()));
  double constants = - (q/2.0) * log2pi;

  arma::vec z = Qinv * (x - mean) ;
  dens = constants - 0.5 * dot(z, z) + sumlogQinv;

  if (!logd) {
    dens = exp(dens);
  }

  return(dens);
}

