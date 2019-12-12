#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

// Adapted from:
// https://github.com/duckmayr/RcppDist
// and
// https://gallery.rcpp.org/articles/dmvnorm_arma/

//' Calculate the multivariate Student t density of a vector
//'
//' @param x numeric vector containing the random numbers
//' @param mean numeric vector
//' @param sigma numeric covariance matrix
//' @param df degrees of freedom
//' @param logd boolean indicating whether the log-transformed density
//'   should be returned
// [[Rcpp::export]]
double dmvt(const arma::vec& x,
            const arma::vec& mean,
            const arma::mat& sigma,
            double df,
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
  double constants = R::lgammafn((df + q) * 0.5) - R::lgammafn(df * 0.5) -
    (q/2.0) * std::log(df * M_PI);

  arma::vec z = Qinv * (x - mean);
  dens = constants - (df + q) * 0.5 * std::log(1 + dot(z, z)/df) + sumlogQinv;

  if (!logd) {
    dens = exp(dens);
  }

  return(dens);
}