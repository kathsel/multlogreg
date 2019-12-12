#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "rtuvn.h"

using namespace Rcpp;


// Implementation adapted from:
// https://github.com/suchitm/tmvn

//' Gibbs sampler for the Truncated Multivariate Normal Distribution
//'
//' Random vector generation for the truncated multivariate normal distribution
//'     using a Gibbs sampler.
//'
//' @param n number of samples to be generated
//' @param mean mean vector
//' @param sigma covariance matrix
//' @param D matrix of linear constraints
//' @param lower vector of lower bounds
//' @param upper vector of upper bounds
//' @param init vector of initial values for the Gibbs sampler.
//' @param burn_in number of burn in samples to use.
//' @param thinning thinning parameter
//'
//' @return matrix of samples with each column being an idependent sample.
//'
//' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
//'     truncated multivariate normal and student-t distributions subject to
//'     linear inequality constraints. Journal of Statistical Theory and
//'     Practice, 9(4), 712-732.
//' @export
// [[Rcpp::export]]
arma::mat rtmvn_gibbs(int n,
                      const arma::vec& mean,
                      const arma::mat& sigma,
                      const arma::mat& D,
                      const arma::vec& lower,
                      const arma::vec& upper,
                      const arma::vec& init,
                      int burn_in = 10,
                      int thinning = 1) {

  // Check of input parameters
  if(mean.n_elem != sigma.n_rows)
    throw Rcpp::exception("Dimensions of mean and sigma must match");

  if(sigma.n_rows != sigma.n_cols)
    throw Rcpp::exception("sigma is not a square matrix");

  if(!sigma.is_sympd())
    throw Rcpp::exception("sigma is not symmetric positive definite");

  if(D.n_cols != mean.n_elem)
    throw Rcpp::exception("Number of columns in D and dimesion of mean must match");

  if(lower.n_elem != upper.n_elem || lower.n_elem != D.n_rows)
    throw Rcpp::exception("Dimension of lower, upper and D must match");

  if(n <= 0 || burn_in <= 0 || thinning <= 0)
    throw Rcpp::exception("n, burn_in, thinning have to be larger than 0");

  // q = Dimension
  int q = sigma.n_cols;
  // Samples (q x n)
  arma::mat x(q, n, arma::fill::zeros);
  arma::vec temp(q, arma::fill::zeros);

  if (q == 1) {

    x.row(0) = rtuvn(n, mean(0), sigma(0, 0), lower(0), upper(0)).t();

  } else {

    arma::vec lb = lower - D * mean;
    arma::vec ub = upper - D * mean;
    arma::mat sigmachol = chol(sigma, "lower");
    arma::mat R = D * sigmachol;
    int p = R.n_cols;

    arma::vec z = inv(trimatl(sigmachol)) * (init - mean);

    double lower_pos, upper_pos, lower_neg, upper_neg;

    for(int i = - burn_in; i < n * thinning; i++) {

      for(int j = 0; j < p; j++) {

        arma::vec rj = R.col(j);

        arma::mat MatRj = R;
        MatRj.shed_col(j);

        arma::vec zj = z;
        zj.shed_row(j);

        arma::vec Rzj = MatRj * zj;

        arma::vec a = lb - Rzj;
        arma::vec b = ub - Rzj;

        arma::uvec pos = find(rj > 0);
        arma::uvec neg = find(rj < 0);

        if(pos.n_rows == 0) {

          lower_pos = R_NegInf;
          upper_pos = R_PosInf;

        } else {

          temp = a(pos) / rj(pos);
          lower_pos = arma::max(temp);

          temp = b(pos) / rj(pos);
          upper_pos = arma::min(temp);
        }

        if(neg.n_rows == 0) {

          lower_neg = R_NegInf;
          upper_neg = R_PosInf;

        } else {

          temp = a(neg) / rj(neg);
          upper_neg = min(temp);

          temp = b(neg) / rj(neg);
          lower_neg = max(temp);
        }

        double lower_j = std::max(lower_pos, lower_neg);
        double upper_j = std::min(upper_pos, upper_neg);

        // univariate sample
        if (fabs(lower_j - upper_j) < 0.0001) {
          z(j) = (lower_j - upper_j)/2.0;
        } else {
          z(j) = rtuvn_single(0, 1, lower_j, upper_j);
        }
      }

      if(i >= 0) {
        if (fmod(i, thinning) == 0) {
          x.col(i) = sigmachol * z + mean;
        }
      }
    }
  }

  return(x);

}