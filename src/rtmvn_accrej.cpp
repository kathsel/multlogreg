#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Generate samples from a truncated multivariate normal distribution
//'
//' Generates a n x q matrix with n samples of a truncated multivariate
//' normal distribution with dimension q using the acceptance-rejection method
//' sampling from a multivariate normal distribution.
//'
//' For the a non-diagonal matrix D, the linear restriction lower <= Dx <= upper
//' hold.
//'
//' @param n number of samples to be generated
//' @param mean mean vector (q x 1)
//' @param sigma covariance matrix (q x q)
//' @param D matrix of linear constraints
//' @param lower vector of lower bounds
//' @param upper vector of upper bounds
//' @param nproposals number of propsal samples per iteration
//' @param verbose should information about sampling be printed?
// [[Rcpp::export]]
arma::mat rtmvn_accrej(int n,
                       const arma::vec& mean,
                       const arma::mat& sigma,
                       const arma::mat& D,
                       const arma::vec& lower,
                       const arma::vec& upper,
                       int nproposals = 0,
                       bool verbose = false) {

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

  if(nproposals > 0 && nproposals < n)
    throw Rcpp::exception("Number of proposals should be at least n");

  // q = Dimension
  int q = mean.n_elem;
  // Samples (n x q)
  arma::mat Y(q, n, arma::fill::zeros);

  // Remaining number of samples to generate
  int rem = n;

  // Number of accepted samples
  int total = 0;
  double alpha = 1.0;
  int noacc = 0;

  // Number of iterations
  int iter = 0;
  int totalproposals = 0;

  // Sample from the multivariate normal distribution and accept those within region
  while(rem > 0) {
    Rcpp::checkUserInterrupt();

    // If the number of proposals is not specified upfront, use 10x number of
    // samples needed in the first run and then adapt to the acceptance rate
    if (nproposals == 0) {
      // Generate n/alpha samples, for very small alpha generate maximum of 100000 proposals
      // but 10 times the number of needed samples
      if (rem/alpha > 100000) {
        nproposals = std::min(10*rem, 100000);
      } else {
        nproposals = std::ceil(std::max(rem/alpha, 10.0));
      }
    }

    arma::mat X = mvnrnd(mean, sigma, nproposals);
    arma::mat Xtransf = D * X;
    int acc = 0;

    for (int i = 0; i < nproposals && total < n; i++) {
      Rcpp::checkUserInterrupt();

      if (all(Xtransf.col(i) >= lower && Xtransf.col(i) <= upper)) {
        Y.col(total) = X.col(i);
        total++;
        acc++;
      }
    }

    if (acc == 0 && total != n) {
      noacc++;
      if (verbose && noacc >= 1000) {
        Rcout << "No samples were accepted for 1000 runs in a row. " <<
          "Number of total propsals: " << totalproposals << "\n";
      }
      // Start over since rem and total have not changed
      iter++;
      totalproposals += nproposals;
      if (totalproposals > pow(10, 6) * nproposals) {
        Rcout << "Acceptance probability is to low (" <<
          total / (pow(10, 6) * nproposals) <<
            ")\nNumber of total propsals: " << totalproposals << "\n" <<
              "Number of accepted values: " << total << "\n";
        return(Y);
      }
      continue;
    }

    if (total != n) {
      alpha = static_cast<double>(acc) / nproposals;
      if (alpha < 0.01) {
        Rcout << "Acceptance rate is very low (" << alpha << ")" << "\n";
      }
    }

    iter++;
    totalproposals += nproposals;
    rem = n - total;
  }

  if (verbose) {
    Rcout << iter << " iterartions with " << totalproposals << " proposals in total" << "\n";
  }

  return(Y);
}