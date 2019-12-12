#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
#include "dmvnrm.h"
#include "dmvt.h"
#include "dmvlogis.h"
#include "rtmvn_gibbs.h"

using namespace Rcpp;

//' Fit the multivariate logistic regression model with Gibbs sampling
//'
//' @param x numeric array containing the covariates for each sample
//' @param y numeric matrix containing the outcomes of interest
//' @param iter integer, the number of iterartions
//' @param beta0 intial values for the coefficients
//' @param r0 initial values for the unique entries of the scale matrix
//' @param phi0 initial values for phis
//' @param mubeta prior mean on beta
//' @param sigbeta prior covariance matrix for beta
//' @param omega variance-covariance matrix for sampling r (Metropolis step)
//' @param burn_in number of burn in samples to use for complete loop
//' @param thinning thinning used
//' @param burn_in_truncated number of burn in samples to use for sampling from
//'   the truncated normal distribution
//'
// [[Rcpp::export]]
List cmvlr(const arma::cube& x,
           const arma::mat& y,
           int iter,
           const arma::vec& beta0,
           const arma::vec& r0,
           const arma::vec& phi0,
           const arma::vec& mubeta,
           const arma::mat& sigbeta,
           const arma::mat& omega,
           int burn_in,
           int thinning,
           int burn_in_truncated = 10) {

  //////////////////////////////////////////////////////////////////////////////
  // CHECK DIMENSIONS OF INPUTS
  //////////////////////////////////////////////////////////////////////////////

  // Number of observations
  // x.n_slices, y.n_rows
  if(y.n_rows != x.n_slices)
    throw Rcpp::exception("Number of rows of y and number of slices of x must match");
  // Number of coefficients
  // x.n_cols, beta0.n_elem, mubeta.n_elem, sigbeta.n_rows, sigbeta.n_cols
  if(x.n_cols != beta0.n_elem || beta0.n_elem != mubeta.n_elem || mubeta.n_elem != sigbeta.n_rows)
    throw Rcpp::exception("Number of columns of x, elements in beta0, mubeta and sigbeta must match");
  // Number of unique correlation coefficients
  //r0.n_elem, omega.n_rows, omega.n_cols, y.n_cols * (y.n_cols - 1) / 2
  if(r0.n_elem != omega.n_rows || omega.n_rows != y.n_cols * (y.n_cols - 1) / 2)
    throw Rcpp::exception("Number of elements in r0, omega, and y must match");
  // Sigbeta, Omega square symmetric positive definte matrix
  if(sigbeta.n_rows != sigbeta.n_cols || !sigbeta.is_sympd())
    throw Rcpp::exception("sigbeta is not a square symmetric positive definite matrix");
  if(omega.n_rows != omega.n_cols || !omega.is_sympd())
    throw Rcpp::exception("sigmabeta or omega are not symmetric positive definite");
  // Check iter, burn_in, thinning and burn_in_truncated
  if(iter <= 0 || burn_in <= 0 || thinning <= 0 || burn_in_truncated <= 0)
    throw Rcpp::exception("iter, burn_in, thinning and burn_in_truncated have to be larger than 0");
  if(any(r0 > 1 || r0 < -1)) {
    throw Rcpp::exception("Absolute values of r0 need to be smaller than 1");
  }

  //////////////////////////////////////////////////////////////////////////////
  // DECLARE AND INITIALIZE VARIABLES AND PARAMETERS
  //////////////////////////////////////////////////////////////////////////////

  // degrees of freedom for the (multivariate) t and logistic distribution
  const double df = 7.3;
  const double sig = std::pow(M_PI, 2.0) * (df - 2.0) / (3.0 * df);

  // Dimensions, number of parameters
  int n = x.n_slices; // observations
  int m = x.n_cols;   // coefficients
  int q = y.n_cols;   // outcomes
  int qstar = q * (q - 1) / 2.0; // unique entries in scale matrix R

  // Objects to return
  //// Racc: total number of acceptations in Metropolis step for updating R
  int Racc = 0;
  //// beta and r: realisations of beta and r after burn-in and thinning
  arma::mat beta(m, iter, arma::fill::zeros);
  arma::mat r(qstar, iter, arma::fill::zeros);
  //// w: weights for correcting the posterior distribution
  arma::vec w(iter, arma::fill::zeros);

  // Auxillary, temporary objects, not returned ...
  // ... for updating z
  //// z: latent variables updated in each iteration,
  ////    samples from the truncated multivariate normal
  arma::mat z(q, n, arma::fill::zeros);
  //// mu: mean matrix updated in each iteration
  arma::mat mu(q, n, arma::fill::zeros);
  //// initialized mu
  for(int i = 0; i < n ; ++i){
    mu.col(i) = x.slice(i) * beta0;
  }
  //// lower and upper bounds for truncated normal
  arma::mat ub(q, n, arma::fill::zeros);
  arma::mat lb(q, n, arma::fill::zeros);

  arma::umat bounds1 = find(y.t() == 1);
  lb.elem(bounds1).fill(0.0);
  ub.elem(bounds1).fill(R_PosInf);

  arma::umat bounds0 = find(y.t() == 0);
  lb.elem(bounds0).fill(R_NegInf);
  ub.elem(bounds0).fill(0.0);

  // ... for updating phi
  //// phi: temporary value of phi, initialized with phi0
  arma::vec phi(phi0);
  //// shape and scale parameter for the gamma distribution
  const double shape = (df + q) / 2.0;
  double scale = 0.0;

  // ... for updating beta
  //// beta_t: temporary value of beta, initialized with beta0
  arma::vec beta_t(beta0);
  //// vectors, matrices for proposal distribution of beta
  arma::mat M1(m, m, arma::fill::zeros);
  arma::vec M2(m, arma::fill::zeros);
  arma::mat sigbetainv(arma::inv_sympd(sigbeta));
  arma::mat sigbeta_t(m, m, arma::fill::zeros);
  arma::vec mubeta_t(m, arma::fill::zeros);

  // ... for updating r
  //// r_t: temporary value of r, initialized with r0
  arma::vec r_t(r0);
  arma::vec rprop(qstar, arma::fill::zeros);   // proposal values
  arma::mat MatRprop(q, q, arma::fill::zeros); // proposal matrix
  double num = 0.0;
  double dem = 0.0;
  //// aux: indices for elements in an upper tringular matrix of size q x q
  arma::mat AUX1(q, q, arma::fill::zeros);
  arma::mat AUX2(q, q, arma::fill::zeros);
  arma::vec indices = arma::linspace<arma::vec>(1, q, q);
  AUX1.each_col() += indices;
  AUX2.each_row() += trans(indices);
  arma::uvec aux = find(AUX2 > AUX1);
  //// MatR: Temporary correlation matrix filled with r_t
  ////      (r_t uniquely determine the correlation matrix)
  //// Fill MatR's upper triangle with the values of r_t (column by column)
  //// using aux and create symmetric matrix with ones on the diagonal
  arma::mat MatR(q, q, arma::fill::zeros);
  MatR.elem(aux) = r_t;
  MatR = arma::symmatu(MatR);
  MatR.diag().ones();
  //// MatRinv: inverse of temporary correlation matrix MatR
  arma::mat MatRinv(arma::inv_sympd(MatR));

  // ... for updating w
  double w_num = 0.0;
  double w_denum = 0.0;

  //////////////////////////////////////////////////////////////////////////////
  // OUTER LOOP
  //////////////////////////////////////////////////////////////////////////////

  for(int t = - (burn_in - 1); t <= (thinning * iter); t++) {
    Rcpp::checkUserInterrupt();

    // Update latent variables z (truncated multivariate normal distribution)

    for(int i = 0; i < n; i++) {
      Rcpp::checkUserInterrupt();

      z.col(i) =
        rtmvn_gibbs(1, mu.col(i), sig * MatR / phi(i), arma::eye(q, q),
                    lb.col(i), ub.col(i), z.col(i), burn_in_truncated, 1);
    }

    // Update phis (gamma distribution)

    for(int i = 0; i < n; i++) {
      Rcpp::checkUserInterrupt();

      scale = 2.0 / (df +
        1/sig * arma::as_scalar((z.col(i) - mu.col(i)).t() * MatRinv * (z.col(i) - mu.col(i))));

      phi(i) = R::rgamma(shape, scale);
    }

    // Update beta (multivariate normal distribution)

    M1.zeros();
    M2.zeros();

    for(int i = 0; i < n; i++) {
      Rcpp::checkUserInterrupt();

      M1 += phi(i) * x.slice(i).t() * MatRinv * x.slice(i);
      M2 += phi(i) * x.slice(i).t() * MatRinv * z.col(i);
    }

    sigbeta_t = arma::inv_sympd(sigbetainv + (1 / sig) * M1);
    mubeta_t  = sigbeta_t * (sigbetainv * mubeta + (1 / sig) * M2);

    beta_t = arma::mvnrnd(mubeta_t, sigbeta_t);

    // Update mu

    for(int i = 0; i < n ; i++) {
      Rcpp::checkUserInterrupt();

      mu.col(i) = x.slice(i) * beta_t;
    }

    // Update r (MatR) (Metropolis step)

    // Proposals from multivariate normal distribution
    rprop = arma::mvnrnd(r_t, omega);

    // Create proposal matrix
    MatRprop.elem(aux) = rprop;
    MatRprop = symmatu(MatRprop);
    MatRprop.diag().ones();

    // If MatRprop matrix is a correlation matrix, accept it with some probability
    // according to the ratio of the values of the complete conditional density for R
    if(all(rprop <= 1) && MatRprop.is_sympd()) {

      num = 0.0;
      dem = 0.0;

      for(int i = 0; i < n; i++) {
        Rcpp::checkUserInterrupt();

        // Sum of the log-densities, use exp then to transform to the desired product
        num += arma::as_scalar(dmvnrm(z.col(i), mu.col(i), sig / phi(i) * MatRprop, true));
        dem += arma::as_scalar(dmvnrm(z.col(i), mu.col(i), sig / phi(i) * MatR, true));
      }

      if(runif(1)[0] < exp(num - dem)){
        r_t = rprop;
        MatR = MatRprop;
        MatRinv = arma::inv_sympd(MatR);
        Racc++;
      }
    }

    //  Apply "burn-in" and "thinning"

    if(t > 0) {
      if(Racc == 0 && t == 1) {
        Rcpp::warning("No scale matrix was accepted during burn-in. Choose a different omega.");
      }
      if(fmod(t, thinning) == 0) {

        beta.col((t - thinning) / thinning) = beta_t;
        r.col((t - thinning) / thinning) = r_t;

        // Calculate weights
        w_num = 0.0;
        w_denum = 0.0;
        for (int j = 0; j < n; j++) {
          Rcpp::checkUserInterrupt();

          w_num += dmvlogis(z.col(j), mu.col(j), MatR, df, true);
          w_denum += dmvt(z.col(j), mu.col(j), sig * MatR, df, true);
        }
        w((t - thinning) / thinning) = exp(w_num - w_denum);
      }
    }
  }
  List result = List::create(_["beta"] = beta,
                             _["r"] = r,
                             _["w"] = w,
                             _["Racc"] = Racc);

  return result;
}
