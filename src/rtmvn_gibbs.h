#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
#include "rtuvn.h"
#ifndef rtmvn_gibbs_H
#define rtmvn_gibbs_H

arma::mat rtmvn_gibbs(int n,
                      const arma::vec& mean,
                      const arma::mat& sigma,
                      const arma::mat& D,
                      const arma::vec& lower,
                      const arma::vec& upper,
                      const arma::vec& init,
                      int burn_in = 10,
                      int thinning = 1);

#endif