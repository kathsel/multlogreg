#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
#ifndef dmvnrm_H
#define dmvnrm_H

double dmvnrm(const arma::vec& x,
              const arma::vec& mean,
              const arma::mat& sigma,
              bool logd = false);

#endif