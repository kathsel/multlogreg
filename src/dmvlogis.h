#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
#ifndef dmvlogis_H
#define dmvlogis_H

double dmvlogis(const arma::vec& x,
                const arma::vec& mean,
                const arma::mat& R,
                double df,
                bool logd = false);

#endif