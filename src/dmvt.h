#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
#ifndef dmvt_H
#define dmvt_H

double dmvt(const arma::vec& x,
            const arma::vec& mean,
            const arma::mat& sigma,
            double df,
            bool logd = false);

#endif