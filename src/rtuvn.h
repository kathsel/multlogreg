#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
#ifndef rtuvn_H
#define rtuvn_H

arma::vec rtuvn(int n,
                double mean,
                double sd,
                double lower,
                double upper);

double rtuvn_single(double mean,
                    double sd,
                    double lower,
                    double upper);

#endif