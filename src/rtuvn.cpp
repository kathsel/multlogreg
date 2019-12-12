#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Implements the ensemble rejection algorithm for the univariate truncated
// normal distribution proposed by
// Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
//   truncated multivariate normal and student-t distributions subject to
//   linear inequality constraints. Journal of Statistical Theory and
//   Practice, 9(4), 712-732.
//
// Implementation adapted from:
// https://github.com/suchitm/tmvn

////////////////////////////////////////////////////////////////////////////////
// Rejection sampling helper functions
////////////////////////////////////////////////////////////////////////////////

// Normal rejection sampling
double norm_rej(double a, double b) {
  // first sample from a random normal
  double x = R::rnorm(0, 1);

  // Accept only x in (a, b)
  while( (x < a) || (x > b) ) {
    x = R::rnorm(0, 1);
  }
  return(x);
}

// Half-normal rejection sampling
double halfnorm_rej(double a, double b) {
  double x = R::rnorm(0, 1);
  double abs_x = fabs(x);

  // Accept only |x| in (a, b)
  while( (abs_x < a) || (abs_x > b) ) {
    x = R::rnorm(0, 1);
    abs_x = fabs(x);
  }
  return(abs_x);
}

// Uniform rejection sampling
double unif_rej(double a, double b) {
  while(true) {
    double x = R::runif(a, b);
    double u = R::runif(0, 1);
    double rho = 0.0;

    // cases for the ratio
    if( (a <= 0) && (0 <= b) ) {
      rho = exp(-0.5 * (x*x));
    }
    if (a > 0) {
      rho = exp(-0.5 * (x*x - a*a));
    }
    if (b < 0) {
      rho = exp(-0.5 * (x*x - b*b));
    }

    // accept step
    if(u <= rho) {
      return(x);
    }
  }
}

// Translated-exponential rejection sampling
double exp_rej(double a, double b) {
  double lambda = (a + sqrt(a*a + 4.0)) / 2.0;
  double x = R::rweibull(1, 1.0/lambda) + a;
  double u = R::runif(0, 1);
  double rho = exp(-0.5 * pow(x - lambda, 2.0));

  // loop to generate the sample
  while((u > rho) || (x > b)) {
    x = R::rweibull(1, 1.0/lambda) + a;
    u = R::runif(0, 1);
    rho = exp(-0.5 * pow(x - lambda, 2.0));
  }

  return(x);
}

////////////////////////////////////////////////////////////////////////////////
// Sample cases for different upper and lower boundaries
////////////////////////////////////////////////////////////////////////////////

// Sample Case 1: [a, \infty)
double sample_case1(double a, double b) {
  double samp;

  if(a <= 0) {
    samp = norm_rej(a, b);
  } else if (a < 0.25696) {
    samp = halfnorm_rej(a, b);
  } else {
    samp = exp_rej(a, b);
  }
  return(samp);
}

// Sample Case 2: [a, b], a <= 0, b >= 0
double sample_case2(double a, double b) {
  double samp;
  double btest = a + sqrt(2.0 * M_PI);

  if(b > btest) {
    samp = norm_rej(a, b);
  } else {
    samp = unif_rej(a, b);
  }
  return(samp);
}

// Sample Case 3: [a, b], a > 0
double sample_case3(double a, double b) {
  double samp;
  double btest = 0;

  if(a < 0.25696) {
    btest = a + sqrt(M_PI / 2.0) * exp(a*a / 2.0);
    if(b <= btest) {
      samp = unif_rej(a, b);
    } else {
      samp = halfnorm_rej(a, b);
    }
  } else {
    btest = a + (2.0 / (a + sqrt(a*a + 4.0))) *
      exp(0.5 + (a*a - a * sqrt(a*a + 4.0)) / 4.0);
    if(b <= btest) {
      samp = unif_rej(a, b);
    } else{
      samp = exp_rej(a, b);
    }
  }
  return(samp);
}

////////////////////////////////////////////////////////////////////////////////
// Sample from the truncate univariate standard normal distribution
////////////////////////////////////////////////////////////////////////////////

// Generate a sample from a truncated univariate standard normal distribution,
// TN(0, 1; a, b).
double sample_tuvsn(double a, double b) {
  double samp;

  if( (a == R_NegInf) || (b == R_PosInf) ) {

    if(b == R_PosInf) {
      samp = sample_case1(a, b);
    } else {
      // Case 4 is symmetric to Case 1
      double temp = sample_case1(-b, -a);
      samp = -temp;
    }

  } else {

    if(a >= 0) {
      samp = sample_case3(a, b);
    } else if (b <= 0) {
      // Case 5 is symmetric to Case 3
      double temp = sample_case3(-b, -a);
      samp = -temp;
    } else{
      samp = sample_case2(a, b);
    }

  }
  return(samp);
}

////////////////////////////////////////////////////////////////////////////////
// Sample from the truncate univariate normal distribution with mean and sd
////////////////////////////////////////////////////////////////////////////////

//' Univariate Truncated Normal Distribution
//'
//' Generates n samples from the truncated normal distribution
//'  TN(mean, sd; lower, upper).
//'
//' @param n number of samples
//' @param mean mean
//' @param sd standard deviation
//' @param lower lower bound
//' @param upper upper bound
//' @return vector of n samples
//'
//' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
//'     truncated multivariate normal and student-t distributions subject to
//'     linear inequality constraints. Journal of Statistical Theory and
//'     Practice, 9(4), 712-732.
//'
//' @examples
//' rtuvn(n = 1, mean = 10, sd = 20, lower = 10, upper = 20)
//' rtuvn(n = 1, mean = 10, sd = 20, lower = 10, upper = Inf)
//'
//' @export
// [[Rcpp::export]]
arma::vec rtuvn(int n,
                double mean,
                double sd,
                double lower,
                double upper) {

  if(sd <= 0)
    throw Rcpp::exception(
        "Error in rtuvn: standard deviation sd has to be larger than zero");

  // transform the boundaries
  double a = 0.0;
  double b = 0.0;

  if(lower > upper) {
    Rcpp::warning(
      "rtuvn_single: lower bound %d larger than upper bound %d\n => switching bounds",
      lower, upper);
    a = (upper - mean) / sd;
    b = (lower - mean) / sd;
  } else {
    a = (lower - mean) / sd;
    b = (upper - mean) / sd;
  }

  // generate sample from TN(0, 1; a, b)
  arma::vec Z(n);
  for(int i = 0; i < n; i++) {
    Z(i) = sample_tuvsn(a, b);
  }

  // transform the data back
  Z = sd * Z + mean;
  return(Z);
}

//' Univariate Truncated Normal Distribution - single sample
//'
//' Generate one sample from the truncated normal distribution
//'  TN(mean, sd; lower, upper).
// [[Rcpp::export]]
double rtuvn_single(double mean,
                    double sd,
                    double lower,
                    double upper) {


  if(sd <= 0)
    throw Rcpp::exception(
        "Error in rtuvn_single: standard deviation sd has to be larger than zero");

  // transform the boundaries
  double a = 0.0;
  double b = 0.0;

  if(lower > upper) {
    Rcpp::warning(
      "rtuvn_single: lower bound %d larger than upper bound %d\n => switching bounds",
      lower, upper);
    a = (upper - mean) / sd;
    b = (lower - mean) / sd;
  } else {
    a = (lower - mean) / sd;
    b = (upper - mean) / sd;
  }

  // generate sample from TN(0, 1; a, b)
  double Z = sample_tuvsn(a, b);

  // transform the data back
  Z = sd * Z + mean;
  return(Z);
}
