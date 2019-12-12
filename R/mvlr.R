#' Fit the multivariate logistic regression model with Gibbs sampling
#'
#' @param coeff named list of coefficients for each outcome,
#'   intercepts are included by default
#' @param data data.frame containing the outcomes and coefficients
#' @param iter integer, the number of iterartions
#' @param betainit intial values for the coefficients
#' @param rinit initial values for the unique entries of the scale matrix
#' @param phiinit initial values for phis
#' @param mubeta prior mean on beta
#' @param sigbeta prior covariance matrix for beta
#' @param omega variance-covariance matrix for sampling r (Metropolis step)
#' @param burn_in number of burn in samples to use for complete loop
#' @param thinning thinning used
#' @param burn_in_truncated number of burn in samples to use for sampling from
#'   the truncated normal distribution
#'
#' @return A list with the samples of the coefficients,
#'   correlation coefficients, weights and the number of accepted values
#' @export
#'
mvlr <- function(coeff, data, iter,
                 beta0, r0, phi0,
                 mubeta, sigbeta, omega,
                 burn_in, thinning, burn_in_truncated = 10) {
  X_list <-
    lapply(names(coeff), function(y) {
      frml <- formula(paste(y, "~", paste(coeff[[y]], collapse = " + ")))
      model.matrix(frml, data = data)
    })

  # Extract dimensions
  q <- length(coeff)
  n <- nrow(data)
  mj <- sapply(X_list, function(x) dim(x)[2])
  m <- sum(mj)

  # Initialize array
  X_array <- array(0, c(q, m, n))

  # Position for starting parameter of each outcome
  pars_m <- cumsum(mj) - mj + 1

  # Fill array with a q x m matrix for each observation j
  for (j in 1:n) {
    for (i in 1:q) {
      X_array[i, pars_m[i]:(pars_m[i] + mj[i] - 1), j] <- X_list[[i]][j, ]
    }
  }

  # Extract outcome columns
  y <- as.matrix(data[, names(coeff)])

  cmvlr(X_array, y, iter,
        beta0, r0, phi0,
        mubeta, sigbeta, omega,
        burn_in, thinning, burn_in_truncated)
}