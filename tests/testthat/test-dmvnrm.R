context("Multivariate normal density calculation")

test_that("the density matches mvtnorm::dmvnorm", {
  set.seed(123)
  sigma1 <- bayesm::rwishart(10, diag(8))$IW
  mean1  <- rnorm(8)
  X1     <- mvtnorm::rmvnorm(1, mean1, sigma1)
  sigma2 <- bayesm::rwishart(10, diag(9))$IW
  mean2  <- rnorm(9)
  X2     <- mvtnorm::rmvnorm(1, mean2, sigma2)
  expect_equal(mvtnorm::dmvnorm(X1, mean1, sigma1, log = FALSE),
               dmvnrm(X1, mean1, sigma1, FALSE))
  expect_equal(mvtnorm::dmvnorm(X1, mean1, sigma1, log = TRUE),
               dmvnrm(X1, mean1, sigma1, TRUE))
  expect_equal(mvtnorm::dmvnorm(X2, mean2, sigma2, log = FALSE),
               dmvnrm(X2, mean2, sigma2, FALSE))
  expect_equal(mvtnorm::dmvnorm(X2, mean2, sigma2, log = TRUE),
               dmvnrm(X2, mean2, sigma2, TRUE))
})

test_that("the density does not accept unequal dimensions", {
  expect_error(dmvnrm(c(0, 1), 0, matrix(1), FALSE))
  expect_error(dmvnrm(c(0, 1), c(0, 0), matrix(1, 0), FALSE))
  expect_error(dmvnrm(0, c(0, 0), rbind(c(1, 0), c(0, 1)), FALSE))
})

test_that("the density does not accept non-symmetric or non-positive definite covariance matrices", {
  expect_error(dmvnrm(c(0, 0), c(0, 0), rbind(c(1, -0.5), c(0, 1)), FALSE))
  expect_error(dmvnrm(c(0, 0), c(0, 0), rbind(c(1, -0.5), c(0, -1)), FALSE))
})

