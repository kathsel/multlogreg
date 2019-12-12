context("Multivariate student t density calculation")

test_that("the density matches mvtnorm::dmvt", {
  set.seed(123)
  sigma1 <- bayesm::rwishart(10, diag(8))$IW
  mean1  <- rnorm(8)
  X1     <- mvtnorm::rmvt(1, delta = mean1, sigma = sigma1, df = 8)
  sigma2 <- bayesm::rwishart(10, diag(9))$IW
  mean2  <- rnorm(9)
  X2     <- mvtnorm::rmvt(1, delta = mean2, sigma = sigma2, df = 7.3)
  expect_equal(mvtnorm::dmvt(X1, mean1, sigma1, df = 8, log = FALSE),
               dmvt(X1, mean1, sigma1, 8, FALSE))
  expect_equal(mvtnorm::dmvt(X1, mean1, sigma1, df = 8, log = TRUE),
               dmvt(X1, mean1, sigma1, 8, TRUE))
  expect_equal(mvtnorm::dmvt(X2, mean2, sigma2, df = 7.3, log = FALSE),
               dmvt(X2, mean2, sigma2, 7.3, FALSE))
  expect_equal(mvtnorm::dmvt(X2, mean2, sigma2, df = 7.3, log = TRUE),
               dmvt(X2, mean2, sigma2, 7.3, TRUE))
})

test_that("the density does not accept unequal dimensions", {
  expect_error(dmvt(c(0, 1), 0, matrix(1), 8, FALSE))
  expect_error(dmvt(c(0, 1), c(0, 0), matrix(1, 0), 8, FALSE))
  expect_error(dmvt(0, c(0, 0), rbind(c(1, 0), c(0, 1)), 8, FALSE))
})

test_that("the density does not accept non-symmetric or non-positive definite covariance matrices", {
  expect_error(dmvt(c(0, 0), c(0, 0), rbind(c(1, -0.5), c(0, 1)), 8, FALSE))
  expect_error(dmvt(c(0, 0), c(0, 0), rbind(c(1, -0.5), c(0, -1)), 8, FALSE))
})

