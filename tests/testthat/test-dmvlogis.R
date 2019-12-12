context("Multivariate logistic density calculation")

test_that("the univariate density matches dlogis", {
  set.seed(123)
  X1    <- rlogis(1)
  mean2 <- runif(1, 0, 20)
  X2    <- rlogis(1, location = mean2)
  expect_equal(stats::dlogis(X1, 0, 1, log = FALSE),
               dmvlogis(X1, 0, matrix(1), df = 7.3, FALSE))
  expect_equal(stats::dlogis(X1, 0, 1, log = TRUE),
               dmvlogis(X1, 0, matrix(1), df = 7.3, TRUE))
  expect_equal(stats::dlogis(X2, mean2, 1, log = FALSE),
               dmvlogis(X2, mean2, matrix(1), df = 7.3, FALSE))
  expect_equal(stats::dlogis(X2, mean2, 1, log = TRUE),
               dmvlogis(X2, mean2, matrix(1), df = 7.3, TRUE))
})

test_that("the univariate density is independent of df", {
  set.seed(478)
  mu  <- runif(1, 0, 20)
  X   <- rlogis(1, location = mu)
  dfs <- sample(1:20, 4)
  expect_equal(dmvlogis(X, mu, matrix(1), df = dfs[[1]], FALSE),
               dmvlogis(X, mu, matrix(1), df = dfs[[2]], FALSE))
  expect_equal(dmvlogis(X, mu, matrix(1), df = dfs[[1]], FALSE),
               dmvlogis(X, mu, matrix(1), df = dfs[[3]], FALSE))
  expect_equal(dmvlogis(X, mu, matrix(1), df = dfs[[1]], FALSE),
               dmvlogis(X, mu, matrix(1), df = dfs[[4]], FALSE))
})

test_that("the multivariate density matches mvtnorm::dmvt, dlogis, dt", {
  set.seed(239)
  # Use a multivariate student t random sample to test the density
  # as we don't have a sampler for the multivariate logistic
  sigma <- bayesm::rwishart(10, diag(8))$IW
  R     <- stats::cov2cor(sigma)
  mu    <- rnorm(8)
  X     <- mvtnorm::rmvt(1, delta = mu, sigma = R, df = 8)

  g       <- stats::qt(exp(X - mu)/(1 + exp(X - mu)), df = 7.3)
  mvtdens <- mvtnorm::dmvt(g, rep(0, 8), R, df = 7.3, log = FALSE)
  ldens   <- mapply(dlogis, X, mu)
  tdens   <- stats::dt(g, df = 7.3)

  expect_equal(mvtdens * prod(ldens / tdens),
               dmvlogis(X, mu, R, df = 7.3, FALSE))
  expect_equal(log(mvtdens) + sum(log(ldens) - log(tdens)),
               dmvlogis(X, mu, R, df = 7.3, TRUE))
})

test_that("the density does not accept unequal dimensions", {
  expect_error(dmvlogis(c(0, 1), 0, matrix(1), 7.3, FALSE))
  expect_error(dmvlogis(c(0, 1), c(0, 0), matrix(1, 0), 7.3, FALSE))
  expect_error(dmvlogis(0, c(0, 0), rbind(c(1, 0), c(0, 1)), 7.3, FALSE))
})

test_that("the density does not accept non-symmetric or non-positive definite correlation matrices", {
  expect_error(dmvlogis(c(0, 0), c(0, 0), rbind(c(1, -0.5), c(0, 1)), 7.3, FALSE))
  expect_error(dmvlogis(c(0, 0), c(0, 0), rbind(c(1, -0.5), c(0, -1)), 7.3, FALSE))
})

test_that("the density does not accept covariance matrices", {
  expect_error(dmvlogis(c(0, 0), c(0, 0), rbind(c(2, 0.5), c(0.5, 2)), 7.3, FALSE))
})
