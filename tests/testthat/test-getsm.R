test_that("the selection refits to the gets coefficients on the GUM's sample", {
  skip_if_not_installed("gets")
  gum <- adl(y ~ x + w, data = adl_sim, order = c(4, 4, 4))
  s <- adl_getsm(gum)
  expect_s3_class(s$fit, "adl")
  expect_equal(unname(coef(s$fit)), unname(coef(s$gets)))
  expect_equal(nobs(s$fit), nobs(gum))
  expect_equal(s$fit$index, gum$index)
  expect_true(all(c("(Intercept)", "x", "w") %in% s$selected))   # never tested
  expect_setequal(c(s$selected, s$dropped), names(coef(gum)))
  out <- capture.output(print(s))
  expect_true(any(grepl("conditional on the search", out)))
})

test_that("single terms: a lag can be kept without the ones before it", {
  skip_if_not_installed("gets")
  set.seed(61)
  n <- 300
  x <- rnorm(n)
  y <- numeric(n)
  e <- rnorm(n, sd = 0.5)
  for (t in 4:n) y[t] <- 0.5 * y[t - 1] + 0.8 * x[t - 3] + e[t]
  d <- stats::ts(cbind(y = y, x = x))
  gum <- adl(y ~ x, data = d, order = c(4, 4))
  s <- adl_getsm(gum, contemporaneous = "search")
  expect_equal(s$fit$lags$x, 3L)
  expect_equal(s$fit$lags$y, 1L)
  # the lag-by-lag search of ADL3 cannot drop lags 1 and 2 alone
  g <- adl_gets(y ~ x, data = d, max_order = 4, contemporaneous = "exclude")
  expect_true(all(1:3 %in% g$fit$lags$x))
})

test_that("kept terms, unlagged terms and the trend are never dropped", {
  skip_if_not_installed("gets")
  set.seed(62)
  n <- 200
  d <- stats::ts(cbind(y = rnorm(n), x = rnorm(n), dum = rbinom(n, 1, 0.1)))
  gum <- adl(y ~ x | dum, data = d, order = c(3, 3), trend = TRUE)
  s <- adl_getsm(gum, keep = "y_lag3")
  expect_true(all(c("dum", "trend", "x", "y_lag3") %in% s$selected))
  expect_true(s$fit$trend)
  expect_equal(s$fit$unlagged, "dum")
  expect_error(adl_getsm(gum, keep = "z_lag1"), "not in the model")
  expect_error(adl_getsm(gum, diagnostics = NA), "TRUE or FALSE")
})

test_that("a search that removes every regressor leaves no ADL", {
  skip_if_not_installed("gets")
  set.seed(63)
  d <- stats::ts(cbind(y = rnorm(150), x = rnorm(150)))
  s <- adl_getsm(adl(y ~ x, data = d, order = c(2, 2)), contemporaneous = "search",
                 alpha = 0.001)
  expect_null(s$fit)
  expect_true(any(grepl("removed every regressor", capture.output(print(s)))))
})
