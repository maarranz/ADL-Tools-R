test_that("order + first_lag builds the lag sets", {
  l <- normalize_lags("y", c("x1", "x2"), order = c(2, 1))
  expect_equal(l, list(y = 1:2, x1 = 0:1, x2 = 0:1))

  l <- normalize_lags("y", c("x1", "x2"), order = c(1, 2, 3), first_lag = c(x2 = 1))
  expect_equal(l, list(y = 1L, x1 = 0:2, x2 = 1:3))

  l <- normalize_lags("y", "x", order = c(0, 0))
  expect_equal(l, list(y = integer(0), x = 0L))
})

test_that("explicit lags allow gaps and any order", {
  l <- normalize_lags("y", c("x1", "x2"),
                      lags = list(x2 = c(4, 1), y = c(1, 4), x1 = 0))
  expect_equal(l, list(y = c(1L, 4L), x1 = 0L, x2 = c(1L, 4L)))
  l <- normalize_lags("y", "x", lags = list(x = 1:2))
  expect_equal(l$y, integer(0))
})

test_that("invalid lag specifications are rejected", {
  expect_error(normalize_lags("y", "x", order = c(1, 0), first_lag = 1), "below")
  expect_error(normalize_lags("y", c("a", "b"), order = c(1, 1, 1, 1)), "length 2 or 3")
  expect_error(normalize_lags("y", "x", lags = list(y = 0, x = 0)), "1 or more")
  expect_error(normalize_lags("y", "x", lags = list(y = 1)), "No lags")
  expect_error(normalize_lags("y", "x", lags = list(x = 0, w = 1)), "unknown")
  expect_error(normalize_lags("y", "x", order = c(1.5, 1)), "whole")
  expect_error(normalize_lags("y", "x"), "either")
})

test_that("design matrix matches dynlm on the same sample", {
  skip_if_not_installed("dynlm")
  suppressPackageStartupMessages(library(dynlm))  # dynlm needs zoo attached
  set.seed(2)
  n <- 60
  x1 <- cumsum(rnorm(n)); x2 <- rnorm(n)
  y <- stats::filter(0.3 * x1 + 0.2 * x2 + rnorm(n), 0.5, method = "recursive")
  m <- stats::ts(cbind(y = as.numeric(y), x1 = x1, x2 = x2), start = c(1990, 1),
                 frequency = 4)
  d <- as_adl_data(m, c("y", "x1", "x2"))
  lags <- normalize_lags("y", c("x1", "x2"), order = c(2, 1, 2), first_lag = c(x2 = 1))
  des <- build_design(d$z, lags)

  expect_equal(des$max_lag, 2L)
  expect_equal(nrow(des$frame), n - 2)
  expect_equal(names(des$frame),
               c("y", "y_lag1", "y_lag2", "x1", "x1_lag1", "x2_lag1", "x2_lag2"))

  ours <- stats::lm(y ~ ., data = des$frame)
  ref <- dynlm::dynlm(y ~ L(y, 1:2) + L(x1, 0:1) + L(x2, 1:2), data = m)
  expect_equal(unname(stats::coef(ours)), unname(stats::coef(ref)))
  expect_equal(stats::nobs(ours), stats::nobs(ref))
})

test_that("max_lag fixes a common sample; trend and unlagged terms", {
  m <- stats::ts(cbind(y = 1:20, x = (1:20)^2, d = rep(0:1, 10)))
  d <- as_adl_data(m, c("y", "x", "d"))
  lags <- normalize_lags("y", "x", order = c(1, 0))
  des <- build_design(d$z, lags, unlagged = "d", trend = TRUE, max_lag = 4)
  expect_equal(nrow(des$frame), 16)
  expect_equal(des$frame$y[1], 5)
  expect_equal(des$frame$y_lag1[1], 4)
  expect_equal(des$frame$trend, 1:16)
  expect_equal(des$terms$role, c("ar", "dl", "unlagged", "trend"))
  expect_error(build_design(d$z, lags, max_lag = 0), "below")
})

test_that("column-name clashes are rejected (audit repro A)", {
  set.seed(7)
  d <- data.frame(x = rnorm(100), y = rnorm(100))
  a <- transform(d, x_lag1 = seq_len(nrow(d)))
  expect_error(suppressMessages(adl(y ~ L(x, 1) + x_lag1, a)),
               "'x_lag1' would be both lag 1 of 'x' and the variable 'x_lag1'")
  expect_error(suppressMessages(adl(y ~ x | x_lag1, a, order = c(0, 1))), "clash")
  # no clash when the lag is not generated
  expect_silent(suppressMessages(adl(y ~ L(x, 2) + x_lag1, a)))
})

test_that("lag values are validated before coercion", {
  expect_error(normalize_lags("y", "x", lags = list(x = 1.9)), "whole")
  expect_error(normalize_lags("y", "x", lags = list(x = Inf)), "whole")
  m <- stats::ts(cbind(y = rnorm(30), x = rnorm(30)))
  d <- as_adl_data(m, c("y", "x"))
  expect_error(build_design(d$z, list(y = 1L, x = 0L), max_lag = c(2, 3)), "single")
})
