sim_data <- function(n = 80) {
  set.seed(3)
  x1 <- cumsum(rnorm(n)); x2 <- rnorm(n); d <- as.numeric(seq_len(n) == 40)
  y <- stats::filter(0.4 * x1 + 0.3 * x2 + 2 * d + rnorm(n), 0.5,
                     method = "recursive")
  stats::ts(cbind(y = as.numeric(y), x1 = x1, x2 = x2, d = d),
            start = c(2000, 1), frequency = 4)
}

test_that("the three entry routes give the same model", {
  m <- sim_data()
  a <- adl(y ~ x1 + x2, data = m, order = c(2, 1))
  b <- adl(y ~ x1 + x2, data = m, lags = list(y = 1:2, x1 = 0:1, x2 = 0:1))
  c <- adl(y ~ L(y, 1:2) + L(x1, 0:1) + L(x2, 0:1), data = m)
  expect_equal(coef(a), coef(b))
  expect_equal(coef(a), coef(c))
  expect_equal(a$label, "ADL(2, 1, 1)")
  expect_equal(a$dynlm, "y ~ L(y, 1:2) + L(x1, 0:1) + L(x2, 0:1)")
})

test_that("estimates match dynlm, with gaps, first_lag and dummies", {
  skip_if_not_installed("dynlm")
  suppressPackageStartupMessages(library(dynlm))
  m <- sim_data()
  a <- adl(y ~ x1 + x2 | d, data = m, order = c(1, 2, 3), first_lag = c(x2 = 1))
  r <- dynlm::dynlm(y ~ L(y, 1) + L(x1, 0:2) + L(x2, 1:3) + d, data = m)
  expect_equal(unname(coef(a)), unname(coef(r)))
  expect_equal(nobs(a), nobs(r))
  expect_match(a$label, "no contemporaneous x2")

  g <- adl(y ~ L(y, c(1, 4)) + x1, data = m)
  rg <- dynlm::dynlm(y ~ L(y, c(1, 4)) + x1, data = m)
  expect_equal(unname(coef(g)), unname(coef(rg)))
  expect_equal(g$label, "ADL (general lag sets)")
})

test_that("max_lag and start/end control the sample", {
  m <- sim_data()
  a <- adl(y ~ x1, data = m, order = c(1, 1), max_lag = 4)
  expect_equal(nobs(a), 76)
  # start is the first estimated observation; lags come from before it
  b <- adl(y ~ x1, data = m, order = c(1, 1), start = c(2005, 1), end = c(2014, 4))
  expect_equal(nobs(b), 40)
  expect_equal(format(b$index[1]), "2005 Q1")
  b3 <- adl(y ~ x1, data = m, order = c(3, 1), start = c(2005, 1), end = c(2014, 4))
  expect_equal(b3$index, b$index)
  skip_if_not_installed("dynlm")
  suppressPackageStartupMessages(library(dynlm))
  ref <- dynlm(y ~ L(y, 1:3) + L(x1, 0:1), data = m, start = c(2005, 1), end = c(2014, 4))
  expect_equal(unname(coef(b3)), unname(coef(ref)))
  # not enough observations before start: estimation starts as early as possible
  expect_message(e <- adl(y ~ x1, data = m, order = c(2, 1), start = c(2000, 1)),
                 "Not enough observations")
  expect_equal(format(e$index[1]), "2000 Q3")
})

test_that("the same model from a data frame with a date column", {
  m <- sim_data()
  df <- data.frame(date = zoo::as.Date(zoo::as.yearqtr(stats::time(m))),
                   zoo::coredata(m))
  expect_equal(coef(adl(y ~ x1, data = df, order = c(1, 1))),
               coef(adl(y ~ x1, data = m, order = c(1, 1))))
})

test_that("methods work", {
  m <- sim_data()
  a <- adl(y ~ x1, data = m, order = c(1, 1), trend = TRUE)
  expect_s3_class(residuals(a), "zoo")
  expect_equal(length(fitted(a)), nobs(a))
  expect_true(is.finite(AIC(a)))
  expect_equal(dim(vcov(a)), c(5L, 5L))
  expect_output(print(a), "ADL\\(1, 1\\)")
  expect_output(print(summary(a)), "R-squared")
  a0 <- adl(y ~ x1, data = m, order = c(0, 0), intercept = FALSE)
  expect_equal(names(coef(a0)), "x1")
})

test_that("formula errors are clear", {
  m <- sim_data()
  expect_error(adl(y ~ L(y, 1) + x1, data = m, order = c(1, 1)), "already contains")
  expect_error(adl(y ~ y + x1, data = m, order = c(1, 1)), "must be written")
  expect_error(adl(y ~ L(y, 0) + x1, data = m), "lag 0")
  expect_error(adl(y ~ log(x1), data = m, order = c(1, 1)), "not understood")
  expect_error(adl(y ~ x1 | x1, data = m, order = c(1, 1)), "both lagged")
  expect_error(adl(~ x1, data = m), "two-sided")
})
