# Reproductions from the fourth external audit (2026-11-18, cc86fb4)

audit_data <- function() {
  set.seed(7)
  data.frame(x = rnorm(100), y = rnorm(100), w = rnorm(100))
}

test_that("F1: ECM standard errors do not depend on the units of x", {
  d <- audit_data()
  e1 <- adl_ecm(adl(y ~ x, data = d, order = c(1, 1)))
  e2 <- adl_ecm(adl(y ~ x, data = transform(d, x = x * 1e10), order = c(1, 1)))
  tab <- e2$coefficients
  expect_true(all(tab$std.error > 0))
  expect_equal(tab$std.error, unname(sqrt(diag(e2$vcov))))
  expect_false(anyNA(tab$p.value))
  # the t-statistics are unit-free
  expect_equal(tab$statistic, e1$coefficients$statistic, tolerance = 1e-6)
  # a genuinely fixed term is still shown as fixed
  e0 <- adl_ecm(adl(y ~ x, data = d, order = c(0, 1)))
  b <- e0$coefficients[e0$coefficients$param == "b", ]
  expect_equal(b$std.error, 0)
  expect_true(is.na(b$p.value))
})

test_that("P3: ECM row labels stay unique", {
  d <- audit_data()
  e <- adl_ecm(adl(y ~ x | x_lag1, data = transform(d, x_lag1 = rnorm(100)),
                   order = c(0, 0)))
  expect_false(anyDuplicated(e$coefficients$term) > 0)
  expect_false(anyDuplicated(rownames(e$vcov)) > 0)
})

test_that("F2: forecast rows one second after the data are kept", {
  d <- audit_data()
  base <- as.POSIXct("2020-01-01", tz = "UTC")
  f <- adl(y ~ x, data = transform(d, t = base + 1:100), order = c(1, 1), time = "t")
  fc <- adl_forecast(f, h = 3, newdata = data.frame(x = 1:3, t = base + 101:103),
                     time = "t")
  expect_equal(fc$forecasts$index, base + 101:103)
})

test_that("F3: regular sub-second grids pass and a missing step still fails", {
  d <- audit_data()
  a <- transform(d, t = as.POSIXct("2020-01-01", tz = "UTC") + (1:100) * 0.1)
  expect_equal(nobs(adl(y ~ x, data = a, order = c(1, 1), time = "t")), 99)
  expect_error(adl(y ~ x, data = a[-50, ], order = c(1, 1), time = "t"),
               "00:00:04.900 and 2020-01-01 00:00:05.100")
})

test_that("F4: fan charts with levels that format differently together", {
  d <- audit_data()
  f <- adl(y ~ x, data = d, order = c(1, 1))
  fc <- adl_forecast(f, newdata = data.frame(x = 1:4), level = c(0.5, 0.95, 0.999))
  expect_true(all(c("upper_50", "upper_95", "upper_99.9") %in% names(fc$forecasts)))
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(fc))
})

test_that("F5: non-syntactic names survive the refits of adl_getsm and adl_outliers", {
  d <- audit_data()
  d[["my x"]] <- d$x
  expect_message(f <- adl(y ~ `my x`, data = d, order = c(1, 1)), "my.x")
  skip_if_not_installed("gets")
  s <- suppressMessages(adl_getsm(f, diagnostics = FALSE))
  expect_s3_class(s$fit, "adl")
  d$y[50] <- d$y[50] + 15
  f2 <- suppressMessages(adl(y ~ `my x`, data = d, order = c(1, 1)))
  o <- suppressMessages(adl_outliers(f2))
  expect_true(50 %in% as.numeric(o$outliers$date))
})

test_that("F6: monthly dates capped at short months are extended correctly", {
  set.seed(8)
  a <- data.frame(y = rnorm(4), x = rnorm(4),
                  t = as.Date(c("2020-10-30", "2020-11-30", "2020-12-30", "2021-01-30")))
  f <- adl(y ~ x, data = a, order = c(0, 0), time = "t")
  fc <- adl_forecast(f, newdata = data.frame(x = 1:2))
  expect_equal(fc$forecasts$index, as.Date(c("2021-02-28", "2021-03-30")))
})

test_that("F7: J-test refuses aliased and weighted source fits", {
  skip_if_not_installed("dynlm")
  suppressPackageStartupMessages(library(dynlm))
  dat <- ts(as.matrix(audit_data()))      # not `d`: dynlm has its own d()
  m <- dynlm(y ~ w + I(2 * w), data = dat)
  expect_error(j_test(adl(y ~ x, data = dat, order = c(0, 0)), m), "aliased")
  mw <- dynlm(y ~ L(y, 1), data = dat, weights = rep(1:2, 50)[-1])
  expect_error(j_test(adl(y ~ x, data = dat, order = c(0, 1)), mw), "weights")
  # rank deficiency that appears only on the common sample
  dd <- audit_data()
  dd$z <- c(rep(0, 99), 1)
  fz <- adl(y ~ x + z, data = dd, order = c(0, 0, 0))
  fl <- adl(y ~ x, data = dd[-100, ], order = c(1, 1))
  expect_error(j_test(fz, fl), "rank deficient")
})

test_that("earlier P2s: tiny CUSUM alpha, tiny outlier samples, max_horizon", {
  d <- audit_data()
  f <- adl(y ~ x, data = d, order = c(1, 1))
  expect_equal(cusum_test(f)$critical, 0.948, tolerance = 1e-3)
  expect_gt(cusum_test(f, alpha = 1e-50)$critical, 5)
  expect_error(adl_outliers(adl(y ~ x, data = d[1:2, ], order = c(0, 0)), types = "step"),
               "Too few observations")
  expect_equal(nrow(adl_outliers(adl(y ~ x, data = d[1:6, ], order = c(0, 0)),
                                 types = "step")$outliers), 0)
  expect_error(adl_dynamics(f, max_horizon = -1), "max_horizon")
})
