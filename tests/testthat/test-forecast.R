# The dynamic recursion written out for an ADL(2, 2) with a trend and a dummy
manual_forecast <- function(fit, z, o, h) {
  b <- coef(fit)
  y <- z[, "y"]
  y[(o + 1):length(y)] <- NA
  t0 <- which(zoo::index(fit$data) == fit$index[1])
  for (t in (o + 1):(o + h)) {
    y[t] <- b[["(Intercept)"]] + b[["y_lag1"]] * y[t - 1] + b[["y_lag2"]] * y[t - 2] +
      b[["x"]] * z[t, "x"] + b[["x_lag1"]] * z[t - 1, "x"] + b[["x_lag2"]] * z[t - 2, "x"] +
      b[["d"]] * z[t, "d"] + b[["trend"]] * (t - t0 + 1)
  }
  unname(y[(o + 1):(o + h)])
}

fc_data <- function(n = 120) {
  set.seed(41)
  x <- as.numeric(stats::arima.sim(list(ar = 0.5), n))
  d <- as.numeric(seq_len(n) %% 10 == 0)
  y <- stats::filter(0.5 * x + 0.2 * c(0, x[-n]) + d + 0.01 * seq_len(n) + rnorm(n),
                     c(0.5, 0.2), method = "recursive")
  stats::ts(cbind(y = as.numeric(y), x = x, d = d), start = c(1990, 1), frequency = 4)
}

test_that("dynamic forecasts follow the recursion, in and beyond the sample", {
  dat <- fc_data()
  fit <- adl(y ~ x | d, data = dat, order = c(2, 2), trend = TRUE)
  z <- zoo::coredata(fit$data)
  # from inside the sample: x from the data, actual values reported
  fc <- adl_forecast(fit, origin = c(2015, 4), h = 10)
  o <- which(abs(as.numeric(zoo::index(fit$data)) - 2015.75) < 1e-8)
  expect_equal(fc$forecasts$forecast, manual_forecast(fit, z, o, 10))
  expect_equal(fc$forecasts$actual, unname(z[(o + 1):(o + 10), "y"]))
  expect_true(fc$in_sample)
  # one step ahead from inside the sample is the fitted value
  k <- which(abs(as.numeric(fit$index) - 2016) < 1e-8)
  expect_equal(fc$forecasts$forecast[1], as.numeric(fitted(fit))[k])

  # beyond the sample: a scenario without dates continues the data
  path <- data.frame(x = c(1, 0.5, 0, 0, 0, 0), d = c(0, 0, 1, 0, 0, 0))
  f2 <- adl_forecast(fit, newdata = path)
  expect_false(f2$in_sample)
  expect_equal(nrow(f2$forecasts), 6)
  expect_equal(as.numeric(f2$forecasts$index[1]), 2020)
  z2 <- rbind(z, cbind(y = NA, as.matrix(path)))
  expect_equal(f2$forecasts$forecast, manual_forecast(fit, z2, nrow(z), 6))
  expect_true(all(is.na(f2$forecasts$actual)))

  # the same scenario with dates, inside a longer data set
  long <- stats::ts(rbind(dat, cbind(y = NA, as.matrix(path))), start = c(1990, 1),
                    frequency = 4)
  f3 <- adl_forecast(fit, newdata = long)
  expect_equal(f3$forecasts$forecast, f2$forecasts$forecast)
  expect_equal(adl_forecast(fit, newdata = path, h = 3)$forecasts$forecast,
               f2$forecasts$forecast[1:3])
})

test_that("psi weights and standard errors", {
  fit <- adl(y ~ x, data = adl_sim, order = c(2, 1))
  a <- coef(fit)[c("y_lag1", "y_lag2")]
  fc <- adl_forecast(fit, origin = c(2017, 4))
  psi <- numeric(8)
  psi[1] <- 1
  psi[2] <- a[1]
  for (j in 3:8) psi[j] <- a[1] * psi[j - 1] + a[2] * psi[j - 2]
  expect_equal(fc$psi, unname(psi))
  s <- summary(fit$lm)$sigma
  expect_equal(fc$forecasts$se, s * sqrt(cumsum(psi^2)))
  expect_equal(fc$forecasts$se[1], s)
  z95 <- qnorm(0.975)
  expect_equal(fc$forecasts$upper_95, fc$forecasts$forecast + z95 * fc$forecasts$se)
  # without lags of y the error does not accumulate
  f0 <- adl_forecast(adl(y ~ x, data = adl_sim, order = c(0, 1)), origin = c(2019, 1))
  expect_equal(f0$psi, c(1, 0, 0))
  expect_equal(length(unique(round(f0$forecasts$se, 12))), 1)
})

test_that("monthly dates are extended at month ends and on fixed days", {
  set.seed(42)
  n <- 60
  for (dates in list(seq(as.Date("2015-02-01"), by = "month", length.out = n + 1)[-1] - 1,
                     seq(as.Date("2015-01-15"), by = "month", length.out = n))) {
    df <- data.frame(date = dates, x = rnorm(n))
    df$y <- 0.5 * df$x + rnorm(n)
    fit <- adl(y ~ x, data = df, order = c(1, 1))
    fc <- adl_forecast(fit, newdata = data.frame(x = c(0, 0, 0)))
    nxt <- seq(dates[n], by = "month", length.out = 4)[-1]
    if (format(dates[1], "%d") != "15") {
      nxt <- seq(dates[n] + 1, by = "month", length.out = 4)[-1] - 1
    }
    expect_equal(fc$forecasts$index, nxt)
  }
})

test_that("argument checks", {
  fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
  expect_error(adl_forecast(fit), "No values of the regressors")
  expect_error(adl_forecast(fit, newdata = data.frame(x = 1:2), h = 3), "covers 2 step")
  expect_error(adl_forecast(fit, newdata = data.frame(z = 1:2)), "lacks: x")
  expect_error(adl_forecast(fit, origin = c(1980, 1)), "origin")
  expect_error(adl_forecast(fit, origin = c(2019, 1), level = 95), "level")
  bad <- data.frame(date = as.Date("2030-01-01") + 0:1, x = 1:2)
  expect_error(adl_forecast(fit, newdata = bad), "same kind")
  # a scenario that skips a quarter
  gap <- stats::ts(cbind(y = NA, x = 1:2), start = c(2020, 2), frequency = 4)
  expect_error(adl_forecast(fit, newdata = gap), "gaps")
})

test_that("printing, fan chart and tidiers", {
  fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
  fc <- adl_forecast(fit, origin = c(2018, 4), level = c(0.9, 0.6))
  out <- capture.output(print(fc))
  expect_true(any(grepl("in-sample", out)))
  expect_true(any(grepl("upper_90", out)))
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(fc))
  expect_silent(plot(adl_forecast(fit, newdata = data.frame(x = rep(0, 5))),
                     n_history = 8))
  expect_equal(nrow(tidy(fc)), 4)
  expect_equal(glance(fc)$h, 4)
})
