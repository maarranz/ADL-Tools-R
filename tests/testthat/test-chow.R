chow_fit <- function() adl(y ~ x, data = adl_sim, order = c(1, 1))

test_that("Chow breakpoint F equals the RSS formula", {
  fit <- chow_fit()
  ch <- chow_break_test(fit, c(2005, 1), dist = "F")
  X <- model.matrix(fit$lm); y <- residuals(fit$lm) + fitted(fit$lm)
  k <- which(fit$index >= zoo::as.yearqtr("2005 Q1"))[1]
  r1 <- lm.fit(X[1:(k - 1), ], y[1:(k - 1)]); r2 <- lm.fit(X[k:nrow(X), ], y[k:nrow(X)])
  rss1 <- sum(r1$residuals^2); rss2 <- sum(r2$residuals^2); rssp <- sum(residuals(fit$lm)^2)
  m <- ncol(X); n <- nrow(X)
  Fc <- ((rssp - rss1 - rss2) / m) / ((rss1 + rss2) / (n - 2 * m))
  expect_equal(ch$statistic, Fc)
  expect_equal(ch$df_denom, n - 2 * m)
  expect_equal(ch$break_date, "2005 Q1")
  expect_equal(chow_break_test(fit, c(2005, 1))$statistic, m * Fc)
  expect_equal(chow_break_test(fit, c(2005, 1), cov_type = "HAC")$cov_type, "HAC")
})

test_that("Chow breakpoint detects a real break", {
  set.seed(40)
  n <- 200
  x <- rnorm(n)
  y <- ifelse(seq_len(n) <= 100, 0.5, 2) * x + rnorm(n)
  fit <- adl(y ~ x, data = stats::ts(cbind(y = y, x = x)), order = c(0, 0))
  expect_true(chow_break_test(fit, 101)$reject_null)
})

test_that("static Chow forecast F equals the classical predictive test", {
  fit <- chow_fit()
  cf <- chow_forecast_test(fit, c(2016, 1), dist = "F")
  X <- model.matrix(fit$lm); y <- residuals(fit$lm) + fitted(fit$lm)
  k <- which(fit$index >= zoo::as.yearqtr("2016 Q1"))[1]
  T1 <- k - 1; T2 <- nrow(X) - T1; m <- ncol(X)
  rss1 <- sum(lm.fit(X[1:T1, ], y[1:T1])$residuals^2)
  rssp <- sum(residuals(fit$lm)^2)
  Fc <- ((rssp - rss1) / T2) / (rss1 / (T1 - m))
  expect_equal(cf$statistic, Fc)
  expect_equal(c(cf$df_num, cf$df_denom), c(T2, T1 - m))
  expect_equal(nrow(cf$forecasts), T2)
  # dummy-variable regression gives the same statistic
  D <- diag(nrow(X))[, k:nrow(X)]
  aug <- lm(y ~ 0 + cbind(X, D))
  r0 <- lm(y[1:T1] ~ 0 + X[1:T1, ])
  expect_equal(cf$statistic,
               ((sum(residuals(fit$lm)^2) - sum(residuals(aug)^2)) / T2) /
                 (sum(residuals(aug)^2) / (T1 - m)))
})

test_that("dynamic Chow forecast: path and Omega for the ADL(1,1)", {
  fit <- chow_fit()
  cf <- chow_forecast_test(fit, c(2016, 1), type = "dynamic")
  X <- model.matrix(fit$lm); y <- residuals(fit$lm) + fitted(fit$lm)
  k <- which(fit$index >= zoo::as.yearqtr("2016 Q1"))[1]
  T1 <- k - 1; T2 <- nrow(X) - T1
  b <- lm.fit(X[1:T1, ], y[1:T1])$coefficients
  s2 <- sum(lm.fit(X[1:T1, ], y[1:T1])$residuals^2) / (T1 - ncol(X))
  # dynamic path: y_lag1 replaced by previous forecast after the first step
  f <- numeric(T2); ylag <- X[k, "y_lag1"]
  for (h in 1:T2) {
    f[h] <- b[1] + b[2] * ylag + b[3] * X[k + h - 1, "x"] + b[4] * X[k + h - 1, "x_lag1"]
    ylag <- f[h]
  }
  expect_equal(cf$forecasts$forecast, f, ignore_attr = TRUE)
  # closed form of ADL4: Cov(e_j, e_l) = s2 a^(l-j) (1 - a^(2j)) / (1 - a^2)
  a <- b[["y_lag1"]]
  Om <- outer(1:T2, 1:T2, function(j, l) {
    jj <- pmin(j, l); d <- abs(l - j)
    s2 * a^d * (1 - a^(2 * jj)) / (1 - a^2)
  })
  e <- y[k:nrow(X)] - f
  expect_equal(cf$statistic, as.numeric(t(e) %*% solve(Om, e)))
  expect_equal(cf$forecasts$se, sqrt(diag(Om)), ignore_attr = TRUE)
  # static and dynamic differ
  expect_false(isTRUE(all.equal(cf$statistic,
                                chow_forecast_test(fit, c(2016, 1))$statistic)))
})

test_that("argument checks", {
  fit <- chow_fit()
  expect_error(chow_break_test(fit, c(1990, 1)), "inside the estimation sample")
  expect_error(chow_break_test(fit, c(2019, 3)), "chow_forecast_test")
  expect_error(chow_forecast_test(fit, c(1990, 4)), "more than")
  expect_output(print(chow_forecast_test(fit, c(2018, 1))), "Chow forecast test \\(static\\)")
})
