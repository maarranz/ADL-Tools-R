
test_that("Wu-Hausman: control-function regression by hand", {
  fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
  # in adl_sim x is AR(1): its further lags are weak instruments
  expect_warning(wh <- wu_hausman_test(fit, instruments = list(y = 2, x = 2)), "weak")
  des <- build_design(fit$data, list(y = 1:2, x = 0:2), max_lag = 2)
  f <- des$frame
  fs <- lm(x ~ y_lag1 + x_lag1 + y_lag2 + x_lag2, data = f)
  aug <- lm(y ~ y_lag1 + x + x_lag1 + residuals(fs), data = f)
  s <- summary(aug)$coefficients
  expect_equal(wh$xi, unname(s[5, 1]))
  expect_equal(wh$statistic, unname(s[5, 3]))
  expect_equal(suppressWarnings(wu_hausman_test(fit, instruments = list(y = 2, x = 2),
                                                dist = "t"))$p.value,
               unname(s[5, 4]))
  expect_equal(wh$nobs, nrow(f))
  expect_lt(wh$first_stage$F, 10)
})

test_that("Wu-Hausman detects simultaneity", {
  set.seed(61)
  n <- 800
  z <- as.numeric(stats::arima.sim(list(ar = 0.7), n))
  e <- rnorm(n)
  x <- z + 0.8 * e + rnorm(n, sd = 0.3)
  y <- numeric(n)
  for (t in 2:n) y[t] <- 0.5 * y[t - 1] + 0.5 * x[t] + e[t]
  # x_t = z_t + 0.8 e_t with z AR(1): lags of x are valid, strong instruments
  fit <- adl(y ~ x, data = stats::ts(cbind(y = y, x = x)), order = c(1, 0))
  wh <- wu_hausman_test(fit, instruments = list(x = 1:2))
  expect_true(wh$reject_null)
  expect_gt(wh$first_stage$F, 10)
})

test_that("weak exogeneity: marginal regression by hand", {
  fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
  we <- weak_exogeneity_test(fit, lags = 2, x_level = FALSE)  # ADL3 as written
  k <- adl_long_run(fit)$estimate[2]
  z <- zoo::coredata(adl_sim)
  N <- nrow(z)
  ec <- z[, "y"] - k * z[, "x"]
  dx <- c(NA, diff(z[, "x"])); dy <- c(NA, diff(z[, "y"]))
  L <- function(v, j) c(rep(NA, j), v[1:(N - j)])
  ref <- lm(dx ~ L(dx, 1) + L(dx, 2) + L(dy, 1) + L(dy, 2) + L(ec, 1))
  s <- summary(ref)$coefficients
  expect_equal(we$kappa, unname(s[6, 1]))
  expect_equal(we$statistic, unname(s[6, 3]))
  expect_equal(we$nobs, nobs(ref))

  # default adds the level x_{t-1}
  wl <- weak_exogeneity_test(fit, lags = 2)
  refl <- lm(dx ~ L(dx, 1) + L(dx, 2) + L(dy, 1) + L(dy, 2) + L(ec, 1) + L(z[, "x"], 1))
  expect_equal(wl$statistic, unname(summary(refl)$coefficients[6, 3]))
  expect_false(wl$reject_null)  # x is exogenous in adl_sim

  st <- weak_exogeneity_test(fit, lags = 2, strong = TRUE, x_level = FALSE)
  expect_equal(st$df, 3)
  V <- vcov(ref); b <- coef(ref)
  sel <- c(4, 5, 6)
  expect_equal(st$statistic, as.numeric(t(b[sel]) %*% solve(V[sel, sel], b[sel])))
})

test_that("with a stationary exogenous x, the level term keeps the size", {
  set.seed(71)
  rej <- function(lev) mean(replicate(60, {
    x <- as.numeric(stats::arima.sim(list(ar = 0.5), 200)); y <- numeric(200); e <- rnorm(200)
    for (t in 2:200) y[t] <- 0.5 * y[t - 1] + 0.4 * x[t] + 0.3 * x[t - 1] + e[t]
    f <- adl(y ~ x, data = stats::ts(cbind(y = y, x = x)), order = c(1, 1))
    weak_exogeneity_test(f, x_level = lev)$reject_null
  }))
  expect_lt(rej(TRUE), 0.2)
  expect_gt(rej(FALSE), 0.3)
})

test_that("weak exogeneity rejected when x error-corrects", {
  set.seed(62)
  n <- 800
  y <- x <- numeric(n); e <- rnorm(n); v <- rnorm(n)
  for (t in 2:n) {
    x[t] <- 0.5 * x[t - 1] + 0.4 * (y[t - 1] - 1.4 * x[t - 1]) + v[t]
    y[t] <- 0.5 * y[t - 1] + 0.4 * x[t] + 0.3 * x[t - 1] + e[t]
  }
  fit <- adl(y ~ x, data = stats::ts(cbind(y = y, x = x)), order = c(1, 1))
  expect_true(weak_exogeneity_test(fit)$reject_null)
})

test_that("argument checks", {
  fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
  expect_error(wu_hausman_test(fit, instruments = list(x = 1)), "already in the model")
  expect_error(wu_hausman_test(fit, instruments = list(y = 0)), "lags >= 1|already")
  expect_error(wu_hausman_test(fit, instruments = c(2, 2)), "named list")
  ds <- adl(y ~ x, data = adl_sim, order = c(1, 1), first_lag = 1)
  expect_error(wu_hausman_test(ds, instruments = list(x = 2)), "lag 0")
  two <- adl(y ~ x + w, data = adl_sim, order = c(1, 1))
  expect_error(weak_exogeneity_test(two), "choose one")
  expect_output(print(weak_exogeneity_test(two, x = "w")), "kappa = 0")
})
