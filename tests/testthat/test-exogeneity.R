
test_that("Wu-Hausman: control-function regression by hand", {
  fit <- adl(y ~ x, data = adl_sim, order = c(1, 1), max_lag = 2)
  # in adl_sim x is AR(1): its further lags are weak instruments
  expect_warning(wh <- wu_hausman_test(fit, instruments = list(y = 2, x = 2)), "weak")
  f <- fit$lm$model
  des <- build_design(fit$data, list(y = 1:2, x = 0:2), max_lag = 2)
  f <- des$frame
  expect_equal(nrow(f), nobs(fit))
  fs <- lm(x ~ y_lag1 + x_lag1 + y_lag2 + x_lag2, data = f)
  aug <- lm(y ~ y_lag1 + x + x_lag1 + residuals(fs), data = f)
  s <- summary(aug)$coefficients
  expect_equal(wh$xi, unname(s[5, 1]))
  expect_equal(wh$statistic, unname(s[5, 3]))
  expect_equal(suppressWarnings(wu_hausman_test(fit, instruments = list(y = 2, x = 2),
                                                dist = "t"))$p.value,
               unname(s[5, 4]))
  expect_equal(wh$nobs, nobs(fit))
  fs0 <- lm(x ~ y_lag1 + x_lag1, data = f)
  expect_equal(wh$first_stage$F, anova(fs0, fs)$F[2])
  expect_lt(wh$first_stage$F, 10)
})

test_that("Wu-Hausman never uses more than the fitted sample (audit repro B)", {
  set.seed(7)
  d <- stats::ts(data.frame(x = rnorm(100), y = rnorm(100)))
  f <- adl(y ~ x, d, order = c(1, 1), max_lag = 10)
  wh <- suppressWarnings(wu_hausman_test(f, instruments = list(x = 2)))
  expect_equal(wh$nobs, nobs(f))
  # instrument lags beyond the fit's initial values: fewer observations, with a message
  f1 <- adl(y ~ x, d, order = c(1, 1))
  expect_message(w1 <- suppressWarnings(wu_hausman_test(f1, instruments = list(x = 3))),
                 "uses 97 of the 99")
  expect_equal(w1$nobs, 97)
  wh_hac <- suppressWarnings(suppressMessages(
    wu_hausman_test(f1, instruments = list(x = 2, y = 2), cov_type = "HAC")))
  expect_true(is.finite(wh_hac$first_stage$hac_F))
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
  wh <- suppressMessages(wu_hausman_test(fit, instruments = list(x = 1:2)))
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

test_that("weak exogeneity: unstable stops; controls, trend and other regressors enter", {
  fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
  bad <- fit
  bad$lm$coefficients[["y_lag1"]] <- 1.05
  expect_error(weak_exogeneity_test(bad), "not dynamically stable")

  d <- adl_sim
  dm <- cbind(d, dmy = as.numeric(seq_len(nrow(d)) %% 4 == 0))
  colnames(dm) <- c("y", "x", "w", "dmy")
  ft <- adl(y ~ x + w | dmy, data = dm, order = c(1, 1), trend = TRUE)
  we <- weak_exogeneity_test(ft, x = "x")
  nm <- names(coef(we$marginal))
  expect_true(all(c("dw_lag1", "w_level_lag1", "dmy_lag1", "trend", "x_lag1") %in% nm))
  # with the levels in, kappa equals the coefficient of y_{t-1} given them
  mf <- we$marginal$model
  ylag <- mf$ec_lag1 + as.numeric(cbind(mf$x_lag1, mf$w_level_lag1) %*%
                                    adl_long_run(ft)$estimate[2:3])
  alt <- lm(mf$dx ~ . , data = cbind(mf[setdiff(names(mf), c("dx", "ec_lag1"))],
                                     ylag = ylag))
  expect_equal(we$kappa, unname(coef(alt)[["ylag"]]))
})

test_that("a control named ec does not replace the equilibrium error (audit 2)", {
  set.seed(41)
  n <- 150
  x <- as.numeric(arima.sim(list(ar = 0.5), n))
  ec <- rnorm(n)
  y <- numeric(n)
  for (t in 2:n) y[t] <- 0.5 * y[t - 1] + 0.4 * x[t] + 0.3 * x[t - 1] + 0.5 * ec[t] + rnorm(1)
  d <- stats::ts(cbind(y = y, x = x, ec = ec))
  f <- adl(y ~ x | ec, data = d, order = c(1, 1))
  w <- weak_exogeneity_test(f)
  nm <- names(stats::coef(w$marginal))
  expect_true(all(c("ec_lag1", "ec_lag1.1") %in% nm))
  # kappa is the coefficient on the equilibrium error, not on the control:
  # the same as with the control renamed
  d2 <- d; colnames(d2)[3] <- "u"
  w2 <- weak_exogeneity_test(adl(y ~ x | u, data = d2, order = c(1, 1)))
  expect_equal(w$kappa, w2$kappa)
  expect_equal(w$p.value, w2$p.value)
  # a regressor called dx, or a control called x, cannot clash either
  d3 <- d; colnames(d3)[3] <- "dx"
  expect_equal(weak_exogeneity_test(adl(y ~ x | dx, data = d3, order = c(1, 1)))$kappa,
               w$kappa)
})

# Supplied equilibrium error (ec =) -------------------------------------------

ec_data <- function(n = 300, seed = 81) {
  set.seed(seed)
  x <- cumsum(rnorm(n))
  y <- 1 + 0.8 * x + as.numeric(stats::arima.sim(list(ar = 0.5), n))
  stats::ts(cbind(y = y, x = x), start = c(1950, 1), frequency = 4)
}

test_that("ec = the ADL's own equilibrium error reproduces the default test", {
  d <- ec_data()
  fit <- adl(y ~ x, data = d, order = c(2, 2))
  k <- adl_long_run(fit)
  ec <- d[, "y"] - k$estimate[k$term == "x"] * d[, "x"]
  for (xl in c(TRUE, FALSE)) for (st in c(FALSE, TRUE)) {
    a <- weak_exogeneity_test(fit, x_level = xl, strong = st)
    b <- weak_exogeneity_test(fit, x_level = xl, strong = st, ec = ec)
    expect_equal(b$statistic, a$statistic)
    expect_equal(b$kappa, a$kappa)
    expect_equal(b$ec_source, "supplied")
    expect_equal(a$ec_source, "adl")
  }
  # a plain vector, one value per row of the data
  v <- weak_exogeneity_test(fit, x_level = FALSE, ec = as.numeric(ec))
  expect_equal(v$statistic, weak_exogeneity_test(fit, x_level = FALSE)$statistic)
})

test_that("ec = a cointegrating residual: matters with x_level = FALSE only", {
  d <- ec_data()
  fit <- adl(y ~ x, data = d, order = c(2, 2))
  z <- stats::ts(residuals(lm(y ~ x, data = as.data.frame(d))),
                 start = c(1950, 1), frequency = 4)
  ad <- weak_exogeneity_test(fit, x_level = FALSE)
  sz <- weak_exogeneity_test(fit, x_level = FALSE, ec = z)
  expect_false(isTRUE(all.equal(sz$kappa, ad$kappa)))
  # with the level of x in the regression, any y - b x gives the same kappa
  expect_equal(weak_exogeneity_test(fit, ec = z)$kappa,
               weak_exogeneity_test(fit)$kappa)
  # direct regression for the change in x
  z1 <- stats::lag(z, -1)
  m <- stats::ts.intersect(dx = diff(d[, "x"]), dx1 = stats::lag(diff(d[, "x"]), -1),
                           dy1 = stats::lag(diff(d[, "y"]), -1), z1 = z1)
  m <- as.data.frame(m)          # the fit's sample: from 1950(3)
  ref <- lm(dx ~ dx1 + dy1 + z1, data = m)
  expect_equal(sz$kappa, unname(coef(ref)["z1"]))
})

test_that("ec = is matched by date, and works for an unstable fit", {
  d <- ec_data()
  fit <- adl(y ~ x, data = d, order = c(1, 1))
  z <- stats::ts(residuals(lm(y ~ x, data = as.data.frame(d))),
                 start = c(1950, 1), frequency = 4)
  # a shorter series covers part of the sample: fewer observations
  zs <- stats::window(z, start = c(1960, 1))
  w <- weak_exogeneity_test(fit, x_level = FALSE, ec = zs)
  expect_lt(w$nobs, weak_exogeneity_test(fit, x_level = FALSE, ec = z)$nobs)
  fit$lm$coefficients[["y_lag1"]] <- 1.02
  expect_error(weak_exogeneity_test(fit), "not dynamically stable")
  expect_s3_class(weak_exogeneity_test(fit, x_level = FALSE, ec = z), "adltools_test")
})

test_that("ec = input checks", {
  d <- ec_data()
  fit <- adl(y ~ x, data = d, order = c(1, 1))
  expect_error(weak_exogeneity_test(fit, ec = 1:10), "300 rows")
  expect_error(weak_exogeneity_test(fit, ec = cbind(a = 1:300, b = 1:300)), "numeric vector")
  expect_error(weak_exogeneity_test(fit, ec = zoo::zoo(rnorm(300), as.Date("2000-01-01") + 0:299)),
               "same kind")
  expect_error(weak_exogeneity_test(fit, ec = stats::ts(rnorm(20), start = c(2100, 1), frequency = 4)),
               "no values on the estimation sample")
  expect_error(weak_exogeneity_test(fit, ec = "a"), "numeric vector")
})
