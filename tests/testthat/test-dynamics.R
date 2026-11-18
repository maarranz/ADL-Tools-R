dyn_fit <- function() {
  set.seed(8)
  n <- 300
  x1 <- cumsum(rnorm(n)); x2 <- rnorm(n)
  y <- stats::filter(0.5 * x1 + 0.3 * x2 + rnorm(n), c(0.6, 0.2), method = "recursive")
  d <- stats::ts(cbind(y = as.numeric(y), x1 = x1, x2 = x2))
  adl(y ~ x1 + x2, data = d, order = c(2, 1, 2), first_lag = c(x2 = 1))
}

# Reference: the recursion written out directly
ref_effects <- function(alpha, beta, H) {
  d <- numeric(H + 1)
  for (h in 0:H) {
    bh <- if (h + 1 <= length(beta)) beta[h + 1] else 0
    ar <- sum(vapply(seq_along(alpha), function(i)
      if (h - i >= 0) alpha[i] * d[h - i + 1] else 0, numeric(1)))
    d[h + 1] <- bh + ar
  }
  d
}

test_that("multipliers follow the recursion and add up to the long run", {
  fit <- dyn_fit()
  b <- coef(fit)
  dy <- adl_dynamics(fit, horizon = 60)
  a <- b[c("y_lag1", "y_lag2")]
  e1 <- dy$effects[dy$effects$variable == "x1", ]
  expect_equal(e1$effect, ref_effects(a, b[c("x1", "x1_lag1")], 60), ignore_attr = TRUE)
  e2 <- dy$effects[dy$effects$variable == "x2", ]
  # x2 has no lag 0: beta_0 = 0
  expect_equal(e2$effect[1], 0)
  expect_equal(e2$effect, ref_effects(a, c(0, b[c("x2_lag1", "x2_lag2")]), 60),
               ignore_attr = TRUE)
  lr <- adl_long_run(fit)
  expect_equal(dy$summary$long_run, lr$estimate[lr$term %in% c("x1", "x2")])
  long <- adl_dynamics(fit, x = "x1", horizon = 400)
  expect_equal(tail(long$effects$cumulative, 1), dy$summary$long_run[1], tolerance = 1e-8)
})

test_that("automatic horizon stops at the tolerance", {
  fit <- dyn_fit()
  dy <- adl_dynamics(fit, tolerance = 0.95)
  for (v in c("x1", "x2")) {
    e <- dy$effects[dy$effects$variable == v, ]
    lr <- dy$summary$long_run[dy$summary$variable == v]
    k <- nrow(e)
    expect_gte(e$cumulative[k] / lr, 0.95)
    expect_lt(e$cumulative[k - 1] / lr, 0.95)
    expect_equal(e$pct_cumulative[k], 100 * e$cumulative[k] / lr)
  }
})

test_that("ADL(1,1) closed forms: mean lag, median lag, multipliers", {
  # y = a y_{-1} + b0 x + b1 x_{-1}: d_0 = b0, d_h = a^(h-1) (a b0 + b1)
  set.seed(9)
  n <- 400
  x <- rnorm(n)
  y <- stats::filter(0.4 * x + 0.3 * c(0, x[-n]) + rnorm(n, sd = 0.1), 0.7,
                     method = "recursive")
  fit <- adl(y ~ x, data = stats::ts(cbind(y = as.numeric(y), x = x)), order = c(1, 1))
  b <- coef(fit)
  a <- b[["y_lag1"]]; b0 <- b[["x"]]; b1 <- b[["x_lag1"]]
  dy <- adl_dynamics(fit, horizon = 10)
  expect_equal(dy$effects$effect, c(b0, a^(0:9) * (a * b0 + b1)))
  expect_equal(dy$summary$mean_lag, b1 / (b0 + b1) + a / (1 - a))
  share <- dy$effects$cumulative / dy$summary$long_run
  expect_equal(dy$summary$median_lag, which(share >= 0.5)[1] - 1L)
})

test_that("roots and stability", {
  fit <- dyn_fit()
  dy <- adl_dynamics(fit, horizon = 5)
  b <- coef(fit)
  r <- polyroot(c(1, -b[["y_lag1"]], -b[["y_lag2"]]))
  expect_equal(sort(dy$roots$modulus), sort(Mod(r)))
  expect_true(dy$stable)
  expect_equal(dy$roots$inverse_modulus, 1 / dy$roots$modulus)

  fit$lm$coefficients[["y_lag1"]] <- 1.1
  expect_warning(du <- adl_dynamics(fit, horizon = 5), "not dynamically stable")
  expect_false(du$stable)
  expect_true(all(is.na(du$summary$long_run)))
  expect_equal(nrow(du$effects), 12)
})

test_that("models without lags of y: multipliers are the betas", {
  set.seed(10)
  d <- stats::ts(cbind(y = rnorm(100), x = rnorm(100)))
  fit <- adl(y ~ x, data = d, order = c(0, 2))
  dy <- adl_dynamics(fit, horizon = 4)
  expect_equal(dy$effects$effect, c(unname(coef(fit)[c("x", "x_lag1", "x_lag2")]), 0, 0))
  expect_equal(nrow(dy$roots), 0)
})

test_that("printing, plotting and argument checks", {
  fit <- dyn_fit()
  dy <- adl_dynamics(fit, x = "x1")
  expect_equal(unique(dy$effects$variable), "x1")
  out <- capture.output(print(dy))
  expect_true(any(grepl("Dynamically stable: yes", out)))
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(dy))
  expect_silent(plot(dy, type = "cumulative"))
  expect_silent(plot(dy, type = "roots"))
  expect_error(adl_dynamics(fit, x = "z"), "Not lagged regressors")
  expect_error(adl_dynamics(fit, tolerance = 1), "tolerance")
  expect_error(adl_dynamics(fit, horizon = -1), "horizon")
})

test_that("automatic horizon does not stop at an early overshoot (audit repro E)", {
  set.seed(7)
  d <- data.frame(x = rnorm(100), y = rnorm(100))
  e <- transform(d, y = 2 * x - c(0, head(x, -1)) + rnorm(100, sd = .01))
  f <- adl(y ~ x, stats::ts(e[c("y", "x")]), order = c(0, 1))
  dy <- adl_dynamics(f)
  expect_equal(dy$summary$horizon, 1)
  expect_equal(nrow(dy$effects), 2)
  # with an AR part the horizon covers the distributed lags and settles
  f2 <- adl(y ~ x, data = adl_sim, order = c(1, 3))
  d2 <- adl_dynamics(f2)
  expect_gte(d2$summary$horizon, 3)
  k <- d2$summary$long_run
  last_M <- tail(d2$effects$cumulative, 1)
  expect_lte(abs(last_M - k), 0.01 * abs(k))
  # and it stays there for the next AR cycle
  longer <- adl_dynamics(f2, horizon = d2$summary$horizon + 1)
  expect_lte(abs(tail(longer$effects$cumulative, 1) - k), 0.01 * abs(k))
})

test_that("unstable model: multipliers over 20 lags, NA summaries", {
  fit <- dyn_fit()
  fit$lm$coefficients[["y_lag1"]] <- 1.1
  expect_message(expect_warning(du <- adl_dynamics(fit), "not dynamically stable"),
                 "up to lag 20")
  expect_equal(du$summary$horizon, c(20L, 20L))
  expect_true(all(is.na(du$summary$long_run)))
  expect_true(all(is.na(du$effects$pct_cumulative)))
})

test_that("a short max_horizon warns instead of failing (audit 2)", {
  fit <- adl(y ~ x, data = adl_sim, order = c(1, 3))
  expect_warning(dy <- adl_dynamics(fit, max_horizon = 2), "did not settle")
  expect_equal(max(dy$effects$lag), 2)
  expect_warning(adl_dynamics(fit, max_horizon = 0), "did not settle")
})

test_that("delta-method standard errors match a numerical Jacobian", {
  fit <- dyn_fit()
  dy <- adl_dynamics(fit, horizon = 12, nsim = 0)
  b <- coef(fit)
  V <- vcov(fit)
  f <- function(th) {
    f2 <- fit
    f2$lm$coefficients[] <- th
    dd <- adl_dynamics(f2, horizon = 12, se = FALSE)
    c(dd$effects$effect, dd$effects$cumulative, dd$summary$long_run,
      dd$summary$mean_lag)
  }
  J <- vapply(seq_along(b), function(k) {
    e <- numeric(length(b)); e[k] <- 1e-6
    (f(b + e) - f(b - e)) / 2e-6
  }, numeric(length(f(b))))
  se_num <- sqrt(rowSums((J %*% V) * J))
  se <- c(dy$effects$se_effect, dy$effects$se_cumulative,
          dy$summary$se_long_run, dy$summary$se_mean_lag)
  expect_equal(se, se_num, tolerance = 1e-6)
  # the long run is the one of adl_long_run(), with either covariance
  for (ct in c("nonrobust", "HAC")) {
    lr <- adl_long_run(fit, cov_type = ct)
    dh <- adl_dynamics(fit, horizon = 3, cov_type = ct, nsim = 0)
    expect_equal(dh$summary$se_long_run, lr$std.error[lr$term %in% c("x1", "x2")])
  }
})

test_that("ADL(1,1) closed forms for the standard errors", {
  # d_0 = b0: se is the coefficient's; d_1 = a b0 + b1
  fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
  b <- coef(fit)
  V <- vcov(fit)
  dy <- adl_dynamics(fit, horizon = 2, nsim = 0)
  expect_equal(dy$effects$se_effect[1], sqrt(V["x", "x"]), ignore_attr = TRUE)
  g <- c(0, b[["x"]], b[["y_lag1"]], 1)[match(names(b),
         c("(Intercept)", "y_lag1", "x", "x_lag1"))]
  expect_equal(dy$effects$se_effect[2], sqrt(drop(t(g) %*% V %*% g)))
})

test_that("median-lag interval by simulation", {
  fit <- adl(y ~ x, data = adl_sim, order = c(1, 3))
  set.seed(3)
  d1 <- adl_dynamics(fit, nsim = 300)
  set.seed(3)
  d2 <- adl_dynamics(fit, nsim = 300)
  expect_identical(d1$summary, d2$summary)
  s <- d1$summary
  expect_true(is.integer(s$median_lag_low))
  expect_lte(s$median_lag_low, s$median_lag_high)
  expect_gte(s$median_lag_undefined, 0)
  expect_lte(s$median_lag_undefined, 1)
  d0 <- adl_dynamics(fit, nsim = 0)
  expect_true(is.na(d0$summary$median_lag_low))
  expect_false(anyNA(d0$summary$se_mean_lag))
})

test_that("se = FALSE keeps the point-estimate output", {
  fit <- dyn_fit()
  dy <- adl_dynamics(fit, horizon = 5, se = FALSE)
  expect_named(dy$effects, c("variable", "lag", "effect", "cumulative",
                             "pct_effect", "pct_cumulative"))
  expect_named(dy$summary, c("variable", "long_run", "mean_lag", "median_lag",
                             "horizon"))
  expect_true(any(grepl("no standard errors", capture.output(print(dy)))))
  expect_error(adl_dynamics(fit, se = FALSE, maxlags = 2), "only used")
  expect_error(adl_dynamics(fit, nsim = -1), "nsim")
  expect_error(adl_dynamics(fit, alpha = 2), "alpha")
})

test_that("unstable model: effects keep their standard errors, summaries are NA", {
  fit <- dyn_fit()
  fit$lm$coefficients[["y_lag1"]] <- 1.1
  expect_warning(du <- adl_dynamics(fit, horizon = 5), "not dynamically stable")
  expect_false(anyNA(du$effects$se_effect))
  expect_true(all(is.na(du$summary$se_long_run)))
  expect_true(all(is.na(du$summary$se_mean_lag)))
  expect_true(all(is.na(du$summary$median_lag_low)))
})

test_that("confidence bands are drawn and can be switched off", {
  fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
  dy <- adl_dynamics(fit, nsim = 50)
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(dy))
  expect_silent(plot(dy, type = "cumulative", ci = FALSE))
  expect_true(any(grepl("delta method", capture.output(print(dy)))))
})
