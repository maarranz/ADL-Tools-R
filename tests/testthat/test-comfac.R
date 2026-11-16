sim_comfac <- function(n = 400, b1 = NULL, seed = 50) {
  # static regression with AR(1) errors: COMFAC holds when b1 is NULL
  set.seed(seed)
  x <- as.numeric(stats::arima.sim(list(ar = 0.5), n))
  u <- as.numeric(stats::filter(rnorm(n), 0.6, method = "recursive"))
  y <- 1 + 0.8 * x + u
  if (!is.null(b1)) {
    y <- numeric(n); e <- rnorm(n)
    for (t in 2:n) y[t] <- 1 + 0.6 * y[t - 1] + 0.8 * x[t] + b1 * x[t - 1] + e[t]
  }
  stats::ts(cbind(y = y, x = x))
}

test_that("ADL(1,1): restricted fit equals nls and Wald follows ADL2", {
  d <- sim_comfac()
  fit <- adl(y ~ x, data = d, order = c(1, 1))
  lr <- comfac_test(fit)
  df <- fit$lm$model
  ref <- nls(y ~ c0 + w * y_lag1 + b0 * x - w * b0 * x_lag1, data = df,
             start = list(c0 = 0.4, w = 0.5, b0 = 0.7))
  expect_equal(lr$rss_restricted, sum(residuals(ref)^2), tolerance = 1e-7)
  expect_equal(lr$omega, unname(coef(ref)["w"]), tolerance = 1e-5)
  expect_equal(lr$statistic, nobs(fit) * log(lr$rss_restricted / lr$rss_unrestricted))
  expect_equal(lr$df, 1)
  expect_false(lr$reject_null)

  w <- comfac_test(fit, method = "Wald")
  b <- coef(fit); V <- vcov(fit)
  g <- b[["y_lag1"]] * b[["x"]] + b[["x_lag1"]]
  G <- c(0, b[["x"]], b[["y_lag1"]], 1)
  expect_equal(w$statistic, g^2 / as.numeric(t(G) %*% V %*% G))
  expect_false(w$reject_null)
  expect_equal(comfac_test(fit, method = "Wald", cov_type = "HAC")$cov_type, "HAC")
})

test_that("COMFAC is rejected when it does not hold", {
  d <- sim_comfac(b1 = 0.4)  # alpha_1 * beta_0 + beta_1 = 0.88, far from 0
  fit <- adl(y ~ x, data = d, order = c(1, 1))
  expect_true(comfac_test(fit)$reject_null)
  expect_true(comfac_test(fit, method = "Wald")$reject_null)
})

test_that("order 2 with two regressors: r * n restrictions and a valid fit", {
  set.seed(51)
  n <- 500
  x1 <- rnorm(n); x2 <- rnorm(n)
  u <- as.numeric(stats::filter(rnorm(n), c(0.5, 0.3), method = "recursive"))
  y <- 0.5 * x1 + 0.3 * c(0, x1[-n]) - 0.4 * x2 + u
  d <- stats::ts(cbind(y = y, x1 = x1, x2 = x2))
  # x1 enters at lags 0 and 1, so with r = 2 the ADL needs q = 3: in an
  # ADL(2,2) delta(L) has degree 0 and the common factor rightly fails
  expect_true(comfac_test(adl(y ~ x1 + x2, data = d, order = c(2, 2)), r = 2)$reject_null)
  fit <- adl(y ~ x1 + x2, data = d, order = c(3, 3))
  t2 <- comfac_test(fit, r = 2)
  expect_equal(t2$df, 4)
  expect_gte(t2$rss_restricted, t2$rss_unrestricted)
  expect_false(t2$reject_null)
  expect_equal(t2$omega, c(0.5, 0.3), tolerance = 0.35)
  # the restricted coefficients reproduce the restricted RSS ...
  X <- model.matrix(fit$lm); yv <- residuals(fit$lm) + fitted(fit$lm)
  map <- comfac_map(fit, 2, colnames(X))
  rss <- function(th) sum((yv - X %*% map$expand(th))^2)
  expect_equal(rss(t2$restricted$estimate), t2$rss_restricted)
  # ... and no full (unconcentrated) NLS from random starts does better
  set.seed(52)
  best <- min(vapply(1:10, function(i) {
    optim(rnorm(length(map$labels), 0, 0.5), rss, method = "BFGS",
          control = list(maxit = 5000, reltol = 1e-14))$value
  }, numeric(1)))
  expect_gte(best, t2$rss_restricted - 1e-6)
})

test_that("polynomial map: C(L) = omega(L) phi(L), B(L) = omega(L) delta(L)", {
  d <- sim_comfac()
  fit <- adl(y ~ x, data = d, order = c(2, 2))
  map <- comfac_map(fit, 1, names(coef(fit)))
  expect_equal(map$labels, c("c", "phi_1", "delta_x_0", "delta_x_1", "omega_1"))
  # phi = 0.3, delta = (2, 1), omega = 0.5:
  # C = (1 - .5L)(1 - .3L) = 1 - .8L + .15L^2 -> alpha = (.8, -.15)
  # B = (1 - .5L)(2 + L) = 2 + 0L - .5L^2
  full <- map$expand(c(1, 0.3, 2, 1, 0.5))
  expect_equal(full, c(1, 0.8, -0.15, 2, 0, -0.5))
})

test_that("argument checks", {
  d <- sim_comfac()
  expect_error(comfac_test(adl(y ~ x, data = d, order = c(1, 1)), r = 2), "p >= 2")
  expect_error(comfac_test(adl(y ~ x, data = d, order = c(1, 0))), "consecutive lags")
  expect_error(comfac_test(adl(y ~ x, data = d, order = c(2, 2)), method = "Wald"),
               "use method")
  expect_output(print(comfac_test(adl(y ~ x, data = d, order = c(1, 1)))),
                "alpha_1 \\* beta_0 \\+ beta_1 = 0")
})
