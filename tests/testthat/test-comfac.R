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
  expect_output(print(comfac_test(adl(y ~ x, data = d, order = c(1, 1)))),
                "alpha_1 \\* beta_0 \\+ beta_1 = 0")
})

test_that("COMFAC: Wald needs no optimizer; RNG state is preserved (audit)", {
  d <- sim_comfac()
  fit <- adl(y ~ x, data = d, order = c(1, 1))
  w <- comfac_test(fit, method = "Wald")
  expect_null(w$restricted)
  set.seed(123); a <- runif(1)
  set.seed(123)
  f3 <- adl(y ~ x, data = d, order = c(3, 3))
  invisible(comfac_test(f3, r = 3))
  expect_equal(runif(1), a)
})

# Sargan's determinantal Wald form -------------------------------------------

pmul <- function(a, b) {
  o <- numeric(length(a) + length(b) - 1)
  for (i in seq_along(a)) o[i:(i + length(b) - 1)] <- o[i:(i + length(b) - 1)] + a[i] * b
  o
}

test_that("Wald: the ADL(1,1) formula of ADL2, with one and two regressors", {
  set.seed(3)
  n <- 300
  x1 <- rnorm(n); x2 <- rnorm(n)
  y <- as.numeric(stats::filter(0.6 * x1 + 0.2 * c(0, x1[-n]) + 0.3 * x2 + rnorm(n),
                                0.5, method = "recursive"))
  d <- stats::ts(cbind(y = y, x1 = x1, x2 = x2))
  for (f in list(adl(y ~ x1, data = d, order = c(1, 1)),
                 adl(y ~ x1 + x2, data = d, order = c(1, 1, 1)))) {
    b <- coef(f); V <- vcov(f); a1 <- b[["y_lag1"]]
    g <- vapply(f$x, function(v) a1 * b[[v]] + b[[paste0(v, "_lag1")]], 1)
    G <- t(vapply(f$x, function(v) {
      r <- stats::setNames(numeric(length(b)), names(b))
      r["y_lag1"] <- b[[v]]; r[v] <- a1; r[paste0(v, "_lag1")] <- 1
      r
    }, numeric(length(b))))
    W <- drop(t(g) %*% solve(G %*% V %*% t(G), g))
    w <- comfac_test(f, method = "Wald")
    expect_equal(w$statistic, W)
    expect_equal(w$df, length(f$x))
  }
})

test_that("Wald: one regressor, two polynomials, r = 1 is the Sylvester resultant", {
  d <- sim_comfac()
  fit <- adl(y ~ x, data = d, order = c(2, 2))
  b <- coef(fit)
  C <- c(1, -b[["y_lag1"]], -b[["y_lag2"]])
  B <- c(b[["x"]], b[["x_lag1"]], b[["x_lag2"]])
  S <- rbind(c(C, 0), c(B, 0), c(0, C), c(0, B))
  w <- comfac_test(fit, method = "Wald")
  expect_equal(abs(w$constraints), abs(det(S)))
})

test_that("Wald: constraints vanish at a common factor, Jacobian of full rank", {
  set.seed(5)
  N <- 120
  base <- stats::ts(cbind(y = rnorm(N), x1 = rnorm(N), x2 = rnorm(N)))
  specs <- list(
    list(lags = list(y = 1:3, x1 = 0:2, x2 = 0:3), r = 1),      # type 1, unequal
    list(lags = list(y = 1:3, x1 = 0:2, x2 = 0:3), r = 2),
    list(lags = list(y = 1:3, x1 = 1:2, x2 = 0:1), r = 1),      # type 2, dead start
    list(lags = list(y = 1:2, x1 = 0:2, x2 = 0:2), r = 2))
  for (sp in specs) {
    fit <- adl(y ~ x1 + x2, data = base, lags = sp$lags)
    b <- coef(fit)
    r <- sp$r
    st <- sargan_structure(fit, r, names(b))
    om <- 1
    for (h in seq_len(r)) om <- pmul(om, c(1, -runif(1, -0.8, 0.8)))
    p <- length(sp$lags$y)
    C <- pmul(om, c(1, -runif(p - r, -0.4, 0.4)))
    bb <- b
    bb[paste0("y_lag", seq_len(p))] <- -C[-1]
    for (v in c("x1", "x2")) {
      l <- sp$lags[[v]]
      bb[ifelse(l == 0, v, paste0(v, "_lag", l))] <- pmul(om, rnorm(length(l) - r))
    }
    con <- sargan_constraints(bb, st)
    expect_length(con$phi, 2 * r)
    expect_lt(max(abs(con$phi)), 1e-12)
    expect_equal(qr(con$G)$rank, 2 * r)
    # exact derivatives
    Gn <- vapply(seq_along(b), function(j) {
      e <- numeric(length(b)); e[j] <- 1e-6
      (sargan_constraints(bb + e, st)$phi - sargan_constraints(bb - e, st)$phi) / 2e-6
    }, numeric(2 * r))
    expect_equal(con$G, matrix(Gn, nrow = 2 * r), tolerance = 1e-6)
  }
})

test_that("Wald: general orders, HAC, and the order of the restrictions", {
  d <- sim_comfac()
  fit <- adl(y ~ x, data = d, order = c(3, 3))
  for (r in 1:3) {
    w <- comfac_test(fit, r = r, method = "Wald")
    expect_equal(w$df, r)
    expect_length(w$constraints, r)
  }
  expect_equal(comfac_test(fit, r = 2, method = "Wald", cov_type = "HAC")$cov_type, "HAC")
})
