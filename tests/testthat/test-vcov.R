vc_fit <- function(n = 120) {
  set.seed(4)
  x <- cumsum(rnorm(n))
  y <- stats::filter(0.5 * x + rnorm(n), 0.6, method = "recursive")
  adl(y ~ x, data = stats::ts(cbind(y = as.numeric(y), x = x)), order = c(1, 1))
}

test_that("nonrobust covariance is the OLS one", {
  fit <- vc_fit()
  V <- adl_vcov(fit)
  expect_equal(unclass(V)[, ], stats::vcov(fit$lm), ignore_attr = TRUE)
  expect_true(is.na(attr(V, "maxlags")))
  expect_equal(rownames(V), names(coef(fit)))
})

test_that("HAC matches sandwich::NeweyWest with the forecastingtools settings", {
  fit <- vc_fit()
  n <- nobs(fit)
  V <- adl_vcov(fit, "HAC")
  L <- floor(4 * (n / 100)^(2 / 9))
  expect_equal(attr(V, "maxlags"), L)
  ref <- sandwich::NeweyWest(fit$lm, lag = L, prewhite = FALSE, adjust = FALSE)
  expect_equal(unclass(V)[, ], ref, ignore_attr = TRUE)

  # Hand computation: (X'X)^-1 S (X'X)^-1, Bartlett weights, divisor n
  X <- model.matrix(fit); u <- residuals(fit$lm)
  g <- X * u
  S <- crossprod(g)
  for (j in seq_len(L)) {
    G <- crossprod(g[-seq_len(j), , drop = FALSE], g[seq_len(n - j), , drop = FALSE])
    S <- S + (1 - j / (L + 1)) * (G + t(G))
  }
  B <- solve(crossprod(X))
  expect_equal(unclass(V)[, ], B %*% S %*% B, ignore_attr = TRUE)

  expect_equal(attr(adl_vcov(fit, "HAC", maxlags = 2), "maxlags"), 2L)
})

test_that("maxlags is validated", {
  fit <- vc_fit()
  expect_error(adl_vcov(fit, maxlags = 2), "only used with HAC")
  expect_error(adl_vcov(fit, "HAC", maxlags = -1), "non-negative")
  expect_error(adl_vcov(fit, "HAC", maxlags = 500), "smaller than")
  expect_error(adl_vcov(lm(1 ~ 1)), "adl\\(\\)")
})

test_that("test results print in the suite format", {
  w <- wald_inference(6, 2, 50, "nonrobust")
  expect_equal(w$statistic, 3)
  expect_equal(w$p.value, stats::pf(3, 2, 50, lower.tail = FALSE))
  expect_equal(wald_inference(6, 2, 50, "HAC")$distribution, "chi2")

  r <- new_adl_test(c(list(method = "demo", null_hypothesis = "b = 0", alpha = 0.05,
                           reject_null = w$p.value < 0.05, nobs = 52,
                           cov_type = "nonrobust"), w), "Demo test")
  expect_s3_class(r, "adltools_test")
  out <- capture.output(print(r))
  expect_match(out[1], "Demo test")
  expect_true(any(grepl("\\[F\\(2, 50\\)\\]", out)))
  expect_true(any(grepl("conventional OLS", out)))
})
