tidy_fit <- function() {
  set.seed(14)
  n <- 150
  x <- cumsum(rnorm(n))
  y <- stats::filter(0.5 * x + rnorm(n), 0.6, method = "recursive")
  adl(y ~ x, data = stats::ts(cbind(y = as.numeric(y), x = x)), order = c(1, 1))
}

test_that("tidy and glance for fits", {
  fit <- tidy_fit()
  td <- tidy(fit)
  expect_s3_class(td, "data.frame")
  expect_equal(td$estimate, unname(coef(fit)))
  expect_equal(td$param, c("alpha_0", "alpha_1", "beta_0", "beta_1"))
  expect_null(attr(td, "cov_type"))
  # conf.level and pass-through arguments
  td90 <- tidy(fit, conf.level = 0.9, dist = "t")
  ci <- stats::confint(fit$lm, level = 0.9)
  expect_equal(td90$conf.low, unname(ci[, 1]))
  expect_equal(tidy(fit, cov_type = "HAC")$std.error,
               unname(sqrt(diag(adl_vcov(fit, "HAC")))))

  gl <- glance(fit)
  expect_equal(nrow(gl), 1)
  expect_equal(gl$nobs, nobs(fit))
  expect_equal(gl$AIC, AIC(fit))
  expect_equal(gl$df, 3)
  expect_equal(gl$statistic, adl_joint_test(fit, terms = c("y_lag1", "x", "x_lag1"))$statistic)
})

test_that("tidy for tests, Granger and dynamics", {
  fit <- tidy_fit()
  tt <- tidy(adl_joint_test(fit, variables = "x", cov_type = "HAC"))
  expect_equal(tt$distribution, "chi2")
  expect_equal(tt$cov_type, "HAC")

  set.seed(15)
  d <- stats::ts(cbind(y = rnorm(200), x = rnorm(200)))
  g <- granger_test("y", "x", p = 2, data = d)
  expect_equal(nrow(tidy(g)), 4)
  gg <- glance(g)
  expect_equal(gg$symbol, g$classification$symbol)

  dy <- adl_dynamics(fit, horizon = 5)
  expect_equal(nrow(tidy(dy)), 6)
  expect_true(glance(dy)$stable)
})
