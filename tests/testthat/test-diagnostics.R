diag_fit <- function() adl(y ~ x, data = adl_sim, order = c(1, 1))

test_that("Breusch-Godfrey matches lmtest::bgtest", {
  skip_if_not_installed("lmtest")
  fit <- diag_fit()
  for (s in c(1, 4)) {
    b <- bg_test(fit, order = s)
    ref <- lmtest::bgtest(fit$lm, order = s, type = "Chisq", fill = 0)
    expect_equal(b$statistic, unname(ref$statistic))
    expect_equal(b$p.value, unname(ref$p.value))
    bf <- bg_test(fit, order = s, dist = "F")
    reff <- lmtest::bgtest(fit$lm, order = s, type = "F", fill = 0)
    expect_equal(bf$statistic, unname(reff$statistic))
    expect_equal(bf$df_denom, unname(reff$parameter[2]))
  }
})

test_that("BG detects serial correlation that is there", {
  set.seed(30)
  n <- 300
  x <- rnorm(n)
  e <- as.numeric(stats::filter(rnorm(n), 0.6, method = "recursive"))
  d <- stats::ts(cbind(y = 1 + 0.5 * x + e, x = x))
  fit <- adl(y ~ x, data = d, order = c(0, 0))
  expect_true(bg_test(fit, order = 1)$reject_null)
})

test_that("ARCH-LM is nR^2 of the squared-residual regression", {
  fit <- diag_fit()
  a <- arch_test(fit, lags = 2)
  e2 <- residuals(fit$lm)^2
  n <- length(e2)
  aux <- lm(e2[3:n] ~ e2[2:(n - 1)] + e2[1:(n - 2)])
  expect_equal(a$statistic, (n - 2) * summary(aux)$r.squared)
  expect_equal(a$df, 2)
})

test_that("RESET matches lmtest::resettest (F) and takes HAC", {
  skip_if_not_installed("lmtest")
  fit <- diag_fit()
  r <- reset_test(fit, dist = "F")
  ref <- lmtest::resettest(fit$lm, power = 2:3, type = "fitted")
  expect_equal(r$statistic, unname(ref$statistic))
  expect_equal(r$p.value, unname(ref$p.value))
  rc <- reset_test(fit)
  expect_equal(rc$statistic, 2 * r$statistic)
  rh <- reset_test(fit, power = 2, cov_type = "HAC")
  expect_equal(rh$cov_type, "HAC")
  expect_false(is.na(rh$maxlags))
})

test_that("Jarque-Bera from skewness and kurtosis", {
  fit <- diag_fit()
  j <- jb_test(fit)
  u <- residuals(fit$lm); e <- u - mean(u); n <- length(u)
  S <- mean(e^3) / mean(e^2)^1.5; K <- mean(e^4) / mean(e^2)^2
  expect_equal(j$statistic, n / 6 * (S^2 + (K - 3)^2 / 4))
  expect_equal(j$skewness, S)
  set.seed(31)
  d <- stats::ts(cbind(y = rexp(300), x = rnorm(300)))
  expect_true(jb_test(adl(y ~ x, data = d, order = c(0, 0)))$reject_null)
})

test_that("diagnostics battery", {
  fit <- diag_fit()
  dg <- adl_diagnostics(fit, bg_order = 4)
  expect_equal(nrow(dg), 4)
  expect_equal(dg$statistic[1], bg_test(fit, order = 4)$statistic)
  expect_output(print(dg), "Breusch-Godfrey \\(4\\)")
  expect_error(bg_test(fit, order = 0), "integer >= 1")
  expect_error(reset_test(fit, power = 1), ">= 2")
})
