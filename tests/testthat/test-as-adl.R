test_that("dynlm fits convert with identical estimates", {
  skip_if_not_installed("dynlm")
  suppressPackageStartupMessages(library(dynlm))
  m <- dynlm(y ~ L(y, 1:2) + L(x, 0:1) + w, data = adl_sim)
  fit <- as_adl(m)
  expect_s3_class(fit, "adl")
  expect_equal(unname(coef(fit)), unname(coef(m)))
  expect_equal(nobs(fit), nobs(m))

  m0 <- dynlm(y ~ L(y, 1) + L(x, 1) - 1, data = adl_sim, start = c(1995, 1))
  f0 <- as_adl(m0)
  expect_equal(unname(coef(f0)), unname(coef(m0)))
  expect_equal(format(f0$index[1]), format(zoo::as.yearqtr(stats::start(m0$residuals)[1] +
                                                           (stats::start(m0$residuals)[2] - 1) / 4)))
  expect_false(f0$intercept)
})

test_that("trend converts with a rescaled coefficient", {
  skip_if_not_installed("dynlm")
  suppressPackageStartupMessages(library(dynlm))
  m <- dynlm(y ~ L(y, 1) + L(x, 0:1) + trend(y), data = adl_sim)
  fit <- as_adl(m)
  b <- coef(fit); bm <- coef(m)
  expect_equal(unname(b[c("y_lag1", "x", "x_lag1")]), unname(bm[2:4]))
  expect_equal(b[["trend"]] * 4, unname(bm[5]))  # per quarter vs per year
})

test_that("unsupported terms and missing data", {
  skip_if_not_installed("dynlm")
  suppressPackageStartupMessages(library(dynlm))
  m <- dynlm(d(y) ~ L(d(y), 1) + L(x, 0:1), data = adl_sim)
  expect_error(as_adl(m), "Unsupported|left-hand")
  expect_error(as_adl(lm(y ~ x, data = as.data.frame(adl_sim))), "dynlm fit")
})
