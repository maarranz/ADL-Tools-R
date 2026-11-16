# Simulate a bivariate VAR(1) with a chosen causal structure.
# a_xy: effect of x_{t-1} on y_t; a_yx: effect of y_{t-1} on x_t;
# rho: correlation of the innovations (instantaneous causality).
sim_var <- function(a_xy = 0, a_yx = 0, rho = 0, n = 400, seed = 1) {
  set.seed(seed)
  e1 <- rnorm(n); e2 <- rho * e1 + sqrt(1 - rho^2) * rnorm(n)
  y <- x <- numeric(n)
  for (t in 2:n) {
    y[t] <- 0.4 * y[t - 1] + a_xy * x[t - 1] + e1[t]
    x[t] <- 0.4 * x[t - 1] + a_yx * y[t - 1] + e2[t]
  }
  stats::ts(cbind(y = y, x = x))
}

test_that("Granger F equals lmtest::grangertest on the same sample", {
  skip_if_not_installed("lmtest")
  d <- sim_var(a_xy = 0.3)
  g <- granger_test("y", "x", p = 3, data = d)
  ref_xy <- lmtest::grangertest(y ~ x, order = 3, data = d)
  ref_yx <- lmtest::grangertest(x ~ y, order = 3, data = d)
  expect_equal(g$tests$F[1], ref_xy$F[2])
  expect_equal(g$tests$p.value_F[1], ref_xy$`Pr(>F)`[2])
  expect_equal(g$tests$F[2], ref_yx$F[2])
  expect_equal(g$nobs, nrow(d) - 3)
  expect_equal(g$tests$chisq, g$tests$F * g$tests$df_num)
})

test_that("with s = r and OLS the two instantaneous statistics are identical", {
  d <- sim_var(a_xy = 0.3, rho = 0.4)
  g <- granger_test("y", "x", p = 2, data = d)
  expect_equal(g$tests$F[3], g$tests$F[4])
  expect_false(g$classification$ambiguous)
  gh <- granger_test("y", "x", p = 2, data = d, cov_type = "HAC")
  expect_false(isTRUE(all.equal(gh$tests$chisq[3], gh$tests$chisq[4])))
})

test_that("the general function reproduces granger_test with equal orders", {
  d <- sim_var(a_yx = 0.3)
  g1 <- granger_test("y", "x", p = 2, data = d)
  g2 <- granger_test_general("y", "x", data = d,
                             orders = list(y = c(2, 2), x = c(2, 2)))
  expect_equal(g1$tests, g2$tests)
  expect_equal(g1$classification, g2$classification)
})

test_that("models use a common sample and the right regressors", {
  d <- sim_var()
  g <- granger_test_general("y", "x", data = d,
                            orders = list(y = c(1, 3), x = c(2, 1)),
                            inst_orders = list(y = c(1, 1), x = c(4, 1)))
  n <- vapply(g$models, nobs, numeric(1))
  expect_true(all(n == nrow(d) - 4))
  expect_equal(names(coef(g$models$granger_y)),
               c("(Intercept)", "y_lag1", "x_lag1", "x_lag2", "x_lag3"))
  expect_equal(names(coef(g$models$inst_x)),
               c("(Intercept)", "x_lag1", "x_lag2", "x_lag3", "x_lag4", "y", "y_lag1"))
  expect_equal(g$tests$df_num, c(3, 1, 1, 1))
})

test_that("the eight causal types are recovered from simulated VARs", {
  cases <- list(
    list(a_xy = 0,   a_yx = 0,   rho = 0,   type = 1L),
    list(a_xy = 0,   a_yx = 0,   rho = 0.5, type = 2L),
    list(a_xy = 0.4, a_yx = 0,   rho = 0,   type = 3L),
    list(a_xy = 0,   a_yx = 0.4, rho = 0,   type = 4L),
    list(a_xy = 0.4, a_yx = 0,   rho = 0.5, type = 5L),
    list(a_xy = 0,   a_yx = 0.4, rho = 0.5, type = 6L),
    list(a_xy = 0.4, a_yx = 0.4, rho = 0,   type = 7L),
    list(a_xy = 0.4, a_yx = 0.4, rho = 0.5, type = 8L)
  )
  for (cs in cases) {
    d <- sim_var(cs$a_xy, cs$a_yx, cs$rho, n = 1000, seed = 11)
    g <- granger_test("y", "x", p = 1, data = d, alpha = 0.01)
    expect_equal(g$classification$type, cs$type, info = paste("type", cs$type))
  }
})

test_that("symbols and descriptions use the variable names", {
  d <- sim_var(a_xy = 0.4, n = 1000, seed = 11)
  colnames(d) <- c("cons", "inc")
  g <- granger_test("cons", "inc", p = 1, data = d, alpha = 0.01)
  expect_equal(g$classification$symbol, "(inc -> cons)")
  expect_equal(g$classification$description, "inc causes cons; no instantaneous causality")
  out <- capture.output(print(g))
  expect_true(any(grepl("type 3 of 8", out)))
  expect_true(any(grepl("assumed stationary", out)))
})

test_that("ambiguous classification when the instantaneous tests disagree", {
  cl <- classify_granger(c(TRUE, FALSE, TRUE, FALSE), "y", "x")
  expect_true(cl$ambiguous)
  expect_equal(cl$type, c(3L, 5L))
  expect_equal(cl$symbol, c("(x -> y)", "(x => y)"))
})

test_that("series can be passed directly, and dist switches the decisions", {
  d <- sim_var(a_xy = 0.3)
  y <- d[, "y"]; x <- d[, "x"]
  g1 <- granger_test(y, x, p = 2)
  g2 <- granger_test("y", "x", p = 2, data = d)
  expect_equal(g1$tests$F, g2$tests$F)
  g3 <- granger_test(as.numeric(y), as.numeric(x), p = 2)
  expect_equal(g3$tests$F, g2$tests$F)
  expect_equal(g2$dist, "chisq")
  expect_equal(g2$tests$reject_null, g2$tests$p.value_chisq < 0.05)
  g4 <- granger_test("y", "x", p = 2, data = d, dist = "F")
  expect_equal(g4$tests$reject_null, g4$tests$p.value_F < 0.05)
})

test_that("argument errors", {
  d <- sim_var()
  expect_error(granger_test("y", "x", p = 0, data = d), ">= 1")
  expect_error(granger_test("y", "y", p = 1, data = d), "different")
  expect_error(granger_test("y", "x", p = 1), "no `data`")
  expect_error(granger_test_general("y", "x", data = d, orders = list(y = c(1, 0), x = c(1, 1))),
               "r >= 1")
  expect_error(granger_test_general("y", "x", data = d,
                                    orders = list(y = c(1, 1), x = c(1, 1)), select = "AIC"),
               "not yet implemented")
})
