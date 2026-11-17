sim_se <- function(lucas, n = 200, seed = 61) {
  set.seed(seed)
  x <- numeric(n); y <- numeric(n)
  for (t in 2:n) {
    x[t] <- 0.5 * x[t - 1] + 2 * (t > n / 2) + rnorm(1)
    b0 <- if (lucas && t > n / 2) 1.2 else 0.4
    y[t] <- 1 + 0.5 * y[t - 1] + b0 * x[t] + 0.3 * x[t - 1] + rnorm(1)
  }
  stats::ts(cbind(y = y, x = x), start = c(1970, 1), frequency = 4)
}
brk <- c(1995, 1)   # observation 101

test_that("invariant model: marginal break found, invariance not rejected", {
  fit <- adl(y ~ x, data = sim_se(FALSE), order = c(1, 1))
  s <- super_exogeneity_test(fit, steps = brk)
  expect_equal(s$coefficients$reject, c(TRUE, FALSE, FALSE))
  expect_false(s$reject_null)
  expect_equal(s$dummies$name, "S(1995 Q1)")
  expect_equal(s$dummies$type, "step")
})

test_that("Lucas critique: the conditional model breaks with the marginal one", {
  fit <- adl(y ~ x, data = sim_se(TRUE), order = c(1, 1))
  s <- super_exogeneity_test(fit, steps = brk)
  expect_true(s$coefficients$reject[1])
  expect_true(s$reject_null)
})

test_that("statistics match the augmented regressions", {
  d <- sim_se(FALSE)
  fit <- adl(y ~ x, data = d, order = c(1, 1))
  s <- super_exogeneity_test(fit, steps = brk, impulses = list(c(1990, 2), c(1999, 4)))
  expect_equal(nrow(s$dummies), 3)
  expect_equal(s$coefficients$df, c(3, 3, 1))
  # step 2 by hand: ADL + dummies, Wald on the dummies
  X <- model.matrix(fit$lm)
  idx <- as.numeric(fit$index)
  D <- cbind(as.numeric(abs(idx - 1990.25) < 1e-8), as.numeric(abs(idx - 1999.75) < 1e-8),
             as.numeric(idx >= 1995 - 1e-8))
  aug <- lm(model.response(model.frame(fit$lm)) ~ 0 + X + D)
  b <- coef(aug)[5:7]
  w <- as.numeric(t(b) %*% solve(vcov(aug)[5:7, 5:7], b))
  expect_equal(s$statistic, w)
  expect_equal(s$coefficients$statistic[2], w)
  # step 3: kappa from the marginal regression with the dummies
  m <- s$marginal
  expect_equal(s$coefficients$statistic[3],
               unname(coef(m)["ec_lag1"]^2 / vcov(m)["ec_lag1", "ec_lag1"]))
  # F form
  sf <- super_exogeneity_test(fit, steps = brk, impulses = list(c(1990, 2), c(1999, 4)),
                              dist = "F")
  expect_equal(sf$statistic, w / 3)
  expect_equal(sf$df_num, 3)
  expect_equal(sf$distribution, "F")
  expect_equal(super_exogeneity_test(fit, steps = brk, cov_type = "HAC")$cov_type, "HAC")
})

test_that("warnings when the test has nothing to work with or weak exogeneity fails", {
  fit <- adl(y ~ x, data = sim_se(FALSE), order = c(1, 1))
  expect_warning(super_exogeneity_test(fit, steps = brk, alpha = 1e-12),
                 "no change at these dates")
  expect_warning(super_exogeneity_test(fit, steps = brk, alpha = 0.99),
                 "Weak exogeneity is rejected")
})

test_that("dates: forms, errors and collinearity", {
  fit <- adl(y ~ x, data = sim_se(FALSE), order = c(1, 1))
  expect_error(super_exogeneity_test(fit), "impulses")
  expect_error(super_exogeneity_test(fit, steps = c(2030, 1)), "not in")
  expect_error(super_exogeneity_test(fit, impulses = c(2019, 4), steps = c(2019, 4)),
               "collinear")
  expect_error(super_exogeneity_test(fit, steps = list(brk, brk)), "repeated")
  # Date index, several dates as a vector
  d <- sim_se(FALSE)
  dd <- data.frame(t = seq(as.Date("1970-01-01"), by = "quarter", length.out = 200),
                   y = as.numeric(d[, "y"]), x = as.numeric(d[, "x"]))
  fd <- adl(y ~ x, data = dd, time = "t", order = c(1, 1))
  s <- super_exogeneity_test(fd, steps = as.Date("1995-01-01"),
                             impulses = as.Date(c("1980-01-01", "1985-04-01")))
  expect_equal(s$dummies$date, c("1980-01-01", "1985-04-01", "1995-01-01"))
  expect_equal(s$statistic,
               super_exogeneity_test(fit, steps = brk,
                                     impulses = list(c(1980, 1), c(1985, 2)))$statistic)
})
