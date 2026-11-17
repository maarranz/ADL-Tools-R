test_that("critical value follows tsoutliers", {
  expect_equal(outlier_cval(30), 3)
  expect_equal(outlier_cval(50), 3)
  expect_equal(outlier_cval(450), 4)
  expect_equal(outlier_cval(1000), 4)
  expect_equal(outlier_cval(250), 3.5)
  expect_equal(outlier_cval(119), round(3 + 0.0025 * 69, 2))
})

test_that("clean data: nothing found, original model returned", {
  fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
  o <- adl_outliers(fit)
  expect_equal(nrow(o$outliers), 0)
  expect_identical(o$fit, fit)
  expect_length(o$impulses, 0)
  expect_output(print(o), "No impulse or step")
})

test_that("an innovation outlier is found as an impulse at its date", {
  set.seed(71)
  n <- 150
  x <- as.numeric(arima.sim(list(ar = 0.5), n))
  e <- rnorm(n)
  e[90] <- 6
  y <- numeric(n)
  for (t in 2:n) y[t] <- 1 + 0.5 * y[t - 1] + 0.4 * x[t] + 0.3 * x[t - 1] + e[t]
  d <- stats::ts(cbind(y = y, x = x), start = c(1980, 1), frequency = 4)
  o <- adl_outliers(adl(y ~ x, data = d, order = c(1, 1)), types = "impulse")
  expect_true("2002 Q2" %in% o$outliers$date)                  # observation 90
  expect_true(all(o$outliers$type == "impulse"))
  expect_true(all(abs(o$outliers$statistic) > o$cval))
})

test_that("a step is found, and the refit is the ADL with the dummy", {
  set.seed(3)
  n <- 160
  x <- as.numeric(arima.sim(list(ar = 0.5), n))
  y <- numeric(n)
  for (t in 2:n) y[t] <- 1 + 1.5 * (t >= 101) + 0.5 * y[t - 1] + 0.4 * x[t] +
    0.3 * x[t - 1] + rnorm(1)
  d <- stats::ts(cbind(y = y, x = x), start = c(1980, 1), frequency = 4)
  fit <- adl(y ~ x, data = d, order = c(1, 1))
  o <- adl_outliers(fit)
  expect_equal(o$outliers$date, "2005 Q1")
  expect_equal(o$outliers$type, "step")
  expect_equal(o$outliers$name, "S_2005Q1")
  # the same model by hand
  S <- as.numeric(time(d) >= 2005)
  d2 <- stats::ts(cbind(d, S_2005Q1 = S), start = c(1980, 1), frequency = 4)
  colnames(d2) <- c("y", "x", "S_2005Q1")
  ref <- adl(y ~ x | S_2005Q1, data = d2, order = c(1, 1))
  expect_equal(coef(o$fit), coef(ref))
  expect_equal(nobs(o$fit), nobs(fit))
  expect_equal(o$outliers$statistic, unname(summary(ref$lm)$coefficients["S_2005Q1", 3]))
  expect_length(o$steps, 1)
  # impulses only: the step is not a candidate
  expect_false(any(adl_outliers(fit, types = "impulse")$outliers$type == "step"))
})

test_that("an additive outlier in y shows as an impulse; isat agrees", {
  d <- adl_sim
  d[76, "y"] <- d[76, "y"] - 5
  fit <- adl(y ~ x, data = d, order = c(1, 1))
  o <- adl_outliers(fit)
  expect_true("2008 Q4" %in% o$outliers$date)
  skip_if_not_installed("gets")
  oi <- adl_outliers(fit, method = "isat")
  expect_true("2008 Q4" %in% oi$outliers$date)
  expect_equal(oi$method, "isat")
})

test_that("dates feed super_exogeneity_test(); Date index; errors", {
  set.seed(1)
  n <- 120
  x <- y <- numeric(n)
  for (t in 2:n) {
    x[t] <- 0.5 * x[t - 1] + 2 * (t >= 61) + rnorm(1)
    y[t] <- 1 + 0.5 * y[t - 1] + 0.4 * x[t] + 0.3 * x[t - 1] + rnorm(1)
  }
  d <- stats::ts(cbind(y = y, x = x), start = c(1990, 1), frequency = 4)
  ox <- adl_outliers(adl(x ~ y, data = d, order = c(1, 1), first_lag = 1))
  expect_equal(ox$outliers$date, "2005 Q1")
  fit <- adl(y ~ x, data = d, order = c(1, 1))
  s <- super_exogeneity_test(fit, impulses = ox$impulses, steps = ox$steps)
  expect_equal(s$statistic, super_exogeneity_test(fit, steps = c(2005, 1))$statistic)
  # Date index
  dd <- data.frame(t = seq(as.Date("1990-01-01"), by = "quarter", length.out = n),
                   y = y, x = x)
  od <- adl_outliers(adl(x ~ y, data = dd, time = "t", order = c(1, 1), first_lag = 1))
  expect_equal(od$outliers$date, "2005-01-01")
  expect_equal(od$outliers$name, "S_20050101")
  expect_error(adl_outliers(fit, cval = -1), "positive")
  expect_error(adl_outliers(fit, types = "AO"))
})
