stab_fit <- function() adl(y ~ x, data = adl_sim, order = c(1, 1))

test_that("recursive residuals match strucchange and add up to the RSS", {
  fit <- stab_fit()
  rr <- recursive_residuals(fit)
  expect_equal(nrow(rr), nobs(fit) - 4)
  expect_equal(sum(rr$residual^2), sum(residuals(fit$lm)^2))
  skip_if_not_installed("strucchange")
  expect_equal(rr$residual, unname(strucchange::recresid(fit$lm)))
})

test_that("CUSUM: EViews scaling, BDE p-value, strucchange equivalence", {
  fit <- stab_fit()
  cs <- cusum_test(fit)
  expect_equal(cs$critical, 0.948, tolerance = 1e-3)
  expect_equal(cusum_test(fit, alpha = 0.01)$critical, 1.143, tolerance = 1e-3)
  w <- recursive_residuals(fit)$residual
  m <- length(w)
  expect_equal(cs$path$value[m], sum(w) / sqrt(sum(w^2) / m))
  expect_equal(cs$path$upper[m], 3 * cs$critical * sqrt(m))
  skip_if_not_installed("strucchange")
  # strucchange scales by the centred standard deviation of w
  sc <- strucchange::sctest(strucchange::efp(formula(fit$lm), data = model.frame(fit$lm),
                                             type = "Rec-CUSUM"))
  ratio <- sqrt(sum(w^2) / m) / sd(w)
  expect_equal(unname(sc$statistic), cs$statistic * ratio)
  expect_equal(unname(sc$p.value), cusum_pvalue(cs$statistic * ratio), tolerance = 1e-3)
})

test_that("CUSUM-of-squares: bounds, reproducible, RNG untouched", {
  fit <- stab_fit()
  set.seed(1)
  before <- .Random.seed
  cq <- cusumsq_test(fit, nsim = 2000)
  expect_identical(.Random.seed, before)
  expect_identical(cq$p.value, cusumsq_test(fit, nsim = 2000)$p.value)
  expect_equal(tail(cq$path$value, 1), 1)
  expect_gt(cq$critical, 0)
  expect_lt(cq$critical, 0.5)
})

test_that("the tests detect breaks they are designed for", {
  set.seed(51)
  n <- 200
  x <- as.numeric(arima.sim(list(ar = 0.5), n))
  e <- rnorm(n)
  # coefficient break halfway
  y <- numeric(n)
  for (t in 2:n) {
    b0 <- if (t <= 100) 0.4 else 1.6
    y[t] <- 1 + 0.5 * y[t - 1] + b0 * x[t] + 0.3 * x[t - 1] + e[t]
  }
  fb <- adl(y ~ x, data = stats::ts(cbind(y = y, x = x)), order = c(1, 1))
  h <- hansen_test(fb)
  expect_true(h$reject_null)
  expect_true(h$coefficients$reject[h$coefficients$term == "x"])
  # variance break
  y2 <- numeric(n)
  sdv <- ifelse(seq_len(n) <= 100, 1, 3)
  for (t in 2:n) y2[t] <- 1 + 0.5 * y2[t - 1] + 0.4 * x[t] + 0.3 * x[t - 1] + sdv[t] * e[t]
  fv <- adl(y ~ x, data = stats::ts(cbind(y = y2, x = x)), order = c(1, 1))
  cq <- cusumsq_test(fv, nsim = 2000)
  expect_true(cq$reject_null)
  expect_false(is.na(cq$crossing))
  hv <- hansen_test(fv)
  expect_true(hv$coefficients$reject[hv$coefficients$term == "variance"])
  # intercept shift in a distributed-lag model: CUSUM (with a lagged y the
  # shift is partly absorbed by alpha_1, and CUSUM has little power there)
  y3 <- (seq_len(n) > 100) * 1.5 + 0.4 * x + 0.3 * c(0, x[-n]) + e
  f3 <- adl(y ~ x, data = stats::ts(cbind(y = y3, x = x)), order = c(0, 1))
  expect_true(cusum_test(f3)$reject_null)
})

test_that("Hansen: statistics, asymptotic critical values, strucchange", {
  fit <- stab_fit()
  h <- hansen_test(fit)
  expect_equal(nrow(h$coefficients), 5)
  expect_equal(h$df, 5)
  X <- model.matrix(fit$lm); e <- residuals(fit$lm); n <- length(e)
  f1 <- X[, "x"] * e
  expect_equal(h$coefficients$statistic[3], sum(cumsum(f1)^2) / (n * sum(f1^2)))
  # Cramer-von Mises quantiles vs Hansen's (1992) table
  q5 <- function(d) uniroot(function(x) cvm_pvalue(x, d) - 0.05, c(0.1, 10))$root
  expect_equal(q5(1), 0.4614, tolerance = 1e-3)
  expect_equal(q5(2), 0.749, tolerance = 5e-3)
  expect_equal(q5(5), 1.47, tolerance = 5e-3)
  skip_if_not_installed("strucchange")
  g <- strucchange::gefp(fit$lm, fit = NULL)
  ref <- strucchange::sctest(g, functional = strucchange::meanL2BB)
  f <- X * e; S <- apply(f, 2, cumsum)
  expect_equal(unname(ref$statistic), sum(S * t(solve(crossprod(f), t(S)))) / n)
})

test_that("stability plots draw", {
  fit <- stab_fit()
  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(cusum_test(fit)))
  expect_invisible(plot(cusumsq_test(fit, nsim = 500)))
})

test_that("CUSUM of squares: decision agrees with the p-value (audit 3)", {
  fit <- stab_fit()
  expect_error(cusumsq_test(fit, nsim = 1), "too small")
  expect_error(cusumsq_test(fit, nsim = 19), "too small")
  for (ns in c(20, 39, 99, 500)) for (a in c(0.05, 0.1)) {
    r <- cusumsq_test(fit, nsim = ns, alpha = a)
    expect_identical(r$reject_null, r$p.value < a)
    expect_identical(r$reject_null, r$statistic > r$critical)
  }
  # and on a model that does reject
  set.seed(52)
  n <- 200
  x <- rnorm(n)
  y <- 1 + 0.5 * x + ifelse(seq_len(n) <= 100, 1, 4) * rnorm(n)
  fv <- adl(y ~ x, data = stats::ts(cbind(y = y, x = x)), order = c(0, 0))
  for (ns in c(20, 99, 999)) {
    r <- cusumsq_test(fv, nsim = ns)
    expect_identical(r$reject_null, r$p.value < 0.05)
    expect_identical(r$reject_null, r$statistic > r$critical)
  }
})
