sim_adl22 <- function(seed = 5, n = 200) {
  set.seed(seed)
  x <- as.numeric(arima.sim(list(ar = 0.5), n))
  y <- numeric(n)
  e <- rnorm(n)
  for (t in 3:n) y[t] <- 0.6 * y[t - 1] - 0.3 * y[t - 2] + 0.4 * x[t] +
    0.3 * x[t - 1] - 0.25 * x[t - 2] + e[t]
  stats::ts(cbind(y = y, x = x), start = c(1970, 1), frequency = 4)
}

test_that("adl_sim: Gets goes from ADL(4,4) down to the true ADL(1,1)", {
  g <- adl_gets(y ~ x, data = adl_sim)
  expect_s3_class(g, "adl_gets")
  expect_equal(g$max_order, 4)                       # quarterly default
  expect_equal(g$fit$label, "ADL(1, 1)")
  expect_equal(g$path$lag, 4:1)
  expect_equal(g$path$action, c("drop", "drop", "drop", "stop"))
  expect_true(g$general_ok)
  expect_match(g$stop_reason, "jointly significant")
  expect_output(print(g), "conditional on the search")
})

test_that("each step is the Wald test of the highest lags on the common sample", {
  g <- adl_gets(y ~ x, data = adl_sim)
  for (p in 4:1) {
    f <- adl(y ~ x, data = adl_sim, order = c(p, p), max_lag = 4)
    expect_equal(nobs(f), nobs(g$general))
    w <- adl_joint_test(f, terms = c(paste0("y_lag", p), paste0("x_lag", p)))
    expect_equal(g$path$statistic[g$path$lag == p], w$statistic)
  }
  expect_equal(coef(g$fit), coef(adl(y ~ x, data = adl_sim, order = c(1, 1), max_lag = 4)))
  # diagnostics columns are those of the model after dropping
  f3 <- adl(y ~ x, data = adl_sim, order = c(3, 3), max_lag = 4)
  expect_equal(g$path$bg_p[1], bg_test(f3, order = 4)$p.value)
  # F form
  gf <- adl_gets(y ~ x, data = adl_sim, dist = "F")
  expect_equal(gf$path$statistic, g$path$statistic / 2)
})

test_that("stops at the true order, backs up when a diagnostic fails", {
  d <- sim_adl22(seed = 9)
  g <- adl_gets(y ~ x, data = d)
  expect_equal(g$fit$label, "ADL(2, 2)")
  # with Wald tests that never reject, the diagnostics stop the search
  gb <- adl_gets(y ~ x, data = d, alpha = 1e-300)
  expect_match(tail(gb$path$action, 1), "^back up \\(")
  expect_match(gb$stop_reason, "would fail")
  expect_equal(gb$fit$label, "ADL(2, 2)")
})

test_that("contemporaneous: exclude, both; smallest model", {
  ge <- adl_gets(y ~ x, data = adl_sim, contemporaneous = "exclude", alpha = 1e-300,
                 diagnostics = character(0))
  expect_false(0 %in% ge$fit$lags$x)
  expect_equal(ge$fit$lags$x, 1L)
  expect_match(ge$stop_reason, "smallest model")
  expect_false(any(grepl("_p$", names(ge$path))))
  gi <- adl_gets(y ~ x, data = adl_sim, alpha = 1e-300, diagnostics = character(0))
  expect_equal(gi$fit$label, "ADL(0, 0)")
  b <- adl_gets(y ~ x + w, data = adl_sim, contemporaneous = "both")
  expect_s3_class(b, "adl_gets_both")
  expect_equal(b$include$contemporaneous, "include")
  expect_equal(b$exclude$contemporaneous, "exclude")
  expect_output(print(b), "wu_hausman_test")
})

test_that("warnings, defaults and errors", {
  expect_warning(adl_gets(y ~ x, data = adl_sim, diag_alpha = 0.999), "mis-specified")
  a <- stats::ts(as.matrix(adl_sim), start = 1900, frequency = 1)
  expect_equal(adl_gets(y ~ x, data = a)$max_order, 2)
  expect_equal(adl_gets(y ~ x, data = adl_sim, max_order = 2)$general$label, "ADL(2, 2)")
  expect_error(adl_gets(y ~ L(x, 0:1), data = adl_sim), "without L")
  expect_error(adl_gets(y ~ x, data = adl_sim, max_order = 0), ">= 1")
  expect_error(adl_gets(y ~ x, data = adl_sim, diagnostics = "dw"))
})
