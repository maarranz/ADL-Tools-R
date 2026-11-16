ic_data <- function(n = 200) {
  set.seed(12)
  x1 <- cumsum(rnorm(n)); x2 <- rnorm(n)
  y <- stats::filter(0.5 * x1 + 0.3 * c(0, x2[-n]) + rnorm(n), 0.5, method = "recursive")
  stats::ts(cbind(y = as.numeric(y), x1 = x1, x2 = x2))
}

test_that("criteria match AIC/BIC of adl() fits on the common sample", {
  d <- ic_data()
  tab <- adl_ic(y ~ x1, data = d, max_order = c(2, 2))
  # p 0:2, q 0:2 with x_t, q 1:2 without: 3 * (3 + 2) = 15 models
  expect_equal(nrow(tab), 15)
  expect_equal(attr(tab, "nobs"), nrow(d) - 2)
  row <- tab[tab$p == 1 & tab$q_x1 == 2 & tab$x_t_x1, ]
  fit <- adl(y ~ x1, data = d, order = c(1, 2), max_lag = 2)
  n <- nobs(fit)
  ll <- as.numeric(logLik(fit))
  # EViews form: k counts the coefficients only (R's AIC() also counts sigma^2)
  expect_equal(row$k, 5)
  expect_equal(row$aic, (-2 * ll + 2 * 5) / n)
  expect_equal(row$bic, (-2 * ll + log(n) * 5) / n)
  expect_equal(row$bic, (BIC(fit) - log(n)) / n)
  row0 <- tab[tab$p == 1 & tab$q_x1 == 2 & !tab$x_t_x1, ]
  fit0 <- adl(y ~ x1, data = d, order = c(1, 2), first_lag = 1, max_lag = 2)
  expect_equal(row0$bic, (BIC(fit0) - log(n)) / n)
  expect_equal(sum(tab$best_bic), 1)
})

test_that("contemporaneous rule and several regressors", {
  d <- ic_data()
  inc <- adl_ic(y ~ x1 + x2, data = d, max_order = c(1, 2), contemporaneous = "include")
  expect_true(all(inc$x_t_x1 & inc$x_t_x2))
  expect_equal(nrow(inc), 2 * 3)
  mix <- adl_ic(y ~ x1 + x2, data = d, max_order = c(1, 2),
                contemporaneous = c(x2 = "exclude"))
  expect_true(all(!mix$x_t_x2))
  expect_true(all(mix$q_x2 >= 1))
  full <- adl_ic(y ~ x1 + x2, data = d, max_order = c(1, 1, 2), same_q = FALSE,
                 contemporaneous = "include")
  expect_equal(nrow(full), 2 * 2 * 3)
  expect_error(adl_ic(y ~ x1 + x2, data = d, max_order = c(1, 1, 2)), "same q_max")
})

test_that("printing and errors", {
  d <- ic_data()
  tab <- adl_ic(y ~ x1, data = d, max_order = c(2, 2))
  out <- capture.output(print(tab, n = 3))
  expect_true(any(grepl("12 more models", out)))
  expect_true(any(grepl("\\*", out)))
  expect_error(adl_ic(y ~ L(y, 1) + x1, data = d, max_order = c(1, 1)), "without L")
  expect_error(adl_ic(y ~ x1, data = d, max_order = c(1, 1, 1)), "length 2")
  expect_error(adl_ic(y ~ x1, data = d, max_order = c(1, 1), contemporaneous = "no"),
               "contemporaneous")
})

test_that("collinear candidates are dropped and an empty grid is diagnosed (audit)", {
  set.seed(7)
  d <- data.frame(x = rnorm(100), y = rnorm(100))
  d$copy <- d$x
  expect_error(suppressMessages(adl_ic(y ~ x + copy, stats::ts(d), c(1, 1))),
               "None of the candidate models")
  d$z <- d$x + c(0, d$x[-100])  # collinear only when both x and x_lag1 are in
  expect_message(tab <- adl_ic(y ~ x + z, stats::ts(d), c(1, 1), contemporaneous = "include"),
                 "dropped")
  expect_false(anyNA(tab$aic))
  expect_error(adl_ic(y ~ x, stats::ts(d), c(1, 0), contemporaneous = "exclude"),
               "No model fits the grid")
})
