inf_fit <- function(order = c(1, 2), n = 150, ...) {
  set.seed(5)
  x <- cumsum(rnorm(n)); z <- rnorm(n)
  y <- stats::filter(0.4 * x + 0.2 * z + rnorm(n), 0.6, method = "recursive")
  d <- stats::ts(cbind(y = as.numeric(y), x = x, z = z))
  adl(y ~ x + z, data = d, order = order, ...)
}

test_that("coefficient table matches summary.lm and sandwich", {
  fit <- inf_fit()
  tab <- adl_coefficients(fit, dist = "t")
  s <- summary(fit$lm)$coefficients
  expect_equal(tab$estimate, unname(s[, 1]))
  expect_equal(tab$std.error, unname(s[, 2]))
  expect_equal(tab$p.value, unname(s[, 4]))
  expect_equal(tab$variable[tab$term == "x_lag2"], "x")
  expect_equal(tab$lag[tab$term == "x_lag2"], 2L)
  ci <- stats::confint(fit$lm)
  expect_equal(tab$conf.low, unname(ci[, 1]))

  n <- adl_coefficients(fit)
  expect_equal(attr(n, "distribution"), "normal")
  expect_equal(n$p.value, 2 * pnorm(-abs(n$statistic)))
  h <- adl_coefficients(fit, "HAC")
  V <- adl_vcov(fit, "HAC")
  expect_equal(h$std.error, unname(sqrt(diag(V))))
  expect_equal(attr(h, "distribution"), "normal")
  expect_equal(h$p.value, 2 * pnorm(-abs(h$estimate / h$std.error)))
})

test_that("joint Wald test equals the restricted-vs-unrestricted F test", {
  fit <- inf_fit()
  t1 <- adl_joint_test(fit, terms = c("x_lag1", "x_lag2"), dist = "F")
  restricted <- stats::lm(y ~ y_lag1 + x + z + z_lag1 + z_lag2, data = fit$lm$model)
  a <- stats::anova(restricted, fit$lm)
  expect_equal(t1$statistic, a$F[2])
  expect_equal(t1$p.value, a$`Pr(>F)`[2])
  expect_equal(t1$df_num, 2)
  expect_equal(t1$df_denom, stats::df.residual(fit$lm))

  t2 <- adl_joint_test(fit, variables = "x")
  expect_equal(t2$n_restrictions, 3)
  expect_equal(t2$distribution, "chi2")
  t3 <- adl_joint_test(fit, variables = "x", cov_type = "HAC")
  expect_equal(t3$distribution, "chi2")
  expect_equal(t3$df, 3)
  expect_equal(adl_joint_test(fit, variables = "x", dist = "F")$statistic,
               t2$statistic / 3)
})

test_that("joint test with R and r", {
  fit <- inf_fit()
  b <- coef(fit)
  R <- matrix(as.numeric(names(b) %in% c("x", "x_lag1", "x_lag2")), nrow = 1)
  t1 <- adl_joint_test(fit, R = R, r = 0.4)
  V <- vcov(fit)
  d <- sum(b[c("x", "x_lag1", "x_lag2")]) - 0.4
  expect_equal(t1$statistic, d^2 / as.numeric(R %*% V %*% t(R)))
  expect_match(t1$null_hypothesis, "x \\+ x_lag1 \\+ x_lag2 = 0.4")
  expect_output(print(t1), "Wald test")
})

test_that("joint test errors", {
  fit <- inf_fit()
  expect_error(adl_joint_test(fit, terms = "w"), "Not coefficient")
  expect_error(adl_joint_test(fit, variables = "w"), "Not variable")
  expect_error(adl_joint_test(fit), "Give")
  # x alone is x_t; variables = "x" is every lag of x
  expect_equal(adl_joint_test(fit, terms = "x")$n_restrictions, 1)
  expect_error(adl_joint_test(fit, R = matrix(1, 1, 2)), "columns")
  expect_error(adl_joint_test(fit, R = rbind(1:8, 2 * (1:8))), "dependent")
  expect_error(adl_joint_test(fit, variables = "x", alpha = 2), "alpha")
})

test_that("long-run multipliers and delta-method standard errors", {
  fit <- inf_fit()
  lr <- adl_long_run(fit)
  b <- coef(fit)
  C1 <- 1 - b[["y_lag1"]]
  kx <- sum(b[c("x", "x_lag1", "x_lag2")]) / C1
  expect_equal(lr$estimate[lr$term == "x"], kx)
  expect_equal(lr$estimate[lr$term == "(Intercept)"], b[["(Intercept)"]] / C1)
  expect_equal(attr(lr, "C1"), C1)
  expect_true(attr(lr, "stable"))

  # Delta method by numerical gradient
  V <- vcov(fit)
  f <- function(bb) sum(bb[c("x", "x_lag1", "x_lag2")]) / (1 - bb[["y_lag1"]])
  g <- vapply(seq_along(b), function(j) {
    h <- 1e-6; e <- b; e[j] <- e[j] + h; (f(e) - f(b)) / h
  }, numeric(1))
  expect_equal(lr$std.error[lr$term == "x"], sqrt(sum(g * (V %*% g))),
               tolerance = 1e-5)
  expect_equal(nrow(adl_long_run(fit, "HAC")), 3)
})

test_that("unstable models: NA long run, with a warning", {
  fit <- inf_fit()
  fit$lm$coefficients[["y_lag1"]] <- 1.05  # root of C(z) inside the unit circle
  expect_false(is_stable(fit))
  expect_warning(lr <- adl_long_run(fit), "do not exist")
  expect_true(all(is.na(lr$estimate)) && all(is.na(lr$std.error)))
  fit$lm$coefficients[["y_lag1"]] <- 1     # unit root: C(1) = 0
  expect_warning(lr1 <- adl_long_run(fit), "do not exist")
  expect_true(all(is.na(lr1$estimate)))
})

test_that("parameter labels follow the chapters' notation", {
  fit <- inf_fit()
  tab <- adl_coefficients(fit)
  expect_equal(tab$param, c("alpha_0", "alpha_1", "beta_1,0", "beta_1,1", "beta_1,2",
                            "beta_2,0", "beta_2,1", "beta_2,2"))
  set.seed(1)
  d <- stats::ts(cbind(y = rnorm(60), x = rnorm(60), dm = rep(0:1, 30)))
  one <- adl(y ~ x | dm, data = d, order = c(2, 1), trend = TRUE)
  expect_equal(adl_coefficients(one)$param,
               c("alpha_0", "alpha_1", "alpha_2", "beta_0", "beta_1", "", ""))
})

test_that("summary uses z by default and t on request", {
  fit <- inf_fit()
  s <- summary(fit)
  expect_equal(s$dist, "normal")
  expect_equal(s$coefficients$p.value, 2 * pnorm(-abs(s$coefficients$statistic)))
  expect_equal(s$overall$distribution, "chi2")
  out <- capture.output(print(s))
  expect_true(any(grepl("z value", out)))
  expect_true(any(grepl("y_lag1 \\[alpha_1\\]", out)))
  st <- summary(fit, dist = "t")
  expect_equal(st$coefficients$p.value, unname(summary(fit$lm)$coefficients[, 4]))
  expect_equal(st$overall$statistic, unname(summary(fit$lm)$fstatistic[1]))
  sh <- summary(fit, cov_type = "HAC")
  expect_output(print(sh), "HAC \\(Bartlett")
})

test_that("named R columns are matched by name (audit)", {
  set.seed(7)
  d <- stats::ts(data.frame(y = rnorm(100), x = rnorm(100)))
  f <- adl(y ~ x, d, order = c(0, 0))
  t1 <- adl_joint_test(f, R = matrix(c(1, 0), 1, dimnames = list(NULL, c("x", "(Intercept)"))))
  expect_equal(t1$statistic, adl_joint_test(f, terms = "x")$statistic)
  expect_match(t1$null_hypothesis, "^x = 0")
  expect_error(adl_joint_test(f, R = matrix(c(1, 0), 1, dimnames = list(NULL, c("x", "w")))),
               "coefficient names")
})
