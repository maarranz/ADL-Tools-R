# The ECM regression run directly: same coefficients and standard errors
ecm_lm <- function(d, p, q) {
  z <- as.data.frame(d)
  N <- nrow(z)
  L <- function(v, j) c(rep(NA, j), v[seq_len(N - j)])
  D <- function(v) c(NA, diff(v))
  cols <- list(Dy = D(z$y), y_lag1 = L(z$y, 1))
  for (v in names(q)) cols[[paste0(v, "_lag1")]] <- L(z[[v]], 1)
  for (j in seq_len(p - 1)) cols[[paste0("D_y_lag", j)]] <- L(D(z$y), j)
  for (v in names(q)) {
    cols[[paste0("D_", v)]] <- D(z[[v]])
    for (j in seq_len(q[[v]] - 1)) cols[[paste0("D_", v, "_lag", j)]] <- L(D(z[[v]]), j)
  }
  df <- as.data.frame(cols)
  df <- df[-seq_len(max(p, unlist(q))), ]
  stats::lm(Dy ~ ., data = df)
}

test_that("ECM coefficients and standard errors equal the ECM regression", {
  set.seed(21)
  n <- 250
  x1 <- cumsum(rnorm(n)); x2 <- rnorm(n)
  y <- stats::filter(0.4 * x1 + 0.3 * x2 + rnorm(n), c(0.5, 0.2), method = "recursive")
  d <- stats::ts(cbind(y = as.numeric(y), x1 = x1, x2 = x2))
  fit <- adl(y ~ x1 + x2, data = d, order = c(3, 2, 1))
  ref <- ecm_lm(d, 3, list(x1 = 2, x2 = 1))
  for (ct in c("nonrobust", "HAC")) {
    e <- adl_ecm(fit, cov_type = ct)
    tab <- e$coefficients
    b_ref <- coef(ref)
    V_ref <- if (ct == "HAC") {
      sandwich::NeweyWest(ref, lag = attr(adl_vcov(fit, "HAC"), "maxlags"),
                          prewhite = FALSE, adjust = FALSE)
    } else vcov(ref)
    m <- match(names(b_ref), tab$term)
    expect_false(anyNA(m))
    expect_equal(tab$estimate[m], unname(b_ref), tolerance = 1e-8)
    expect_equal(tab$std.error[m], unname(sqrt(diag(V_ref))), tolerance = 1e-8)
  }
  expect_equal(nrow(e$coefficients), length(coef(ref)))
})

test_that("ADL(1,1): the formulas of ADL2", {
  fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
  th <- coef(fit)
  e <- adl_ecm(fit)
  tab <- e$coefficients
  expect_equal(tab$term, c("(Intercept)", "y_lag1", "x_lag1", "D_x"))
  expect_equal(tab$param, c("alpha_0", "b", "beta*", "delta_0"))
  expect_equal(tab$estimate, unname(c(th[1], th["y_lag1"] - 1, th["x"] + th["x_lag1"],
                                      th["x"])))
  # b has the standard error of alpha_1
  expect_equal(tab$std.error[2], sqrt(vcov(fit)["y_lag1", "y_lag1"]))
  lr <- adl_long_run(fit)
  eq <- e$equilibrium
  expect_equal(eq$term, c("b", "k_x"))
  expect_equal(eq$estimate[2], lr$estimate[lr$term == "x"])
  expect_equal(eq$std.error[2], lr$std.error[lr$term == "x"])
  expect_equal(eq$estimate[1], tab$estimate[2])
  # equilibrium error
  z <- zoo::coredata(fit$data)
  expect_equal(as.numeric(e$ec), z[, "y"] - eq$estimate[2] * z[, "x"])
  expect_s3_class(e$ec, "zoo")
})

test_that("missing lags, dead start, no lags of y, controls and trend", {
  set.seed(22)
  n <- 200
  d <- stats::ts(cbind(y = rnorm(n), x = rnorm(n), dum = rbinom(n, 1, 0.1)))
  # dead start: no D_x term; x lags {1, 3}: delta_1 = delta_2 = -beta_3
  fit <- adl(y ~ x | dum, data = d, lags = list(y = 1:2, x = c(1, 3)), trend = TRUE)
  th <- coef(fit)
  e <- adl_ecm(fit)
  tab <- e$coefficients
  expect_false("D_x" %in% tab$term)
  expect_equal(tab$estimate[tab$term == "D_x_lag1"], -unname(th["x_lag3"]))
  expect_equal(tab$estimate[tab$term == "D_x_lag2"], -unname(th["x_lag3"]))
  expect_equal(tab$estimate[tab$term == "x_lag1"], unname(th["x_lag1"] + th["x_lag3"]))
  expect_equal(tab$estimate[tab$term == "D_y_lag1"], -unname(th["y_lag2"]))
  expect_equal(tab$estimate[tab$term == "dum"], unname(th["dum"]))
  expect_equal(tab$estimate[tab$term == "trend"], unname(th["trend"]))
  # no lags of y: b = -1 is fixed
  f0 <- adl(y ~ x, data = d, order = c(0, 1))
  e0 <- adl_ecm(f0)
  b <- e0$coefficients[e0$coefficients$term == "y_lag1", ]
  expect_equal(b$estimate, -1)
  expect_equal(b$std.error, 0)
  expect_true(is.na(b$statistic))
})

test_that("unstable models, printing and tidiers", {
  fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
  fit$lm$coefficients[["y_lag1"]] <- 1.05
  expect_warning(e <- adl_ecm(fit), "not dynamically stable")
  expect_null(e$ec)
  expect_true(is.na(e$equilibrium$estimate[2]))
  expect_equal(e$coefficients$estimate[2], 0.05)

  fit <- adl(y ~ x, data = adl_sim, order = c(2, 2))
  e <- adl_ecm(fit, dist = "t")
  out <- capture.output(print(e))
  expect_true(any(grepl("Equilibrium-error form", out)))
  expect_true(any(grepl("D_y_lag1", out)))
  expect_equal(nrow(tidy(e)), nrow(e$coefficients) + nrow(e$equilibrium))
  expect_equal(glance(e)$b, e$equilibrium$estimate[1])
  expect_error(adl_ecm(fit, maxlags = 2), "maxlags")
  expect_error(adl_ecm(lm(1 ~ 1)), "adl()")
})
