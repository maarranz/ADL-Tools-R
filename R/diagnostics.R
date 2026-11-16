# Residual diagnostics (ADL3) -------------------------------------------------
#
# Breusch-Godfrey and ARCH-LM test *for* serial correlation and conditional
# heteroskedasticity under the null of well-behaved errors, so they always use
# OLS (no HAC option): wrapping them in HAC would presume the answer. RESET
# is an ordinary Wald test on an augmented regression and takes `cov_type`.

#' Breusch-Godfrey LM test for serial correlation
#'
#' Regresses the residuals on all the regressors of the model and on
#' `order` lags of the residuals (pre-sample residuals set to zero) and tests
#' that the coefficients of the lagged residuals are zero. Valid with lagged
#' dependent variables, unlike Durbin-Watson.
#'
#' @param fit An `"adl"` object.
#' @param order Number of lagged residuals \eqn{s}.
#' @param dist `"chisq"` (default): the LM statistic \eqn{nR^2 \sim
#'   \chi^2_s}; `"F"`: the F form of the same auxiliary regression.
#' @param alpha Significance level for the decision.
#' @return An object of class `"adltools_test"`.
#' @references Breusch, T. S. (1978). Testing for autocorrelation in dynamic
#'   linear models. *Australian Economic Papers*, 17(31), 334-355.
#'
#'   Godfrey, L. G. (1978). Testing against general autoregressive and
#'   moving average error models when the regressors include lagged
#'   dependent variables. *Econometrica*, 46(6), 1293-1301.
#' @examples
#' fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
#' bg_test(fit, order = 4)
#' @export
bg_test <- function(fit, order = 2, dist = c("chisq", "F"), alpha = 0.05) {
  check_adl(fit)
  dist <- match.arg(dist)
  check_alpha(alpha)
  check_count(order, "order")
  u <- stats::residuals(fit$lm)
  n <- length(u)
  if (order >= n - length(stats::coef(fit))) {
    stop("`order` is too large for the sample.", call. = FALSE)
  }
  U <- vapply(seq_len(order), function(j) c(rep(0, j), u[seq_len(n - j)]), numeric(n))
  colnames(U) <- paste0("u_lag", seq_len(order))
  X <- stats::model.matrix(fit$lm)
  aux <- stats::lm.fit(cbind(X, U), u)
  if (aux$rank < ncol(X) + order) {
    stop("The auxiliary regression is rank deficient; reduce `order`.", call. = FALSE)
  }
  df_res <- n - ncol(X) - order
  # LM = n R^2 with R^2 = 1 - RSS_u / sum(u^2): the share of the residuals
  # explained by the auxiliary regression. Uncentred, so it is also valid
  # (and non-negative) without an intercept; equal to the centred form when
  # the model has one, since the residuals then have mean zero.
  rss_u <- sum(aux$residuals^2)
  lm_stat <- n * (1 - rss_u / sum(u^2))
  # F form: (RSS_r - RSS_u)/s / (RSS_u / df), with RSS_r = sum(u^2)
  f_stat <- ((sum(u^2) - rss_u) / order) / (rss_u / df_res)
  inf <- if (dist == "chisq") {
    list(statistic = lm_stat, distribution = "chi2", df = order,
         p.value = stats::pchisq(lm_stat, order, lower.tail = FALSE))
  } else {
    list(statistic = f_stat, distribution = "F", df_num = order, df_denom = df_res,
         p.value = stats::pf(f_stat, order, df_res, lower.tail = FALSE))
  }
  new_adl_test(c(list(
    method = "breusch_godfrey",
    null_hypothesis = paste0("no serial correlation up to order ", order),
    alpha = alpha, reject_null = inf$p.value < alpha, nobs = n,
    n_restrictions = order, df_resid = df_res
  ), inf), "Breusch-Godfrey LM test for serial correlation")
}

#' Engle's ARCH-LM test
#'
#' Regresses the squared residuals on a constant and `lags` of their own
#' lags and tests that the lag coefficients are zero with
#' \eqn{nR^2 \sim \chi^2_q}, where \eqn{n} is the number of observations in
#' the auxiliary regression.
#'
#' @param fit An `"adl"` object.
#' @param lags Number of lagged squared residuals \eqn{q}.
#' @param alpha Significance level for the decision.
#' @return An object of class `"adltools_test"`.
#' @references Engle, R. F. (1982). Autoregressive conditional
#'   heteroscedasticity with estimates of the variance of United Kingdom
#'   inflation. *Econometrica*, 50(4), 987-1007.
#' @examples
#' fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
#' arch_test(fit, lags = 4)
#' @export
arch_test <- function(fit, lags = 1, alpha = 0.05) {
  check_adl(fit)
  check_alpha(alpha)
  check_count(lags, "lags")
  e2 <- stats::residuals(fit$lm)^2
  n0 <- length(e2)
  if (lags >= n0 - 2) stop("`lags` is too large for the sample.", call. = FALSE)
  rows <- (lags + 1):n0
  Z <- cbind(1, vapply(seq_len(lags), function(j) e2[rows - j], numeric(length(rows))))
  y <- e2[rows]
  aux <- stats::lm.fit(Z, y)
  n <- length(y)
  r2 <- 1 - sum(aux$residuals^2) / sum((y - mean(y))^2)
  stat <- n * r2
  p <- stats::pchisq(stat, lags, lower.tail = FALSE)
  new_adl_test(list(
    method = "arch_lm",
    null_hypothesis = paste0("no ARCH effects up to order ", lags),
    statistic = stat, distribution = "chi2", df = lags, p.value = p,
    alpha = alpha, reject_null = p < alpha, nobs = n, n_restrictions = lags
  ), "Engle's ARCH-LM test")
}

#' Ramsey's RESET test
#'
#' Adds powers of the fitted values to the regression and tests that their
#' coefficients are zero with a Wald test, with conventional or HAC
#' covariance.
#'
#' @param fit An `"adl"` object.
#' @param power Powers of the fitted values added; `2:3` adds
#'   \eqn{\hat y^2} and \eqn{\hat y^3}.
#' @inheritParams adl_vcov
#' @param dist `"chisq"` (default) or `"F"`.
#' @param alpha Significance level for the decision.
#' @return An object of class `"adltools_test"`.
#' @references Ramsey, J. B. (1969). Tests for specification errors in
#'   classical linear least-squares regression analysis. *Journal of the
#'   Royal Statistical Society, Series B*, 31(2), 350-371.
#' @examples
#' fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
#' reset_test(fit)
#' reset_test(fit, power = 2, cov_type = "HAC")
#' @export
reset_test <- function(fit, power = 2:3, cov_type = c("nonrobust", "HAC"),
                       maxlags = NULL, dist = c("chisq", "F"), alpha = 0.05) {
  check_adl(fit)
  cov_type <- match.arg(cov_type)
  dist <- match.arg(dist)
  check_alpha(alpha)
  check_whole(power, "power")
  if (any(power < 2)) stop("`power` must contain integers >= 2.", call. = FALSE)
  if (anyDuplicated(power)) stop("`power` must not repeat values.", call. = FALSE)
  yhat <- stats::fitted(fit$lm)
  # scale before taking powers, for numerical stability (the test is invariant)
  z <- yhat / max(abs(yhat))
  P <- vapply(power, function(k) z^k, numeric(length(z)))
  colnames(P) <- paste0("fitted_pow", power)
  X <- stats::model.matrix(fit$lm)
  y <- stats::model.response(stats::model.frame(fit$lm))
  aug <- stats::lm(y ~ 0 + cbind(X, P))
  if (aug$rank < ncol(X) + ncol(P)) {
    stop("The powers of the fitted values are collinear with the regressors ",
         "(e.g. a model with only dummies); RESET is not defined here.", call. = FALSE)
  }
  names(aug$coefficients) <- c(colnames(X), colnames(P))
  n <- length(y)
  lags <- hac_lags(maxlags, cov_type, n)
  V <- if (cov_type == "HAC") {
    sandwich::NeweyWest(aug, lag = lags, prewhite = FALSE, adjust = FALSE)
  } else stats::vcov(aug)
  b <- stats::coef(aug)
  q <- length(power)
  R <- cbind(matrix(0, q, ncol(X)), diag(q))
  w <- wald_stat(b, V, R, rep(0, q))
  inf <- wald_inference(w, q, n - length(b), dist)
  new_adl_test(c(list(
    method = "reset",
    null_hypothesis = paste0("coefficients of fitted^", paste(power, collapse = ", fitted^"),
                             " are zero"),
    alpha = alpha, reject_null = inf$p.value < alpha, nobs = n,
    cov_type = cov_type, maxlags = lags, n_restrictions = q,
    df_resid = n - length(b)
  ), inf), "Ramsey's RESET test")
}

#' Jarque-Bera normality test
#'
#' \eqn{JB = \frac{n}{6}\left(S^2 + \frac{(K-3)^2}{4}\right) \sim \chi^2_2},
#' with skewness \eqn{S} and kurtosis \eqn{K} of the residuals (moments
#' with divisor \eqn{n}, as in EViews).
#'
#' @param fit An `"adl"` object.
#' @param alpha Significance level for the decision.
#' @return An object of class `"adltools_test"`, with `skewness` and
#'   `kurtosis` added.
#' @references Jarque, C. M. and Bera, A. K. (1980). Efficient tests for
#'   normality, homoscedasticity and serial independence of regression
#'   residuals. *Economics Letters*, 6(3), 255-259.
#' @examples
#' fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
#' jb_test(fit)
#' @export
jb_test <- function(fit, alpha = 0.05) {
  check_adl(fit)
  check_alpha(alpha)
  u <- stats::residuals(fit$lm)
  n <- length(u)
  e <- u - mean(u)
  m2 <- mean(e^2)
  S <- mean(e^3) / m2^1.5
  K <- mean(e^4) / m2^2
  stat <- n / 6 * (S^2 + (K - 3)^2 / 4)
  p <- stats::pchisq(stat, 2, lower.tail = FALSE)
  new_adl_test(list(
    method = "jarque_bera",
    null_hypothesis = "residuals are normally distributed",
    statistic = stat, distribution = "chi2", df = 2L, p.value = p,
    alpha = alpha, reject_null = p < alpha, nobs = n, n_restrictions = 2L,
    skewness = S, kurtosis = K
  ), "Jarque-Bera normality test")
}

#' Battery of residual diagnostics
#'
#' Runs [bg_test()], [arch_test()], [reset_test()] and [jb_test()] and
#' collects the results in one table. Breusch-Godfrey and ARCH-LM always
#' use OLS; RESET follows `cov_type`.
#'
#' @param fit An `"adl"` object.
#' @param bg_order,arch_lags,reset_power Passed to the individual tests.
#' @inheritParams adl_vcov
#' @param alpha Significance level for the decisions.
#' @return A data frame of class `"adl_diagnostics"`, one row per test, with
#'   the individual test objects in attribute `"tests"`.
#' @examples
#' fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
#' adl_diagnostics(fit)
#' @export
adl_diagnostics <- function(fit, bg_order = 2, arch_lags = 1, reset_power = 2:3,
                            cov_type = c("nonrobust", "HAC"), maxlags = NULL,
                            alpha = 0.05) {
  cov_type <- match.arg(cov_type)
  tests <- list(
    bg = bg_test(fit, order = bg_order, alpha = alpha),
    arch = arch_test(fit, lags = arch_lags, alpha = alpha),
    reset = reset_test(fit, power = reset_power, cov_type = cov_type,
                       maxlags = maxlags, alpha = alpha),
    jb = jb_test(fit, alpha = alpha)
  )
  tab <- data.frame(
    test = c(paste0("Breusch-Godfrey (", bg_order, ")"),
             paste0("ARCH-LM (", arch_lags, ")"),
             paste0("RESET (", paste(reset_power, collapse = ","), ")"),
             "Jarque-Bera"),
    statistic = vapply(tests, `[[`, numeric(1), "statistic"),
    df = vapply(tests, function(t) as.integer(t$df), integer(1)),
    p.value = vapply(tests, `[[`, numeric(1), "p.value"),
    decision = ifelse(vapply(tests, `[[`, logical(1), "reject_null"),
                      "reject", "do not reject"),
    row.names = NULL
  )
  structure(tab, class = c("adl_diagnostics", "data.frame"), tests = tests,
            alpha = alpha, cov_type = cov_type, label = fit$label)
}

#' Print a diagnostics table
#' @param x An `"adl_diagnostics"` object.
#' @param digits Significant digits.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.adl_diagnostics <- function(x, digits = 4, ...) {
  title <- paste0("Residual diagnostics: ", attr(x, "label"))
  cat(title, "\n", strrep("-", nchar(title)), "\n", sep = "")
  d <- as.data.frame(unclass(x))
  d$statistic <- signif(d$statistic, digits)
  d$p.value <- signif(d$p.value, digits)
  names(d)[names(d) == "df"] <- "df (chi2)"
  print(d, row.names = FALSE)
  cat("Chi-square statistics; decisions at ", format(100 * attr(x, "alpha")),
      "%. BG and ARCH with OLS; RESET with ",
      if (attr(x, "cov_type") == "HAC") "HAC" else "conventional", " covariance.\n",
      sep = "")
  invisible(x)
}

check_count <- function(v, what) {
  if (!is.numeric(v) || length(v) != 1 || is.na(v) || v < 1 || v != round(v)) {
    stop("`", what, "` must be a single integer >= 1.", call. = FALSE)
  }
  invisible(v)
}
