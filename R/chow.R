# Chow breakpoint and forecast tests (ADL3, ADL4) ---------------------------

#' Chow breakpoint test
#'
#' Tests that the coefficients are the same before and after a break date.
#' The model is re-estimated with every regressor interacted with a dummy
#' for the second period, and the interaction coefficients are tested with a
#' Wald test. With the conventional covariance and `dist = "F"` this is
#' exactly the classical Chow statistic
#' \deqn{F = \frac{(RSS_p - RSS_1 - RSS_2)/m}{(RSS_1 + RSS_2)/(T - 2m)};}
#' the interaction form also allows HAC standard errors.
#'
#' @param fit An `"adl"` object.
#' @param break_date First observation of the second period, in the units of
#'   the time index (e.g. `c(2008, 3)` for quarterly `ts` data, a date for a
#'   Date index, or a row number for data without dates).
#' @inheritParams adl_vcov
#' @param dist `"chisq"` (default) or `"F"`.
#' @param alpha Significance level for the decision.
#' @return An object of class `"adltools_test"`, with `break_date`, `n1`
#'   and `n2` added.
#' @references Chow, G. C. (1960). Tests of equality between sets of
#'   coefficients in two linear regressions. *Econometrica*, 28(3), 591-605.
#' @examples
#' fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
#' chow_break_test(fit, c(2005, 1))
#' chow_break_test(fit, c(2005, 1), dist = "F")
#' @export
chow_break_test <- function(fit, break_date, cov_type = c("nonrobust", "HAC"),
                            maxlags = NULL, dist = c("chisq", "F"), alpha = 0.05) {
  check_adl(fit)
  cov_type <- match.arg(cov_type)
  dist <- match.arg(dist)
  check_alpha(alpha)
  k <- split_position(fit, break_date, "break_date")
  X <- stats::model.matrix(fit$lm)
  y <- stats::model.response(stats::model.frame(fit$lm))
  n <- length(y)
  m <- ncol(X)
  n1 <- k - 1
  n2 <- n - n1
  if (n1 <= m || n2 <= m) {
    stop("Each period needs more than ", m, " observations (has ", n1, " and ",
         n2, "). With a short second period use chow_forecast_test().",
         call. = FALSE)
  }
  D <- as.numeric(seq_len(n) >= k)
  XD <- X * D
  colnames(XD) <- paste0(colnames(X), ":break")
  aug <- stats::lm(y ~ 0 + cbind(X, XD))
  names(aug$coefficients) <- c(colnames(X), colnames(XD))
  if (anyNA(stats::coef(aug))) {
    stop("A regressor is constant within one period; the test is not defined.",
         call. = FALSE)
  }
  lags <- hac_lags(maxlags, cov_type, n)
  V <- if (cov_type == "HAC") {
    sandwich::NeweyWest(aug, lag = lags, prewhite = FALSE, adjust = FALSE)
  } else stats::vcov(aug)
  R <- cbind(matrix(0, m, m), diag(m))
  w <- wald_stat(stats::coef(aug), V, R, rep(0, m))
  inf <- wald_inference(w, m, n - 2 * m, dist)
  bd <- format(fit$index[k])
  new_adl_test(c(list(
    method = "chow_break",
    null_hypothesis = paste0("no break at ", bd, " (all ", m, " coefficients constant)"),
    alpha = alpha, reject_null = inf$p.value < alpha, nobs = n,
    cov_type = cov_type, maxlags = lags, n_restrictions = m,
    df_resid = n - 2 * m, break_date = bd, n1 = n1, n2 = n2
  ), inf), "Chow breakpoint test")
}

#' Chow forecast (predictive failure) test
#'
#' Estimates the model on the observations before `start` (\eqn{T_1}),
#' forecasts the remaining \eqn{T_2} with those estimates fixed, and tests
#' whether the forecast errors \eqn{\hat e} are too large with
#' \deqn{F = \hat e' \hat\Omega^{-1} \hat e / T_2 \sim F_{T_2,\, T_1 - m},
#' \qquad \chi^2 = \hat e' \hat\Omega^{-1} \hat e \sim \chi^2_{T_2}.}
#'
#' Two constructions, which give different statistics (ADL4):
#'
#' * `type = "static"` (the 1-step Chow forecast test): one-step forecasts
#'   with the actual lagged values. \eqn{\hat\Omega = \hat\sigma_1^2 (I +
#'   X_2 (X_1'X_1)^{-1} X_2')} includes parameter uncertainty, so the F form
#'   is the classical Chow predictive test reported by EViews.
#' * `type = "dynamic"` (the N-step test): a single dynamic forecast path
#'   from \eqn{T_1}, with forecasts replacing the lagged dependent variable
#'   and the actual values of the regressors. \eqn{\hat\Omega} is built
#'   from the \eqn{\psi}-weights of \eqn{1/C(L)}, with
#'   \eqn{\mathrm{Cov}(e_{T_1+j}, e_{T_1+l}) = \hat\sigma_1^2
#'   \sum_{i=0}^{j-1} \psi_i \psi_{i+l-j}} for \eqn{j \le l}, treating the
#'   parameters as known (Hendry, 1979).
#'
#' @param fit An `"adl"` object (estimated on the whole sample).
#' @param start First observation of the forecast period, in the units of
#'   the time index.
#' @param type `"static"` (default) or `"dynamic"`.
#' @param dist `"chisq"` (default) or `"F"`.
#' @param alpha Significance level for the decision.
#' @return An object of class `"adltools_test"`, with `type`, `start`,
#'   `n1`, `n2` and `forecasts` (a data frame with the actual values,
#'   forecasts and errors) added.
#' @references Chow, G. C. (1960). Tests of equality between sets of
#'   coefficients in two linear regressions. *Econometrica*, 28(3), 591-605.
#'
#'   Hendry, D. F. (1979). Predictive failure and econometric modelling in
#'   macroeconomics: the transactions demand for money. In P. Ormerod (ed.),
#'   *Economic Modelling*, 217-242. London: Heinemann.
#' @examples
#' fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
#' chow_forecast_test(fit, c(2016, 1))
#' chow_forecast_test(fit, c(2016, 1), type = "dynamic")
#' @export
chow_forecast_test <- function(fit, start, type = c("static", "dynamic"),
                               dist = c("chisq", "F"), alpha = 0.05) {
  check_adl(fit)
  type <- match.arg(type)
  dist <- match.arg(dist)
  check_alpha(alpha)
  k <- split_position(fit, start, "start")
  X <- stats::model.matrix(fit$lm)
  y <- stats::model.response(stats::model.frame(fit$lm))
  n <- length(y)
  m <- ncol(X)
  T1 <- k - 1
  T2 <- n - T1
  if (T1 <= m) stop("The estimation period needs more than ", m,
                    " observations (has ", T1, ").", call. = FALSE)
  X1 <- X[seq_len(T1), , drop = FALSE]
  X2 <- X[k:n, , drop = FALSE]
  y1 <- y[seq_len(T1)]
  y2 <- y[k:n]
  est <- stats::lm.fit(X1, y1)
  if (est$rank < m) stop("The model cannot be estimated on the first period alone.",
                         call. = FALSE)
  b1 <- est$coefficients
  s2 <- sum(est$residuals^2) / (T1 - m)

  if (type == "static") {
    f <- as.numeric(X2 %*% b1)
    XtX_inv <- solve(crossprod(X1))
    Omega <- s2 * (diag(T2) + X2 %*% XtX_inv %*% t(X2))
  } else {
    terms <- term_info(fit, colnames(X))
    ar_cols <- which(terms$role == "ar")
    ar_lags <- terms$lag[ar_cols]
    f <- numeric(T2)
    for (h in seq_len(T2)) {
      row <- X2[h, ]
      for (j in seq_along(ar_cols)) {
        lag <- ar_lags[j]
        if (h - lag >= 1) row[ar_cols[j]] <- f[h - lag]
      }
      f[h] <- sum(row * b1)
    }
    alpha_poly <- numeric(max(c(0, ar_lags)))
    alpha_poly[ar_lags] <- b1[ar_cols]
    psi <- numeric(T2)
    psi[1] <- 1
    if (T2 > 1) for (j in 2:T2) {
      acc <- 0
      for (i in seq_along(alpha_poly)) if (j - i >= 1) acc <- acc + alpha_poly[i] * psi[j - i]
      psi[j] <- acc
    }
    Omega <- matrix(0, T2, T2)
    for (j in seq_len(T2)) for (l in j:T2) {
      v <- s2 * sum(psi[seq_len(j)] * psi[seq_len(j) + (l - j)])
      Omega[j, l] <- Omega[l, j] <- v
    }
  }
  e <- y2 - f
  w <- as.numeric(t(e) %*% solve(Omega, e))
  inf <- wald_inference(w, T2, T1 - m, dist)
  st <- format(fit$index[k])
  new_adl_test(c(list(
    method = paste0("chow_forecast_", type),
    null_hypothesis = paste0("no predictive failure from ", st, " (", T2,
                             " ", type, " forecasts)"),
    alpha = alpha, reject_null = inf$p.value < alpha, nobs = n,
    n_restrictions = T2, df_resid = T1 - m,
    type = type, start = st, n1 = T1, n2 = T2,
    forecasts = data.frame(index = format(fit$index[k:n]), actual = y2,
                           forecast = f, error = e,
                           se = sqrt(diag(Omega)))
  ), inf), paste0("Chow forecast test (", type, ")"))
}

#' Row position of a date in the estimation sample
#' @noRd
split_position <- function(fit, date, what) {
  idx <- fit$index
  v <- as_index_value(date, idx)
  k <- which(idx >= v)[1]
  if (is.na(k) || k <= 1) {
    stop("`", what, "` must fall inside the estimation sample (",
         format(idx[1]), " - ", format(idx[length(idx)]), "), after its first observation.",
         call. = FALSE)
  }
  k
}
