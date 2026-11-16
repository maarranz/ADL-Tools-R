# Exogeneity tests (ADL3) ------------------------------------------------------

#' Wu-Hausman test of predeterminedness (control-function form)
#'
#' Tests whether a contemporaneous regressor \eqn{x_t} is correlated with
#' the error of the ADL equation (simultaneity), following ADL3:
#'
#' 1. regress \eqn{x_t} on the instruments and on the other regressors of
#'    the ADL, and save the residuals \eqn{\hat v_t};
#' 2. add \eqn{\hat v_t} to the ADL equation;
#' 3. test \eqn{H_0: \xi = 0} on its coefficient.
#'
#' Instruments are further lags of the variables in the model, given as a
#' named list of lags, e.g. `list(y = 2, x = 2)` for \eqn{y_{t-2}} and
#' \eqn{x_{t-2}}. They are valid only if they do not belong in the ADL
#' equation itself, so run this test after checking the lag order (e.g. with
#' [bg_test()]). The model and the auxiliary regressions are re-estimated on
#' the sample where the instruments exist.
#'
#' @param fit An `"adl"` object.
#' @param x The regressor tested; it must enter at lag 0. By default the
#'   only such regressor.
#' @param instruments Named list of lags of model variables used as
#'   instruments; they must not already be in the model.
#' @inheritParams adl_vcov
#' @param dist `"normal"` (default) or `"t"`.
#' @param alpha Significance level for the decision.
#' @return An object of class `"adltools_test"`, with the estimate `xi`,
#'   and `first_stage` (a Wald test that the instruments are jointly zero in
#'   the first-stage regression, and its F form). A warning is given when
#'   the first-stage F is below 10 (weak instruments): the test is then not
#'   reliable evidence either way.
#' @references Wu, D.-M. (1973). Alternative tests of independence between
#'   stochastic regressors and disturbances. *Econometrica*, 41(4), 733-750.
#'
#'   Hausman, J. A. (1978). Specification tests in econometrics.
#'   *Econometrica*, 46(6), 1251-1271.
#' @examples
#' # Simultaneity: x_t depends on the error of the y equation
#' set.seed(1)
#' n <- 400
#' z <- as.numeric(arima.sim(list(ar = 0.7), n))
#' e <- rnorm(n)
#' x <- z + 0.8 * e
#' y <- numeric(n)
#' for (t in 2:n) y[t] <- 0.5 * y[t - 1] + 0.5 * x[t] + e[t]
#' fit <- adl(y ~ x, data = ts(cbind(y = y, x = x)), order = c(1, 0))
#' wu_hausman_test(fit, instruments = list(x = 1:2))
#' @export
wu_hausman_test <- function(fit, x = NULL, instruments,
                            cov_type = c("nonrobust", "HAC"), maxlags = NULL,
                            dist = c("normal", "t"), alpha = 0.05) {
  check_adl(fit)
  cov_type <- match.arg(cov_type)
  dist <- match.arg(dist)
  check_alpha(alpha)
  x <- contemporaneous_regressor(fit, x)
  if (!is.list(instruments) || is.null(names(instruments)) ||
      !all(names(instruments) %in% names(fit$lags))) {
    stop("`instruments` must be a named list of lags of model variables, ",
         "e.g. list(y = 2, x = 2).", call. = FALSE)
  }
  for (v in names(instruments)) {
    check_whole(instruments[[v]], paste0("instruments$", v))
    clash <- intersect(instruments[[v]], fit$lags[[v]])
    if (length(clash)) {
      stop("Lag(s) ", paste(clash, collapse = ", "), " of '", v,
           "' are already in the model and cannot be instruments.", call. = FALSE)
    }
    if (v == fit$y && any(instruments[[v]] < 1)) {
      stop("Instruments from the dependent variable need lags >= 1.", call. = FALSE)
    }
  }
  all_lags <- fit$lags
  for (v in names(instruments)) all_lags[[v]] <- sort(unique(c(all_lags[[v]], instruments[[v]])))
  m <- max(unlist(all_lags))
  des <- build_design(fit$data, all_lags, unlagged = fit$unlagged, trend = fit$trend,
                      max_lag = m)
  model_cols <- fit$terms$column
  inst_cols <- setdiff(des$terms$column, model_cols)
  yv <- des$frame[[des$y_col]]
  Xm <- as.matrix(des$frame[model_cols])
  if (fit$intercept) Xm <- cbind(`(Intercept)` = 1, Xm)
  Z <- as.matrix(des$frame[inst_cols])
  # first stage: x_t on instruments and the other regressors
  W <- cbind(Xm[, colnames(Xm) != x, drop = FALSE], Z)
  fs <- stats::lm(des$frame[[x]] ~ 0 + W)
  names(fs$coefficients) <- colnames(W)
  v_hat <- stats::residuals(fs)
  n <- length(yv)
  lags <- hac_lags(maxlags, cov_type, n)
  vc <- function(m) if (cov_type == "HAC") {
    sandwich::NeweyWest(m, lag = lags, prewhite = FALSE, adjust = FALSE)
  } else stats::vcov(m)
  Vfs <- vc(fs)
  q <- ncol(Z)
  Rfs <- cbind(matrix(0, q, ncol(W) - q), diag(q))
  fs_w <- wald_stat(stats::coef(fs), Vfs, Rfs, rep(0, q))
  first_stage <- list(statistic = fs_w, df = q, F = fs_w / q,
                      p.value = stats::pchisq(fs_w, q, lower.tail = FALSE))
  if (fs_w / q < 10) {
    warning("The instruments look weak (first-stage F = ", format(signif(fs_w / q, 3)),
            " < 10): the test is not reliable evidence either way.", call. = FALSE)
  }
  # augmented ADL
  aug <- stats::lm(yv ~ 0 + cbind(Xm, v_hat))
  names(aug$coefficients) <- c(colnames(Xm), "v_hat")
  b <- stats::coef(aug)[["v_hat"]]
  se <- sqrt(vc(aug)["v_hat", "v_hat"])
  inf <- single_coef_inference(b, se, n - length(stats::coef(aug)), dist)
  new_adl_test(c(list(
    method = "wu_hausman",
    null_hypothesis = paste0(x, "_t is predetermined (xi = 0)"),
    alpha = alpha, reject_null = inf$p.value < alpha, nobs = n,
    cov_type = cov_type, maxlags = lags, n_restrictions = 1L,
    xi = b, first_stage = first_stage,
    instruments = inst_cols
  ), inf), "Wu-Hausman test of predeterminedness")
}

#' Weak and strong exogeneity test (regression for the change in x)
#'
#' Following ADL3 (and Urbain, 1992), computes the lagged equilibrium error
#' \eqn{\widehat{ec}_{t-1} = y_{t-1} - \sum_i \hat k_i x_{i,t-1}} from the
#' long-run multipliers of the ADL, and estimates the marginal equation
#' \deqn{\Delta x_t = \lambda_0 + \sum_{j=1}^{s} \lambda_{1j} \Delta x_{t-j}
#' + \sum_{j=1}^{s} \lambda_{2j} \Delta y_{t-j} + \kappa\, \widehat{ec}_{t-1}
#' + \nu_t.}
#'
#' * Weak exogeneity: \eqn{H_0: \kappa = 0}. A significant \eqn{\hat\kappa}
#'   means that \eqn{x} responds to the equilibrium error, so the long-run
#'   parameters also appear in its marginal model.
#' * Strong exogeneity (`strong = TRUE`): weak exogeneity and no Granger
#'   causality from \eqn{y} to \eqn{x}, tested jointly as \eqn{\kappa = 0}
#'   and \eqn{\lambda_{2j} = 0} for all \eqn{j}.
#'
#' By default (`x_level = TRUE`) the regression also includes the level
#' \eqn{x_{t-1}}. Without it, when \eqn{x} is stationary its own mean
#' reversion (\eqn{\Delta x_t} depends on \eqn{x_{t-1}}) is picked up by
#' \eqn{\widehat{ec}_{t-1}}, which contains \eqn{x_{t-1}}, and the test
#' rejects far too often (about 60% at a nominal 5% in simulations with an
#' exogenous AR(1) regressor). With the level term the size is correct both
#' for stationary variables and for cointegrated I(1) ones, where
#' \eqn{\kappa} multiplies a stationary combination and keeps its normal
#' limit. `x_level = FALSE` gives the regression exactly as written in ADL3,
#' appropriate only for cointegrated I(1) variables.
#'
#' Choose `lags` so that \eqn{\nu_t} is free of serial correlation (the
#' marginal regression is returned in `$marginal` for checking).
#'
#' @param fit An `"adl"` object.
#' @param x The regressor tested; by default the only regressor.
#' @param lags Number of lags \eqn{s} of \eqn{\Delta x} and \eqn{\Delta y}.
#' @param strong Test strong exogeneity instead of weak?
#' @param x_level Include \eqn{x_{t-1}} in the marginal regression (see
#'   Details).
#' @inheritParams adl_vcov
#' @param dist For the weak test, `"normal"` (default) or `"t"`; for the
#'   strong test, `"chisq"` (default) or `"F"`.
#' @param alpha Significance level for the decision.
#' @return An object of class `"adltools_test"`, with `kappa` and the
#'   marginal regression (`marginal`, an `lm` object).
#' @references Urbain, J.-P. (1992). On weak exogeneity in error correction
#'   models. *Oxford Bulletin of Economics and Statistics*, 54(2), 187-207.
#'
#'   Engle, R. F., Hendry, D. F. and Richard, J.-F. (1983). Exogeneity.
#'   *Econometrica*, 51(2), 277-304.
#' @examples
#' fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
#' weak_exogeneity_test(fit)
#' weak_exogeneity_test(fit, strong = TRUE)
#' @export
weak_exogeneity_test <- function(fit, x = NULL, lags = 1, strong = FALSE,
                                 x_level = TRUE,
                                 cov_type = c("nonrobust", "HAC"),
                                 maxlags = NULL, dist = NULL, alpha = 0.05) {
  check_adl(fit)
  cov_type <- match.arg(cov_type)
  check_alpha(alpha)
  check_count(lags, "lags")
  if (is.null(x)) {
    if (length(fit$x) != 1) stop("Several regressors: choose one with `x`.", call. = FALSE)
    x <- fit$x
  }
  if (!x %in% fit$x) stop("'", x, "' is not a regressor of the model.", call. = FALSE)
  dist <- dist %||% if (strong) "chisq" else "normal"
  dist <- match.arg(dist, if (strong) c("chisq", "F") else c("normal", "t"))

  lr <- suppressWarnings(adl_long_run(fit))
  k <- stats::setNames(lr$estimate, lr$term)
  z <- zoo::coredata(fit$data)
  N <- nrow(z)
  ec <- z[, fit$y] - as.numeric(z[, fit$x, drop = FALSE] %*% k[fit$x])
  dx <- c(NA, diff(z[, x]))
  dy <- c(NA, diff(z[, fit$y]))
  lagv <- function(v, j) c(rep(NA, j), v[seq_len(N - j)])
  cols <- list(dx = dx)
  for (j in seq_len(lags)) cols[[paste0("dx_lag", j)]] <- lagv(dx, j)
  for (j in seq_len(lags)) cols[[paste0("dy_lag", j)]] <- lagv(dy, j)
  cols$ec_lag1 <- lagv(ec, 1)
  if (x_level) cols$x_lag1 <- lagv(z[, x], 1)
  df <- as.data.frame(cols)
  # same observations as the ADL fit where possible
  keep <- stats::complete.cases(df) & zoo::index(fit$data) %in% fit$index
  df <- df[keep, , drop = FALSE]
  marg <- stats::lm(dx ~ ., data = df)
  n <- nrow(df)
  lg <- hac_lags(maxlags, cov_type, n)
  V <- if (cov_type == "HAC") {
    sandwich::NeweyWest(marg, lag = lg, prewhite = FALSE, adjust = FALSE)
  } else stats::vcov(marg)
  b <- stats::coef(marg)
  dfr <- stats::df.residual(marg)
  if (!strong) {
    inf <- single_coef_inference(b[["ec_lag1"]], sqrt(V["ec_lag1", "ec_lag1"]), dfr, dist)
    h0 <- paste0(x, " is weakly exogenous (kappa = 0)")
    q <- 1L
    title <- "Weak exogeneity test (regression for the change in x)"
  } else {
    terms <- c("ec_lag1", paste0("dy_lag", seq_len(lags)))
    q <- length(terms)
    R <- matrix(0, q, length(b))
    R[cbind(seq_len(q), match(terms, names(b)))] <- 1
    w <- wald_stat(b, V, R, rep(0, q))
    inf <- wald_inference(w, q, dfr, dist)
    h0 <- paste0(x, " is strongly exogenous (kappa = 0 and no lagged change in ",
                 fit$y, ")")
    title <- "Strong exogeneity test (regression for the change in x)"
  }
  new_adl_test(c(list(
    method = if (strong) "strong_exogeneity" else "weak_exogeneity",
    null_hypothesis = h0,
    alpha = alpha, reject_null = inf$p.value < alpha, nobs = n,
    cov_type = cov_type, maxlags = lg, n_restrictions = q, df_resid = dfr,
    kappa = b[["ec_lag1"]], marginal = marg
  ), inf), title)
}

#' The regressor entering at lag 0 (default: the only one)
#' @noRd
contemporaneous_regressor <- function(fit, x) {
  at0 <- fit$x[vapply(fit$x, function(v) 0L %in% fit$lags[[v]], logical(1))]
  if (is.null(x)) {
    if (length(at0) != 1) {
      stop("Choose the regressor with `x` (", length(at0),
           " regressors enter at lag 0).", call. = FALSE)
    }
    return(at0)
  }
  if (!x %in% at0) stop("'", x, "' does not enter the model at lag 0.", call. = FALSE)
  x
}

#' z or t inference for one coefficient
#' @noRd
single_coef_inference <- function(b, se, df, dist) {
  stat <- b / se
  if (dist == "t") {
    list(statistic = stat, distribution = "t", df = df,
         p.value = 2 * stats::pt(abs(stat), df, lower.tail = FALSE))
  } else {
    list(statistic = stat, distribution = "normal", df = NA_integer_,
         p.value = 2 * stats::pnorm(abs(stat), lower.tail = FALSE))
  }
}
