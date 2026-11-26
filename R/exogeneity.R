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
#' [bg_test()]). The test uses the fitted sample, never more: if the
#' instrument lags need observations from before it, the first fitted
#' observations are dropped with a message (refit with `max_lag` equal to the
#' largest instrument lag to keep the same sample).
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
#'   the `sample` used (the fitted sample, or the part of it where the
#'   instruments exist), and
#'   `first_stage`: the conventional F test that the instruments are
#'   jointly zero (`F`, `df`, `p.value`) and, with HAC, the robust Wald
#'   statistic and its value divided by the number of instruments
#'   (`hac_wald`, `hac_F`). A warning is given when the conventional F is
#'   below the rule-of-thumb value 10 (weak instruments): the test is then
#'   not reliable evidence either way.
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
  m_inst <- max(unlist(all_lags))
  # build on all the data, then keep exactly the fitted sample
  des <- build_design(fit$data, all_lags, unlagged = fit$unlagged, trend = fit$trend,
                      max_lag = m_inst)
  keep <- des$index %in% fit$index
  if (sum(keep) < length(fit$index)) {
    message("The instrument lags need ", m_inst, " earlier observations: the test ",
            "uses ", sum(keep), " of the ", length(fit$index), " fitted observations (",
            format(des$index[keep][1]), " - ", format(des$index[keep][sum(keep)]),
            "). Refit with `max_lag = ", m_inst, "` to test on the same sample.")
  }
  if (sum(keep) <= length(fit$terms$column) + 3) {
    stop("Too few observations left for the test.", call. = FALSE)
  }
  frame <- des$frame[keep, , drop = FALSE]
  if (fit$trend) frame$trend <- seq_len(nrow(frame))
  model_cols <- fit$terms$column
  x_col <- fit$terms$column[fit$terms$variable == x & fit$terms$lag == 0 &
                              fit$terms$role == "dl"]
  inst_cols <- setdiff(des$terms$column, model_cols)
  yv <- frame[[des$y_col]]
  Xm <- as.matrix(frame[model_cols])
  if (fit$intercept) Xm <- cbind(`(Intercept)` = 1, Xm)
  Z <- as.matrix(frame[inst_cols])
  # first stage: x_t on the other regressors and the instruments
  W0 <- Xm[, colnames(Xm) != x_col, drop = FALSE]
  W <- cbind(W0, Z)
  xv <- frame[[x_col]]
  fs <- stats::lm(xv ~ 0 + W)
  if (fs$rank < ncol(W)) {
    stop("The instruments are collinear with each other or with the regressors ",
         "of the model.", call. = FALSE)
  }
  names(fs$coefficients) <- colnames(W)
  v_hat <- stats::residuals(fs)
  n <- length(yv)
  lags <- hac_lags(maxlags, cov_type, n)
  vc <- function(m) if (cov_type == "HAC") {
    sandwich::NeweyWest(m, lag = lags, prewhite = FALSE, adjust = FALSE)
  } else stats::vcov(m)
  q <- ncol(Z)
  # conventional first-stage F (RSS form) and, with HAC, the robust Wald / q
  rss1 <- sum(v_hat^2)
  rss0 <- sum(stats::lm.fit(W0, xv)$residuals^2)
  F_conv <- ((rss0 - rss1) / q) / (rss1 / (n - ncol(W)))
  Rfs <- cbind(matrix(0, q, ncol(W0)), diag(q))
  w_rob <- if (cov_type == "HAC") wald_stat(stats::coef(fs), vc(fs), Rfs, rep(0, q))
  first_stage <- list(F = F_conv, df = c(q, n - ncol(W)),
                      p.value = stats::pf(F_conv, q, n - ncol(W), lower.tail = FALSE),
                      hac_wald = w_rob,
                      hac_F = if (!is.null(w_rob)) w_rob / q)
  if (F_conv < 10) {
    warning("The instruments look weak (first-stage F = ", format(signif(F_conv, 3)),
            ", below the rule-of-thumb 10): the test is not reliable evidence ",
            "either way.", call. = FALSE)
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
    instruments = inst_cols,
    sample = c(format(des$index[keep][1]), format(des$index[keep][sum(keep)]))
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
#' The marginal equation conditions on: lags of \eqn{\Delta x} and
#' \eqn{\Delta y}; with several regressors, lags of the changes of the
#' others (and, with `x_level`, their levels at \eqn{t-1}); the controls
#' entered after `|` (at \eqn{t-1}); and the trend, if the model has one.
#' Because \eqn{\widehat{ec}_{t-1}} is linear, including the controls and
#' trend there is equivalent to building their long-run terms into the
#' equilibrium error. With `x_level = TRUE`, \eqn{\kappa} is the coefficient
#' of \eqn{y_{t-1}} given the regressor levels, so it does not depend on the
#' estimated long-run multipliers. The strong test adds the lagged changes
#' of \eqn{y} to the null, within this same conditioning set.
#'
#' **A supplied equilibrium error** (`ec =`), for instance the residual
#' \eqn{\hat z_t} of a cointegrating regression, replaces the one built from
#' the ADL's long-run multipliers; give it at time \eqn{t} (it is lagged
#' here), as a `ts` or `zoo` series matched by date, or as a vector with one
#' value per row of the fit's data. The model then need not be stable.
#' Note that with `x_level = TRUE`, any \eqn{ec} that is a linear
#' combination of \eqn{y} and the model's regressors (plus constant or
#' trend) gives the same \eqn{\hat\kappa} and test, since \eqn{\kappa} is
#' then the coefficient of \eqn{y_{t-1}} given the levels; a supplied
#' \eqn{\hat z} makes a difference with `x_level = FALSE`, the ADL3 form for
#' cointegrated I(1) variables. \eqn{\hat z} is a generated regressor:
#' conventional inference on \eqn{\kappa} is justified when it comes from a
#' super-consistent cointegrating regression, not in general.
#'
#' Without `ec`, the test needs a stable model (there is no equilibrium error otherwise)
#' and an adequate marginal equation: choose `lags` so that \eqn{\nu_t} is
#' free of serial correlation (the marginal regression is returned in
#' `$marginal` for checking). It addresses the part of weak exogeneity the
#' data can speak to (whether \eqn{x} responds to the equilibrium error);
#' variation-freeness of the parameters has to be argued from the economics.
#'
#' @param fit An `"adl"` object.
#' @param x The regressor tested; by default the only regressor.
#' @param lags Number of lags \eqn{s} of the changes in the marginal
#'   equation.
#' @param strong Test strong exogeneity instead of weak?
#' @param x_level Include \eqn{x_{t-1}} in the marginal regression (see
#'   Details).
#' @param ec Optional equilibrium error at time \eqn{t} (a `ts`/`zoo`
#'   series, or a vector with one value per row of the fit's data) to use
#'   instead of the one implied by the long-run multipliers (see Details).
#' @inheritParams adl_vcov
#' @param dist For the weak test, `"normal"` (default) or `"t"`; for the
#'   strong test, `"chisq"` (default) or `"F"`.
#' @param alpha Significance level for the decision.
#' @return An object of class `"adltools_test"`, with `kappa`, the
#'   marginal regression (`marginal`, an `lm` object) and `ec_source`
#'   (`"adl"` or `"supplied"`).
#' @references Urbain, J.-P. (1992). On weak exogeneity in error correction
#'   models. *Oxford Bulletin of Economics and Statistics*, 54(2), 187-207.
#'
#'   Engle, R. F., Hendry, D. F. and Richard, J.-F. (1983). Exogeneity.
#'   *Econometrica*, 51(2), 277-304.
#' @seealso [super_exogeneity_test()], [wu_hausman_test()]
#' @examples
#' fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
#' weak_exogeneity_test(fit)
#' weak_exogeneity_test(fit, strong = TRUE)
#' @export
weak_exogeneity_test <- function(fit, x = NULL, lags = 1, strong = FALSE,
                                 x_level = TRUE, ec = NULL,
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

  df <- marginal_frame(fit, x, lags, x_level, ec)
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
    kappa = b[["ec_lag1"]], marginal = marg,
    ec_source = if (is.null(ec)) "adl" else "supplied"
  ), inf), title)
}

#' Data for the regression for the change in x (weak, strong and super
#' exogeneity tests). Returns a data frame with `dx` first and the index of
#' its rows in attribute "index".
#' @noRd
marginal_frame <- function(fit, x, lags, x_level, ec = NULL) {
  z <- zoo::coredata(fit$data)
  N <- nrow(z)
  if (is.null(ec)) {
    if (!is_stable(fit)) {
      stop("The model is not dynamically stable, so there is no long-run ",
           "relation and no equilibrium error to test against (supply one ",
           "with `ec =`).", call. = FALSE)
    }
    lr <- adl_long_run(fit)
    k <- stats::setNames(lr$estimate, lr$term)
    # equilibrium error from y and the distributed-lag regressors; the
    # intercept, controls and trend enter the marginal equation directly,
    # which gives the same kappa as building their long-run terms into ec
    ec <- z[, fit$y] - as.numeric(z[, fit$x, drop = FALSE] %*% k[fit$x])
  } else {
    ec <- align_ec(ec, fit)
  }
  lagv <- function(v, j) c(rep(NA, j), v[seq_len(N - j)])
  dv <- function(v) c(NA, diff(z[, v]))
  # the test's own columns come first and keep their names; columns from
  # other regressors and controls are renamed if they would clash (a control
  # called `ec` must not replace the equilibrium error)
  cols <- list(dx = dv(x))
  others <- setdiff(fit$x, x)
  for (j in seq_len(lags)) cols[[paste0("dx_lag", j)]] <- lagv(dv(x), j)
  for (j in seq_len(lags)) cols[[paste0("dy_lag", j)]] <- lagv(dv(fit$y), j)
  cols$ec_lag1 <- lagv(ec, 1)
  if (x_level) cols$x_lag1 <- lagv(z[, x], 1)
  extra <- list()
  for (v in others) for (j in seq_len(lags)) {
    extra <- c(extra, stats::setNames(list(lagv(dv(v), j)), paste0("d", v, "_lag", j)))
  }
  if (x_level) for (v in others) {
    extra <- c(extra, stats::setNames(list(lagv(z[, v], 1)), paste0(v, "_level_lag1")))
  }
  for (v in fit$unlagged) {
    extra <- c(extra, stats::setNames(list(lagv(z[, v], 1)), paste0(v, "_lag1")))
  }
  if (fit$trend) extra <- c(extra, list(trend = seq_len(N)))
  cols <- c(cols, extra)
  names(cols) <- make.unique(names(cols), sep = ".")
  df <- as.data.frame(cols, check.names = FALSE)
  # same observations as the ADL fit where possible
  keep <- stats::complete.cases(df) & zoo::index(fit$data) %in% fit$index
  df <- df[keep, , drop = FALSE]
  if (nrow(df) <= ncol(df)) {
    stop("Too few observations for the marginal regression (", nrow(df),
         ").", call. = FALSE)
  }
  attr(df, "index") <- zoo::index(fit$data)[keep]
  df
}

#' A supplied equilibrium error, as one value per row of the fit's data
#'
#' A `ts` or `zoo` series is matched by date (dates outside its range are
#' NA); a plain numeric vector must have one value per row of `fit$data`.
#' @noRd
align_ec <- function(ec, fit) {
  idx <- zoo::index(fit$data)
  if (stats::is.ts(ec)) ec <- zoo::as.zoo(ec)
  if (inherits(ec, "zoo")) {
    if (NCOL(ec) != 1) stop("`ec` must be a single series.", call. = FALSE)
    eidx <- standard_index(zoo::index(ec))
    num <- function(i) class(i)[1] %in% c("integer", "numeric")
    if (!identical(class(eidx), class(idx)) && !(num(eidx) && num(idx))) {
      stop("The time index of `ec` is not of the same kind as that of the ",
           "data.", call. = FALSE)
    }
    vals <- as.numeric(zoo::coredata(ec))
    pos <- vapply(seq_along(idx), function(i) index_position(idx[i], eidx), 1L)
    out <- vals[pos]
  } else if (is.numeric(ec) && is.null(dim(ec))) {
    if (length(ec) != length(idx)) {
      stop("`ec` has ", length(ec), " values but the data of the fit have ",
           length(idx), " rows: give it as a ts or zoo series to match it ",
           "by date.", call. = FALSE)
    }
    out <- as.numeric(ec)
  } else {
    stop("`ec` must be a numeric vector, or a ts or zoo series.", call. = FALSE)
  }
  if (!any(is.finite(out[idx %in% fit$index]))) {
    stop("`ec` has no values on the estimation sample of the fit.", call. = FALSE)
  }
  out
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
