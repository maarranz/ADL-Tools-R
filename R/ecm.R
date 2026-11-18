# ECM view of an estimated ADL (ADL2) -----------------------------------------
#
# The ECM is a reparametrization: every ECM coefficient is a linear function
# of the ADL coefficients (plus a constant, -1, for the adjustment
# coefficient), so ECM = A theta + c and its covariance is A V A' exactly.
# Nothing is re-estimated. With
#   C(L) y_t = alpha_0 + sum_i B_i(L) x_it + controls + e_t,
# the identities
#   y_t     = y_{t-1} + Dy_t,
#   C(L)y_t = C(1) y_{t-1} + Dy_t - sum_{j=1}^{p-1} gamma_j Dy_{t-j}  (rearranged)
#   B(L)x_t = B(1) x_{t-1} + beta_0 Dx_t - sum_{j=1}^{q-1} (sum_{k>j} beta_k) Dx_{t-j}
# give
#   Dy_t = alpha_0 + b y_{t-1} + sum_i B_i(1) x_{i,t-1}
#          + sum_{j=1}^{p-1} gamma_j Dy_{t-j} + sum_i sum_{j=0}^{q_i-1} delta_ij Dx_{i,t-j}
#          + controls + e_t,
# with b = -C(1) = sum alpha - 1, gamma_j = -sum_{k>j} alpha_k,
# delta_i0 = beta_i0 and delta_ij = -sum_{k>j} beta_ik.

#' Error correction (ECM) form of an ADL fit
#'
#' Rewrites an estimated ADL model in error correction form, as in ADL2. The
#' ECM is a reparametrization of the same model, not a new estimate: the
#' residuals, the fit and every test are unchanged, and only the way the
#' coefficients are read changes. With one regressor and the general orders,
#' \deqn{\Delta y_t = \alpha_0 + b\, y_{t-1} + \beta^\star x_{t-1} +
#'   \sum_{j=1}^{p-1} \gamma_j \Delta y_{t-j} +
#'   \sum_{j=0}^{q-1} \delta_j \Delta x_{t-j} + \epsilon_t,}
#' where \eqn{b = \alpha_1 + \dots + \alpha_p - 1 = -C(1)},
#' \eqn{\beta^\star = B(1)}, \eqn{\gamma_j = -\sum_{k > j} \alpha_k},
#' \eqn{\delta_0 = \beta_0} and \eqn{\delta_j = -\sum_{k > j} \beta_k}. For
#' the ADL(1,1) these are the \eqn{\alpha^\star_1 = \alpha_1 - 1} and
#' \eqn{\beta^\star_1 = \beta_0 + \beta_1} of ADL2. Each regressor gets its
#' own \eqn{\beta^\star_i} and \eqn{\delta_{ij}}; unlagged terms and the trend
#' keep their coefficients.
#'
#' Every ECM coefficient is a linear function of the ADL coefficients, so its
#' standard error is exact given the covariance of the ADL (conventional or
#' HAC): it is the one the ECM regression itself would report. A lag that is
#' not in the ADL is a zero restriction, and the ECM terms are reported with
#' the restrictions it implies (e.g. without lag 0 of \eqn{x}, the
#' \eqn{\Delta x_t} term is absent). Without lags of \eqn{y}, \eqn{b = -1}
#' is not estimated and is shown with a zero standard error.
#'
#' The **equilibrium-error form** of ADL2,
#' \deqn{\Delta y_t = \alpha_0 + b\,(y_{t-1} - k_1 x_{1,t-1} - \dots -
#'   k_n x_{n,t-1}) + \dots,}
#' is reported separately: the speed of adjustment \eqn{b} and the long-run
#' multipliers \eqn{k_i = -\beta^\star_i / b = B_i(1)/C(1)}, with
#' delta-method standard errors (those of [adl_long_run()]). The
#' equilibrium error \eqn{y_t - \sum_i k_i x_{it}} is returned as a series,
#' defined as in [weak_exogeneity_test()].
#'
#' Reading \eqn{b} and \eqn{k_i} off one equation requires a predetermined
#' and weakly exogenous \eqn{x} (ADL1, ADL3), and a stable model. With
#' integrated series, the long-run relation is meaningful only if the
#' variables are cointegrated, which this function does not test; the usual
#' t-statistic of \eqn{b} does not then have a normal distribution.
#'
#' @param fit An `"adl"` object.
#' @inheritParams adl_vcov
#' @param alpha Significance level for the confidence intervals.
#' @param dist `"normal"` (default) or `"t"`, as in [adl_coefficients()].
#' @return An object of class `"adl_ecm"`: a list with `coefficients` (one row
#'   per ECM term: `term`, `param`, `estimate`, `std.error`, `statistic`,
#'   `p.value`, `conf.low`, `conf.high`), `equilibrium` (the same columns
#'   for \eqn{b} and each \eqn{k_i}), `ec` (the equilibrium error, a `zoo`
#'   series over the data of the fit; `NULL` for an unstable model), `vcov`
#'   (covariance of the ECM coefficients), `A` and `c` (the map
#'   ECM = A theta + c from the ADL coefficients theta), `stable`, `y`,
#'   `x`, and the settings `cov_type`, `maxlags`, `alpha` and `dist`.
#' @seealso [adl_long_run()], [adl_dynamics()], [weak_exogeneity_test()].
#' @examples
#' fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
#' ecm <- adl_ecm(fit)
#' ecm
#' plot(ecm$ec, ylab = "Equilibrium error")
#'
#' # several lags: the ECM has lagged differences
#' adl_ecm(adl(y ~ x, data = adl_sim, order = c(2, 2)), cov_type = "HAC")
#' @export
adl_ecm <- function(fit, cov_type = c("nonrobust", "HAC"), maxlags = NULL,
                    alpha = 0.05, dist = c("normal", "t")) {
  check_adl(fit)
  cov_type <- match.arg(cov_type)
  dist <- match.arg(dist)
  check_alpha(alpha)
  theta <- stats::coef(fit)
  V <- adl_vcov(fit, cov_type, maxlags)
  info <- term_info(fit, names(theta))
  K <- length(theta)
  several <- length(fit$x) > 1

  rows <- list()
  add <- function(term, param, a, c0 = 0) {
    rows[[length(rows) + 1]] <<- list(term = term, param = param, a = a, c = c0)
  }
  unit <- function(pos) { a <- numeric(K); a[pos] <- 1; a }

  if (fit$intercept) add("(Intercept)", "alpha_0", unit(which(names(theta) == "(Intercept)")))

  ar_pos <- poly_positions(info, fit$y, from = 1)
  p <- length(ar_pos)
  a <- numeric(K)
  a[ar_pos[!is.na(ar_pos)]] <- 1
  add(paste0(fit$y, "_lag1"), "b", a, -1)

  dl_pos <- lapply(stats::setNames(fit$x, fit$x), function(v) poly_positions(info, v, from = 0))
  for (i in seq_along(fit$x)) {
    v <- fit$x[i]
    a <- numeric(K)
    a[dl_pos[[v]][!is.na(dl_pos[[v]])]] <- 1
    add(paste0(v, "_lag1"), if (several) paste0("beta*_", i) else "beta*", a)
  }
  for (j in seq_len(max(0, p - 1))) {
    a <- numeric(K)
    later <- ar_pos[seq.int(j + 1, p)]
    a[later[!is.na(later)]] <- -1
    add(paste0("D_", fit$y, "_lag", j), paste0("gamma_", j), a)
  }
  for (i in seq_along(fit$x)) {
    v <- fit$x[i]
    pos <- dl_pos[[v]]
    q <- length(pos) - 1
    lab <- function(j) if (several) paste0("delta_", i, ",", j) else paste0("delta_", j)
    if (!is.na(pos[1])) add(paste0("D_", v), lab(0), unit(pos[1]))
    for (j in seq_len(max(0, q - 1))) {
      a <- numeric(K)
      later <- pos[seq.int(j + 2, q + 1)]
      later <- later[!is.na(later)]
      if (!length(later)) next
      a[later] <- -1
      add(paste0("D_", v, "_lag", j), lab(j), a)
    }
  }
  other <- which(!info$role %in% c("ar", "dl", "intercept"))
  for (k in other) add(names(theta)[k], "", unit(k))

  A <- do.call(rbind, lapply(rows, `[[`, "a"))
  colnames(A) <- names(theta)
  cvec <- vapply(rows, `[[`, numeric(1), "c")
  terms <- vapply(rows, `[[`, character(1), "term")
  rownames(A) <- names(cvec) <- terms
  est <- drop(A %*% theta) + cvec
  Vecm <- A %*% V %*% t(A)
  se <- sqrt(pmax(diag(Vecm), 0))
  tab <- ecm_table(terms, est, se, dist, df_resid(fit), alpha)
  tab <- cbind(tab[1], param = vapply(rows, `[[`, character(1), "param"), tab[-1])

  # equilibrium-error form: b and the long-run multipliers
  stable <- is_stable(fit)
  lr <- suppressWarnings(adl_long_run(fit, cov_type = cov_type, maxlags = maxlags,
                                      alpha = alpha, dist = dist))
  lr <- lr[lr$term %in% fit$x, , drop = FALSE]
  lr$term <- paste0("k_", lr$term)
  b_row <- tab[tab$param == "b", setdiff(names(tab), "param"), drop = FALSE]
  b_row$term <- "b"
  equil <- rbind(b_row, lr)
  equil <- cbind(equil[1],
                 param = c("b", if (several) paste0("k_", seq_along(fit$x)) else "k_1"),
                 equil[-1])
  rownames(equil) <- NULL
  if (!stable) {
    warning("The model is not dynamically stable: there is no long-run ",
            "relation, so the k_i and the equilibrium error are NA/NULL. The ",
            "ECM coefficients are still a valid reparametrization.", call. = FALSE)
  }

  ec <- NULL
  if (stable) {
    k <- stats::setNames(lr$estimate, fit$x)
    z <- zoo::coredata(fit$data)
    ec <- zoo::zoo(z[, fit$y] - as.numeric(z[, fit$x, drop = FALSE] %*% k),
                   zoo::index(fit$data))
  }

  structure(list(
    coefficients = tab,
    equilibrium = equil,
    ec = ec,
    vcov = Vecm,
    A = A,
    c = cvec,
    stable = stable,
    y = fit$y,
    x = fit$x,
    unlagged = fit$unlagged,
    label = fit$label,
    nobs = stats::nobs(fit),
    cov_type = cov_type,
    maxlags = attr(V, "maxlags"),
    alpha = alpha,
    dist = dist
  ), class = "adl_ecm")
}

#' coef_table() with fixed (zero-variance) terms
#' @noRd
ecm_table <- function(terms, b, se, dist, df, alpha) {
  tab <- coef_table(terms, b, se, dist, df, alpha)
  fixed <- se <= sqrt(.Machine$double.eps) * pmax(1, abs(b))
  tab$std.error[fixed] <- 0
  tab$statistic[fixed] <- tab$p.value[fixed] <- NA_real_
  tab$conf.low[fixed] <- tab$conf.high[fixed] <- b[fixed]
  tab
}

#' Print an ECM view
#'
#' @param x An `"adl_ecm"` object.
#' @param digits Significant digits.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.adl_ecm <- function(x, digits = 4, ...) {
  title <- paste0("Error correction form of the ", x$label)
  cat(title, "\n", strrep("-", nchar(title)), "\n", sep = "")
  cat("Dependent variable: D_", x$y, " (first difference of ", x$y, ")\n", sep = "")
  cat("Same model as the ADL: only the parametrization changes.\n\n")
  fmt <- function(tab) {
    tab[] <- lapply(tab, function(col) if (is.double(col)) signif(col, digits) else col)
    tab
  }
  print(fmt(x$coefficients), row.names = FALSE)
  cat("\nEquilibrium-error form: D_", x$y, " = ... + b (", x$y, "_lag1 - ",
      paste0("k_", x$x, " ", x$x, "_lag1", collapse = " - "), ") + ...\n", sep = "")
  print(fmt(x$equilibrium), row.names = FALSE)
  if (!x$stable) cat("\nNot dynamically stable: no long-run relation.\n")
  cat("\nStandard errors: ", x$cov_type,
      if (identical(x$cov_type, "HAC")) paste0(" (bandwidth ", x$maxlags, ")"),
      "; exact for the ECM coefficients, delta method for the k's. ",
      "Distribution: ", x$dist, ".\n", sep = "")
  cat("b and k assume a weakly exogenous x; with integrated series, a ",
      "cointegrated one.\n", sep = "")
  invisible(x)
}
