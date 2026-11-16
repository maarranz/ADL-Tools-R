# Inference on ADL coefficients ---------------------------------------------
#
# Distributions: the asymptotic normal and chi-square by default, with either
# covariance. In a dynamic regression the exact t and F results do not hold
# even with nonrobust errors; `dist` gives the t/F versions on request (the
# figures EViews and summary.lm() report).

#' Coefficient table with conventional or HAC standard errors
#'
#' @param fit An `"adl"` object.
#' @inheritParams adl_vcov
#' @param alpha Significance level for the confidence intervals.
#' @param dist `"normal"` (default) or `"t"`: the reference distribution for
#'   p-values and intervals. The regression is dynamic, so both are
#'   asymptotic approximations; `"t"` reproduces `summary.lm()` and EViews.
#' @return A data frame with one row per coefficient: `term` (column name),
#'   `param` (the parameter in the chapters' notation: `alpha_0`, `alpha_i`,
#'   `beta_j` with one regressor, `beta_i,j` with several; empty for
#'   unlagged terms and the trend), `variable`, `lag`, `role`, `estimate`,
#'   `std.error`, `statistic`, `p.value`,
#'   `conf.low`, `conf.high`. Attributes `cov_type`, `maxlags` and
#'   `distribution` record how it was computed.
#' @examples
#' set.seed(1)
#' x <- cumsum(rnorm(100))
#' y <- stats::filter(0.5 * x + rnorm(100), 0.6, method = "recursive")
#' fit <- adl(y ~ x, data = ts(cbind(y = as.numeric(y), x = x)), order = c(1, 1))
#' adl_coefficients(fit)
#' adl_coefficients(fit, cov_type = "HAC")
#' @export
adl_coefficients <- function(fit, cov_type = c("nonrobust", "HAC"),
                             maxlags = NULL, alpha = 0.05,
                             dist = c("normal", "t")) {
  check_adl(fit)
  cov_type <- match.arg(cov_type)
  dist <- match.arg(dist)
  check_alpha(alpha)
  V <- adl_vcov(fit, cov_type, maxlags)
  b <- stats::coef(fit)
  se <- sqrt(diag(V))
  tab <- coef_table(names(b), b, se, dist, df_resid(fit), alpha)
  info <- term_info(fit, names(b))
  tab <- cbind(tab[1], param = param_labels(fit, info), info, tab[-1])
  structure(tab, cov_type = cov_type, maxlags = attr(V, "maxlags"),
            distribution = dist)
}

#' Wald test of linear restrictions on the coefficients
#'
#' Tests \eqn{H_0: R\beta = r}. The restrictions can be given in three ways:
#'
#' * `terms = c("x_lag1", "x_lag2")` sets those coefficients to zero (the
#'   column of a regressor at lag 0 has the plain name, e.g. `"x"`);
#' * `variables = "x"` sets all the coefficients of `x` (every lag,
#'   including lag 0) to zero;
#' * `R` and `r` give general linear restrictions, with the columns of `R`
#'   in the order of `coef(fit)`.
#'
#' `terms` and `variables` can be combined.
#'
#' @param fit An `"adl"` object.
#' @param terms Coefficient names set to zero under the null.
#' @param variables Variable names whose coefficients (all lags) are zero
#'   under the null.
#' @param R,r Restriction matrix and right-hand side (default zero), used
#'   when `terms` and `variables` are `NULL`.
#' @inheritParams adl_vcov
#' @param alpha Significance level for the decision.
#' @param dist `"chisq"` (default) or `"F"`. The regression is dynamic, so
#'   both are asymptotic approximations; `"F"` (the Wald statistic divided
#'   by the number of restrictions) reproduces EViews and `anova()`.
#' @return An object of class `"adltools_test"`.
#' @examples
#' set.seed(1)
#' x <- cumsum(rnorm(100))
#' y <- stats::filter(0.5 * x + rnorm(100), 0.6, method = "recursive")
#' fit <- adl(y ~ x, data = ts(cbind(y = as.numeric(y), x = x)), order = c(1, 2))
#' adl_joint_test(fit, terms = c("x_lag1", "x_lag2"))
#' adl_joint_test(fit, variables = "x", cov_type = "HAC")
#' # beta_0 + beta_1 + beta_2 = 0.5
#' R <- matrix(c(0, 0, 1, 1, 1), nrow = 1)
#' adl_joint_test(fit, R = R, r = 0.5)
#' @export
adl_joint_test <- function(fit, terms = NULL, variables = NULL, R = NULL, r = NULL,
                           cov_type = c("nonrobust", "HAC"), maxlags = NULL,
                           alpha = 0.05, dist = c("chisq", "F")) {
  check_adl(fit)
  cov_type <- match.arg(cov_type)
  check_alpha(alpha)
  dist <- match.arg(dist)
  b <- stats::coef(fit)
  k <- length(b)

  if (!is.null(terms) || !is.null(variables)) {
    if (!is.null(R)) stop("Give `terms`/`variables` or `R`, not both.", call. = FALSE)
    cols <- resolve_terms(fit, terms, variables)
    R <- matrix(0, length(cols), k)
    R[cbind(seq_along(cols), match(cols, names(b)))] <- 1
    r <- rep(0, length(cols))
    h0 <- paste(cols, collapse = " = ")
    h0 <- paste0(h0, " = 0")
  } else {
    if (is.null(R)) stop("Give `terms`, `variables` or `R`.", call. = FALSE)
    R <- if (is.matrix(R)) R else matrix(R, nrow = 1)
    if (ncol(R) != k) stop("`R` must have ", k, " columns (one per coefficient).",
                           call. = FALSE)
    r <- if (is.null(r)) rep(0, nrow(R)) else r
    if (length(r) != nrow(R)) stop("`r` must have one value per row of `R`.",
                                   call. = FALSE)
    if (qr(R)$rank < nrow(R)) stop("The restrictions in `R` are linearly dependent.",
                                   call. = FALSE)
    h0 <- restriction_text(R, r, names(b))
  }

  V <- adl_vcov(fit, cov_type, maxlags)
  q <- nrow(R)
  w <- wald_stat(b, V, R, r)
  inf <- wald_inference(w, q, df_resid(fit), dist)
  new_adl_test(c(list(
    method = "wald",
    null_hypothesis = h0,
    alpha = alpha,
    reject_null = inf$p.value < alpha,
    nobs = stats::nobs(fit),
    cov_type = cov_type,
    maxlags = attr(V, "maxlags"),
    n_restrictions = q,
    df_resid = df_resid(fit)
  ), inf), "Wald test of linear restrictions")
}

#' Long-run multipliers
#'
#' Computes the long-run multiplier of each regressor, \eqn{k_i = B_i(1) /
#' C(1)}, and the long-run intercept \eqn{k_0 = \alpha_0 / C(1)}, where
#' \eqn{C(1) = 1 - \sum \alpha_i} and \eqn{B_i(1)} is the sum of the lag
#' coefficients of regressor i. Standard errors by the delta method.
#'
#' The long run exists only when the model is dynamically stable (all roots
#' of \eqn{C(z)} outside the unit circle); otherwise the multipliers are
#' `NA`, with a warning. The delta-method standard errors assume a stable
#' model with \eqn{C(1)} well away from zero: near a unit root the ratio is
#' weakly identified and symmetric intervals can be misleading. For an
#' unlagged dummy, the reported value is the long-run effect of keeping the
#' dummy at one permanently; for a one-off pulse dummy that is a
#' hypothetical, not the cumulated effect of the pulse.
#'
#' @param fit An `"adl"` object.
#' @inheritParams adl_vcov
#' @param alpha Significance level for the confidence intervals.
#' @param dist `"normal"` (default) or `"t"`. The delta method is asymptotic
#'   in any case.
#' @return A data frame with `term`, `estimate`, `std.error`, `statistic`
#'   (for \eqn{H_0: k = 0}), `p.value`, `conf.low`, `conf.high`, and
#'   attributes `cov_type`, `maxlags`, `distribution`, `C1` and `stable`.
#' @examples
#' set.seed(1)
#' x <- cumsum(rnorm(100))
#' y <- stats::filter(0.5 * x + rnorm(100), 0.6, method = "recursive")
#' fit <- adl(y ~ x, data = ts(cbind(y = as.numeric(y), x = x)), order = c(1, 1))
#' adl_long_run(fit)
#' @export
adl_long_run <- function(fit, cov_type = c("nonrobust", "HAC"), maxlags = NULL,
                         alpha = 0.05, dist = c("normal", "t")) {
  check_adl(fit)
  cov_type <- match.arg(cov_type)
  dist <- match.arg(dist)
  check_alpha(alpha)
  b <- stats::coef(fit)
  V <- adl_vcov(fit, cov_type, maxlags)
  info <- term_info(fit, names(b))

  ar <- info$role == "ar"
  C1 <- 1 - sum(b[ar])
  stable <- is_stable(fit)
  if (!stable) {
    warning("The model is not dynamically stable (a root of C(z) is on or ",
            "inside the unit circle): the long-run multipliers do not exist and ",
            "are reported as NA.", call. = FALSE)
  }
  rows <- list()
  if (fit$intercept) rows[["(Intercept)"]] <- names(b) == "(Intercept)"
  for (v in c(fit$x, fit$unlagged)) rows[[v]] <- info$variable == v & info$role != "ar"
  est <- se <- numeric(length(rows))
  for (i in seq_along(rows)) {
    sel <- rows[[i]]
    B1 <- sum(b[sel])
    est[i] <- B1 / C1
    g <- numeric(length(b))
    g[sel] <- 1 / C1
    g[ar] <- B1 / C1^2
    se[i] <- sqrt(as.numeric(t(g) %*% V %*% g))
  }
  if (!stable) est[] <- se[] <- NA_real_
  tab <- coef_table(names(rows), est, se, dist, df_resid(fit), alpha)
  structure(tab, cov_type = cov_type, maxlags = attr(V, "maxlags"),
            distribution = dist,
            C1 = C1, stable = stable)
}

# Helpers -------------------------------------------------------------------

df_resid <- function(fit) stats::df.residual(fit$lm)

check_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a number between 0 and 1.", call. = FALSE)
  }
  invisible(alpha)
}

#' Estimates, t/z statistics, p-values and intervals
#' @noRd
coef_table <- function(terms, b, se, dist, df, alpha) {
  stat <- b / se
  if (dist == "t") {
    p <- 2 * stats::pt(abs(stat), df, lower.tail = FALSE)
    crit <- stats::qt(1 - alpha / 2, df)
  } else {
    p <- 2 * stats::pnorm(abs(stat), lower.tail = FALSE)
    crit <- stats::qnorm(1 - alpha / 2)
  }
  data.frame(term = terms, estimate = unname(b), std.error = unname(se),
             statistic = unname(stat), p.value = unname(p),
             conf.low = unname(b - crit * se), conf.high = unname(b + crit * se),
             row.names = NULL, stringsAsFactors = FALSE)
}

#' Variable, lag and role of each coefficient
#' @noRd
term_info <- function(fit, terms) {
  m <- match(terms, fit$terms$column)
  data.frame(
    variable = ifelse(is.na(m), terms, fit$terms$variable[m]),
    lag = ifelse(is.na(m), NA_integer_, fit$terms$lag[m]),
    role = ifelse(is.na(m), ifelse(terms == "(Intercept)", "intercept", "other"),
                  fit$terms$role[m]),
    stringsAsFactors = FALSE
  )
}

#' Parameter labels in the notation of the chapters
#' @noRd
param_labels <- function(fit, info) {
  several <- length(fit$x) > 1
  vapply(seq_len(nrow(info)), function(k) {
    switch(info$role[k],
      intercept = "alpha_0",
      ar = paste0("alpha_", info$lag[k]),
      dl = if (several) paste0("beta_", match(info$variable[k], fit$x), ",", info$lag[k])
           else paste0("beta_", info$lag[k]),
      "")
  }, "")
}

#' Expand variable names to their coefficient columns
#' @noRd
resolve_terms <- function(fit, terms = NULL, variables = NULL) {
  b <- names(stats::coef(fit))
  bad <- setdiff(terms, b)
  if (length(bad)) {
    stop("Not coefficient(s) of the model: ", paste(bad, collapse = ", "),
         ". Coefficients: ", paste(b, collapse = ", "), call. = FALSE)
  }
  vars <- unique(fit$terms$variable)
  bad <- setdiff(variables, vars)
  if (length(bad)) {
    stop("Not variable(s) of the model: ", paste(bad, collapse = ", "),
         ". Variables: ", paste(vars, collapse = ", "), call. = FALSE)
  }
  from_vars <- fit$terms$column[fit$terms$variable %in% variables]
  b[b %in% c(terms, from_vars)]
}

#' Readable form of R b = r
#' @noRd
restriction_text <- function(R, r, names) {
  rows <- vapply(seq_len(nrow(R)), function(i) {
    nz <- which(R[i, ] != 0)
    lhs <- paste0(ifelse(R[i, nz] == 1, "", paste0(format(R[i, nz]), "*")),
                  names[nz], collapse = " + ")
    paste(lhs, "=", format(r[i]))
  }, "")
  paste(rows, collapse = "; ")
}

#' Are all roots of C(z) = 1 - a_1 z - ... - a_p z^p outside the unit circle?
#' @noRd
is_stable <- function(fit) {
  lags <- fit$lags[[1]]
  if (length(lags) == 0) return(TRUE)
  b <- stats::coef(fit)
  a <- numeric(max(lags))
  info <- term_info(fit, names(b))
  a[info$lag[info$role == "ar"]] <- b[info$role == "ar"]
  roots <- polyroot(c(1, -a))
  all(Mod(roots) > 1 + sqrt(.Machine$double.eps))
}
