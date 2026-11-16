# Estimation ----------------------------------------------------------------

#' Estimate an autoregressive distributed lag model
#'
#' Fits an ADL model by OLS. The model can be written in three ways:
#'
#' * a formula with variables only, `y ~ x1 + x2`, plus `order` (and
#'   optionally `first_lag`) for the lag lengths;
#' * a formula with variables only plus an explicit `lags` list;
#' * a `dynlm`-style formula with the lags written in it,
#'   `y ~ L(y, 1:2) + L(x1, 0:2) + L(x2, 1)`.
#'
#' Variables after `|` enter without lags, e.g. dummies:
#' `y ~ x | d2008 + d2020`.
#'
#' @param formula Model formula (see Details).
#' @param data A `ts`/`mts`, `zoo`/`xts`, data frame, tibble or tsibble.
#' @param order Lag orders: `c(p, q)` with the same `q` for every regressor,
#'   or `c(p, q_1, ..., q_n)`. Lags of `y` run from 1 to `p`; lags of each
#'   regressor from `first_lag` to `q`.
#' @param first_lag 0 (the default) includes the regressor at time t; 1
#'   leaves it out. One value for all regressors, or a named vector such as
#'   `c(x2 = 1)` (unnamed regressors keep 0).
#' @param lags Named list of lag sets, e.g. `list(y = 1:2, x1 = 0:1,
#'   x2 = c(1, 4))`. Overrides `order`.
#' @param intercept Include an intercept?
#' @param trend Include a linear trend?
#' @param time Name of the time column when `data` is a data frame; detected
#'   automatically when there is a single date column.
#' @param start,end Estimation sample, in the units of the time index (e.g.
#'   `c(1990, 1)` for quarterly `ts` data, or `"1990-01-01"` for dates). The
#'   first `max_lag` observations of this sample are used as initial values.
#' @param max_lag Number of initial observations to set aside. Defaults to
#'   the largest lag; set it larger to estimate several models on a common
#'   sample.
#' @return An object of class `"adl"`; see [print.adl()] and
#'   [summary.adl()]. Standard methods (`coef()`, `vcov()`, `residuals()`,
#'   `fitted()`, `nobs()`, `logLik()`, `AIC()`, `BIC()`) work on it, and the
#'   underlying `lm` fit is in `$lm`.
#' @examples
#' set.seed(1)
#' n <- 100
#' x <- cumsum(rnorm(n))
#' y <- stats::filter(0.5 * x + rnorm(n), 0.6, method = "recursive")
#' d <- ts(cbind(y = as.numeric(y), x = x), start = c(2000, 1), frequency = 4)
#'
#' adl(y ~ x, data = d, order = c(1, 1))
#' adl(y ~ L(y, 1) + L(x, 0:1), data = d)          # the same model
#' adl(y ~ x, data = d, order = c(1, 1), first_lag = 1)  # without x_t
#' @export
adl <- function(formula, data, order = NULL, first_lag = 0, lags = NULL,
                intercept = TRUE, trend = FALSE, time = NULL,
                start = NULL, end = NULL, max_lag = NULL) {
  cl <- match.call()
  spec <- parse_adl_formula(formula)

  if (!is.null(spec$lags)) {
    if (!is.null(order) || !is.null(lags)) {
      stop("The formula already contains L() terms; do not also give ",
           "`order` or `lags`.", call. = FALSE)
    }
    lags <- spec$lags
  }
  lag_sets <- normalize_lags(spec$y, spec$x, order = order,
                             first_lag = first_lag, lags = lags)

  d <- as_adl_data(data, c(spec$y, spec$x, spec$unlagged), time = time,
                   start = start, end = end)
  des <- build_design(d$z, lag_sets, unlagged = spec$unlagged, trend = trend,
                      max_lag = max_lag)

  rhs <- setdiff(names(des$frame), des$y_col)
  if (length(rhs) == 0 && !intercept) stop("The model has no regressors.", call. = FALSE)
  f <- stats::reformulate(if (length(rhs)) rhs else "1", response = des$y_col,
                          intercept = intercept)
  frame <- des$frame
  fit <- stats::lm(f, data = frame)
  if (any(is.na(stats::coef(fit)))) {
    bad <- names(stats::coef(fit))[is.na(stats::coef(fit))]
    stop("Perfect collinearity; cannot estimate: ", paste(bad, collapse = ", "),
         call. = FALSE)
  }

  structure(list(
    lm = fit,
    call = cl,
    y = spec$y,
    x = spec$x,
    unlagged = spec$unlagged,
    lags = lag_sets,
    intercept = intercept,
    trend = trend,
    terms = des$terms,
    index = des$index,
    max_lag = des$max_lag,
    data = d$z,
    label = adl_label(lag_sets),
    dynlm = dynlm_formula(lag_sets, spec$unlagged, trend, intercept)
  ), class = "adl")
}

# Methods -------------------------------------------------------------------

#' Print and summarize an ADL fit
#'
#' `summary()` gives the coefficient table, with the parameter of each
#' coefficient in the notation of the ATSE chapters (\eqn{\alpha_0} for the
#' intercept, \eqn{\alpha_i} for lag i of the dependent variable,
#' \eqn{\beta_j} for lag j of a single regressor, \eqn{\beta_{ij}} for lag
#' j of regressor i), and a Wald test that all coefficients except the
#' intercept are zero. As in [adl_coefficients()], p-values use the normal
#' and chi-square distributions by default; `dist = "t"` reproduces
#' `summary.lm()`.
#'
#' @param x,object An `"adl"` object (`x` is a `"summary.adl"` object in its
#'   print method).
#' @param digits Number of significant digits.
#' @inheritParams adl_vcov
#' @param dist `"normal"` (default) or `"t"`.
#' @param ... Not used.
#' @return `print()` returns `x` invisibly; `summary()` returns an object of
#'   class `"summary.adl"`.
#' @export
print.adl <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(x$label, "\n", sep = "")
  cat("Dependent variable: ", x$y, "\n", sep = "")
  cat("Sample: ", sample_text(x), "\n", sep = "")
  cat("dynlm equivalent: ", x$dynlm, "\n\n", sep = "")
  print.default(format(stats::coef(x), digits = digits), print.gap = 2L,
                quote = FALSE)
  invisible(x)
}

#' @rdname print.adl
#' @export
summary.adl <- function(object, cov_type = c("nonrobust", "HAC"),
                        maxlags = NULL, dist = c("normal", "t"), ...) {
  cov_type <- match.arg(cov_type)
  dist <- match.arg(dist)
  tab <- adl_coefficients(object, cov_type = cov_type, maxlags = maxlags,
                          dist = dist)
  slopes <- setdiff(names(stats::coef(object)), "(Intercept)")
  overall <- if (length(slopes)) {
    adl_joint_test(object, terms = slopes, cov_type = cov_type,
                   maxlags = maxlags, dist = if (dist == "t") "F" else "chisq")
  }
  s <- summary(object$lm)
  structure(list(object = object, coefficients = tab, overall = overall,
                 sigma = s$sigma, df_resid = s$df[2], r.squared = s$r.squared,
                 adj.r.squared = s$adj.r.squared, cov_type = cov_type,
                 maxlags = attr(tab, "maxlags"), dist = dist),
            class = "summary.adl")
}

#' @rdname print.adl
#' @export
print.summary.adl <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  o <- x$object
  cat(o$label, " estimated by OLS\n", sep = "")
  cat("Dependent variable: ", o$y, "\n", sep = "")
  cat("Sample: ", sample_text(o), "\n", sep = "")
  cat("dynlm equivalent: ", o$dynlm, "\n", sep = "")
  cov <- if (x$cov_type == "HAC") paste0("HAC (Bartlett, ", x$maxlags, " lags)")
         else "conventional OLS"
  cat("Covariance: ", cov, "\n\n", sep = "")
  tab <- x$coefficients
  m <- cbind(tab$estimate, tab$std.error, tab$statistic, tab$p.value)
  zt <- if (x$dist == "t") "t" else "z"
  colnames(m) <- c("Estimate", "Std. Error", paste(zt, "value"),
                   paste0("Pr(>|", zt, "|)"))
  rownames(m) <- ifelse(tab$param == "", tab$term,
                        paste0(tab$term, " [", tab$param, "]"))
  stats::printCoefmat(m, digits = digits, signif.legend = TRUE)
  cat("\nResidual standard error: ", format(x$sigma, digits = digits),
      " on ", x$df_resid, " degrees of freedom\n", sep = "")
  cat("R-squared: ", formatC(x$r.squared, digits = digits),
      ",  Adjusted R-squared: ", formatC(x$adj.r.squared, digits = digits), "\n",
      sep = "")
  if (!is.null(x$overall)) {
    w <- x$overall
    lab <- if (w$distribution == "F") paste0("F(", w$df_num, ", ", w$df_denom, ")")
           else paste0("Chi-square(", w$df, ")")
    cat("Wald test, all slopes zero: ", lab, " = ", format(signif(w$statistic, digits)),
        ",  p-value: ", format.pval(w$p.value, digits = digits), "\n", sep = "")
  }
  invisible(x)
}

sample_text <- function(o) {
  paste0(format(o$index[1]), " - ", format(o$index[length(o$index)]),
         " (", length(o$index), " observations)")
}

#' @export
coef.adl <- function(object, ...) stats::coef(object$lm)
#' @export
vcov.adl <- function(object, ...) stats::vcov(object$lm)
#' @export
residuals.adl <- function(object, ...) {
  zoo::zoo(stats::residuals(object$lm), object$index)
}
#' @export
fitted.adl <- function(object, ...) {
  zoo::zoo(stats::fitted(object$lm), object$index)
}
#' @export
nobs.adl <- function(object, ...) stats::nobs(object$lm)
#' @export
logLik.adl <- function(object, ...) stats::logLik(object$lm)
#' @export
model.matrix.adl <- function(object, ...) stats::model.matrix(object$lm)
#' @export
formula.adl <- function(x, ...) stats::formula(x$lm)
