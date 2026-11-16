# broom-style methods ---------------------------------------------------------
#
# tidy() and glance() generics come from the lightweight 'generics' package
# (the one broom uses), so they work whether or not broom is loaded.

#' @importFrom generics tidy glance
#' @export
generics::tidy

#' @export
generics::glance

#' Tidy and glance methods
#'
#' `tidy()` returns one row per coefficient, test or lag as a plain data
#' frame; `glance()` returns a one-row summary of a fit or of a set of
#' Granger tests. The methods use the same conventions as the rest of the
#' package (normal/chi-square inference by default).
#'
#' @param x An object from adltools.
#' @param conf.int Ignored; confidence limits are always included.
#' @param conf.level Confidence level for the intervals.
#' @param ... Passed on: for `"adl"` fits, `cov_type`, `maxlags` and `dist`
#'   (see [adl_coefficients()]).
#' @return A data frame.
#' @name adltools_tidiers
#' @examples
#' set.seed(1)
#' x <- cumsum(rnorm(100))
#' y <- stats::filter(0.5 * x + rnorm(100), 0.6, method = "recursive")
#' fit <- adl(y ~ x, data = ts(cbind(y = as.numeric(y), x = x)), order = c(1, 1))
#' tidy(fit)
#' glance(fit)
#' tidy(adl_joint_test(fit, variables = "x"))
NULL

#' @rdname adltools_tidiers
#' @export
tidy.adl <- function(x, conf.int = TRUE, conf.level = 0.95, ...) {
  tab <- adl_coefficients(x, alpha = 1 - conf.level, ...)
  out <- as.data.frame(tab)
  attr(out, "cov_type") <- attr(out, "maxlags") <- attr(out, "distribution") <- NULL
  out
}

#' @rdname adltools_tidiers
#' @export
glance.adl <- function(x, ...) {
  s <- summary(x$lm)
  slopes <- setdiff(names(stats::coef(x)), "(Intercept)")
  w <- if (length(slopes)) adl_joint_test(x, terms = slopes)
  n <- stats::nobs(x)
  k <- length(stats::coef(x))
  ll <- as.numeric(stats::logLik(x))
  data.frame(
    model = x$label,
    r.squared = s$r.squared,
    adj.r.squared = s$adj.r.squared,
    sigma = s$sigma,
    statistic = if (is.null(w)) NA_real_ else w$statistic,
    p.value = if (is.null(w)) NA_real_ else w$p.value,
    df = length(slopes),
    logLik = ll,
    AIC = stats::AIC(x),
    BIC = stats::BIC(x),
    aic_eviews = (-2 * ll + 2 * k) / n,
    bic_eviews = (-2 * ll + log(n) * k) / n,
    nobs = n,
    df.residual = stats::df.residual(x$lm),
    start = format(x$index[1]),
    end = format(x$index[length(x$index)])
  )
}

#' @rdname adltools_tidiers
#' @export
tidy.adltools_test <- function(x, ...) {
  data.frame(
    method = x$method,
    null_hypothesis = x$null_hypothesis,
    statistic = x$statistic,
    distribution = x$distribution,
    df = x$df, df_num = x$df_num, df_denom = x$df_denom,
    p.value = x$p.value,
    reject_null = x$reject_null,
    alpha = x$alpha,
    cov_type = x$cov_type %||% NA_character_,
    maxlags = x$maxlags,
    nobs = x$nobs
  )
}

#' @rdname adltools_tidiers
#' @export
tidy.granger_test <- function(x, ...) {
  x$tests
}

#' @rdname adltools_tidiers
#' @export
glance.granger_test <- function(x, ...) {
  cl <- x$classification
  data.frame(
    y = x$y, x = x$x,
    type = paste(cl$type, collapse = " or "),
    symbol = paste(cl$symbol, collapse = " or "),
    description = paste(cl$description, collapse = " or "),
    ambiguous = cl$ambiguous,
    dist = x$dist, alpha = x$alpha, cov_type = x$cov_type, maxlags = x$maxlags,
    nobs = x$nobs, start = x$sample[1], end = x$sample[2]
  )
}

#' @rdname adltools_tidiers
#' @export
tidy.adl_dynamics <- function(x, ...) {
  x$effects
}

#' @rdname adltools_tidiers
#' @export
glance.adl_dynamics <- function(x, ...) {
  out <- x$summary
  out$stable <- x$stable
  out$max_inverse_root <- if (nrow(x$roots)) max(x$roots$inverse_modulus) else 0
  out
}
