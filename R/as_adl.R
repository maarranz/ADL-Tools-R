# Existing dynlm fits ---------------------------------------------------------

#' Convert a dynlm fit into an adl fit
#'
#' Re-estimates a model fitted with `dynlm::dynlm()` as an `"adl"` object,
#' so that the tests and tools of this package can be applied to it. The
#' formula may contain variables, `L(variable, lags)` terms and `trend()`;
#' other `dynlm` operators (`d()`, `season()`, transformations) are not
#' supported. The data are taken from the original call unless `data` is
#' given.
#'
#' The estimates are the same, except that `adl()` counts the trend in
#' observations (1, 2, ...) while `dynlm::trend()` counts in years (1/4 per
#' quarter), so only the trend coefficient (and the intercept) differ in
#' scale.
#'
#' @param x A `"dynlm"` object.
#' @param data The data used in the original call, when they cannot be found
#'   from it.
#' @param ... Passed to [adl()] (for example `time`).
#' @return An object of class `"adl"`.
#' @examples
#' if (requireNamespace("dynlm", quietly = TRUE)) {
#'   library(dynlm)
#'   m <- dynlm(y ~ L(y, 1) + L(x, 0:1), data = adl_sim)
#'   fit <- as_adl(m)
#'   summary(fit)
#' }
#' @export
as_adl <- function(x, data = NULL, ...) {
  if (!inherits(x, "dynlm")) stop("`x` must be a dynlm fit.", call. = FALSE)
  f <- stats::formula(x)
  env <- environment(f) %||% parent.frame()
  if (is.null(data)) {
    if (is.null(x$call$data)) {
      stop("The dynlm call has no `data` argument; supply `data`.", call. = FALSE)
    }
    data <- tryCatch(eval(x$call$data, env), error = function(e) NULL)
    if (is.null(data)) {
      stop("Could not find the data of the dynlm call; supply `data`.", call. = FALSE)
    }
  }
  tt <- stats::terms(f)
  intercept <- attr(tt, "intercept") == 1
  labs <- attr(tt, "term.labels")
  is_trend <- grepl("^trend\\(", labs)
  bad <- labs[!is_trend & !grepl("^L\\(", labs) & !vapply(labs, function(l)
    is.symbol(str2lang(l)), logical(1))]
  if (length(bad)) {
    stop("Unsupported dynlm term(s): ", paste(bad, collapse = ", "),
         ". Only variables, L() and trend() can be converted.", call. = FALSE)
  }
  rhs <- labs[!is_trend]
  if (!length(rhs)) stop("The model has no regressors besides the trend.", call. = FALSE)
  newf <- stats::as.formula(paste(deparse(f[[2]]), "~", paste(rhs, collapse = " + ")),
                            env = env)
  start <- if (!is.null(x$call$start)) eval(x$call$start, env)
  end <- if (!is.null(x$call$end)) eval(x$call$end, env)
  adl(newf, data = data, intercept = intercept, trend = any(is_trend),
      start = start, end = end, ...)
}
