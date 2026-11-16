# Test results ----------------------------------------------------------------
#
# Every test returns a list of class "adltools_test" with the same fields as
# forecastingtools ("forecastingtools_test"), so results from the suite look
# and extract alike: method, null_hypothesis, statistic, distribution
# ("t", "normal", "F", "chi2"), p.value, alpha, reject_null, nobs, cov_type,
# maxlags, coefficients (data frame or NULL), n_restrictions, df_resid, df,
# df_num, df_denom. Tests can add their own fields.

#' Build a test result
#' @noRd
new_adl_test <- function(x, title) {
  defaults <- list(alternative = NULL, cov_type = NULL, maxlags = NA_integer_,
                   coefficients = NULL, n_restrictions = NA_integer_,
                   df_resid = NA_integer_, df = NA_integer_,
                   df_num = NA_integer_, df_denom = NA_integer_)
  for (nm in names(defaults)) if (!nm %in% names(x)) x[nm] <- list(defaults[[nm]])
  structure(x, class = "adltools_test", title = title)
}

#' p-value and degrees of freedom for a Wald statistic
#'
#' `dist = "chisq"` uses the Wald statistic with a chi-square(q) reference;
#' `dist = "F"` divides it by q and uses F(q, df_resid).
#' @param w Wald statistic (chi-square form) for `q` restrictions.
#' @noRd
wald_inference <- function(w, q, df_resid, dist = "chisq") {
  if (dist == "F") {
    stat <- w / q
    list(statistic = stat, distribution = "F", df_num = q, df_denom = df_resid,
         df = NA_integer_,
         p.value = stats::pf(stat, q, df_resid, lower.tail = FALSE))
  } else {
    list(statistic = w, distribution = "chi2", df = q, df_num = NA_integer_,
         df_denom = NA_integer_,
         p.value = stats::pchisq(w, q, lower.tail = FALSE))
  }
}

#' Wald quadratic form for H0: R b = r
#' @noRd
wald_stat <- function(b, V, R, r) {
  d <- R %*% b - r
  M <- R %*% V %*% t(R)
  if (anyNA(M) || qr(M)$rank < nrow(M)) {
    stop("The covariance of the restrictions is singular: the restrictions ",
         "are not identified (collinear regressors, or too few degrees of ",
         "freedom).", call. = FALSE)
  }
  as.numeric(t(d) %*% solve(M, d))
}

fmt_num <- function(x, digits) {
  if (is.null(x) || (length(x) == 1L && is.na(x))) return("-")
  if (is.numeric(x) && !is.integer(x)) format(signif(x, digits)) else format(x)
}

#' Print a test result
#'
#' Prints the hypothesis, statistic, p-value and decision, followed by the
#' coefficient table when there is one. Every element can also be extracted
#' directly, e.g. `result$p.value`.
#'
#' @param x A test result from adltools.
#' @param digits Significant digits.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.adltools_test <- function(x, digits = 4, ...) {
  title <- attr(x, "title")
  cat(title, "\n", strrep("-", nchar(title)), "\n", sep = "")
  dist <- x$distribution
  df <- if (identical(dist, "F")) paste0("(", x$df_num, ", ", x$df_denom, ")")
        else if (!is.null(x$df) && !is.na(x$df)) paste0("(", x$df, ")") else ""
  label <- switch(dist, normal = "z", t = "t", chi2 = "Chi-square", F = "F", dist)
  cat("H0:          ", x$null_hypothesis, "\n", sep = "")
  if (!is.null(x$alternative)) cat("Alternative: ", x$alternative, "\n", sep = "")
  cat("Statistic:   ", fmt_num(x$statistic, digits), "  [", label, df, "]\n", sep = "")
  cat("p-value:     ", fmt_num(x$p.value, digits), "\n", sep = "")
  cat("Decision:    ", if (isTRUE(x$reject_null)) "reject H0" else "do not reject H0",
      " at ", format(100 * x$alpha), "%\n", sep = "")
  if (!is.null(x$cov_type)) {
    cov <- if (x$cov_type == "HAC") paste0("HAC (Bartlett, ", x$maxlags, " lags)")
           else "conventional OLS"
    cat("Covariance:  ", cov, "\n", sep = "")
  }
  cat("N:           ", x$nobs, "\n", sep = "")
  if (is.data.frame(x$coefficients) && nrow(x$coefficients)) {
    cat("\n")
    tab <- x$coefficients
    num <- vapply(tab, is.numeric, logical(1))
    tab[num] <- lapply(tab[num], signif, digits = digits)
    print(tab, row.names = FALSE)
  }
  invisible(x)
}
