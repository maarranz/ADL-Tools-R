# Multi-path term-by-term Gets through the gets package (tier 2) ------------
#
# The general ADL's design matrix is handed to gets::arx() and searched with
# gets::getsm() (Hoover-Perez / Autometrics-style multi-path search: every
# insignificant term starts a path; diagnostics and backtesting against the
# GUM checked along each path; terminals compared by an information
# criterion). The selected terms are refitted as an adl on the GUM's sample.

#' Multi-path general-to-specific selection with the gets package
#'
#' Term-by-term, multi-path general-to-specific search of an ADL, through
#' [gets::getsm()] (Pretis, Reade and Sucarrat, 2018). Where [adl_gets()]
#' drops whole lags (\eqn{y_{t-p}} and every \eqn{x_{i,t-p}} together), as
#' in ADL3, this search tests single terms, so the selected model may keep,
#' say, \eqn{x_{t-3}} without \eqn{x_{t-2}}. Every insignificant term of the
#' general model starts a deletion path; along each path the residual
#' diagnostics are checked and the simplification is tested against the
#' general model (encompassing); the terminal models of the paths are
#' compared with an information criterion.
#'
#' The general unrestricted model (GUM) is an ADL fit, for instance
#' `adl_gets(...)$general` or `adl(y ~ x, data, order = c(4, 4))`, and it
#' should pass its diagnostics: simplifying a mis-specified model is not a
#' valid search (ADL3). The selected model is refitted with [adl()] on the
#' sample of the GUM (lags with gaps are allowed), so every function of the
#' package can be applied to it.
#'
#' Always kept (never tested): the intercept, unlagged terms and the trend,
#' which are modelling decisions; and, with `contemporaneous = "include"`
#' (default), every \eqn{x_t}, which can be significant only because of
#' simultaneity bias (ADL1); compare with `"search"`, or use
#' [wu_hausman_test()]. `keep` adds other terms.
#'
#' The standard errors of the selected model are conditional on the search:
#' they do not account for the terms having been chosen with the same data.
#'
#' @param fit The general model, an `"adl"` object.
#' @param alpha Significance level of the t-tests of the search (`t.pval`
#'   of [gets::getsm()]), also used for the test against the GUM.
#' @param contemporaneous `"include"` (default) keeps every \eqn{x_t};
#'   `"search"` lets the search drop them.
#' @param keep Further coefficient names (columns of the fit, e.g.
#'   `"y_lag4"`) that the search may not drop.
#' @param diagnostics Check the residual diagnostics of [gets::getsm()]
#'   along each path (Ljung-Box tests for autocorrelation and ARCH).
#' @param diag_alpha Significance level of those diagnostics.
#' @param ic Information criterion that chooses among terminal models:
#'   `"bic"` (default), `"aic"` or `"hq"`.
#' @inheritParams adl_vcov
#' @param ... Further arguments to [gets::getsm()].
#' @return An object of class `"adl_getsm"`: a list with `fit` (the selected
#'   model as an `"adl"` fit, or `NULL` if the search removed every
#'   regressor), `general` (the GUM), `selected` and `dropped` (coefficient
#'   names), `kept` (the terms never tested), `terminals` (the terminal
#'   models and their criteria), `paths` (the deletion paths, as column
#'   names), `gets` (the [gets::getsm()] object) and the settings.
#' @references Pretis, F., Reade, J. J. and Sucarrat, G. (2018). Automated
#'   general-to-specific (GETS) regression modeling and indicator saturation
#'   for outliers and structural breaks. *Journal of Statistical Software*,
#'   86(3), 1-44.
#'
#'   Hoover, K. D. and Perez, S. J. (1999). Data mining reconsidered:
#'   encompassing and the general-to-specific approach to specification
#'   search. *Econometrics Journal*, 2(2), 167-191.
#' @seealso [adl_gets()] (lag by lag, as in ADL3), [adl_ic()].
#' @examples
#' if (requireNamespace("gets", quietly = TRUE)) {
#'   gum <- adl(y ~ x + w, data = adl_sim, order = c(4, 4, 4))
#'   sel <- adl_getsm(gum)
#'   sel
#'   summary(sel$fit)
#' }
#' @export
adl_getsm <- function(fit, alpha = 0.05,
                      contemporaneous = c("include", "search"),
                      keep = NULL, diagnostics = TRUE, diag_alpha = 0.025,
                      ic = c("bic", "aic", "hq"),
                      cov_type = c("nonrobust", "HAC"), ...) {
  check_adl(fit)
  if (!requireNamespace("gets", quietly = TRUE)) {
    stop("adl_getsm() needs the 'gets' package: install.packages(\"gets\").",
         call. = FALSE)
  }
  check_alpha(alpha)
  check_alpha(diag_alpha)
  contemporaneous <- match.arg(contemporaneous)
  ic <- match.arg(ic)
  cov_type <- match.arg(cov_type)
  if (!is.logical(diagnostics) || length(diagnostics) != 1 || is.na(diagnostics)) {
    stop("`diagnostics` must be TRUE or FALSE.", call. = FALSE)
  }

  X <- stats::model.matrix(fit$lm)
  y <- as.numeric(stats::model.response(stats::model.frame(fit$lm)))
  b <- stats::coef(fit)
  info <- term_info(fit, colnames(X))
  bad <- setdiff(keep, colnames(X))
  if (length(bad)) {
    stop("`keep` names terms that are not in the model: ",
         paste(bad, collapse = ", "), ".", call. = FALSE)
  }
  always <- info$role %in% c("intercept", "unlagged", "trend") |
    (contemporaneous == "include" & info$role == "dl" & info$lag == 0) |
    colnames(X) %in% keep
  mx <- X[, info$role != "intercept", drop = FALSE]
  # positions in gets' mean equation: mconst first, then the mxreg columns
  mx_names <- colnames(mx)
  pos <- if (fit$intercept) c("(Intercept)", mx_names) else mx_names
  keep_idx <- which(pos %in% colnames(X)[always])

  vt <- if (cov_type == "HAC") "newey-west" else "ordinary"
  a <- gets::arx(y, mc = fit$intercept, mxreg = if (ncol(mx)) mx, vcov.type = vt,
                 plot = FALSE)
  diag_list <- if (diagnostics) list(lag = NULL, pval = diag_alpha) else NULL
  g <- gets::getsm(a, t.pval = alpha, keep = if (length(keep_idx)) keep_idx,
                   ar.LjungB = diag_list, arch.LjungB = diag_list,
                   info.method = c(bic = "sc", aic = "aic", hq = "hq")[[ic]],
                   print.searchinfo = FALSE, plot = FALSE, ...)
  if (is.null(g$specific.spec)) {
    stop("gets::getsm() did not return a specific model (the general model ",
         "fails its diagnostics?). Check adl_diagnostics() on the general ",
         "model, or set diagnostics = FALSE.", call. = FALSE)
  }
  selected <- pos[g$specific.spec]
  dropped <- setdiff(colnames(X), selected)

  # refit the selection as an adl on the GUM's sample
  sel_info <- info[match(selected, colnames(X)), ]
  lagsel <- function(v, role) sort(sel_info$lag[sel_info$variable == v & sel_info$role == role])
  xs <- fit$x[vapply(fit$x, function(v) length(lagsel(v, "dl")) > 0, logical(1))]
  unl <- fit$unlagged[fit$unlagged %in% sel_info$variable[sel_info$role == "unlagged"]]
  new_fit <- NULL
  if (length(xs)) {
    lags <- c(stats::setNames(list(lagsel(fit$y, "ar")), fit$y),
              stats::setNames(lapply(xs, lagsel, role = "dl"), xs))
    f <- adl_formula(fit$y, xs, unl)
    new_fit <- suppressMessages(adl(f, data = fit$data, lags = lags,
                                    intercept = fit$intercept,
                                    trend = fit$trend && "trend" %in% selected,
                                    max_lag = fit$max_lag))
  }

  paths <- lapply(g$paths, function(p) pos[p])
  structure(list(
    fit = new_fit,
    general = fit,
    selected = selected,
    dropped = dropped,
    kept = colnames(X)[always],
    terminals = g$terminals.results,
    paths = paths,
    gets = g,
    alpha = alpha, diag_alpha = diag_alpha, diagnostics = diagnostics,
    ic = ic, cov_type = cov_type, contemporaneous = contemporaneous
  ), class = "adl_getsm")
}

#' Print a multi-path Gets selection
#'
#' @param x An `"adl_getsm"` object.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.adl_getsm <- function(x, ...) {
  title <- "Multi-path general-to-specific search (gets::getsm)"
  cat(title, "\n", strrep("-", nchar(title)), "\n", sep = "")
  cat("General model: ", x$general$label, " on ", sample_text(x$general),
      " (", length(stats::coef(x$general)), " terms)\n", sep = "")
  cat("Search:        t-tests at ", format(100 * x$alpha), "%, ",
      if (x$cov_type == "HAC") "HAC" else "conventional", " covariance; ",
      if (x$diagnostics) paste0("diagnostics at ", format(100 * x$diag_alpha), "%")
      else "no diagnostics", "; terminals by ", toupper(x$ic), "\n", sep = "")
  cat("Never tested:  ", if (length(x$kept)) paste(x$kept, collapse = ", ") else "-",
      "\n", sep = "")
  cat("Paths:         ", length(x$paths), "; terminal models: ",
      max(1L, NROW(x$terminals)), "\n\n", sep = "")
  cat("Selected:      ", paste(x$selected, collapse = ", "), "\n", sep = "")
  cat("Dropped:       ", if (length(x$dropped)) paste(x$dropped, collapse = ", ")
      else "-", "\n", sep = "")
  if (is.null(x$fit)) {
    cat("\nThe search removed every regressor: no ADL is left (see $gets).\n")
  } else {
    cat("Model:         ", x$fit$label, "\n", sep = "")
  }
  cat("Standard errors of the selected model are conditional on the search.\n")
  invisible(x)
}
