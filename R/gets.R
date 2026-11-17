# General-to-specific lag selection (ADL3) -----------------------------------

#' General-to-specific (Gets) lag selection
#'
#' The general-to-specific procedure of ADL3. Start from a deliberately
#' generous ADL\eqn{(p, p)} (the same \eqn{p} for \eqn{y} and every
#' regressor) and check that it passes the diagnostics: there is no point in
#' simplifying a model that is mis-specified to start with. Then, at each
#' step, test whether the highest-order lags, \eqn{y_{t-p}} and every
#' \eqn{x_{i,t-p}}, are jointly zero (a Wald test). If not rejected, drop to
#' ADL\eqn{(p-1, p-1)} and repeat; stop at the first rejection.
#'
#' Every model in the search is estimated on the same sample (the one of the
#' general model), so the steps are nested tests. At every step the
#' diagnostics are re-run on the simplified model; if it fails one that the
#' current model passes, the simplification is not accepted and the search
#' stops there ("back up", ADL3).
#'
#' The contemporaneous \eqn{x_t} is never dropped by a significance test:
#' with simultaneity \eqn{\hat\beta_0} can be significant because of the bias
#' (Chapter 1). Whether it belongs is a modelling decision, set with
#' `contemporaneous`: `"include"` (default) keeps it in every model,
#' `"exclude"` searches among dead-start models, and `"both"` runs both
#' searches so they can be compared (with [wu_hausman_test()], for instance).
#'
#' The standard errors of the selected model are conditional on the search:
#' they do not account for the lags having been chosen with the same data.
#' Information criteria over a grid of orders are in [adl_ic()].
#'
#' @param formula A formula with the variables only, e.g. `y ~ x1 + x2 | d`
#'   (unlagged regressors after `|` are kept in every model).
#' @param data The data (any type accepted by [adl()]).
#' @param max_order The lag \eqn{p} of the general model. By default 4 for
#'   quarterly, 12 for monthly and 2 for other data (the rule of thumb in
#'   ADL3).
#' @param contemporaneous `"include"` (default), `"exclude"` or `"both"`;
#'   see Details.
#' @param diagnostics Diagnostics checked at every step: any of `"bg"`
#'   (Breusch-Godfrey), `"arch"`, `"reset"` and `"jb"` (Jarque-Bera).
#'   Normality is not required for consistency, so `"jb"` is not in the
#'   default. `character(0)` skips the checks.
#' @param bg_order Order of the Breusch-Godfrey test; by default `max_order`.
#' @param alpha Significance level of the Wald tests of the search.
#' @param diag_alpha Significance level of the diagnostics.
#' @inheritParams adl_vcov
#' @param dist `"chisq"` (default) or `"F"`, for the Wald tests.
#' @param intercept,trend,time,start,end As in [adl()]. `start` is the first
#'   observation of the common sample.
#' @return An object of class `"adl_gets"`: a list with `fit` (the selected
#'   model), `general` (the general model), `path` (a data frame with one row
#'   per step: the model, the lag tested, the Wald test that all terms at that
#'   lag are zero, the action (`drop`, `stop`, or `back up` with the failing
#'   diagnostic) and the diagnostic p-values of the model without them), `stop_reason`, `general_ok` (whether
#'   the general model passes the diagnostics) and the settings. With
#'   `contemporaneous = "both"`, a list of two such objects (`include`,
#'   `exclude`) of class `"adl_gets_both"`.
#' @references Hendry, D. F. (1995). *Dynamic Econometrics*. Oxford
#'   University Press.
#' @seealso [adl_ic()], [adl_diagnostics()], [adl_joint_test()]
#' @examples
#' g <- adl_gets(y ~ x, data = adl_sim)
#' g
#' summary(g$fit)
#' adl_gets(y ~ x + w, data = adl_sim, contemporaneous = "both")
#' @export
adl_gets <- function(formula, data, max_order = NULL,
                     contemporaneous = c("include", "exclude", "both"),
                     diagnostics = c("bg", "arch", "reset"), bg_order = NULL,
                     alpha = 0.05, diag_alpha = 0.05,
                     cov_type = c("nonrobust", "HAC"), maxlags = NULL,
                     dist = c("chisq", "F"), intercept = TRUE, trend = FALSE,
                     time = NULL, start = NULL, end = NULL) {
  contemporaneous <- match.arg(contemporaneous)
  cov_type <- match.arg(cov_type)
  dist <- match.arg(dist)
  diagnostics <- if (length(diagnostics)) {
    match.arg(diagnostics, c("bg", "arch", "reset", "jb"), several.ok = TRUE)
  } else character(0)
  check_alpha(alpha)
  check_alpha(diag_alpha)
  spec <- parse_adl_formula(formula)
  if (!is.null(spec$lags)) {
    stop("Write the formula without L() terms; the search sets the lags.",
         call. = FALSE)
  }
  if (!length(spec$x)) stop("The formula needs at least one regressor.", call. = FALSE)
  if (is.null(max_order)) {
    z <- to_zoo(data, time = time)
    freq <- tryCatch(stats::frequency(z), error = function(e) 1)
    max_order <- if (freq == 4) 4L else if (freq == 12) 12L else 2L
  }
  check_whole(max_order, "max_order")
  if (length(max_order) != 1 || max_order < 1) {
    stop("`max_order` must be a single integer >= 1.", call. = FALSE)
  }
  if (is.null(bg_order)) bg_order <- max_order
  args <- list(formula = formula, data = data, max_order = max_order,
               diagnostics = diagnostics, bg_order = bg_order, alpha = alpha,
               diag_alpha = diag_alpha, cov_type = cov_type, maxlags = maxlags,
               dist = dist, intercept = intercept, trend = trend, time = time,
               start = start, end = end)
  if (contemporaneous == "both") {
    return(structure(list(
      include = do.call(gets_search, c(args, list(first_lag = 0L))),
      exclude = do.call(gets_search, c(args, list(first_lag = 1L)))
    ), class = "adl_gets_both"))
  }
  do.call(gets_search, c(args, list(first_lag = if (contemporaneous == "include") 0L else 1L)))
}

#' One Gets search for a given first lag of the regressors
#' @noRd
gets_search <- function(formula, data, max_order, first_lag, diagnostics,
                        bg_order, alpha, diag_alpha, cov_type, maxlags, dist,
                        intercept, trend, time, start, end) {
  spec <- parse_adl_formula(formula)
  m <- max_order
  fit_p <- function(p) {
    lags <- c(list(seq_len(p)), rep(list(first_lag:p), length(spec$x)))
    names(lags) <- c(spec$y, spec$x)
    adl(formula, data = data, lags = lags, intercept = intercept, trend = trend,
        time = time, start = start, end = end, max_lag = m)
  }
  diag_p <- function(fit) {
    if (!length(diagnostics)) return(stats::setNames(numeric(0), character(0)))
    out <- c(
      bg = if ("bg" %in% diagnostics) bg_test(fit, order = bg_order)$p.value,
      arch = if ("arch" %in% diagnostics) arch_test(fit)$p.value,
      reset = if ("reset" %in% diagnostics)
        tryCatch(reset_test(fit)$p.value, error = function(e) NA_real_),
      jb = if ("jb" %in% diagnostics) jb_test(fit)$p.value
    )
    out[diagnostics]
  }
  passes <- function(pv) all(is.na(pv) | pv >= diag_alpha)

  min_p <- if (first_lag == 0L) 0L else 1L
  general <- fit_p(m)
  gen_diag <- diag_p(general)
  general_ok <- passes(gen_diag)
  if (!general_ok) {
    bad <- names(gen_diag)[!is.na(gen_diag) & gen_diag < diag_alpha]
    warning("The general ADL(", m, ", ", m, ") fails ",
            paste(toupper(bad), collapse = ", "), " at ", 100 * diag_alpha,
            "%: simplifying a mis-specified model is not a valid search (ADL3). ",
            "Consider a larger `max_order`, other regressors or dummies.",
            call. = FALSE)
  }

  rows <- list()
  current <- general
  cur_diag <- gen_diag
  p <- m
  stop_reason <- NULL
  while (is.null(stop_reason)) {
    if (p <= min_p) {
      stop_reason <- paste0("smallest model reached (",
                            if (first_lag == 0L) "static regression" else
                              "one lag of each regressor", ")")
      break
    }
    info <- term_info(current, names(stats::coef(current)))
    top <- names(stats::coef(current))[!is.na(info$lag) & info$lag == p &
                                         info$role %in% c("ar", "dl")]
    w <- adl_joint_test(current, terms = top, cov_type = cov_type,
                        maxlags = maxlags, alpha = alpha, dist = dist)
    row <- data.frame(model = current$label, lag = p,
                      statistic = w$statistic, df = length(top), p.value = w$p.value)
    if (w$p.value < alpha) {
      row$action <- "stop"
      rows[[length(rows) + 1]] <- add_diag(row, stats::setNames(
        rep(NA_real_, length(diagnostics)), diagnostics))
      stop_reason <- paste0("the lags at ", p, " are jointly significant (p = ",
                            format(signif(w$p.value, 3)), ")")
      break
    }
    simpler <- fit_p(p - 1)
    sd <- diag_p(simpler)
    newly_bad <- !is.na(sd) & sd < diag_alpha & (is.na(cur_diag) | cur_diag >= diag_alpha)
    if (any(newly_bad)) {
      row$action <- paste0("back up (", paste(toupper(names(sd)[newly_bad]),
                                              collapse = ", "), ")")
      rows[[length(rows) + 1]] <- add_diag(row, sd)
      stop_reason <- paste0("dropping the lags at ", p, " would fail ",
                            paste(toupper(names(sd)[newly_bad]), collapse = ", "))
      break
    }
    row$action <- "drop"
    rows[[length(rows) + 1]] <- add_diag(row, sd)
    current <- simpler
    cur_diag <- sd
    p <- p - 1
  }
  path <- if (length(rows)) do.call(rbind, rows) else NULL
  structure(list(
    fit = current, general = general, path = path, stop_reason = stop_reason,
    general_ok = general_ok, general_diagnostics = gen_diag,
    final_diagnostics = cur_diag, max_order = m,
    contemporaneous = if (first_lag == 0L) "include" else "exclude",
    alpha = alpha, diag_alpha = diag_alpha, cov_type = cov_type, dist = dist
  ), class = "adl_gets")
}

#' @export
print.adl_gets <- function(x, digits = 4, ...) {
  title <- paste0("General-to-specific search (x_t ",
                  if (x$contemporaneous == "include") "included" else "excluded", ")")
  cat(title, "\n", strrep("-", nchar(title)), "\n", sep = "")
  cat("General model: ", x$general$label, " on ", sample_text(x$general), "\n", sep = "")
  if (length(x$general_diagnostics)) {
    cat("Diagnostics:   ", gets_diag_text(x$general_diagnostics, digits),
        if (x$general_ok) "" else "  <- FAILS", "\n", sep = "")
  }
  cat("Wald tests:    ", if (x$dist == "F") "F" else "chi-square", ", ",
      if (x$cov_type == "HAC") "HAC" else "conventional", " covariance, ",
      format(100 * x$alpha), "%; diagnostics at ", format(100 * x$diag_alpha), "%\n\n",
      sep = "")
  if (!is.null(x$path)) {
    tab <- x$path
    num <- vapply(tab, is.numeric, logical(1))
    tab[num] <- lapply(tab[num], signif, digits = digits)
    print(tab, row.names = FALSE)
    cat("\n")
  }
  cat("Selected:      ", x$fit$label, "\n", sep = "")
  cat("Stopped:       ", x$stop_reason, "\n", sep = "")
  cat("Standard errors of the selected model are conditional on the search.\n")
  invisible(x)
}

#' @export
print.adl_gets_both <- function(x, ...) {
  print(x$include, ...)
  cat("\n")
  print(x$exclude, ...)
  cat("\nx_t is not chosen by significance: compare the two, e.g. with",
      "wu_hausman_test().\n")
  invisible(x)
}

#' Diagnostic p-values as columns of a path row
#' @noRd
add_diag <- function(row, pv) {
  for (nm in names(pv)) row[[paste0(nm, "_p")]] <- unname(pv[nm])
  row
}

gets_diag_text <- function(pv, digits) {
  paste(paste0(toupper(names(pv)), " p = ", format(signif(pv, 3))), collapse = ", ")
}
