# Information criteria over a grid of lag orders ----------------------------

#' Information criteria for ADL models on a common sample
#'
#' Estimates every ADL model in a grid of lag orders on the same sample and
#' reports AIC, BIC (SC) and Hannan-Quinn, so the orders can be compared on
#' equal terms. The table is a guide for choosing orders; it does not select
#' a model by itself.
#'
#' The grid covers \eqn{p = 0, \dots, p_{max}} lags of the dependent
#' variable and, for each regressor, lags from `first_lag` up to
#' \eqn{q = 0, \dots, q_{max}}. Following `DESIGN.md` (section 3.7), the
#' contemporaneous term is never dropped by the search alone: with
#' `contemporaneous = "both"` (the default) every order is estimated with
#' and without \eqn{x_t}, and both versions are listed side by side.
#'
#' The criteria are the per-observation forms reported by EViews,
#' \eqn{-2\ell/n + k\,c(n)/n} with \eqn{c(n) = 2}, \eqn{\log n} and
#' \eqn{2 \log\log n}, where \eqn{\ell} is the Gaussian log-likelihood and
#' \eqn{k} the number of coefficients (as in EViews, \eqn{\sigma^2} is
#' not counted; R's `AIC()` and `BIC()` count it, which shifts every model
#' by the same amount and leaves the ranking unchanged).
#'
#' @param formula A formula with the variables only, e.g. `y ~ x1 + x2 | d`
#'   (no `L()` terms).
#' @inheritParams adl
#' @param max_order `c(p_max, q_max)`, with the same `q_max` for every
#'   regressor, or `c(p_max, q_max_1, ..., q_max_n)`.
#' @param contemporaneous `"both"`, `"include"` or `"exclude"`: whether the
#'   regressors enter at lag 0. One value for all regressors or one per
#'   regressor (named or in order).
#' @param same_q If `TRUE` (default with several regressors) all regressors
#'   share the same \eqn{q} in each model; if `FALSE` the full grid over
#'   every regressor's \eqn{q} is estimated.
#' @return A data frame of class `"adl_ic"` with one row per model: `p`, a
#'   `q_` and an `x_t_` column per regressor (lags and whether lag 0 is
#'   included), `k`, `aic`, `bic`, `hq`, and `best_aic`, `best_bic`,
#'   `best_hq` flags. Attribute `nobs` holds the common sample size.
#' @examples
#' set.seed(1)
#' n <- 200
#' x <- cumsum(rnorm(n))
#' y <- stats::filter(0.5 * x + rnorm(n), 0.6, method = "recursive")
#' d <- ts(cbind(y = as.numeric(y), x = x))
#' adl_ic(y ~ x, data = d, max_order = c(3, 3))
#' @export
adl_ic <- function(formula, data, max_order, contemporaneous = "both",
                   same_q = TRUE, intercept = TRUE, trend = FALSE, time = NULL,
                   start = NULL, end = NULL) {
  spec <- parse_adl_formula(formula)
  if (!is.null(spec$lags)) {
    stop("Write the formula without L() terms; the grid sets the lags.",
         call. = FALSE)
  }
  x <- spec$x
  n_x <- length(x)
  check_whole(max_order, "max_order")
  if (any(max_order < 0)) stop("`max_order` must be non-negative.", call. = FALSE)
  if (length(max_order) == 2) {
    qmax <- rep(max_order[2], n_x)
  } else if (length(max_order) == n_x + 1) {
    qmax <- max_order[-1]
  } else {
    stop("`max_order` must have length 2 or ", n_x + 1, ".", call. = FALSE)
  }
  pmax <- max_order[1]
  cont <- recycle_named_chr(contemporaneous, x)
  if (!all(cont %in% c("both", "include", "exclude"))) {
    stop("`contemporaneous` must be \"both\", \"include\" or \"exclude\".",
         call. = FALSE)
  }
  if (same_q && length(unique(qmax)) > 1) {
    stop("With `same_q = TRUE` all regressors need the same q_max.", call. = FALSE)
  }

  # first-lag options per regressor
  fl_opts <- lapply(cont, function(c) switch(c, both = c(0L, 1L), include = 0L,
                                             exclude = 1L))
  # grid of (p, q_1..q_n, fl_1..fl_n)
  if (same_q) {
    qgrid <- data.frame(q = 0:qmax[1])
    qgrid <- qgrid[rep(seq_len(nrow(qgrid)), 1), , drop = FALSE]
    qgrid <- as.data.frame(matrix(rep(qgrid$q, n_x), ncol = n_x))
  } else {
    qgrid <- expand.grid(lapply(qmax, function(m) 0:m), KEEP.OUT.ATTRS = FALSE)
  }
  names(qgrid) <- paste0("q_", x)
  flgrid <- expand.grid(fl_opts, KEEP.OUT.ATTRS = FALSE)
  names(flgrid) <- paste0("fl_", x)
  grid <- merge(merge(data.frame(p = 0:pmax), qgrid, by = NULL), flgrid, by = NULL)
  # a regressor without lag 0 needs q >= 1
  ok <- rep(TRUE, nrow(grid))
  for (v in x) ok <- ok & !(grid[[paste0("fl_", v)]] == 1 & grid[[paste0("q_", v)]] < 1)
  grid <- grid[ok, , drop = FALSE]

  m <- max(pmax, qmax)
  d <- as_adl_data(data, c(spec$y, x, spec$unlagged), time = time,
                   start = start, end = end, pre = m)
  rows <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    lags <- c(list(seq_len(g$p)),
              lapply(x, function(v) g[[paste0("fl_", v)]]:g[[paste0("q_", v)]]))
    names(lags) <- c(spec$y, x)
    lags <- normalize_lags(spec$y, x, lags = lags)
    des <- build_design(d$z, lags, unlagged = spec$unlagged, trend = trend,
                        max_lag = m)
    rhs <- setdiff(names(des$frame), des$y_col)
    f <- stats::reformulate(if (length(rhs)) rhs else "1", response = des$y_col,
                            intercept = intercept)
    fit <- stats::lm(f, data = des$frame)
    n <- stats::nobs(fit)
    k <- length(stats::coef(fit))
    ll <- as.numeric(stats::logLik(fit))
    rows[[i]] <- data.frame(k = k, aic = (-2 * ll + 2 * k) / n,
                            bic = (-2 * ll + log(n) * k) / n,
                            hq = (-2 * ll + 2 * log(log(n)) * k) / n)
  }
  res <- cbind(grid, do.call(rbind, rows))
  for (v in x) {
    res[[paste0("x_t_", v)]] <- res[[paste0("fl_", v)]] == 0
    res[[paste0("fl_", v)]] <- NULL
  }
  res <- res[c("p", paste0("q_", x), paste0("x_t_", x), "k", "aic", "bic", "hq")]
  res$best_aic <- res$aic == min(res$aic)
  res$best_bic <- res$bic == min(res$bic)
  res$best_hq <- res$hq == min(res$hq)
  rownames(res) <- NULL
  structure(res, class = c("adl_ic", "data.frame"), nobs = n, y = spec$y,
            contemporaneous = cont)
}

recycle_named_chr <- function(v, x) {
  if (length(v) == 1 && is.null(names(v))) return(stats::setNames(rep(v, length(x)), x))
  if (!is.null(names(v))) {
    unknown <- setdiff(names(v), x)
    if (length(unknown)) stop("`contemporaneous` names unknown regressor(s): ",
                              paste(unknown, collapse = ", "), call. = FALSE)
    out <- stats::setNames(rep("both", length(x)), x)
    out[names(v)] <- v
    return(out)
  }
  if (length(v) != length(x)) stop("`contemporaneous` must have length 1 or ",
                                   length(x), ".", call. = FALSE)
  stats::setNames(v, x)
}

#' Print an information-criteria table
#'
#' Prints the models ordered by the chosen criterion, marking the best
#' model under each criterion.
#'
#' @param x An `"adl_ic"` object.
#' @param sort_by `"bic"` (default), `"aic"` or `"hq"`.
#' @param n Number of models shown; `NULL` shows them all.
#' @param digits Significant digits.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.adl_ic <- function(x, sort_by = c("bic", "aic", "hq"), n = 15, digits = 5, ...) {
  sort_by <- match.arg(sort_by)
  title <- "Information criteria (per observation; smaller is better)"
  cat(title, "\n", strrep("-", nchar(title)), "\n", sep = "")
  cat("Dependent variable: ", attr(x, "y"), ";  common sample of ", attr(x, "nobs"),
      " observations;  ", nrow(x), " models\n\n", sep = "")
  d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  d <- d[order(d[[sort_by]]), , drop = FALSE]
  mark <- function(b) ifelse(b, "*", "")
  d$aic <- paste0(formatC(d$aic, digits = digits, format = "fg"), mark(d$best_aic))
  d$bic <- paste0(formatC(d$bic, digits = digits, format = "fg"), mark(d$best_bic))
  d$hq <- paste0(formatC(d$hq, digits = digits, format = "fg"), mark(d$best_hq))
  d$best_aic <- d$best_bic <- d$best_hq <- NULL
  total <- nrow(d)
  if (!is.null(n) && total > n) d <- d[seq_len(n), , drop = FALSE]
  print(d, row.names = FALSE)
  if (nrow(d) < total) cat("... ", total - nrow(d), " more models\n", sep = "")
  cat("* best model under each criterion. Sorted by ", toupper(sort_by), ".\n", sep = "")
  invisible(x)
}
