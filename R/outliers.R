# Outlier and break detection in the ADL equation ---------------------------

#' Outlier and break detection
#'
#' Searches for impulse and step dummies in the ADL equation, keeps those
#' that are significant, and refits the model with them as unlagged
#' regressors.
#'
#' **What the dummies mean.** A dummy in the ADL equation passes through
#' \eqn{1/C(L)}, like any shock to the equation:
#'
#' * An **impulse** \eqn{I_t(\tau)} (1 at \eqn{\tau}, 0 elsewhere) is an
#'   *innovation outlier*: an unusual shock at \eqn{\tau} whose effect on
#'   \eqn{y} dies out as the dynamic multipliers do.
#' * A **step** \eqn{S_t(\tau)} (0 before \eqn{\tau}, 1 from \eqn{\tau} on)
#'   shifts the intercept: \eqn{y} moves gradually to a new level, and the
#'   long-run shift is the step coefficient divided by \eqn{C(1) = 1 - \sum
#'   \alpha_i}.
#'
#' Outliers defined on \eqn{y} itself look different in the equation. A
#' single wrong observation (an additive outlier, AO) enters both
#' \eqn{y_\tau} and the lag \eqn{y_{\tau-1}} of the next equation: in an
#' ADL(1, q) it shows as an impulse at \eqn{\tau} followed by one of
#' opposite sign and size \eqn{\alpha_1} times as large at \eqn{\tau+1}. An
#' immediate level shift in \eqn{y} (LS) shows as a step plus an impulse at
#' the same date.
#'
#' **Search.** With `method = "search"` (default), a Chen-Liu-type search
#' adapted to OLS: at each step, the \eqn{t}-statistic of every candidate
#' dummy (an impulse at each date, a step at each date after the first and
#' before the last) is computed given the model and the dummies already
#' found; the largest is added if \eqn{|t| > } `cval`, and the search goes on
#' until none is. Then dummies whose \eqn{|t|} has fallen below `cval` are
#' removed one at a time, and forward and backward passes alternate until
#' nothing changes.
#'
#' With `method = "isat"`, the search is impulse and step indicator
#' saturation by [gets::isat()] (Hendry, Johansen and Santos, 2008; Castle,
#' Doornik, Hendry and Pretis, 2015), with the regressors of the ADL kept and
#' the indicators selected at the significance level that corresponds to
#' `cval` (two-sided normal). It needs the `gets` package. `isat()` chooses
#' among its terminal models by an information criterion, so it can keep a
#' dummy whose \eqn{|t|} is slightly below `cval`.
#'
#' **Breaks in x.** To find the dates for [super_exogeneity_test()], run
#' the search on a model for \eqn{x}, e.g. `adl(x ~ y, order = c(1, 1),
#' first_lag = 1)` (\eqn{x} on its own lags and lagged \eqn{y}), and pass
#' `$impulses` and `$steps` on.
#'
#' **Critical value.** By default as in `tsoutliers` (Chen and Liu, 1993;
#' López-de-Lacalle): 3 for \eqn{T \le 50}, 4 for \eqn{T \ge 450}, and
#' \eqn{3 + 0.0025(T - 50)} in between. It grows with \eqn{T} because more
#' dates are searched: with a fixed value, a long sample would produce
#' spurious outliers.
#'
#' Detected dates are not tested hypotheses chosen in advance: the usual
#' p-values of the dummies in the refitted model overstate their
#' significance. And a dummy is a description, not an explanation: each date
#' found deserves a look at what happened then.
#'
#' @param fit An `"adl"` object.
#' @param types `"impulse"`, `"step"` or both (default).
#' @param cval Critical value for \eqn{|t|}; `NULL` uses the rule above.
#' @param method `"search"` (default) or `"isat"`.
#' @param max_iter Maximum number of forward-backward passes (search only).
#' @return An object of class `"adl_outliers"`: a list with `outliers` (data
#'   frame with the date, type, dummy name, coefficient and \eqn{t}-statistic
#'   in the refitted model), `fit` (the ADL refitted with the dummies),
#'   `dummies` (zoo matrix of the dummies on the data index), `impulses` and
#'   `steps` (lists of dates, ready for [super_exogeneity_test()]), `cval`,
#'   `method` and `original` (the model without dummies).
#' @references Chen, C. and Liu, L.-M. (1993). Joint estimation of model
#'   parameters and outlier effects in time series. *Journal of the American
#'   Statistical Association*, 88(421), 284-297.
#'
#'   Hendry, D. F., Johansen, S. and Santos, C. (2008). Automatic selection
#'   of indicators in a fully saturated regression. *Computational
#'   Statistics*, 23(2), 317-335.
#'
#'   Castle, J. L., Doornik, J. A., Hendry, D. F. and Pretis, F. (2015).
#'   Detecting location shifts during model selection by step-indicator
#'   saturation. *Econometrics*, 3(2), 240-264.
#' @seealso [super_exogeneity_test()], [jb_test()], [hansen_test()]
#' @examples
#' fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
#' adl_outliers(fit)
#'
#' # a wrong observation (additive outlier) in 2008 Q4
#' d <- adl_sim
#' d[76, "y"] <- d[76, "y"] - 5
#' out <- adl_outliers(adl(y ~ x, data = d, order = c(1, 1)))
#' out
#' summary(out$fit)
#' @export
adl_outliers <- function(fit, types = c("impulse", "step"), cval = NULL,
                         method = c("search", "isat"), max_iter = 20) {
  check_adl(fit)
  types <- match.arg(types, c("impulse", "step"), several.ok = TRUE)
  method <- match.arg(method)
  check_count(max_iter, "max_iter")
  X <- stats::model.matrix(fit$lm)
  y <- as.numeric(stats::model.response(stats::model.frame(fit$lm)))
  n <- length(y)
  if (is.null(cval)) cval <- outlier_cval(n)
  if (!is.numeric(cval) || length(cval) != 1 || !is.finite(cval) || cval <= 0) {
    stop("`cval` must be a single positive number.", call. = FALSE)
  }
  found <- if (method == "search") {
    outlier_search(y, X, types, cval, max_iter)
  } else {
    outlier_isat(y, X, types, cval)
  }
  build_outlier_result(fit, found, cval, method)
}

#' Critical value of tsoutliers
#' @noRd
outlier_cval <- function(n) {
  if (n <= 50) 3 else if (n >= 450) 4 else round(3 + 0.0025 * (n - 50), 2)
}

#' Candidate dummies: impulses at every date, steps at dates 2..n-1
#' @noRd
outlier_candidates <- function(n, types) {
  out <- list()
  if ("impulse" %in% types) {
    out$impulse <- data.frame(type = "impulse", pos = seq_len(n))
  }
  if ("step" %in% types && n > 2) {
    out$step <- data.frame(type = "step", pos = 2:(n - 1))
  }
  do.call(rbind, unname(out))
}

dummy_col <- function(type, pos, n) {
  if (type == "impulse") as.numeric(seq_len(n) == pos) else as.numeric(seq_len(n) >= pos)
}

#' Forward-backward search on t-statistics
#' @noRd
outlier_search <- function(y, X, types, cval, max_iter) {
  n <- length(y)
  cand <- outlier_candidates(n, types)
  Dall <- vapply(seq_len(nrow(cand)), function(i) dummy_col(cand$type[i], cand$pos[i], n),
                 numeric(n))
  sel <- integer(0)
  tstats <- function(sel) {
    Z <- cbind(X, Dall[, sel, drop = FALSE])
    q <- qr(Z)
    e <- qr.resid(q, y)
    Dt <- qr.resid(q, Dall)                      # candidates net of the model
    dd <- colSums(Dt^2)
    de <- colSums(Dt * e)
    df <- n - q$rank - 1
    rss <- sum(e^2) - de^2 / dd
    t <- (de / dd) / sqrt(pmax(rss, 0) / df / dd)
    t[dd < 1e-8 * max(1, max(dd)) | seq_along(t) %in% sel] <- 0
    t
  }
  current_t <- function(sel) {
    Z <- cbind(X, Dall[, sel, drop = FALSE])
    f <- stats::lm.fit(Z, y)
    if (f$rank < ncol(Z)) return(rep(0, length(sel)))
    s2 <- sum(f$residuals^2) / (n - ncol(Z))
    V <- s2 * chol2inv(qr.R(f$qr))
    piv <- f$qr$pivot
    b <- f$coefficients
    se <- sqrt(diag(V))[order(piv)]
    idx <- ncol(X) + seq_along(sel)
    b[idx] / se[idx]
  }
  for (it in seq_len(max_iter)) {
    changed <- FALSE
    repeat {                                      # forward
      if (length(sel) >= n - ncol(X) - 2) break
      t <- tstats(sel)
      j <- which.max(abs(t))
      if (abs(t[j]) <= cval) break
      sel <- c(sel, j)
      changed <- TRUE
    }
    repeat {                                      # backward
      if (!length(sel)) break
      ct <- current_t(sel)
      j <- which.min(abs(ct))
      if (abs(ct[j]) > cval) break
      sel <- sel[-j]
      changed <- TRUE
    }
    if (!changed) break
  }
  cand[sel, , drop = FALSE]
}

#' Indicator saturation with gets::isat()
#' @noRd
outlier_isat <- function(y, X, types, cval) {
  if (!requireNamespace("gets", quietly = TRUE)) {
    stop("method = \"isat\" needs the 'gets' package: install.packages(\"gets\").",
         call. = FALSE)
  }
  colnames(X) <- paste0("reg", seq_len(ncol(X)))
  res <- gets::isat(y, mc = FALSE, mxreg = X, iis = "impulse" %in% types,
                    sis = "step" %in% types, t.pval = 2 * stats::pnorm(-cval),
                    print.searchinfo = FALSE, plot = FALSE)
  nm <- res$ISnames
  if (is.null(nm) || !length(nm)) {
    return(data.frame(type = character(0), pos = integer(0)))
  }
  data.frame(type = ifelse(startsWith(nm, "iis"), "impulse", "step"),
             pos = as.integer(sub("^(iis|sis)", "", nm)))
}

#' Refit with the dummies and assemble the result
#' @noRd
build_outlier_result <- function(fit, found, cval, method) {
  found <- found[order(found$pos, found$type), , drop = FALSE]
  idx <- fit$index
  dat <- fit$data
  didx <- zoo::index(dat)
  tag <- function(type, pos) {
    paste0(if (type == "impulse") "I_" else "S_",
           gsub("[^0-9A-Za-z]", "", format(idx[pos])))
  }
  names_d <- if (nrow(found)) mapply(tag, found$type, found$pos, USE.NAMES = FALSE) else character(0)
  names_d <- make.unique(c(colnames(dat), names_d), sep = "_")[ncol(dat) + seq_along(names_d)]
  D <- vapply(seq_len(nrow(found)), function(i) {
    v <- as.numeric(as_index_value(idx[found$pos[i]], didx))
    pos <- which(abs(as.numeric(didx) - v) < 1e-8 * max(1, abs(v)))
    if (found$type[i] == "impulse") as.numeric(seq_along(didx) == pos) else
      as.numeric(seq_along(didx) >= pos)
  }, numeric(length(didx)))
  D <- matrix(D, nrow = length(didx), dimnames = list(NULL, names_d))
  dummies <- zoo::zoo(D, didx)

  refit <- fit
  tab <- data.frame(date = character(0), type = character(0), name = character(0),
                    estimate = numeric(0), statistic = numeric(0))
  if (nrow(found)) {
    newdata <- zoo::zoo(cbind(zoo::coredata(dat), D), didx)
    rhs_x <- if (length(fit$x)) paste(fit$x, collapse = " + ") else "1"
    f <- stats::as.formula(paste(fit$y, "~", rhs_x, "|",
                                 paste(c(fit$unlagged, names_d), collapse = " + ")))
    refit <- adl(f, data = newdata, lags = fit$lags, intercept = fit$intercept,
                 trend = fit$trend, max_lag = fit$max_lag)
    if (stats::nobs(refit) != stats::nobs(fit)) {
      stop("Internal error: the refit changed the sample.", call. = FALSE)
    }
    cf <- summary(refit$lm)$coefficients
    tab <- data.frame(date = format(idx[found$pos]), type = found$type,
                      name = names_d, estimate = unname(cf[names_d, 1]),
                      statistic = unname(cf[names_d, 3]))
  }
  dates <- lapply(found$pos, function(p) idx[p])
  structure(list(
    outliers = tab, fit = refit, original = fit, dummies = dummies,
    impulses = dates[found$type == "impulse"], steps = dates[found$type == "step"],
    cval = cval, method = method, nobs = stats::nobs(fit)
  ), class = "adl_outliers")
}

#' @export
print.adl_outliers <- function(x, digits = 4, ...) {
  title <- "Outliers and breaks in the ADL equation"
  cat(title, "\n", strrep("-", nchar(title)), "\n", sep = "")
  cat("Method:      ", if (x$method == "search") "forward-backward search on t" else
    "indicator saturation (gets::isat)", "\n", sep = "")
  cat("Critical |t|: ", format(x$cval), " (T = ", x$nobs, ")\n", sep = "")
  if (!nrow(x$outliers)) {
    cat("\nNo impulse or step dummy exceeds the critical value.\n")
  } else {
    cat("\n")
    tab <- x$outliers
    tab$estimate <- signif(tab$estimate, digits)
    tab$statistic <- signif(tab$statistic, digits)
    print(tab, row.names = FALSE)
    cat("\nModel refitted with the dummies in $fit. The dummies' p-values in it",
        "\noverstate their significance: the dates were searched for.\n")
  }
  invisible(x)
}
