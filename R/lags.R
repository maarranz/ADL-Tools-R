# Lag sets and design matrix ------------------------------------------------
#
# A model is described internally by a named list of integer lag sets, one
# element per variable: the dependent variable (lags >= 1, possibly empty)
# and each regressor (lags >= 0, at least one). The three entry routes in
# DESIGN.md section 3.3 (order + first_lag, explicit lags, formula) all end
# here.

#' Normalize the lag specification to a named list of lag sets
#'
#' @param y Name of the dependent variable.
#' @param x Names of the regressors.
#' @param order Either `c(p, q)` (q recycled for every regressor) or
#'   `c(p, q_1, ..., q_n)`.
#' @param first_lag First lag of each regressor: 0 includes x_t, 1 excludes
#'   it. A single value, or one per regressor (named or in order).
#' @param lags Explicit named list of lag sets; overrides `order`.
#' @return Named list of sorted unique integer vectors, `y` first.
#' @noRd
normalize_lags <- function(y, x, order = NULL, first_lag = 0, lags = NULL) {
  if (length(x) == 0) stop("At least one regressor is needed.", call. = FALSE)
  if (anyDuplicated(c(y, x))) stop("Variable names must be distinct.", call. = FALSE)

  if (!is.null(lags)) {
    if (!is.null(order)) {
      warning("Both `order` and `lags` supplied; `order` is ignored.", call. = FALSE)
    }
    if (!is.list(lags) || is.null(names(lags)) || any(names(lags) == "")) {
      stop("`lags` must be a named list, e.g. list(y = 1:2, x = 0:1).", call. = FALSE)
    }
    unknown <- setdiff(names(lags), c(y, x))
    if (length(unknown) > 0) {
      stop("`lags` refers to unknown variable(s): ",
           paste(unknown, collapse = ", "), call. = FALSE)
    }
    no_lags <- setdiff(x, names(lags))
    if (length(no_lags) > 0) {
      stop("No lags given for regressor(s): ", paste(no_lags, collapse = ", "),
           call. = FALSE)
    }
    out <- c(list(lags[[y]] %||% integer(0)), lags[x])
    names(out) <- c(y, x)
  } else {
    if (is.null(order)) stop("Supply either `order` or `lags`.", call. = FALSE)
    check_whole(order, "order")
    n <- length(x)
    if (length(order) == 2) {
      order <- c(order[1], rep(order[2], n))
    } else if (length(order) != n + 1) {
      stop("`order` must have length 2 or ", n + 1,
           " (p and one q per regressor).", call. = FALSE)
    }
    fl <- recycle_named(first_lag, x, "first_lag")
    check_whole(fl, "first_lag")
    if (any(!fl %in% c(0, 1))) stop("`first_lag` must be 0 or 1.", call. = FALSE)
    q <- order[-1]
    if (any(q < fl)) {
      bad <- x[q < fl]
      stop("Maximum lag below `first_lag` for: ", paste(bad, collapse = ", "),
           ". With first_lag = 1 the order must be at least 1.", call. = FALSE)
    }
    out <- c(list(seq_len(order[1])), Map(seq, fl, q))
    names(out) <- c(y, x)
  }

  out <- lapply(out, function(l) sort(unique(as.integer(l))))
  for (v in names(out)) check_whole(out[[v]], paste0("lags for '", v, "'"))
  if (any(out[[y]] < 1)) {
    stop("Lags of the dependent variable must be 1 or more.", call. = FALSE)
  }
  for (v in x) {
    if (length(out[[v]]) == 0) stop("Regressor '", v, "' has no lags.", call. = FALSE)
    if (any(out[[v]] < 0)) stop("Lags of '", v, "' must be 0 or more.", call. = FALSE)
  }
  out
}

#' Build the estimation data frame from a zoo matrix and lag sets
#'
#' Columns are named `y`, `y_lag1`, `x`, `x_lag1`, ... (lag 0 keeps the
#' plain name). All rows start after `max_lag`, so models with different lag
#' sets estimated with the same `max_lag` share the sample.
#'
#' @param z Numeric zoo matrix from as_adl_data().
#' @param lags Named list from normalize_lags(); first element is y.
#' @param unlagged Names of variables entered as they are (dummies, etc.).
#' @param trend Add a linear trend column `trend` (1, 2, ...)?
#' @param max_lag Observations to drop at the start; at least the largest lag.
#' @return A list with `frame`, `index`, `terms`, `y_col`, `max_lag`.
#' @noRd
build_design <- function(z, lags, unlagged = character(0), trend = FALSE,
                         max_lag = NULL) {
  vars <- names(lags)
  m <- max(c(0L, unlist(lags)))
  if (!is.null(max_lag)) {
    check_whole(max_lag, "max_lag")
    if (max_lag < m) stop("`max_lag` (", max_lag, ") is below the largest lag (",
                          m, ").", call. = FALSE)
    m <- as.integer(max_lag)
  }
  n <- nrow(z)
  if (n - m < 1) stop("Not enough observations for the lags requested.", call. = FALSE)
  rows <- (m + 1):n
  zz <- zoo::coredata(z)

  safe <- stats::setNames(make.names(c(vars, unlagged), unique = TRUE),
                          c(vars, unlagged))
  if (any(safe != names(safe))) {
    message("Non-syntactic names changed: ",
            paste0(names(safe)[safe != names(safe)], " -> ",
                   safe[safe != names(safe)], collapse = ", "))
  }

  col_name <- function(v, k) if (k == 0) safe[[v]] else paste0(safe[[v]], "_lag", k)
  cols <- list()
  terms <- list()
  y <- vars[1]
  cols[[safe[[y]]]] <- zz[rows, y]
  for (v in vars) {
    for (k in lags[[v]]) {
      nm <- col_name(v, k)
      cols[[nm]] <- zz[rows - k, v]
      terms[[length(terms) + 1]] <- data.frame(variable = v, lag = k, column = nm,
                                               role = if (v == y) "ar" else "dl")
    }
  }
  for (v in unlagged) {
    cols[[safe[[v]]]] <- zz[rows, v]
    terms[[length(terms) + 1]] <- data.frame(variable = v, lag = 0L,
                                             column = safe[[v]], role = "unlagged")
  }
  if (trend) {
    cols[["trend"]] <- seq_along(rows)
    terms[[length(terms) + 1]] <- data.frame(variable = "trend", lag = 0L,
                                             column = "trend", role = "trend")
  }

  frame <- as.data.frame(cols, check.names = FALSE)
  terms <- if (length(terms)) do.call(rbind, terms) else
    data.frame(variable = character(), lag = integer(), column = character(),
               role = character())
  list(frame = frame, index = zoo::index(z)[rows], terms = terms,
       y_col = safe[[y]], max_lag = m)
}

# Small helpers -------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

check_whole <- function(v, what) {
  if (!is.numeric(v) || anyNA(v) || any(v != round(v))) {
    stop("`", what, "` must contain whole numbers.", call. = FALSE)
  }
  invisible(v)
}

#' Recycle a scalar or match a named/ordered vector to `x`
#' @noRd
recycle_named <- function(v, x, what) {
  if (length(v) == 1 && is.null(names(v))) return(stats::setNames(rep(v, length(x)), x))
  if (!is.null(names(v))) {
    unknown <- setdiff(names(v), x)
    if (length(unknown)) stop("`", what, "` names unknown regressor(s): ",
                              paste(unknown, collapse = ", "), call. = FALSE)
    out <- stats::setNames(rep(0, length(x)), x)
    out[names(v)] <- v
    return(out)
  }
  if (length(v) != length(x)) stop("`", what, "` must have length 1 or ",
                                   length(x), ".", call. = FALSE)
  stats::setNames(v, x)
}
