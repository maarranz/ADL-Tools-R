# Input layer ---------------------------------------------------------------
#
# Every user-facing function converts its data with as_adl_data(), so the
# rest of the package only ever sees one structure: a numeric zoo matrix with
# a sorted, unique time index and no missing values.

#' Convert supported inputs to a numeric zoo matrix
#'
#' Accepts `ts`/`mts`, `zoo`/`xts`, data frames and tibbles (with a time
#' column or in row order) and tsibbles. Selects `vars`, applies the sample
#' window, trims leading and trailing missing values and stops on internal
#' missing values.
#'
#' @param data Input data.
#' @param vars Character vector of variable names to keep.
#' @param time Name of the time column for data frames; `NULL` to detect it.
#' @param start,end Optional limits of the estimation sample, in the units
#'   of the time index: `start` is the first observation estimated.
#' @param pre Number of observations kept before `start` for the initial
#'   lags (the largest lag of the model).
#' @return A list with `z` (zoo matrix), `index`, `n` and `trimmed`.
#' @noRd
as_adl_data <- function(data, vars, time = NULL, start = NULL, end = NULL,
                        pre = 0) {
  z <- to_zoo(data, time = time)

  missing_vars <- setdiff(vars, colnames(z))
  if (length(missing_vars) > 0) {
    stop("Variable(s) not found in the data: ",
         paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  z <- z[, vars, drop = FALSE]

  if (!is.null(start) || !is.null(end)) {
    # start/end delimit the estimation sample (as in dynlm and EViews): the
    # `pre` observations before `start` are kept as initial lags
    idx <- zoo::index(z)
    first <- 1L
    last <- length(idx)
    if (!is.null(start)) {
      s <- as_index_value(start, idx)
      pos <- which(idx >= s)[1]
      if (is.na(pos)) stop("`start` is after the end of the data.", call. = FALSE)
      first <- pos - pre
      if (first < 1) {
        message("Not enough observations before `start` for ", pre,
                " lag(s); the estimation sample starts at ",
                format(idx[min(pre + 1, length(idx))]), ".")
        first <- 1L
      }
    }
    if (!is.null(end)) {
      e <- as_index_value(end, idx)
      last <- max(c(0L, which(idx <= e)))
    }
    if (last < first) stop("No observations between `start` and `end`.", call. = FALSE)
    z <- z[first:last, , drop = FALSE]
  }

  ok <- stats::complete.cases(zoo::coredata(z))
  if (!any(ok)) stop("No complete observations in the data.", call. = FALSE)
  first <- min(which(ok))
  last <- max(which(ok))
  if (!all(ok[first:last])) {
    bad <- zoo::index(z)[first:last][!ok[first:last]]
    stop("Missing values inside the sample (at ", format(bad[1]),
         if (length(bad) > 1) paste0(" and ", length(bad) - 1, " more"),
         "). Missing values would break the lag structure; ",
         "fill or remove them before estimating.", call. = FALSE)
  }
  trimmed <- (first - 1) + (nrow(z) - last)
  if (trimmed > 0) {
    message("Trimmed ", trimmed,
            " observation(s) with missing values at the start or end of the sample.")
    z <- z[first:last, , drop = FALSE]
  }

  # gaps matter only inside the observations used (estimation sample and
  # its initial lags), not elsewhere in the data
  check_regular(zoo::index(z))
  storage.mode(z) <- "double"
  list(z = z, index = zoo::index(z), n = nrow(z), trimmed = trimmed)
}

#' Convert any supported input to a zoo matrix with column names
#' @noRd
to_zoo <- function(data, time = NULL) {
  if (inherits(data, "tbl_ts")) {
    if (!requireNamespace("tsibble", quietly = TRUE)) {
      stop("Package 'tsibble' is needed to read tsibble objects.", call. = FALSE)
    }
    if (length(tsibble::key_vars(data)) > 0) {
      stop("tsibbles with key variables (several series per variable) are ",
           "not supported; filter to a single series first.", call. = FALSE)
    }
    time <- tsibble::index_var(data)
    data <- as.data.frame(data)
  }

  if (inherits(data, "zoo")) {
    z <- zoo::as.zoo(data)
  } else if (stats::is.ts(data)) {
    z <- zoo::as.zoo(data)
  } else if (is.data.frame(data)) {
    z <- df_to_zoo(as.data.frame(data), time = time)
  } else {
    stop("Unsupported data type: ", paste(class(data), collapse = "/"),
         ". Use a ts, zoo, xts, data frame, tibble or tsibble.", call. = FALSE)
  }

  if (is.null(dim(z))) {
    stop("A single unnamed series was supplied; pass all the variables ",
         "together (multivariate ts, zoo/xts with columns, or a data frame).",
         call. = FALSE)
  }
  if (is.null(colnames(z)) || anyNA(colnames(z))) {
    stop("The data must have column names.", call. = FALSE)
  }
  idx <- zoo::index(z)
  if (anyDuplicated(idx)) stop("The time index has duplicated values.", call. = FALSE)
  z
}

#' Row of an index value: the one row within half the smallest step (exact
#' for dates and times, robust to floating point for ts/yearqtr), or NA
#' @noRd
index_position <- function(v, idx) {
  num <- as.numeric(idx)
  val <- as.numeric(v)
  steps <- diff(sort(unique(num)))
  tol <- if (length(steps) && any(steps > 0)) 1e-6 * min(steps[steps > 0]) else
    1e-8 * max(1, abs(val))
  k <- which(abs(num - val) <= tol)
  if (length(k) == 1) k else NA_integer_
}

#' Stop if a time index has gaps
#'
#' Lags are taken by row, so a missing period would make the "lag" of the
#' next observation two periods old. Accepted: equal steps (ts, integer or
#' daily/weekly dates), or consecutive months, quarters or years for dates
#' and yearmon/yearqtr indexes. For dates, every observation must also sit
#' at the same point of its period (the same day of the month, or the last
#' day of the month; for quarters and years also the same month of the
#' period), so 1 Jan, 25 Feb, 1 Mar is not read as three consecutive months.
#' Irregular data are not supported.
#' @noRd
check_regular <- function(idx) {
  if (length(idx) < 3) return(invisible(TRUE))
  # position of the first irregular step, or 0 if all steps are equal
  bad_step <- function(v, step = NULL) {
    d <- diff(v)
    step <- step %||% d[1]
    tol <- 1e-8 * max(1, abs(step))
    k <- which(abs(d - step) >= tol)
    if (length(k)) k[1] else 0L
  }
  k <- if (inherits(idx, "yearqtr")) {
    bad_step(as.numeric(idx) * 4, 1)
  } else if (inherits(idx, "yearmon")) {
    bad_step(as.numeric(idx) * 12, 1)
  } else if (inherits(idx, "POSIXt")) {
    # times: equal steps in seconds (minutes, hours, days); or, when every
    # observation is at the same clock time, consecutive calendar periods
    secs <- as.numeric(idx)
    clock <- format(idx, "%H:%M:%OS6")
    if (bad_step(secs) == 0L) 0L else if (length(unique(clock)) == 1) {
      dd <- as.Date(format(idx, "%Y-%m-%d"))
      # daily or weekly across daylight-saving changes, or calendar periods
      if (bad_step(as.numeric(dd)) == 0L) 0L else calendar_bad_step(dd)
    } else bad_step(secs, stats::median(diff(secs)))
  } else if (inherits(idx, "Date")) {
    dd <- as.Date(idx)
    if (bad_step(as.numeric(dd)) == 0L) 0L else calendar_bad_step(dd)
  } else if (is.numeric(idx)) {
    bad_step(as.numeric(idx))
  } else {
    0L
  }
  if (k > 0L) {
    stop("The time index has gaps or is irregular (between ", format(idx[k]),
         " and ", format(idx[k + 1]), "). Lags are taken period by period, so ",
         "a missing period would be skipped. Insert the missing periods or ",
         "estimate on a stretch without gaps; irregular data are not supported.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Dates as consecutive months, quarters or years
#'
#' Returns 0 if the dates are consecutive periods of one of these
#' frequencies, each at the same point of its period; otherwise the position
#' of the first offending step, for the frequency that fits longest.
#' @noRd
calendar_bad_step <- function(dd) {
  n <- length(dd)
  mon <- as.integer(format(dd, "%m"))
  mday <- as.integer(format(dd, "%d"))
  ym <- as.integer(format(dd, "%Y")) * 12L + mon - 1L          # month number
  yr <- as.integer(format(dd, "%Y"))
  next1 <- as.Date(paste(yr + (mon == 12L), mon %% 12L + 1L, 1, sep = "-"))
  dim <- as.integer(format(next1 - 1, "%d"))                   # days in month
  eom <- mday == dim
  # same day of the month: either the same day, or the month end, or a day
  # capped by short months (the 30th, and the 28th/29th in February)
  same_day <- function(i) {
    D <- mday[i[1]]
    Dmax <- max(mday[i])
    all(mday[i] == D) || all(eom[i]) || all(mday[i] == pmin(Dmax, dim[i]))
  }
  first_bad <- function(m) {
    for (k in seq_len(n - 1)) {
      if (ym[k + 1] - ym[k] != m || !same_day(seq_len(k + 1))) return(k)
    }
    0L
  }
  ks <- vapply(c(1L, 3L, 12L), first_bad, integer(1))
  if (any(ks == 0L)) 0L else max(ks)
}

#' Convert tsibble index classes to their zoo equivalents
#'
#' `yearquarter` -> `yearqtr`, `yearmonth` -> `yearmon`, `yearweek` -> `Date`
#' (first day of the week). Other classes are left as they are.
#' @noRd
standard_index <- function(idx) {
  if (inherits(idx, "yearquarter")) return(zoo::as.yearqtr(as.Date(idx)))
  if (inherits(idx, "yearmonth")) return(zoo::as.yearmon(as.Date(idx)))
  if (inherits(idx, "yearweek")) return(as.Date(idx))
  idx
}

#' Translate a user sample limit to the class of the index
#'
#' Accepts ts-style `c(year, period)` for yearqtr, yearmon and numeric
#' indexes, character dates for Date indexes, and plain values otherwise.
#' @noRd
as_index_value <- function(v, idx) {
  if (inherits(idx, "yearqtr")) {
    if (is.numeric(v) && length(v) == 2) return(zoo::as.yearqtr(v[1] + (v[2] - 1) / 4))
    return(zoo::as.yearqtr(v))
  }
  if (inherits(idx, "yearmon")) {
    if (is.numeric(v) && length(v) == 2) return(zoo::as.yearmon(v[1] + (v[2] - 1) / 12))
    return(zoo::as.yearmon(v))
  }
  if (inherits(idx, "Date")) return(as.Date(v))
  if (is.numeric(idx) && is.numeric(v) && length(v) == 2) {
    freq <- stats::frequency(zoo::zoo(seq_along(idx), idx))
    return(v[1] + (v[2] - 1) / freq)
  }
  if (length(v) != 1) stop("`start`/`end` must be a single value for this index.", call. = FALSE)
  v
}

#' Data frame to zoo, using a time column or row order
#' @noRd
df_to_zoo <- function(df, time = NULL) {
  if (is.null(time)) {
    is_time <- vapply(df, function(col) {
      inherits(col, c("Date", "POSIXt", "yearmon", "yearqtr", "yearweek",
                      "yearquarter", "yearmonth"))
    }, logical(1))
    if (sum(is_time) == 1) {
      time <- names(df)[is_time]
    } else if (sum(is_time) > 1) {
      stop("Several time columns found (", paste(names(df)[is_time], collapse = ", "),
           "); choose one with `time =`.", call. = FALSE)
    }
  }

  if (is.null(time)) {
    message("No time column found; using row order as the time index.")
    idx <- seq_len(nrow(df))
    vals <- df
  } else {
    if (!time %in% names(df)) {
      stop("Time column '", time, "' not found in the data.", call. = FALSE)
    }
    idx <- standard_index(df[[time]])
    vals <- df[setdiff(names(df), time)]
  }

  num <- vapply(vals, is.numeric, logical(1))
  vals <- vals[num]
  if (ncol(vals) == 0) stop("No numeric columns in the data.", call. = FALSE)

  ord <- order(idx)
  mat <- as.matrix(vals)[ord, , drop = FALSE]
  rownames(mat) <- NULL
  zoo::zoo(mat, idx[ord])
}
