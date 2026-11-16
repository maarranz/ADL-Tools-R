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
#' @param start,end Optional sample limits, in the units of the time index.
#' @return A list with `z` (zoo matrix), `index`, `n` and `trimmed`.
#' @noRd
as_adl_data <- function(data, vars, time = NULL, start = NULL, end = NULL) {
  z <- to_zoo(data, time = time)

  missing_vars <- setdiff(vars, colnames(z))
  if (length(missing_vars) > 0) {
    stop("Variable(s) not found in the data: ",
         paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  z <- z[, vars, drop = FALSE]

  if (!is.null(start) || !is.null(end)) {
    idx <- zoo::index(z)
    s <- if (is.null(start)) NULL else as_index_value(start, idx)
    e <- if (is.null(end)) NULL else as_index_value(end, idx)
    z <- stats::window(z, start = s, end = e)
    if (nrow(z) == 0) stop("No observations between `start` and `end`.", call. = FALSE)
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
    freq <- zoo::frequency(zoo::zoo(seq_along(idx), idx))
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
    idx <- df[[time]]
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
