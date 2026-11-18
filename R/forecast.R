# Dynamic forecasts from an ADL (ADL4) ----------------------------------------
#
# Iterated (dynamic) forecasts: the lags of y are replaced by their
# forecasts, the regressors follow a known path (scenario conditioning, as in
# ADL4). With parameters treated as known, the h-step error is
# sum_{i<h} psi_i e_{T+h-i}, with psi the weights of 1/C(L), so its variance
# is sigma^2 sum_{i<h} psi_i^2. Parameter uncertainty is not included.

#' Dynamic forecasts with fan charts
#'
#' Iterated forecasts of \eqn{y} from an estimated ADL, as in ADL4: from the
#' forecast origin \eqn{T}, each step uses the forecasts of the earlier steps
#' in place of the unknown lags of \eqn{y}, and the regressors follow a
#' **known path** (the scenario for \eqn{x}). The path is taken from the data
#' for periods inside the data, and from `newdata` beyond its end.
#'
#' The forecast error at horizon \eqn{h} is
#' \deqn{y_{T+h} - \hat y_{T+h|T} = \sum_{i=0}^{h-1} \psi_i \epsilon_{T+h-i},}
#' where the \eqn{\psi_i} are the weights of \eqn{1/C(L)}
#' (\eqn{\psi_0 = 1}, \eqn{\psi_i = \sum_j \alpha_j \psi_{i-j}}), so its
#' variance is \eqn{\sigma^2 \sum_{i<h} \psi_i^2}; \eqn{\sigma} is the
#' standard error of the regression. The intervals are normal,
#' \eqn{\hat y \pm z\, se}, one for each `level`: together they draw the fan
#' chart.
#'
#' Two sources of uncertainty are **not** included, so the intervals are too
#' narrow: the estimation of the coefficients (treated as known), and the
#' path of \eqn{x} (treated as known; if \eqn{x} must itself be forecast,
#' its forecast errors add to these). The first vanishes as the sample
#' grows; the second does not.
#'
#' With `origin` inside the estimation sample the forecasts are an
#' in-sample dynamic simulation: the coefficients were estimated with the
#' observations being forecast, so this is not an out-of-sample evaluation
#' (re-estimate up to the origin with `end =` in [adl()] for that). The
#' actual values, when the data have them, are reported alongside.
#'
#' @param fit An `"adl"` object.
#' @param h Number of steps ahead. By default, as many as the path of the
#'   regressors allows (from the data and `newdata`).
#' @param newdata Values of the regressors (and unlagged terms) after the
#'   end of the data: any data type accepted by [adl()]. With a time index of
#'   the same kind as the data, its rows after the end of the data are used
#'   (so the full data set, extended with the scenario, can be passed);
#'   without one (a data frame with no time column), its rows are taken as
#'   the periods that follow the data, in order. Their dates follow the
#'   calendar when the data are consecutive months, quarters or years (a
#'   month end stays a month end, 31 December stays 31 December, a day
#'   capped in short months stays capped), and a fixed number of days
#'   otherwise (daily or weekly data); for any other convention give
#'   `newdata` with dates. If it contains the dependent
#'   variable, those values are reported as the actual outcomes.
#' @param origin Forecast origin: the last observation used, in the units
#'   of the time index (e.g. `c(2016, 4)` for quarterly data). By default
#'   the end of the estimation sample.
#' @param level Confidence levels of the intervals (default 50%, 80% and
#'   95%).
#' @param time Name of the time column of `newdata`, if it is a data frame
#'   and the column cannot be detected.
#' @return An object of class `"adl_forecast"`: a list with `forecasts` (one
#'   row per step: `index`, `h`, `forecast`, `se`, `actual` (`NA` where
#'   unknown), and `lower_XX`, `upper_XX` for each level), `history` (a
#'   `zoo` series of \eqn{y} up to the origin), `psi` (the weights
#'   \eqn{\psi_0, \dots, \psi_{h-1}}), `sigma`, `origin`, `level`, `y`,
#'   `in_sample` (`TRUE` when the forecast periods overlap the estimation
#'   sample) and `label`.
#' @seealso [chow_forecast_test()] for the test of predictive failure based
#'   on the same dynamic forecasts.
#' @examples
#' fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
#'
#' # x held at its last value for eight quarters (a scenario)
#' x_path <- data.frame(x = rep(tail(adl_sim[, "x"], 1), 8))
#' fc <- adl_forecast(fit, newdata = x_path)
#' fc
#' plot(fc)
#'
#' # a dynamic simulation from inside the sample, with the actual values
#' fc2 <- adl_forecast(fit, origin = c(2016, 4))
#' plot(fc2)
#' @export
adl_forecast <- function(fit, h = NULL, newdata = NULL, origin = NULL,
                         level = c(0.5, 0.8, 0.95), time = NULL) {
  check_adl(fit)
  if (!is.numeric(level) || !length(level) || anyNA(level) || any(level <= 0) ||
      any(level >= 1)) {
    stop("`level` must contain numbers strictly between 0 and 1.", call. = FALSE)
  }
  level <- sort(unique(level))
  tags <- level_tag(level)
  if (anyDuplicated(tags)) {
    stop("`level` contains values too close to tell apart (",
         paste(level[duplicated(tags) | duplicated(tags, fromLast = TRUE)],
               collapse = ", "), "): give distinct levels.", call. = FALSE)
  }
  vars <- c(fit$y, fit$x, fit$unlagged)
  z <- fit$data[, vars, drop = FALSE]
  idx <- zoo::index(z)

  # forecast origin: a row of the data, from the first estimation row on
  first_est <- index_position(fit$index[1], idx)
  last_est <- index_position(fit$index[length(fit$index)], idx)
  o <- if (is.null(origin)) last_est else {
    k <- index_position(as_index_value(origin, idx), idx)
    if (is.na(k) || k < first_est) {
      stop("`origin` must be a date of the data from the start of the ",
           "estimation sample (", format(fit$index[1]), ") on.", call. = FALSE)
    }
    k
  }

  # extend the data with the scenario in newdata
  if (!is.null(newdata)) z <- extend_with_newdata(z, newdata, vars, fit$y, time)
  idx <- zoo::index(z)
  zz <- zoo::coredata(z)
  N <- nrow(zz)

  # periods after the origin with a complete path of the regressors
  needed <- setdiff(vars, fit$y)
  avail <- 0L
  if (o < N) {
    ok <- stats::complete.cases(zz[(o + 1):N, needed, drop = FALSE])
    avail <- if (all(ok)) length(ok) else which(!ok)[1] - 1L
  }
  if (is.null(h)) {
    if (avail == 0L) {
      stop("No values of the regressors after the origin: give their future ",
           "path in `newdata` (or set `h` with `newdata`).", call. = FALSE)
    }
    h <- avail
  } else {
    check_whole(h, "h")
    if (length(h) != 1 || h < 1) stop("`h` must be a single integer >= 1.", call. = FALSE)
    if (h > avail) {
      stop("The path of the regressors covers ", avail, " step(s) after the ",
           "origin, fewer than h = ", h, ". Give the missing values (",
           paste(needed, collapse = ", "), ") in `newdata`.", call. = FALSE)
    }
  }

  b <- stats::coef(fit)
  info <- term_info(fit, names(b))
  cols <- fit$terms
  # forecast recursion, row by row of the design
  ypath <- zz[, fit$y]
  ypath[seq.int(o + 1, N)] <- NA_real_
  f <- numeric(h)
  for (s in seq_len(h)) {
    t <- o + s
    xrow <- numeric(length(b))
    names(xrow) <- names(b)
    for (j in seq_along(b)) {
      nm <- names(b)[j]
      if (nm == "(Intercept)") { xrow[j] <- 1; next }
      role <- info$role[j]
      v <- info$variable[j]
      lag <- info$lag[j]
      xrow[j] <- switch(role,
        ar = ypath[t - lag],
        dl = zz[t - lag, v],
        unlagged = zz[t, v],
        trend = t - first_est + 1,
        stop("Unexpected term '", nm, "'.", call. = FALSE))
    }
    f[s] <- sum(xrow * b)
    ypath[t] <- f[s]
  }

  # psi weights of 1/C(L) and the forecast standard errors
  alpha <- poly_coefs(b, info, fit$y, from = 1)
  psi <- multipliers(alpha, 1, h - 1)
  sigma <- summary(fit$lm)$sigma
  se <- sigma * sqrt(cumsum(psi^2))

  rows <- (o + 1):(o + h)
  actual <- if (fit$y %in% colnames(zz)) zz[rows, fit$y] else rep(NA_real_, h)
  out <- data.frame(index = idx[rows], h = seq_len(h), forecast = f, se = se,
                    actual = unname(actual))
  for (lv in level) {
    zc <- stats::qnorm(1 - (1 - lv) / 2)
    tag <- level_tag(lv)
    out[[paste0("lower_", tag)]] <- f - zc * se
    out[[paste0("upper_", tag)]] <- f + zc * se
  }

  structure(list(
    forecasts = out,
    history = zoo::zoo(zz[seq_len(o), fit$y], idx[seq_len(o)]),
    psi = psi,
    sigma = sigma,
    origin = idx[o],
    level = level,
    y = fit$y,
    in_sample = o < last_est,
    label = fit$label
  ), class = "adl_forecast")
}

#' Append the rows of newdata that follow the data
#' @noRd
extend_with_newdata <- function(z, newdata, vars, y, time) {
  idx <- zoo::index(z)
  undated <- is.data.frame(newdata) && !inherits(newdata, "tbl_ts") &&
    is.null(time) && !any(vapply(newdata, function(col) inherits(col, c(
      "Date", "POSIXt", "yearmon", "yearqtr", "yearweek", "yearquarter",
      "yearmonth")), logical(1)))
  nz <- if (undated) {
    num <- vapply(newdata, is.numeric, logical(1))
    zoo::zoo(as.matrix(newdata[num]), next_index(idx, nrow(newdata)))
  } else {
    to_zoo(newdata, time = time)
  }
  missing_vars <- setdiff(setdiff(vars, y), colnames(nz))
  if (length(missing_vars)) {
    stop("`newdata` lacks: ", paste(missing_vars, collapse = ", "), ".", call. = FALSE)
  }
  nidx <- zoo::index(nz)
  num <- function(i) class(i)[1] %in% c("integer", "numeric")
  if (!identical(class(nidx), class(idx)) && !(num(nidx) && num(idx))) {
    stop("The time index of `newdata` is not of the same kind as that of the ",
         "data; give its rows in order without a time column instead.", call. = FALSE)
  }
  # rows after the end of the data: compare with half the smallest spacing
  # of the data, never with an absolute or epoch-relative tolerance
  num_idx <- as.numeric(idx)
  steps <- diff(num_idx)
  half <- if (length(steps)) min(steps[steps > 0]) / 2 else 0.5
  keep <- as.numeric(nidx) > num_idx[length(num_idx)] + half
  if (!any(keep)) return(z)
  nz <- nz[keep, , drop = FALSE]
  cols <- intersect(vars, colnames(nz))
  add <- matrix(NA_real_, nrow(nz), length(vars), dimnames = list(NULL, vars))
  add[, cols] <- as.matrix(zoo::coredata(nz[, cols, drop = FALSE]))
  ext <- zoo::zoo(rbind(as.matrix(zoo::coredata(z)), add),
                  c(idx, zoo::index(nz)))
  # the scenario must continue the data without gaps
  check_regular(zoo::index(ext))
  ext
}

#' The k periods after the last value of a regular index
#' @noRd
next_index <- function(idx, k) {
  last <- idx[length(idx)]
  if (inherits(idx, "yearqtr")) return(last + seq_len(k) / 4)
  if (inherits(idx, "yearmon")) return(last + seq_len(k) / 12)
  if (length(idx) < 2) {
    stop("Cannot extend a time index with a single observation; give ",
         "`newdata` with dates.", call. = FALSE)
  }
  if (inherits(idx, "Date")) {
    # calendar convention first: consecutive months, quarters or years (as
    # accepted by the validation) are extended in calendar periods, so 31
    # December stays 31 December across leap years and a capped day stays
    # capped; only otherwise a fixed number of days (daily, weekly data)
    mon <- as.integer(format(idx, "%Y")) * 12L + as.integer(format(idx, "%m"))
    msteps <- unique(diff(mon))
    if (length(msteps) == 1 && msteps %in% c(1L, 3L, 12L) &&
        calendar_bad_step(as.Date(idx)) == 0L) {
      ym <- zoo::as.yearmon(last) + msteps * seq_len(k) / 12
      month_end <- all(zoo::as.Date(zoo::as.yearmon(idx), frac = 1) == idx)
      if (month_end) return(zoo::as.Date(ym, frac = 1))
      # the day of the month, capped at the end of short months (the 30th
      # is the 28th in February), as accepted when the data were validated
      anchor <- max(as.integer(format(idx, "%d")))
      mlen <- as.integer(format(zoo::as.Date(ym, frac = 1), "%d"))
      return(zoo::as.Date(ym) + (pmin(anchor, mlen) - 1L))
    }
    d <- as.numeric(diff(idx))
    if (length(unique(d)) == 1) return(last + d[1] * seq_len(k))
    stop("Cannot tell how to extend these dates (neither calendar months, ",
         "quarters or years nor a fixed number of days): give `newdata` ",
         "with dates.", call. = FALSE)
  }
  if (inherits(idx, "POSIXt")) {
    stop("Give `newdata` with times for data indexed by date-times.", call. = FALSE)
  }
  if (is.numeric(idx)) {
    step <- stats::median(diff(as.numeric(idx)))
    return(as.numeric(last) + step * seq_len(k))
  }
  stop("Cannot extend this kind of time index; give `newdata` with dates.",
       call. = FALSE)
}

# Methods -------------------------------------------------------------------

#' Print and plot dynamic forecasts
#'
#' `plot()` draws the fan chart: the history of \eqn{y}, the forecasts and
#' one shaded band per confidence level (darker for the narrower ones), with
#' the actual values as points when they are known.
#'
#' @param x An `"adl_forecast"` object.
#' @param digits Significant digits.
#' @param n_history Number of observations of the history drawn before the
#'   origin; by default three times the horizon, at least 12.
#' @param main,xlab,ylab Graphical labels.
#' @param ... Further arguments to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
print.adl_forecast <- function(x, digits = 4, ...) {
  title <- paste0("Dynamic forecasts of ", x$y, " from the ", x$label)
  cat(title, "\n", strrep("-", nchar(title)), "\n", sep = "")
  cat("Origin: ", format(x$origin), "; ", nrow(x$forecasts), " step(s) ahead; ",
      "sigma = ", format(signif(x$sigma, digits)), "\n", sep = "")
  if (x$in_sample) {
    cat("The forecast periods are inside the estimation sample (in-sample ",
        "dynamic simulation).\n", sep = "")
  }
  cat("\n")
  tab <- x$forecasts
  tab$index <- format(tab$index)
  if (all(is.na(tab$actual))) tab$actual <- NULL
  tab[] <- lapply(tab, function(col) if (is.double(col)) signif(col, digits) else col)
  print(tab, row.names = FALSE)
  cat("\nRegressors on the given path; intervals from the error variance ",
      "only (no parameter or x uncertainty).\n", sep = "")
  invisible(x)
}

#' @rdname print.adl_forecast
#' @export
plot.adl_forecast <- function(x, n_history = NULL, main = NULL, xlab = "",
                              ylab = x$y, ...) {
  fc <- x$forecasts
  H <- nrow(fc)
  n_history <- n_history %||% max(12L, 3L * H)
  nh <- length(x$history)
  hist <- x$history[seq.int(max(1L, nh - n_history + 1L), nh)]
  th <- zoo::index(hist)
  tf <- fc$index
  tt <- c(th, tf)
  levs <- sort(x$level, decreasing = TRUE)
  lows <- unlist(fc[paste0("lower_", level_tag(levs))])
  highs <- unlist(fc[paste0("upper_", level_tag(levs))])
  ylim <- range(c(zoo::coredata(hist), fc$forecast, lows, highs, fc$actual),
                na.rm = TRUE)
  main <- main %||% paste0("Dynamic forecasts of ", x$y, " (origin ",
                           format(x$origin), ")")
  graphics::plot(tt, c(zoo::coredata(hist), rep(NA, H)), type = "l", ylim = ylim,
                 main = main, xlab = xlab, ylab = ylab, ...)
  # bands joined to the last observation, widest first
  last_t <- as.numeric(th[length(th)])
  last_y <- as.numeric(hist[length(hist)])
  greys <- grDevices_grey(length(levs))
  for (i in seq_along(levs)) {
    tag <- level_tag(levs[i])
    lo <- c(last_y, fc[[paste0("lower_", tag)]])
    hi <- c(last_y, fc[[paste0("upper_", tag)]])
    xs <- c(last_t, as.numeric(tf))
    graphics::polygon(c(xs, rev(xs)), c(lo, rev(hi)), col = greys[i], border = NA)
  }
  graphics::lines(c(last_t, as.numeric(tf)), c(last_y, fc$forecast), lwd = 2)
  if (any(!is.na(fc$actual))) graphics::points(tf, fc$actual, pch = 19, cex = 0.6)
  graphics::abline(v = last_t, lty = 3, col = "grey40")
  graphics::legend("topleft", bty = "n", fill = greys,
                   legend = paste0(level_tag(levs), "%"), border = NA,
                   cex = 0.8)
  invisible(x)
}

#' Column tag of a confidence level (one per level, so a vector of levels
#' gets exactly the tags the columns were built with: 0.95 -> "95",
#' 0.999 -> "99.9")
#' @noRd
level_tag <- function(lv) vapply(lv, function(v) format(100 * v, digits = 15), "")

#' Shades of grey, light (widest band) to dark (narrowest)
#' @noRd
grDevices_grey <- function(k) {
  vals <- seq(0.88, 0.62, length.out = k)
  vapply(vals, function(v) {
    hx <- sprintf("%02X", round(255 * v))
    paste0("#", hx, hx, hx)
  }, character(1))
}
