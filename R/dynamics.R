# Dynamic effects of an estimated ADL model ----------------------------------
#
# Same definitions as the adl-theory Python library, applied to the estimated
# coefficients: the recursion psi_h = beta_h + sum_i alpha_i psi_{h-i}, the
# automatic horizon at `tolerance` of the long-run effect, the analytical
# mean lag and the median lag (first lag with half the long-run effect).

#' Dynamic multipliers, mean and median lag, and stability
#'
#' Computes, for each regressor, the dynamic multipliers \eqn{d_h} (the
#' effect on \eqn{y_{t+h}} of a one-off unit change in \eqn{x_t}), the
#' cumulative multipliers \eqn{M_h} (the effect of a permanent change),
#' each as a percentage of the long-run multiplier, and the mean and median
#' lags. Also reports the roots of \eqn{C(z) = 1 - \alpha_1 z - \dots -
#' \alpha_p z^p} and their inverses: the model is dynamically stable when
#' all inverse roots lie inside the unit circle.
#'
#' The multipliers follow the recursion
#' \deqn{d_h = \beta_h + \sum_{i=1}^{p} \alpha_i d_{h-i},}
#' with \eqn{\beta_h = 0} beyond the last lag. The mean lag is
#' \eqn{B'(1)/B(1) + \sum i\alpha_i / C(1)}, the first moment of the
#' multipliers weighted by \eqn{d_h/k}; when the multipliers change sign it
#' is a signed moment and can be negative or exceed the horizon, so read it
#' as an average lag only for responses of one sign. The median lag is the
#' first lag at which the cumulative multiplier reaches half the long-run
#' multiplier (the first crossing, for responses that overshoot).
#'
#' For an unstable model the multipliers are still computed for a finite
#' horizon, showing how the response grows, but the long-run multiplier,
#' the percentages and the mean and median lags are `NA`.
#'
#' These are point estimates; their sampling uncertainty is not reported.
#'
#' @param fit An `"adl"` object.
#' @param x Regressors to analyse; by default all the lagged regressors.
#' @param horizon Last lag computed. `NULL` (default) chooses it
#'   automatically: at least the last distributed lag, and then the first lag
#'   from which the cumulative multiplier stays within `1 - tolerance` of the
#'   long-run multiplier for a full AR cycle (at most `max_horizon`). For an
#'   unstable model, or one with a zero long run, the default is 20 lags.
#' @param tolerance Share of the long-run multiplier that sets the automatic
#'   horizon.
#' @param max_horizon Upper limit for the automatic horizon.
#' @return An object of class `"adl_dynamics"`: a list with `effects` (one
#'   row per regressor and lag: `variable`, `lag`, `effect`, `cumulative`,
#'   `pct_effect`, `pct_cumulative`), `summary` (one row per regressor:
#'   `long_run`, `mean_lag`, `median_lag`, `horizon`), `roots` (roots,
#'   inverse roots and their moduli) and `stable`.
#' @examples
#' set.seed(1)
#' x <- cumsum(rnorm(200))
#' y <- stats::filter(0.5 * x + rnorm(200), 0.8, method = "recursive")
#' fit <- adl(y ~ x, data = ts(cbind(y = as.numeric(y), x = x)), order = c(1, 1))
#' dyn <- adl_dynamics(fit)
#' dyn
#' plot(dyn)
#' plot(dyn, type = "roots")
#' @export
adl_dynamics <- function(fit, x = NULL, horizon = NULL, tolerance = 0.99,
                         max_horizon = 500) {
  check_adl(fit)
  if (!is.numeric(tolerance) || length(tolerance) != 1 || tolerance <= 0 ||
      tolerance >= 1) {
    stop("`tolerance` must lie strictly between 0 and 1.", call. = FALSE)
  }
  check_whole(max_horizon, "max_horizon")
  if (!is.null(horizon)) {
    check_whole(horizon, "horizon")
    if (length(horizon) != 1 || horizon < 0) {
      stop("`horizon` must be a single integer >= 0.", call. = FALSE)
    }
  }
  x <- x %||% fit$x
  bad <- setdiff(x, fit$x)
  if (length(bad)) {
    stop("Not lagged regressors of the model: ", paste(bad, collapse = ", "),
         call. = FALSE)
  }

  b <- stats::coef(fit)
  info <- term_info(fit, names(b))
  alpha <- poly_coefs(b, info, fit$y, from = 1)
  roots <- ar_roots(alpha)
  stable <- roots$stable
  if (!stable) {
    warning("The estimated model is not dynamically stable; multipliers are ",
            "computed for a finite horizon, but the long run, the percentages ",
            "and the mean and median lags are not defined.", call. = FALSE)
  }
  C1 <- 1 - sum(alpha)
  default_horizon <- 20L
  if (!stable && is.null(horizon)) {
    message("Unstable model: multipliers shown up to lag ", default_horizon,
            " (set `horizon` for more).")
  }

  effects <- list()
  summ <- list()
  for (v in x) {
    beta <- poly_coefs(b, info, v, from = 0)
    lr <- if (stable && abs(C1) > sqrt(.Machine$double.eps)) sum(beta) / C1 else NA_real_
    has_lr <- is.finite(lr) && abs(lr) > sqrt(.Machine$double.eps)
    # multipliers up to the largest horizon that may be needed
    H <- if (!is.null(horizon)) horizon else if (has_lr) max_horizon else
      max(default_horizon, length(beta) - 1)
    d <- numeric(H + 1)
    for (h in 0:H) {
      bh <- if (h + 1 <= length(beta)) beta[h + 1] else 0
      ar <- 0
      for (i in seq_along(alpha)) if (h - i >= 0) ar <- ar + alpha[i] * d[h - i + 1]
      d[h + 1] <- bh + ar
    }
    if (is.null(horizon) && has_lr) {
      # first lag, at or after the last distributed-lag term, from which the
      # cumulative multiplier stays within (1 - tolerance) of the long run for
      # a full AR cycle (p further lags): an early overshoot does not stop it
      Mfull <- cumsum(d)
      close <- abs(lr - Mfull) <= (1 - tolerance) * abs(lr)
      p_ar <- length(alpha)
      first_ok <- NA_integer_
      for (h in seq(max(0L, length(beta) - 1L), H - p_ar, length.out =
                    max(0L, H - p_ar - max(0L, length(beta) - 1L) + 1L))) {
        if (all(close[(h + 1):(h + 1 + p_ar)])) { first_ok <- h; break }
      }
      if (is.na(first_ok)) {
        warning("The cumulative multiplier of '", v, "' did not settle within ",
                format(100 * (1 - tolerance)), "% of the long-run multiplier by lag ",
                max_horizon, ".", call. = FALSE)
        first_ok <- H
      }
      d <- d[seq_len(first_ok + 1)]
    }
    M <- cumsum(d)
    pct <- if (has_lr) 100 * d / lr else rep(NA_real_, length(d))
    pctM <- if (has_lr) 100 * M / lr else rep(NA_real_, length(d))

    mean_lag <- NA_real_
    if (stable && abs(sum(beta)) > sqrt(.Machine$double.eps)) {
      mean_lag <- sum((seq_along(beta) - 1) * beta) / sum(beta) +
        sum(seq_along(alpha) * alpha) / C1
    }
    median_lag <- NA_integer_
    if (has_lr) {
      hit <- which(M / lr >= 0.5)
      if (length(hit)) median_lag <- as.integer(hit[1] - 1)
    }

    effects[[v]] <- data.frame(variable = v, lag = seq_along(d) - 1L,
                               effect = d, cumulative = M,
                               pct_effect = pct, pct_cumulative = pctM)
    summ[[v]] <- data.frame(variable = v, long_run = lr, mean_lag = mean_lag,
                            median_lag = median_lag,
                            horizon = length(d) - 1L)
  }

  structure(list(
    effects = do.call(rbind, unname(effects)),
    summary = do.call(rbind, unname(summ)),
    roots = roots$table,
    stable = stable,
    y = fit$y,
    tolerance = if (is.null(horizon)) tolerance else NA_real_
  ), class = "adl_dynamics")
}

#' Coefficients of one lag polynomial as a dense vector
#'
#' Element k (1-based) is the coefficient at lag `from + k - 1`; lags not in
#' the model are zero.
#' @noRd
poly_coefs <- function(b, info, variable, from) {
  sel <- info$variable == variable & info$role %in% c("ar", "dl")
  if (!any(sel)) return(numeric(0))
  lags <- info$lag[sel]
  out <- numeric(max(lags) - from + 1)
  out[lags - from + 1] <- b[sel]
  out
}

#' Roots and inverse roots of C(z) = 1 - a_1 z - ... - a_p z^p
#' @noRd
ar_roots <- function(alpha) {
  if (length(alpha) == 0 || all(alpha == 0)) {
    return(list(stable = TRUE,
                table = data.frame(root = complex(0), modulus = numeric(0),
                                   inverse_root = complex(0),
                                   inverse_modulus = numeric(0))))
  }
  # drop trailing zeros so polyroot gets the right degree
  last <- max(which(alpha != 0))
  r <- polyroot(c(1, -alpha[seq_len(last)]))
  tab <- data.frame(root = r, modulus = Mod(r), inverse_root = 1 / r,
                    inverse_modulus = 1 / Mod(r))
  tab <- tab[order(-tab$inverse_modulus), , drop = FALSE]
  rownames(tab) <- NULL
  list(stable = all(Mod(r) > 1 + sqrt(.Machine$double.eps)), table = tab)
}

# Methods -------------------------------------------------------------------

#' Print and plot dynamic effects
#'
#' @param x An `"adl_dynamics"` object.
#' @param digits Significant digits.
#' @param n Number of lags shown per regressor; `NULL` shows them all.
#' @param type `"effects"` (dynamic multipliers), `"cumulative"` (cumulative
#'   multipliers with the long-run multiplier as a dashed line) or
#'   `"roots"` (inverse roots and the unit circle).
#' @param ... Further graphical arguments.
#' @return `x`, invisibly.
#' @export
print.adl_dynamics <- function(x, digits = 4, n = 12, ...) {
  title <- "Dynamic effects of the estimated ADL model"
  cat(title, "\n", strrep("-", nchar(title)), "\n", sep = "")
  cat("Dynamically stable: ", if (x$stable) "yes" else "NO", sep = "")
  if (nrow(x$roots)) {
    cat("  (largest inverse root modulus: ",
        format(signif(max(x$roots$inverse_modulus), digits)), ")", sep = "")
  }
  cat("\n\n")
  s <- x$summary
  s[] <- lapply(s, function(col) if (is.double(col)) signif(col, digits) else col)
  print(s, row.names = FALSE)
  for (v in unique(x$effects$variable)) {
    e <- x$effects[x$effects$variable == v, -1]
    total <- nrow(e)
    if (!is.null(n) && total > n) e <- e[seq_len(n), ]
    e[] <- lapply(e, function(col) if (is.double(col)) signif(col, digits) else col)
    cat("\nMultipliers of ", v, " on ", x$y, ":\n", sep = "")
    print(e, row.names = FALSE)
    if (nrow(e) < total) cat("... ", total - nrow(e), " more lags (see $effects)\n", sep = "")
  }
  cat("\nPoint estimates; no standard errors.\n")
  invisible(x)
}

#' @rdname print.adl_dynamics
#' @export
plot.adl_dynamics <- function(x, type = c("effects", "cumulative", "roots"), ...) {
  type <- match.arg(type)
  if (type == "roots") {
    th <- seq(0, 2 * pi, length.out = 361)
    ir <- x$roots$inverse_root
    lim <- max(1.05, Mod(ir) * 1.1)
    graphics::plot(cos(th), sin(th), type = "l", asp = 1, xlim = c(-lim, lim),
                   ylim = c(-lim, lim), xlab = "Real", ylab = "Imaginary",
                   main = "Inverse roots of C(L)", ...)
    graphics::abline(h = 0, v = 0, col = "grey")
    if (length(ir)) graphics::points(Re(ir), Im(ir), pch = 19,
                                     col = ifelse(Mod(ir) < 1, "black", "red"))
    return(invisible(x))
  }
  vars <- unique(x$effects$variable)
  old <- graphics::par(mfrow = c(1, length(vars)))
  on.exit(graphics::par(old))
  for (v in vars) {
    e <- x$effects[x$effects$variable == v, ]
    yv <- if (type == "effects") e$effect else e$cumulative
    lr <- x$summary$long_run[x$summary$variable == v]
    rng <- range(c(0, yv, if (type == "cumulative" && is.finite(lr)) lr))
    graphics::plot(e$lag, yv, type = "h", lwd = 2, ylim = rng, xlab = "Lag",
                   ylab = if (type == "effects") "Dynamic multiplier"
                          else "Cumulative multiplier",
                   main = paste0(v, " on ", x$y), ...)
    graphics::points(e$lag, yv, pch = 19, cex = 0.6)
    graphics::abline(h = 0, col = "grey")
    if (type == "cumulative" && is.finite(lr)) graphics::abline(h = lr, lty = 2)
  }
  invisible(x)
}
