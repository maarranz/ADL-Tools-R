# Dynamic effects of an estimated ADL model ----------------------------------
#
# Same definitions as the adl-theory Python library, applied to the estimated
# coefficients: the recursion psi_h = beta_h + sum_i alpha_i psi_{h-i}, the
# automatic horizon at `tolerance` of the long-run effect, the analytical
# mean lag and the median lag (first lag with half the long-run effect).
#
# Standard errors: delta method for the multipliers (the derivatives follow
# the same recursion), the long run and the mean lag; the median lag is an
# integer (a first crossing), not a smooth function of the coefficients, so
# its interval comes from draws of the coefficients from their asymptotic
# normal distribution.

#' Dynamic multipliers, mean and median lag, and stability
#'
#' Computes, for each regressor, the dynamic multipliers \eqn{d_h} (the
#' effect on \eqn{y_{t+h}} of a one-off unit change in \eqn{x_t}), the
#' cumulative multipliers \eqn{M_h} (the effect of a permanent change),
#' each as a percentage of the long-run multiplier, and the mean and median
#' lags, with their standard errors. Also reports the roots of
#' \eqn{C(z) = 1 - \alpha_1 z - \dots - \alpha_p z^p} and their inverses:
#' the model is dynamically stable when all inverse roots lie inside the
#' unit circle.
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
#' **Standard errors** (`se = TRUE`). The multipliers, the cumulative
#' multipliers, the long-run multiplier and the mean lag are smooth
#' functions of the coefficients, and their standard errors come from the
#' delta method with the covariance of [adl_vcov()]. The derivatives of
#' \eqn{d_h} follow the same recursion as the multipliers:
#' \deqn{\frac{\partial d_h}{\partial \beta_j} = 1\{h = j\} +
#'   \sum_i \alpha_i \frac{\partial d_{h-i}}{\partial \beta_j}, \qquad
#'   \frac{\partial d_h}{\partial \alpha_k} = d_{h-k} +
#'   \sum_i \alpha_i \frac{\partial d_{h-i}}{\partial \alpha_k}.}
#' The long-run standard error is the one of [adl_long_run()].
#'
#' The median lag is an integer (the first lag past one half), so it has no
#' derivative and the delta method does not apply. Its interval is
#' obtained by simulation: `nsim` coefficient vectors are drawn from the
#' asymptotic normal distribution \eqn{N(\hat\theta, \hat V)}, the median lag
#' is recomputed for each, and the interval is formed by the
#' \eqn{\alpha/2} and \eqn{1 - \alpha/2} quantiles of the draws. Draws for
#' which the median lag is not defined (an unstable model, a zero long run,
#' or no crossing within `max_horizon` lags) are dropped, and their share is
#' reported as `median_lag_undefined`; a large share means the interval
#' should not be trusted. The draws use R's random-number generator: call
#' [set.seed()] first for reproducible intervals.
#'
#' All this inference is asymptotic and conditional on the lag orders, and
#' assumes a stable model with \eqn{C(1)} well away from zero: near a unit
#' root the long run and the mean lag are weakly identified and symmetric
#' intervals can be misleading (the median-lag interval then shows a large
#' undefined share or a very long upper limit). The same holds when the
#' long-run multiplier is not significantly different from zero: the
#' percentages and the median lag are then shares of a quantity that may be
#' zero, and their intervals are not informative.
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
#' @param max_horizon Upper limit for the automatic horizon, and for the
#'   search for the median lag in the simulated draws.
#' @param se Compute standard errors and intervals (default `TRUE`).
#' @inheritParams adl_vcov
#' @param alpha Significance level for the intervals.
#' @param dist `"normal"` (default) or `"t"`: the reference distribution for
#'   the delta-method intervals. Both are asymptotic approximations.
#' @param nsim Number of simulated coefficient vectors for the median-lag
#'   interval; `0` skips the simulation (the interval is then `NA`).
#' @return An object of class `"adl_dynamics"`: a list with `effects` (one
#'   row per regressor and lag: `variable`, `lag`, `effect`, `cumulative`,
#'   `pct_effect`, `pct_cumulative`, and with `se = TRUE` also `se_effect`
#'   and `se_cumulative`), `summary` (one row per regressor: `long_run`,
#'   `mean_lag`, `median_lag`, `horizon`, and with `se = TRUE` also
#'   `se_long_run`, `se_mean_lag`, `median_lag_low`, `median_lag_high` and
#'   `median_lag_undefined`), `roots` (roots, inverse roots and their
#'   moduli), `stable`, and the settings of the inference (`se`, `cov_type`,
#'   `maxlags`, `alpha`, `dist`, `crit`, `nsim`).
#' @examples
#' set.seed(1)
#' x <- cumsum(rnorm(200))
#' y <- stats::filter(0.5 * x + rnorm(200), 0.8, method = "recursive")
#' fit <- adl(y ~ x, data = ts(cbind(y = as.numeric(y), x = x)), order = c(1, 1))
#' dyn <- adl_dynamics(fit)
#' dyn
#' plot(dyn)
#' plot(dyn, type = "cumulative")
#' plot(dyn, type = "roots")
#' adl_dynamics(fit, cov_type = "HAC", nsim = 500)
#' @export
adl_dynamics <- function(fit, x = NULL, horizon = NULL, tolerance = 0.99,
                         max_horizon = 500, se = TRUE,
                         cov_type = c("nonrobust", "HAC"), maxlags = NULL,
                         alpha = 0.05, dist = c("normal", "t"), nsim = 1000) {
  check_adl(fit)
  cov_type <- match.arg(cov_type)
  dist <- match.arg(dist)
  if (!is.numeric(tolerance) || length(tolerance) != 1 || tolerance <= 0 ||
      tolerance >= 1) {
    stop("`tolerance` must lie strictly between 0 and 1.", call. = FALSE)
  }
  check_whole(max_horizon, "max_horizon")
  if (length(max_horizon) != 1 || max_horizon < 0) {
    stop("`max_horizon` must be a single integer >= 0.", call. = FALSE)
  }
  if (!is.null(horizon)) {
    check_whole(horizon, "horizon")
    if (length(horizon) != 1 || horizon < 0) {
      stop("`horizon` must be a single integer >= 0.", call. = FALSE)
    }
  }
  if (!is.logical(se) || length(se) != 1 || is.na(se)) {
    stop("`se` must be TRUE or FALSE.", call. = FALSE)
  }
  check_alpha(alpha)
  check_whole(nsim, "nsim")
  if (length(nsim) != 1 || nsim < 0) {
    stop("`nsim` must be a single integer >= 0.", call. = FALSE)
  }
  x <- x %||% fit$x
  bad <- setdiff(x, fit$x)
  if (length(bad)) {
    stop("Not lagged regressors of the model: ", paste(bad, collapse = ", "),
         call. = FALSE)
  }

  b <- stats::coef(fit)
  info <- term_info(fit, names(b))
  ar_coef <- poly_coefs(b, info, fit$y, from = 1)
  ar_pos <- poly_positions(info, fit$y, from = 1)
  roots <- ar_roots(ar_coef)
  stable <- roots$stable
  if (!stable) {
    warning("The estimated model is not dynamically stable; multipliers are ",
            "computed for a finite horizon, but the long run, the percentages ",
            "and the mean and median lags are not defined.", call. = FALSE)
  }
  C1 <- 1 - sum(ar_coef)
  default_horizon <- 20L
  if (!stable && is.null(horizon)) {
    message("Unstable model: multipliers shown up to lag ", default_horizon,
            " (set `horizon` for more).")
  }

  V <- maxlags_used <- NULL
  if (se) {
    V <- adl_vcov(fit, cov_type, maxlags)
    maxlags_used <- attr(V, "maxlags")
    if (anyNA(b)) {
      stop("Standard errors need all coefficients estimated (NA in the fit).",
           call. = FALSE)
    }
  } else if (!is.null(maxlags)) {
    stop("`maxlags` is only used with `se = TRUE`.", call. = FALSE)
  }
  crit <- if (dist == "t") stats::qt(1 - alpha / 2, df_resid(fit)) else
    stats::qnorm(1 - alpha / 2)
  draws <- if (se && nsim > 0) draw_coefficients(b, V, nsim) else NULL

  effects <- list()
  summ <- list()
  for (v in x) {
    beta <- poly_coefs(b, info, v, from = 0)
    beta_pos <- poly_positions(info, v, from = 0)
    lr <- if (stable && abs(C1) > sqrt(.Machine$double.eps)) sum(beta) / C1 else NA_real_
    has_lr <- is.finite(lr) && abs(lr) > sqrt(.Machine$double.eps)
    # multipliers up to the largest horizon that may be needed
    H <- if (!is.null(horizon)) horizon else if (has_lr) max_horizon else
      max(default_horizon, length(beta) - 1)
    d <- multipliers(ar_coef, beta, H)
    if (is.null(horizon) && has_lr) {
      # first lag, at or after the last distributed-lag term, from which the
      # cumulative multiplier stays within (1 - tolerance) of the long run for
      # a full AR cycle (p further lags): an early overshoot does not stop it
      Mfull <- cumsum(d)
      close <- abs(lr - Mfull) <= (1 - tolerance) * abs(lr)
      p_ar <- length(ar_coef)
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
    has_mean <- stable && abs(sum(beta)) > sqrt(.Machine$double.eps)
    if (has_mean) {
      mean_lag <- sum((seq_along(beta) - 1) * beta) / sum(beta) +
        sum(seq_along(ar_coef) * ar_coef) / C1
    }
    median_lag <- if (has_lr) first_half(M, lr) else NA_integer_

    eff <- data.frame(variable = v, lag = seq_along(d) - 1L,
                      effect = d, cumulative = M,
                      pct_effect = pct, pct_cumulative = pctM)
    sm <- data.frame(variable = v, long_run = lr, mean_lag = mean_lag,
                     median_lag = median_lag, horizon = length(d) - 1L)

    if (se) {
      G <- multiplier_gradient(ar_coef, ar_pos, beta_pos, d, length(b))
      Gc <- apply(G, 2, cumsum)
      if (!is.matrix(Gc)) Gc <- matrix(Gc, nrow = 1)
      eff$se_effect <- quad_se(G, V)
      eff$se_cumulative <- quad_se(Gc, V)
      eff <- eff[c("variable", "lag", "effect", "se_effect", "cumulative",
                   "se_cumulative", "pct_effect", "pct_cumulative")]

      se_lr <- se_ml <- NA_real_
      if (is.finite(lr)) {
        g <- numeric(length(b))
        g[beta_pos[!is.na(beta_pos)]] <- 1 / C1
        g[ar_pos[!is.na(ar_pos)]] <- sum(beta) / C1^2
        se_lr <- quad_se(matrix(g, nrow = 1), V)
      }
      if (has_mean) {
        B1 <- sum(beta)
        Bp <- sum((seq_along(beta) - 1) * beta)
        Ap <- sum(seq_along(ar_coef) * ar_coef)
        g <- numeric(length(b))
        j <- seq_along(beta) - 1
        ok <- !is.na(beta_pos)
        g[beta_pos[ok]] <- j[ok] / B1 - Bp / B1^2
        k <- seq_along(ar_coef)
        ok <- !is.na(ar_pos)
        g[ar_pos[ok]] <- k[ok] / C1 + Ap / C1^2
        se_ml <- quad_se(matrix(g, nrow = 1), V)
      }
      med_low <- med_high <- undefined <- NA_real_
      if (!is.null(draws) && has_lr) {
        sim <- simulate_median_lag(draws, ar_pos, beta_pos, max_horizon)
        undefined <- mean(is.na(sim))
        if (undefined < 1) {
          q <- stats::quantile(sim, c(alpha / 2, 1 - alpha / 2), type = 1,
                               na.rm = TRUE, names = FALSE)
          med_low <- q[1]
          med_high <- q[2]
        }
      }
      sm <- data.frame(variable = v, long_run = lr, se_long_run = se_lr,
                       mean_lag = mean_lag, se_mean_lag = se_ml,
                       median_lag = median_lag,
                       median_lag_low = as.integer(med_low),
                       median_lag_high = as.integer(med_high),
                       median_lag_undefined = undefined,
                       horizon = length(d) - 1L)
    }
    effects[[v]] <- eff
    summ[[v]] <- sm
  }

  structure(list(
    effects = do.call(rbind, unname(effects)),
    summary = do.call(rbind, unname(summ)),
    roots = roots$table,
    stable = stable,
    y = fit$y,
    tolerance = if (is.null(horizon)) tolerance else NA_real_,
    se = se,
    cov_type = if (se) cov_type else NA_character_,
    maxlags = if (se) maxlags_used else NA_integer_,
    alpha = alpha,
    dist = dist,
    crit = crit,
    nsim = if (se) nsim else 0L
  ), class = "adl_dynamics")
}

#' Multipliers d_0, ..., d_H by the recursion
#' @noRd
multipliers <- function(ar_coef, beta, H) {
  bb <- numeric(H + 1)
  n <- min(length(beta), H + 1)
  bb[seq_len(n)] <- beta[seq_len(n)]
  if (length(ar_coef) == 0 || all(ar_coef == 0)) return(bb)
  if (length(ar_coef) < length(bb)) {
    return(as.numeric(stats::filter(bb, ar_coef, method = "recursive")))
  }
  d <- bb
  for (h in seq_along(d)) {
    for (i in seq_along(ar_coef)) if (h - i >= 1) d[h] <- d[h] + ar_coef[i] * d[h - i]
  }
  d
}

#' First lag at which the cumulative multiplier reaches half the long run
#' @noRd
first_half <- function(M, lr) {
  hit <- which(M / lr >= 0.5)
  if (length(hit)) as.integer(hit[1] - 1) else NA_integer_
}

#' Jacobian of d_0, ..., d_H with respect to the coefficient vector
#'
#' Rows are lags, columns the coefficients of the fit; `ar_pos` and
#' `beta_pos` give the column of each element of the dense lag polynomials
#' (NA where the lag is not in the model, a coefficient fixed at zero).
#' @noRd
multiplier_gradient <- function(ar_coef, ar_pos, beta_pos, d, K) {
  H1 <- length(d)
  G <- matrix(0, H1, K)
  for (h in seq_len(H1)) {
    g <- numeric(K)
    if (h <= length(beta_pos) && !is.na(beta_pos[h])) g[beta_pos[h]] <- 1
    for (i in seq_along(ar_coef)) {
      if (h - i >= 1) {
        g <- g + ar_coef[i] * G[h - i, ]
        if (!is.na(ar_pos[i])) g[ar_pos[i]] <- g[ar_pos[i]] + d[h - i]
      }
    }
    G[h, ] <- g
  }
  G
}

#' Delta-method standard errors: sqrt(diag(G V G'))
#' @noRd
quad_se <- function(G, V) {
  sqrt(pmax(rowSums((G %*% V) * G), 0))
}

#' Draws from N(b, V), one per row (V may be singular)
#' @noRd
draw_coefficients <- function(b, V, nsim) {
  e <- eigen((V + t(V)) / 2, symmetric = TRUE)
  root <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), length(b))
  z <- matrix(stats::rnorm(nsim * length(b)), nsim)
  sweep(z %*% t(root), 2, b, "+")
}

#' Median lag of each simulated coefficient vector (NA where undefined)
#' @noRd
simulate_median_lag <- function(draws, ar_pos, beta_pos, max_horizon) {
  out <- rep(NA_integer_, nrow(draws))
  for (s in seq_len(nrow(draws))) {
    th <- draws[s, ]
    a <- ifelse(is.na(ar_pos), 0, th[ar_pos])
    be <- ifelse(is.na(beta_pos), 0, th[beta_pos])
    if (length(a) && !ar_roots(a)$stable) next
    C1 <- 1 - sum(a)
    if (abs(C1) <= sqrt(.Machine$double.eps)) next
    lr <- sum(be) / C1
    if (!is.finite(lr) || abs(lr) <= sqrt(.Machine$double.eps)) next
    out[s] <- first_half(cumsum(multipliers(a, be, max_horizon)), lr)
  }
  out
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

#' Positions in the coefficient vector of a dense lag polynomial (NA = absent)
#' @noRd
poly_positions <- function(info, variable, from) {
  sel <- which(info$variable == variable & info$role %in% c("ar", "dl"))
  if (!length(sel)) return(integer(0))
  lags <- info$lag[sel]
  out <- rep(NA_integer_, max(lags) - from + 1)
  out[lags - from + 1] <- sel
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
#' @param ci Draw the pointwise confidence band (the interval of each
#'   multiplier, at the level and with the distribution set in
#'   [adl_dynamics()]) when standard errors are available.
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
  if (isTRUE(x$se)) {
    cat("\nStandard errors: delta method, ", x$cov_type,
        if (identical(x$cov_type, "HAC")) paste0(" (bandwidth ", x$maxlags, ")"),
        " covariance. Median-lag interval (", format(100 * (1 - x$alpha)),
        "%): quantiles of the median lag over ",
        if (x$nsim > 0) paste0(x$nsim, " draws of the coefficients")
        else "no draws (nsim = 0)", ".\n", sep = "")
  } else {
    cat("\nPoint estimates; no standard errors.\n")
  }
  invisible(x)
}

#' @rdname print.adl_dynamics
#' @export
plot.adl_dynamics <- function(x, type = c("effects", "cumulative", "roots"),
                              ci = TRUE, ...) {
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
    sv <- if (type == "effects") e$se_effect else e$se_cumulative
    band <- isTRUE(ci) && !is.null(sv)
    lo <- if (band) yv - x$crit * sv
    hi <- if (band) yv + x$crit * sv
    lr <- x$summary$long_run[x$summary$variable == v]
    rng <- range(c(0, yv, lo, hi, if (type == "cumulative" && is.finite(lr)) lr),
                 na.rm = TRUE)
    graphics::plot(e$lag, yv, type = "h", lwd = 2, ylim = rng, xlab = "Lag",
                   ylab = if (type == "effects") "Dynamic multiplier"
                          else "Cumulative multiplier",
                   main = paste0(v, " on ", x$y), ...)
    if (band) {
      graphics::polygon(c(e$lag, rev(e$lag)), c(lo, rev(hi)),
                        col = "grey85", border = NA)
      graphics::lines(e$lag, yv, type = "h", lwd = 2)
    }
    graphics::points(e$lag, yv, pch = 19, cex = 0.6)
    graphics::abline(h = 0, col = "grey")
    if (type == "cumulative" && is.finite(lr)) graphics::abline(h = lr, lty = 2)
  }
  invisible(x)
}
