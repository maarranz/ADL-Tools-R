# Common factor (COMFAC) restrictions, generalized (ADL2) --------------------
#
# Order-r common factor: C(L) and every B_i(L) share a degree-r factor
# omega(L) = 1 - omega_1 L - ... - omega_r L^r, so that
#   C(L) = omega(L) phi(L),   B_i(L) = omega(L) delta_i(L),
# and the ADL reduces to phi(L) y_t = c + sum_i delta_i(L) x_it + u_t with
# AR(r) errors omega(L) u_t = eps_t. That imposes r * n restrictions.
#
# LR: the restricted model is fitted by nonlinear least squares on the same
# design matrix (same sample) as the unrestricted ADL, its coefficients being
# the polynomial products above; LR = T log(RSS_r / RSS_u) ~ chi2(r n).

#' Common factor (COMFAC) test
#'
#' Tests whether the lag polynomials of the ADL model share a common factor
#' of order `r`, i.e. whether the model is a regression with fewer lags and
#' AR(`r`) errors. For the ADL(1,1),
#' \deqn{y_t = \alpha_0 + \alpha_1 y_{t-1} + \beta_0 x_t + \beta_1 x_{t-1}
#' + \epsilon_t,}
#' the order-1 restriction is \eqn{\alpha_1 \beta_0 + \beta_1 = 0}: the
#' model is a static regression \eqn{y_t = \mu + \beta_0 x_t + u_t} with
#' AR(1) errors (ADL2).
#'
#' In general, with \eqn{C(L) = 1 - \alpha_1 L - \dots - \alpha_p L^p} and
#' regressor polynomials \eqn{B_i(L)}, the null is \eqn{C(L) =
#' \omega(L)\phi(L)} and \eqn{B_i(L) = \omega(L)\delta_i(L)} for a common
#' \eqn{\omega(L)} of degree \eqn{r}: \eqn{r \cdot n} restrictions with
#' \eqn{n} regressors. It requires \eqn{p \ge r} and at least \eqn{r + 1}
#' consecutive lags of each regressor.
#'
#' * `method = "LR"` (default): the restricted model is fitted by
#'   nonlinear least squares on the same sample, and
#'   \eqn{LR = T \log(RSS_r/RSS_u) \sim \chi^2_{rn}}. It does not depend on
#'   how the restriction is written and handles complex common roots.
#' * `method = "Wald"`: available for \eqn{r = 1} and models with one lag of
#'   \eqn{y} and lags 0 and 1 of each regressor, testing
#'   \eqn{\alpha_1 \beta_{i0} + \beta_{i1} = 0} with the delta method, as in
#'   ADL2. Wald tests of nonlinear restrictions depend on how the
#'   restriction is written (Gregory and Veall, 1985); the LR test does not.
#'
#' Test \eqn{r = 1, 2, \dots} in sequence (Sargan, 1980). Intercepts are
#' left free; models with unlagged terms or a trend are not supported yet.
#'
#' @param fit An `"adl"` object.
#' @param r Order of the common factor.
#' @param method `"LR"` (default) or `"Wald"`.
#' @param cov_type,maxlags Covariance for the Wald test (see [adl_vcov()]);
#'   not used by the LR test.
#' @param alpha Significance level for the decision.
#' @return An object of class `"adltools_test"`. For LR it also holds
#'   `omega` (the AR coefficients of the common factor), `restricted` (a
#'   data frame of the restricted model's coefficients: \eqn{c}, \eqn{\phi},
#'   \eqn{\delta_i}, \eqn{\omega}), `rss_restricted` and
#'   `rss_unrestricted`; the Wald form does not estimate the restricted
#'   model.
#' @references Sargan, J. D. (1980). Some tests of dynamic specification for
#'   a single equation. *Econometrica*, 48(4), 879-897.
#'
#'   Hendry, D. F. and Mizon, G. E. (1978). Serial correlation as a
#'   convenient simplification, not a nuisance: a comment on a study of the
#'   demand for money by the Bank of England. *Economic Journal*, 88(351),
#'   549-563.
#'
#'   Gregory, A. W. and Veall, M. R. (1985). Formulating Wald tests of
#'   nonlinear restrictions. *Econometrica*, 53(6), 1465-1468.
#' @examples
#' fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
#' comfac_test(fit)
#' comfac_test(fit, method = "Wald")
#' fit2 <- adl(y ~ x, data = adl_sim, order = c(2, 2))
#' comfac_test(fit2, r = 2)
#' @export
comfac_test <- function(fit, r = 1, method = c("LR", "Wald"),
                        cov_type = c("nonrobust", "HAC"), maxlags = NULL,
                        alpha = 0.05) {
  check_adl(fit)
  method <- match.arg(method)
  cov_type <- match.arg(cov_type)
  check_alpha(alpha)
  check_count(r, "r")
  if (length(fit$unlagged) || fit$trend) {
    stop("comfac_test() does not support unlagged terms or a trend yet.",
         call. = FALSE)
  }
  lags <- fit$lags
  y <- fit$y
  p <- length(lags[[y]])
  if (p == 0 || !identical(lags[[y]], seq_len(p))) {
    stop("The lags of the dependent variable must be 1, ..., p.", call. = FALSE)
  }
  if (p < r) stop("COMFAC of order ", r, " needs p >= ", r, " (p = ", p, ").",
                  call. = FALSE)
  for (v in fit$x) {
    l <- lags[[v]]
    if (!all(diff(l) == 1) || length(l) < r + 1) {
      stop("Regressor '", v, "' needs at least ", r + 1,
           " consecutive lags for a common factor of order ", r, ".", call. = FALSE)
    }
  }
  n_x <- length(fit$x)
  q_restr <- r * n_x
  X <- stats::model.matrix(fit$lm)
  yv <- stats::model.response(stats::model.frame(fit$lm))
  T <- length(yv)
  b <- stats::coef(fit)

  col_of <- function(v, k) {
    fit$terms$column[fit$terms$variable == v & fit$terms$lag == k &
                       fit$terms$role %in% c("ar", "dl")]
  }

  if (method == "Wald") {
    if (r != 1 || p != 1 || any(vapply(fit$x, function(v) !identical(lags[[v]], 0:1),
                                       logical(1)))) {
      stop("The Wald form is available for r = 1 in models with one lag of y and ",
           "lags 0 and 1 of each regressor; use method = \"LR\".", call. = FALSE)
    }
    V <- adl_vcov(fit, cov_type, maxlags)
    ya <- col_of(y, 1)
    a1 <- b[[ya]]
    g <- numeric(n_x)
    G <- matrix(0, n_x, length(b), dimnames = list(NULL, names(b)))
    for (i in seq_along(fit$x)) {
      c0 <- col_of(fit$x[i], 0); c1 <- col_of(fit$x[i], 1)
      g[i] <- a1 * b[[c0]] + b[[c1]]
      G[i, ya] <- b[[c0]]
      G[i, c0] <- a1
      G[i, c1] <- 1
    }
    stat <- wald_stat(g, G %*% V %*% t(G), diag(n_x), rep(0, n_x))
    p_val <- stats::pchisq(stat, n_x, lower.tail = FALSE)
    inf <- list(statistic = stat, distribution = "chi2", df = n_x, p.value = p_val)
    extra <- list()
    lags_hac <- attr(V, "maxlags")
    ctype <- cov_type
    omega <- NULL
    restricted <- NULL
  } else {
    # Restricted model, by least squares on the same sample: given omega,
    # filtering y and every regressor by omega(L) leaves a model linear in
    # (c, phi, delta), so the RSS is concentrated on omega (r values) and
    # minimized by a grid search followed by local refinement.
    map <- comfac_map(fit, r, colnames(X))
    conc <- comfac_concentrated(fit, r, X, yv)
    om_hat <- comfac_minimize(conc$rss, r)
    ols <- conc$fit(om_hat)
    theta <- numeric(length(map$labels))
    theta[map$omega] <- om_hat
    theta[-map$omega] <- ols$coefficients
    rss_r <- sum(ols$residuals^2)
    rss_u <- sum(stats::residuals(fit$lm)^2)
    omega <- theta[map$omega]
    restricted <- data.frame(parameter = map$labels, estimate = unname(theta))
    stat <- T * log(rss_r / rss_u)
    p_val <- stats::pchisq(stat, q_restr, lower.tail = FALSE)
    inf <- list(statistic = stat, distribution = "chi2", df = q_restr,
                p.value = p_val)
    extra <- list(rss_restricted = rss_r, rss_unrestricted = rss_u)
    lags_hac <- NA_integer_
    ctype <- NULL
  }
  h0 <- if (r == 1 && p == 1 && n_x == 1 && identical(lags[[fit$x]], 0:1)) {
    paste0("common factor of order 1 (alpha_1 * beta_0 + beta_1 = 0)")
  } else paste0("common factor of order ", r, " (", q_restr, " restrictions)")
  new_adl_test(c(list(
    method = paste0("comfac_", tolower(method)),
    null_hypothesis = h0,
    alpha = alpha, reject_null = inf$p.value < alpha, nobs = T,
    cov_type = ctype, maxlags = lags_hac, n_restrictions = q_restr,
    r = r, omega = if (!is.null(omega)) unname(omega), restricted = restricted
  ), inf, extra), paste0("COMFAC test, order ", r, " (", method, ")"))
}

#' Parameter map for the restricted (common factor) model
#'
#' theta = (c, phi_1..phi_{p-r}, delta_i coefficients, omega_1..omega_r).
#' expand(theta) returns the full ADL coefficient vector in the order of the
#' design matrix columns.
#' @noRd
comfac_map <- function(fit, r, cols) {
  y <- fit$y
  p <- length(fit$lags[[y]])
  info <- term_info(fit, cols)
  labels <- character(0)
  idx <- list()
  pos <- 0
  add <- function(lab) {
    labels <<- c(labels, lab)
    pos <<- pos + length(lab)
    (pos - length(lab) + 1):pos
  }
  if (fit$intercept) idx$c <- add("c")
  if (p > r) idx$phi <- add(paste0("phi_", seq_len(p - r)))
  for (v in fit$x) {
    l <- fit$lags[[v]]
    fl <- l[1]
    nd <- length(l) - r
    idx[[paste0("delta_", v)]] <- add(paste0("delta_", v, "_", fl + seq_len(nd) - 1))
  }
  idx$omega <- add(paste0("omega_", seq_len(r)))

  expand <- function(theta) {
    out <- numeric(length(cols))
    om <- c(1, -theta[idx$omega])                 # omega(L)
    ph <- c(1, if (!is.null(idx$phi)) -theta[idx$phi])  # phi(L)
    Cp <- poly_mult(om, ph)                       # C(L) = 1 - a_1 L - ...
    a <- -Cp[-1]
    for (j in seq_len(p)) out[info$role == "ar" & info$lag == j] <- a[j]
    if (fit$intercept) out[info$role == "intercept"] <- theta[idx$c]
    for (v in fit$x) {
      l <- fit$lags[[v]]
      d <- theta[idx[[paste0("delta_", v)]]]
      B <- poly_mult(om, d)                        # coefficients at lags fl, fl+1, ...
      for (k in seq_along(B)) out[info$variable == v & info$role == "dl" &
                                  info$lag == l[1] + k - 1] <- B[k]
    }
    out
  }
  list(expand = expand, labels = labels, omega = idx$omega, idx = idx)
}

poly_mult <- function(a, b) {
  out <- numeric(length(a) + length(b) - 1)
  for (i in seq_along(a)) out[i:(i + length(b) - 1)] <- out[i:(i + length(b) - 1)] + a[i] * b
  out
}

#' Concentrated least squares for the common factor model
#'
#' For a given omega, returns the OLS fit of omega(L) y_t on the filtered
#' regressors, with coefficients in the order (c, phi, delta_i).
#' @noRd
comfac_concentrated <- function(fit, r, X, yv) {
  info <- term_info(fit, colnames(X))
  y <- fit$y
  p <- length(fit$lags[[y]])
  # column of y at lag k (k = 0 is the response)
  ycol <- function(k) if (k == 0) yv else X[, info$role == "ar" & info$lag == k]
  xcol <- function(v, k) X[, info$variable == v & info$role == "dl" & info$lag == k]
  filt <- function(get, k, om) {
    out <- get(k)
    for (j in seq_along(om)) out <- out - om[j] * get(k + j)
    out
  }
  build <- function(om) {
    z <- filt(ycol, 0, om)
    cols <- list()
    if (fit$intercept) cols[["c"]] <- rep(1, length(yv))
    if (p > r) for (k in seq_len(p - r)) cols[[paste0("phi_", k)]] <- filt(ycol, k, om)
    for (v in fit$x) {
      l <- fit$lags[[v]]
      for (k in l[1]:(l[length(l)] - r)) {
        cols[[paste0("delta_", v, "_", k)]] <- filt(function(kk) xcol(v, kk), k, om)
      }
    }
    list(z = z, Z = do.call(cbind, cols))
  }
  fit_om <- function(om) {
    b <- build(om)
    stats::lm.fit(b$Z, b$z)
  }
  list(fit = fit_om, rss = function(om) sum(fit_om(om)$residuals^2))
}

#' Global minimum of the concentrated RSS over omega
#'
#' Grid over the stationary region (r = 1: (-0.99, 0.99); r = 2: the
#' stationarity triangle; r > 2: random stationary starts), then local
#' refinement without bounds from the best few grid points.
#' @noRd
comfac_minimize <- function(rss, r) {
  if (r == 1) {
    g <- seq(-0.99, 0.99, by = 0.01)
    starts <- matrix(g, ncol = 1)
  } else if (r == 2) {
    g <- seq(-1.98, 1.98, by = 0.06)
    h <- seq(-0.98, 0.98, by = 0.04)
    starts <- as.matrix(expand.grid(w1 = g, w2 = h))
    ok <- starts[, 2] + starts[, 1] < 1 & starts[, 2] - starts[, 1] < 1
    starts <- starts[ok, , drop = FALSE]
  } else {
    # random stationary starts, without disturbing the caller's RNG state
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(1)
    starts <- t(replicate(2000, {
      roots <- stats::runif(r, -0.95, 0.95)
      -poly_from_roots(roots)[-1]
    }))
  }
  vals <- apply(starts, 1, function(om) rss(om))
  best <- starts[order(vals)[seq_len(min(5, nrow(starts)))], , drop = FALSE]
  sols <- lapply(seq_len(nrow(best)), function(i) {
    if (r == 1) {
      o <- stats::optim(best[i, ], rss, method = "BFGS",
                        control = list(reltol = 1e-14, maxit = 1000))
    } else {
      o <- stats::optim(best[i, ], rss, method = "Nelder-Mead",
                        control = list(reltol = 1e-14, maxit = 5000))
      o <- stats::optim(o$par, rss, method = "BFGS",
                        control = list(reltol = 1e-14, maxit = 1000))
    }
    o
  })
  conv <- vapply(sols, function(o) o$convergence == 0, logical(1))
  if (!any(conv)) {
    warning("The optimizer did not report convergence for the restricted model; ",
            "the LR statistic may be inaccurate.", call. = FALSE)
    conv[] <- TRUE
  }
  vals <- vapply(sols, `[[`, numeric(1), "value")
  vals[!conv] <- Inf
  sols[[which.min(vals)]]$par
}

#' Coefficients of prod_j (1 - r_j L), as c(1, -w_1, ..., -w_r)
#' @noRd
poly_from_roots <- function(roots) {
  out <- 1
  for (rt in roots) out <- poly_mult(out, c(1, -rt))
  out
}
