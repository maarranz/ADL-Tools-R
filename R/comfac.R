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
#' * `method = "Wald"`: Sargan's (1980) determinantal form, for any order
#'   \eqn{r} and any number of regressors, with the unrestricted ADL only
#'   (no restricted model is estimated). The polynomials \eqn{C(z)} and
#'   \eqn{B_i(z)} (without leading zeros) have a common factor of degree
#'   \eqn{r} if and only if a matrix of their shifted coefficients has a
#'   given rank; Sargan's elimination algorithm (sections 3-4 of the paper:
#'   the matrices \eqn{\Psi(k)} for his Type 1 cases and
#'   \eqn{\Psi_b} for Type 2, with "vanishing rows") turns that into
#'   exactly \eqn{r n} determinants that must be zero. They are tested
#'   jointly by the delta method, with exact derivatives (sums of
#'   cofactors). For the ADL(1,1) the single determinant is
#'   \eqn{\alpha_1 \beta_0 + \beta_1}, the restriction of ADL2; for one
#'   regressor and \eqn{r = 1} it is the resultant of \eqn{C} and \eqn{B}.
#'   The determinants are returned in `constraints`.
#'
#' The Wald test does not need the restricted model, which may have
#' several local optima when \eqn{r} is too small or the common roots are
#' complex (Sargan, 1980, section 5). But Wald tests of nonlinear restrictions
#' depend on how the restrictions are written (Gregory and Veall, 1985): the
#' statistic here follows Sargan's construction (polynomials ordered by
#' degree, \eqn{C} first among equals), and another program may write the
#' same hypothesis differently. In small samples it tends to under-reject
#' (about 3% at the 5% level with T = 200 in the simulations recorded in
#' DESIGN.md; close to 5% with T = 1000).
#'
#' **Choosing r.** Sargan (1980, section 5) shows that if the true order is
#' \eqn{\bar r}, the Wald test is asymptotically valid for
#' \eqn{r \ge \bar r}, but for \eqn{r < \bar r} the constraints and their
#' derivatives vanish together, the statistic is not
#' \eqn{\chi^2}, and it is biased downward. Differences between successive
#' Wald statistics are therefore not \eqn{\chi^2} either. Compute the tests
#' for \eqn{r = 1, 2, \dots} and read them together, as in Sargan's
#' Table I, rather than as a sequence of independent tests.
#'
#' Intercepts are left free; models with unlagged terms or a trend are not
#' supported yet.
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
#'   model and holds `constraints` (the values of Sargan's determinants).
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
#'
#' # Sargan's Wald tests for r = 1, 2, read together
#' sapply(1:2, function(r) comfac_test(fit2, r = r, method = "Wald")$statistic)
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
    V <- adl_vcov(fit, cov_type, maxlags)
    st <- sargan_structure(fit, r, names(b))
    con <- sargan_constraints(b, st)
    stat <- wald_stat(con$phi, con$G %*% V %*% t(con$G), diag(length(con$phi)),
                      rep(0, length(con$phi)))
    p_val <- stats::pchisq(stat, q_restr, lower.tail = FALSE)
    inf <- list(statistic = stat, distribution = "chi2", df = q_restr, p.value = p_val)
    extra <- list(constraints = con$phi)
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

# Sargan's determinantal Wald form ---------------------------------------------
#
# Sargan (1980, sections 3-4 and appendix B). The polynomials psi_j(z) --
# C(z) and each B_i(z) without its leading zeros -- of degrees
# f = f_1 >= f_2 >= ... >= f_n have a common factor of degree r iff a matrix
# of shifted coefficient rows has a given rank. The rows are:
#   Type 1 (sum_j (f_j - r) > (n - 1)(f - r)): Psi(k) of eq. (8), i.e. Psi_0
#     (psi_1, and psi_j shifted 0..d_j places, d_j = f - f_j) followed by k
#     blocks of the right-aligned Psi_a shifted 1..k places, with k the
#     smallest integer such that (n - 1) k >= sum_j (f_j - r) - (n - 1)(f - r),
#     keeping the first f + k + n - r rows;
#   Type 2 (otherwise): the rows of Psi_0 reordered as Psi_b (by decreasing
#     number of right zeros, psi_1 before psi_2 ... within a group), with
#     "surplus vanishing rows" dropped.
# Elimination (section 4): rows are taken in order; a row reduced by the
# pivots so far that has only r structurally nonzero elements is a
# "vanishing row" and is not used as a pivot; its r elements are the
# (algorithmic) constraints. Multiplied by the pivots, they are the
# determinantal constraints used here: the determinant of the pivot rows plus
# the vanishing row, on the pivot columns plus the element's column. There
# are r(n - 1) of them. Every entry is 0, 1 or +/- one coefficient, so the
# derivatives are exact sums of cofactors.

#' Sargan's rows and the determinants that define the constraints
#' @noRd
sargan_structure <- function(fit, r, cols) {
  info <- term_info(fit, cols)
  polys <- list()
  p <- length(fit$lags[[fit$y]])
  ar <- vapply(seq_len(p), function(j) which(info$role == "ar" & info$lag == j), 1L)
  polys[[fit$y]] <- list(idx = c(-1L, ar), sign = c(1, rep(-1, p)))   # -1: the constant 1
  for (v in fit$x) {
    l <- fit$lags[[v]]
    dl <- vapply(l, function(j) which(info$variable == v & info$role == "dl" &
                                        info$lag == j), 1L)
    polys[[v]] <- list(idx = dl, sign = rep(1, length(l)))
  }
  deg <- vapply(polys, function(q) length(q$idx) - 1L, 1L)
  ord <- order(-deg, seq_along(deg))
  polys <- polys[ord]
  deg <- unname(deg[ord])
  n <- length(polys)
  f <- deg[1]
  d <- f - deg
  # a row: polynomial j placed after `off` zeros
  row_list <- function(j, off) list(j = j, off = off)
  psi0 <- list(row_list(1L, 0L))
  for (j in seq_len(n)[-1]) for (o in 0:d[j]) psi0[[length(psi0) + 1]] <- row_list(j, o)
  type1 <- sum(deg - r) > (n - 1) * (f - r)
  if (type1) {
    need <- sum(deg - r) - (n - 1) * (f - r)
    k <- max(0L, as.integer(ceiling(need / (n - 1))))
    rows <- psi0
    for (s in seq_len(k)) for (j in seq_len(n)) rows[[length(rows) + 1]] <- row_list(j, d[j] + s)
    rows <- rows[seq_len(f + k + n - r)]
    ncol <- f + k + 1
  } else {
    k <- 0L
    rz <- vapply(psi0, function(w) f - deg[w$j] - w$off, 1L)
    jj <- vapply(psi0, `[[`, 1L, "j")
    rows <- psi0[order(-rz, jj)]
    ncol <- f + 1
  }
  # structural elimination: decide pivots and vanishing rows, tracking the
  # columns where each reduced row can be nonzero (a reduction by a pivot
  # brings in the pivot's columns)
  piv <- list()                                    # row, column, support
  dets <- list()
  vanished <- list()
  for (i in seq_along(rows)) {
    w <- rows[[i]]
    # surplus: a right-translated copy of a vanishing row (type 2)
    if (!type1 && any(vapply(vanished, function(v) v$j == w$j && v$off < w$off,
                             logical(1)))) next
    sup <- w$off + seq_len(deg[w$j] + 1)
    used <- 0L
    frozen <- length(sup) <= r
    for (t in seq_along(piv)) {
      if (frozen) break
      used <- t
      if (piv[[t]]$col %in% sup) sup <- setdiff(union(sup, piv[[t]]$sup), piv[[t]]$col)
      frozen <- length(sup) <= r
    }
    if (!frozen) {
      piv[[length(piv) + 1]] <- list(row = i, col = min(sup), sup = sort(sup))
      next
    }
    if (length(sup) != r) {
      stop("Sargan's construction failed (a row with ", length(sup),
           " elements); use method = \"LR\".", call. = FALSE)
    }
    prow <- vapply(piv[seq_len(used)], `[[`, 1L, "row")
    pcol <- vapply(piv[seq_len(used)], `[[`, 1L, "col")
    for (cc in sort(sup)) {
      dets[[length(dets) + 1]] <- list(rows = c(prow, i), cols = c(pcol, cc))
    }
    vanished[[length(vanished) + 1]] <- w
  }
  if (length(dets) != r * (n - 1)) {
    stop("Sargan's construction gave ", length(dets), " constraints instead of ",
         r * (n - 1), "; use method = \"LR\".", call. = FALSE)
  }
  idx <- matrix(0L, length(rows), ncol)
  sgn <- matrix(0, length(rows), ncol)
  for (i in seq_along(rows)) {
    w <- rows[[i]]
    q <- polys[[w$j]]
    pos <- w$off + seq_along(q$idx)
    idx[i, pos] <- q$idx
    sgn[i, pos] <- q$sign
  }
  list(idx = idx, sgn = sgn, dets = dets, n = n, r = r, f = f, k = k,
       type = if (type1) 1L else 2L, order = names(polys))
}

#' Constraint values and their exact Jacobian
#' @noRd
sargan_constraints <- function(b, st) {
  M <- matrix(0, nrow(st$idx), ncol(st$idx))
  M[st$idx == -1L] <- 1
  pos <- st$idx > 0L
  M[pos] <- st$sgn[pos] * b[st$idx[pos]]
  phi <- numeric(length(st$dets))
  G <- matrix(0, length(st$dets), length(b))
  for (h in seq_along(st$dets)) {
    rr <- st$dets[[h]]$rows
    cl <- st$dets[[h]]$cols
    A <- M[rr, cl, drop = FALSE]
    phi[h] <- det(A)
    I <- st$idx[rr, cl, drop = FALSE]
    S <- st$sgn[rr, cl, drop = FALSE]
    for (a in seq_len(nrow(A))) for (e in seq_len(ncol(A))) {
      j <- I[a, e]
      if (j > 0L) {
        minor <- A[-a, -e, drop = FALSE]
        cof <- (-1)^(a + e) * (if (length(minor)) det(minor) else 1)
        G[h, j] <- G[h, j] + S[a, e] * cof
      }
    }
  }
  list(phi = phi, G = G)
}
