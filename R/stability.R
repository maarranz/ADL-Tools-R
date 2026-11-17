# Recursive stability tests and Hansen's (1992) test (ADL3) -----------------
#
# CUSUM and CUSUM-of-squares (Brown, Durbin and Evans, 1975) are built from
# recursive residuals; Hansen's test from the score contributions of the
# full-sample OLS estimates. None of them takes a HAC option (ADL3).

#' Recursive residuals
#'
#' The standardized one-step-ahead prediction errors of Brown, Durbin and
#' Evans (1975): for each \eqn{t > k}, the model is estimated on the first
#' \eqn{t-1} observations and
#' \deqn{w_t = \frac{y_t - x_t'\hat b_{t-1}}{\sqrt{1 + x_t'(X_{t-1}'X_{t-1})^{-1}x_t}}.}
#' Under the classical assumptions they are independent \eqn{N(0, \sigma^2)},
#' and their sum of squares equals the residual sum of squares of the
#' full-sample fit.
#'
#' @param fit An `"adl"` object.
#' @return A data frame with the time index and the recursive residual, one
#'   row per observation after the first \eqn{k} (\eqn{k} = number of
#'   coefficients). If the first observations do not identify the model
#'   (e.g. a dummy that is still zero), the recursion starts later.
#' @references Brown, R. L., Durbin, J. and Evans, J. M. (1975). Techniques
#'   for testing the constancy of regression relationships over time.
#'   *Journal of the Royal Statistical Society B*, 37(2), 149-192.
#' @examples
#' fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
#' head(recursive_residuals(fit))
#' @export
recursive_residuals <- function(fit) {
  check_adl(fit)
  X <- stats::model.matrix(fit$lm)
  y <- stats::model.response(stats::model.frame(fit$lm))
  n <- nrow(X)
  k <- ncol(X)
  w <- rep(NA_real_, n)
  for (t in seq_len(n)[-seq_len(k)]) {
    X1 <- X[seq_len(t - 1), , drop = FALSE]
    q <- qr(X1)
    if (q$rank < k) next
    b <- qr.coef(q, y[seq_len(t - 1)])
    xt <- X[t, ]
    h <- sum(xt * solve(crossprod(X1), xt))
    w[t] <- (y[t] - sum(xt * b)) / sqrt(1 + h)
  }
  ok <- !is.na(w)
  if (sum(ok) < 2) stop("Too few observations for recursive residuals.", call. = FALSE)
  first <- which(ok)[1]
  data.frame(index = fit$index[first:n], residual = w[first:n])
}

#' CUSUM test
#'
#' Brown, Durbin and Evans' (1975) test on the cumulative sum of the recursive
#' residuals, scaled by the standard error of the regression,
#' \deqn{W_r = \frac{1}{\hat\sigma}\sum_{j=1}^{r} w_{k+j}, \qquad r = 1, \ldots, m = T-k.}
#' A stable model keeps \eqn{W_r} inside the bounds
#' \eqn{\pm a(\sqrt{m} + 2r/\sqrt{m})}, the lines EViews draws; drift
#' outside them signals instability, typically a gradual change in the
#' coefficients. The statistic is \eqn{\max_r |W_r| / (\sqrt m + 2r/\sqrt m)},
#' compared with \eqn{a} (0.948 at 5%), with the asymptotic boundary-crossing
#' probability as p-value.
#'
#' @inheritParams recursive_residuals
#' @param alpha Significance level for the bounds and the decision.
#' @return An object of class `"adltools_test"` (and `"adltools_stability"`,
#'   which has a `plot()` method), with `path` (a data frame with the index,
#'   \eqn{W_r} and the bounds), `crossing` (first date outside the bounds, or
#'   `NA`) and `critical` (\eqn{a}) added.
#' @references Brown, R. L., Durbin, J. and Evans, J. M. (1975). Techniques
#'   for testing the constancy of regression relationships over time.
#'   *Journal of the Royal Statistical Society B*, 37(2), 149-192.
#' @seealso [cusumsq_test()], [hansen_test()], [chow_break_test()]
#' @examples
#' fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
#' cs <- cusum_test(fit)
#' cs
#' plot(cs)
#' @export
cusum_test <- function(fit, alpha = 0.05) {
  check_alpha(alpha)
  rr <- recursive_residuals(fit)
  w <- rr$residual
  m <- length(w)
  s <- sqrt(sum(w^2) / m)
  W <- cumsum(w) / s
  r <- seq_len(m)
  scale <- sqrt(m) + 2 * r / sqrt(m)
  stat <- max(abs(W) / scale)
  a <- stats::uniroot(function(a) cusum_pvalue(a) - alpha, c(0.1, 5), tol = 1e-10)$root
  out <- abs(W) > a * scale
  res <- new_adl_test(list(
    method = "cusum",
    null_hypothesis = "constant coefficients (CUSUM of recursive residuals)",
    statistic = stat, distribution = "cusum", p.value = cusum_pvalue(stat),
    alpha = alpha, reject_null = stat > a, nobs = nrow(stats::model.matrix(fit$lm)),
    critical = a,
    crossing = if (any(out)) format(rr$index[which(out)[1]]) else NA_character_,
    path = data.frame(index = rr$index, value = W, lower = -a * scale,
                      upper = a * scale)
  ), "CUSUM test")
  class(res) <- c("adltools_stability", class(res))
  res
}

#' Asymptotic p-value of the CUSUM statistic (Brown, Durbin and Evans, 1975)
#' @noRd
cusum_pvalue <- function(a) {
  min(1, 2 * (1 - stats::pnorm(3 * a) + exp(-4 * a^2) * stats::pnorm(a)))
}

#' CUSUM-of-squares test
#'
#' Brown, Durbin and Evans' (1975) test on the cumulative sum of the squared
#' recursive residuals,
#' \deqn{S_r = \sum_{j=1}^{r} w_{k+j}^2 \Big/ \sum_{j=1}^{m} w_{k+j}^2,}
#' which rises from 0 to 1 and, for a stable model, stays within
#' \eqn{\pm c_0} of the line \eqn{r/m}. It is more sensitive than CUSUM to a
#' change in the error variance or a sharp break.
#'
#' The statistic is \eqn{\max_r |S_r - r/m|}. Its critical value \eqn{c_0}
#' and p-value come from the exact distribution under independent normal
#' recursive residuals, by simulation (`nsim` draws, with a fixed seed so
#' results are reproducible; the caller's random numbers are not affected).
#' EViews takes \eqn{c_0} from Durbin's (1969) table, an approximation to
#' the same distribution, so the bounds can differ slightly.
#'
#' @inheritParams cusum_test
#' @param nsim Number of simulated samples for the critical value and the
#'   p-value.
#' @return An object of class `"adltools_test"` (and `"adltools_stability"`),
#'   with `path`, `crossing` and `critical` (\eqn{c_0}) added.
#' @references Brown, R. L., Durbin, J. and Evans, J. M. (1975). Techniques
#'   for testing the constancy of regression relationships over time.
#'   *Journal of the Royal Statistical Society B*, 37(2), 149-192.
#'
#'   Durbin, J. (1969). Tests for serial correlation in regression analysis
#'   based on the periodogram of least-squares residuals. *Biometrika*,
#'   56(1), 1-15.
#' @seealso [cusum_test()], [hansen_test()]
#' @examples
#' fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
#' cq <- cusumsq_test(fit)
#' cq
#' plot(cq)
#' @export
cusumsq_test <- function(fit, alpha = 0.05, nsim = 10000) {
  check_alpha(alpha)
  check_count(nsim, "nsim")
  rr <- recursive_residuals(fit)
  w2 <- rr$residual^2
  m <- length(w2)
  S <- cumsum(w2) / sum(w2)
  line <- seq_len(m) / m
  stat <- max(abs(S - line))
  sim <- with_seed(20261117, vapply(seq_len(nsim), function(i) {
    e2 <- stats::rnorm(m)^2
    max(abs(cumsum(e2) / sum(e2) - line))
  }, numeric(1)))
  c0 <- unname(stats::quantile(sim, 1 - alpha, type = 1))
  out <- abs(S - line) > c0
  res <- new_adl_test(list(
    method = "cusumsq",
    null_hypothesis = "constant coefficients and variance (CUSUM of squares)",
    statistic = stat, distribution = "cusumsq",
    p.value = (1 + sum(sim >= stat)) / (nsim + 1),
    alpha = alpha, reject_null = stat > c0, nobs = nrow(stats::model.matrix(fit$lm)),
    critical = c0,
    crossing = if (any(out)) format(rr$index[which(out)[1]]) else NA_character_,
    path = data.frame(index = rr$index, value = S, lower = line - c0,
                      upper = line + c0, expected = line)
  ), "CUSUM-of-squares test")
  class(res) <- c("adltools_stability", class(res))
  res
}

#' Run an expression with a fixed seed, restoring the caller's RNG state
#' @noRd
with_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
  on.exit({
    if (is.null(old)) {
      rm(".Random.seed", envir = globalenv())
    } else {
      assign(".Random.seed", old, envir = globalenv())
    }
  })
  set.seed(seed)
  expr
}

#' Hansen's (1992) parameter stability test
#'
#' Tests the constancy of each coefficient and of the error variance, and of
#' all of them jointly, against instability of unknown timing, without a
#' break date and without re-estimating the model. With OLS residuals
#' \eqn{\hat\epsilon_t} and \eqn{\hat\sigma^2 = RSS/T}, the score
#' contributions are \eqn{f_{it} = x_{it}\hat\epsilon_t} for each
#' coefficient and \eqn{f_t = \hat\epsilon_t^2 - \hat\sigma^2} for the
#' variance; with \eqn{S_{it} = \sum_{j \le t} f_{ij}},
#' \deqn{L_i = \frac{1}{T V_i}\sum_{t=1}^{T} S_{it}^2, \quad V_i = \sum_t f_{it}^2,
#' \qquad L_c = \frac{1}{T}\sum_{t=1}^{T} S_t' V^{-1} S_t, \quad V = \sum_t f_t f_t'.}
#' Large values reject constancy. Each \eqn{L_i} tells which parameter is
#' unstable, which CUSUM and Chow tests cannot.
#'
#' p-values are from the asymptotic distribution, \eqn{\int_0^1 B(r)'B(r)\,dr}
#' for a Brownian bridge \eqn{B} of dimension 1 (individual) or \eqn{k+1}
#' (joint). Hansen (1992) tabulates simulated critical values (0.470 at 5%
#' for one parameter), very close to the asymptotic ones used here (0.461).
#' The asymptotic theory assumes stationary regressors.
#'
#' @inheritParams cusum_test
#' @return An object of class `"adltools_test"` with the joint test, and in
#'   `coefficients` a data frame with \eqn{L_i} and its p-value for each
#'   coefficient and for the variance.
#' @references Hansen, B. E. (1992). Testing for parameter instability in
#'   linear models. *Journal of Policy Modeling*, 14(4), 517-533.
#'
#'   Nyblom, J. (1989). Testing for the constancy of parameters over time.
#'   *Journal of the American Statistical Association*, 84(405), 223-230.
#' @seealso [cusum_test()], [chow_break_test()]
#' @examples
#' fit <- adl(y ~ x, data = adl_sim, order = c(1, 1))
#' hansen_test(fit)
#' @export
hansen_test <- function(fit, alpha = 0.05) {
  check_adl(fit)
  check_alpha(alpha)
  X <- stats::model.matrix(fit$lm)
  e <- stats::residuals(fit$lm)
  n <- length(e)
  s2 <- sum(e^2) / n
  f <- cbind(X * e, e^2 - s2)
  colnames(f) <- c(colnames(X), "variance")
  S <- apply(f, 2, cumsum)
  V <- crossprod(f)
  Li <- colSums(S^2) / (n * diag(V))
  if (qr(V)$rank < ncol(V)) {
    stop("The score covariance is singular; the joint test is not defined.",
         call. = FALSE)
  }
  Lc <- sum(S * t(solve(V, t(S)))) / n
  d <- ncol(f)
  pi <- vapply(Li, cvm_pvalue, numeric(1), d = 1)
  pc <- cvm_pvalue(Lc, d)
  param <- c(param_labels(fit, term_info(fit, colnames(X))), "sigma^2")
  new_adl_test(list(
    method = "hansen",
    null_hypothesis = paste0("all ", d, " parameters constant (",
                             d - 1, " coefficients and the variance)"),
    statistic = Lc, distribution = "hansen", df = d, p.value = pc,
    alpha = alpha, reject_null = pc < alpha, nobs = n,
    coefficients = data.frame(term = colnames(f), param = param,
                              statistic = unname(Li), p.value = unname(pi),
                              reject = unname(pi < alpha))
  ), "Hansen (1992) parameter stability test")
}

#' P(Q > x) for Q = sum_j chi2_d / (j pi)^2, the Cramer-von Mises
#' distribution of int B'B for a d-dimensional Brownian bridge, by Imhof's
#' (1961) inversion of the characteristic function. The first J terms are
#' exact; the tail, whose variance is negligible, enters through its mean.
#' @noRd
cvm_pvalue <- function(x, d, J = 200) {
  lam <- 1 / ((seq_len(J) * pi)^2)
  tail_mean <- d * (1 / 6 - sum(lam))           # sum_j 1/(j pi)^2 = 1/6
  xs <- x - tail_mean
  if (xs <= 0) return(1)
  integrand <- function(u) {
    vapply(u, function(v) {
      th <- 0.5 * d * sum(atan(lam * v)) - 0.5 * xs * v
      rho <- exp(0.25 * d * sum(log1p((lam * v)^2)))
      sin(th) / (v * rho)
    }, numeric(1))
  }
  val <- stats::integrate(integrand, 0, Inf, subdivisions = 2000L,
                          rel.tol = 1e-8)$value
  min(1, max(0, 0.5 + val / pi))
}

#' Plot a CUSUM or CUSUM-of-squares test
#'
#' Draws the statistic over time with its bounds at the chosen significance
#' level (and, for CUSUM of squares, the expected line).
#'
#' @param x A result of [cusum_test()] or [cusumsq_test()].
#' @param ... Passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
plot.adltools_stability <- function(x, ...) {
  p <- x$path
  idx <- p$index
  if (!inherits(idx, c("Date", "POSIXt"))) idx <- as.numeric(idx)
  sq <- identical(x$method, "cusumsq")
  rng <- range(c(p$value, p$lower, p$upper))
  graphics::plot(idx, p$value, type = "l", lwd = 2, ylim = rng,
                 xlab = "", ylab = if (sq) "CUSUM of squares" else "CUSUM",
                 main = paste0(attr(x, "title"), " (", format(100 * x$alpha),
                               "% bounds)"), ...)
  graphics::lines(idx, p$lower, lty = 2, col = "red")
  graphics::lines(idx, p$upper, lty = 2, col = "red")
  if (sq) graphics::lines(idx, p$expected, lty = 3) else graphics::abline(h = 0, lty = 3)
  invisible(x)
}
