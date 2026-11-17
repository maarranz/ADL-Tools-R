# Super exogeneity (ADL3; Engle and Hendry, 1993) ----------------------------

#' Super exogeneity test
#'
#' Super exogeneity is weak exogeneity plus invariance: the parameters of the
#' conditional (ADL) model do not change when the process generating
#' \eqn{x} changes, which is what policy analysis needs (the econometric
#' answer to the Lucas critique). Following the outline in ADL3 (Engle and
#' Hendry, 1993), the test needs dates at which the marginal process of
#' \eqn{x} did change, summarized by impulse dummies (1 at the date, 0
#' elsewhere) and step dummies (0 before the date, 1 from it on):
#'
#' 1. **Breaks in the marginal model.** The dummies are added to the
#'    regression for \eqn{\Delta x_t} of [weak_exogeneity_test()] and tested
#'    jointly. They should be significant: if the marginal process did not
#'    change at these dates, the test has nothing to work with.
#' 2. **Invariance of the conditional model.** The same dummies are added to
#'    the ADL equation and tested jointly. Under super exogeneity they are
#'    insignificant: the conditional model stays constant across the episodes
#'    that changed the marginal model. This is the main result.
#' 3. **Weak exogeneity.** \eqn{\kappa = 0} in the marginal regression with
#'    the dummies (super exogeneity requires it).
#'
#' Marginal breaks with an invariant conditional model and a weakly exogenous
#' \eqn{x} support super exogeneity. If the conditional model breaks at the
#' same dates, invariance, and any policy simulation built on the model, is
#' not credible. If the marginal model shows no break, not rejecting
#' invariance is uninformative, and the function warns.
#'
#' The dates are a modelling choice: find them with the stability tests on
#' the marginal model (e.g. [hansen_test()] or [cusum_test()] on an ADL for
#' \eqn{x}, or [chow_break_test()]), from the history of the policy, or from
#' the residuals of the marginal regression. Dates chosen by looking at the
#' conditional model would bias the test towards invariance.
#'
#' @param fit An `"adl"` object.
#' @param impulses,steps Dates of the impulse and step dummies, in the units
#'   of the time index: one date (e.g. `c(2008, 4)` for quarterly data, or a
#'   `Date`), or several as a list (`list(c(2008, 4), c(2009, 1))`) or as a
#'   vector of dates. At least one dummy is needed.
#' @param x The regressor whose marginal process is examined; by default the
#'   only regressor.
#' @param lags,x_level The marginal regression, as in
#'   [weak_exogeneity_test()].
#' @inheritParams adl_vcov
#' @param dist `"chisq"` (default) or `"F"`.
#' @param alpha Significance level for the decisions.
#' @return An object of class `"adltools_test"` whose main result is the
#'   invariance test (step 2), with a `coefficients` table of the three steps
#'   and `marginal` and `conditional` (the augmented regressions, `lm`
#'   objects) and `dummies` (their names and dates) added.
#' @references Engle, R. F. and Hendry, D. F. (1993). Testing super
#'   exogeneity and invariance in regression models. *Journal of
#'   Econometrics*, 56(1-2), 119-139.
#'
#'   Engle, R. F., Hendry, D. F. and Richard, J.-F. (1983). Exogeneity.
#'   *Econometrica*, 51(2), 277-304.
#'
#'   Lucas, R. E. (1976). Econometric policy evaluation: a critique.
#'   *Carnegie-Rochester Conference Series on Public Policy*, 1, 19-46.
#' @seealso [weak_exogeneity_test()], [hansen_test()], [chow_break_test()]
#' @examples
#' # x shifts its mean in 2005 Q1; y given x is unchanged
#' set.seed(1)
#' n <- 120
#' x <- numeric(n); y <- numeric(n)
#' for (t in 2:n) {
#'   x[t] <- 0.5 * x[t - 1] + 2 * (t >= 61) + rnorm(1)
#'   y[t] <- 1 + 0.5 * y[t - 1] + 0.4 * x[t] + 0.3 * x[t - 1] + rnorm(1)
#' }
#' d <- ts(cbind(y = y, x = x), start = c(1990, 1), frequency = 4)
#' fit <- adl(y ~ x, data = d, order = c(1, 1))
#' super_exogeneity_test(fit, steps = c(2005, 1))
#' @export
super_exogeneity_test <- function(fit, impulses = NULL, steps = NULL, x = NULL,
                                  lags = 1, x_level = TRUE,
                                  cov_type = c("nonrobust", "HAC"),
                                  maxlags = NULL, dist = c("chisq", "F"),
                                  alpha = 0.05) {
  check_adl(fit)
  cov_type <- match.arg(cov_type)
  dist <- match.arg(dist)
  check_alpha(alpha)
  check_count(lags, "lags")
  if (is.null(x)) {
    if (length(fit$x) != 1) stop("Several regressors: choose one with `x`.", call. = FALSE)
    x <- fit$x
  }
  if (!x %in% fit$x) stop("'", x, "' is not a regressor of the model.", call. = FALSE)
  imp <- split_dates(impulses, fit$index)
  stp <- split_dates(steps, fit$index)
  if (length(imp) + length(stp) == 0) {
    stop("Give the dates of the breaks in the marginal process of x with ",
         "`impulses` and/or `steps`.", call. = FALSE)
  }

  # 1. marginal model with the dummies
  mf <- marginal_frame(fit, x, lags, x_level)
  Dm <- break_dummies(imp, stp, attr(mf, "index"), "the marginal regression")
  dnames <- colnames(Dm)
  # syntactic internal names (lm would backtick "S(2005 Q1)")
  dcol_m <- paste0("dummy.", seq_along(dnames))
  colnames(Dm) <- dcol_m
  mdat <- cbind(mf, Dm)
  marg <- stats::lm(dx ~ ., data = mdat)
  # 2. conditional model with the dummies
  X <- stats::model.matrix(fit$lm)
  yv <- stats::model.response(stats::model.frame(fit$lm))
  Dc <- break_dummies(imp, stp, fit$index, "the estimation sample")
  dcol_c <- dcol_m
  Xd <- cbind(X, as.matrix(Dc))
  colnames(Xd) <- c(colnames(X), dcol_c)
  cond <- stats::lm(yv ~ 0 + Xd)
  names(cond$coefficients) <- colnames(Xd)

  test_cols <- function(model, cols) {
    b <- stats::coef(model)
    if (anyNA(b[cols])) {
      stop("A dummy is collinear with the other regressors (e.g. a step at ",
           "the first observation, or two dummies that coincide in the sample).",
           call. = FALSE)
    }
    n <- stats::nobs(model)
    lg <- hac_lags(maxlags, cov_type, n)
    V <- if (cov_type == "HAC") {
      sandwich::NeweyWest(model, lag = lg, prewhite = FALSE, adjust = FALSE)
    } else stats::vcov(model)
    q <- length(cols)
    R <- matrix(0, q, length(b))
    R[cbind(seq_len(q), match(cols, names(b)))] <- 1
    w <- wald_stat(b, V, R, rep(0, q))
    c(wald_inference(w, q, stats::df.residual(model), dist), list(maxlags = lg))
  }
  t_marg <- test_cols(marg, dcol_m)
  t_cond <- test_cols(cond, dcol_c)
  t_weak <- test_cols(marg, "ec_lag1")

  steps_tab <- data.frame(
    step = c("1. breaks in the marginal model of x",
             "2. invariance of the conditional model",
             "3. weak exogeneity (kappa = 0)"),
    wanted = c("reject", "do not reject", "do not reject"),
    statistic = c(t_marg$statistic, t_cond$statistic, t_weak$statistic),
    df = c(length(dnames), length(dnames), 1L),
    p.value = c(t_marg$p.value, t_cond$p.value, t_weak$p.value)
  )
  steps_tab$reject <- steps_tab$p.value < alpha

  if (!steps_tab$reject[1]) {
    warning("The dummies are not significant in the marginal model of ", x,
            " (p = ", format(signif(t_marg$p.value, 3)), "): the marginal ",
            "process shows no change at these dates, so not rejecting ",
            "invariance says nothing about super exogeneity.", call. = FALSE)
  }
  if (steps_tab$reject[3]) {
    warning("Weak exogeneity is rejected (p = ", format(signif(t_weak$p.value, 3)),
            "), and super exogeneity requires it.", call. = FALSE)
  }

  new_adl_test(c(list(
    method = "super_exogeneity",
    null_hypothesis = paste0("the ADL for ", fit$y, " is invariant to the breaks ",
                             "in the marginal process of ", x, " (",
                             length(dnames), if (length(dnames) == 1) " dummy)" else " dummies)"),
    alpha = alpha, reject_null = t_cond$p.value < alpha, nobs = stats::nobs(cond),
    cov_type = cov_type, maxlags = t_cond$maxlags,
    n_restrictions = length(dnames), df_resid = stats::df.residual(cond),
    coefficients = steps_tab,
    dummies = data.frame(name = dnames,
                         type = rep(c("impulse", "step"), c(length(imp), length(stp))),
                         date = vapply(c(imp, stp), function(v)
                           format(as_index_value(v, fit$index)), character(1))),
    marginal = marg, conditional = cond
  ), t_cond[setdiff(names(t_cond), "maxlags")]), "Super exogeneity test")
}

#' Dates given as one date, a list, or a vector of dates
#' @noRd
split_dates <- function(d, idx) {
  if (is.null(d)) return(list())
  if (is.list(d)) return(d)
  periodic <- inherits(idx, c("yearqtr", "yearmon")) ||
    (is.numeric(idx) && stats::frequency(zoo::zoo(seq_along(idx), idx)) > 1)
  if (is.numeric(d) && length(d) == 2 && periodic) return(list(d))
  as.list(d)
}

#' Impulse and step dummies on an index
#' @noRd
break_dummies <- function(imp, stp, idx, where) {
  num <- as.numeric(idx)
  pos <- function(v) {
    val <- as.numeric(as_index_value(v, idx))
    k <- which(abs(num - val) < 1e-8 * max(1, abs(val)))
    if (!length(k)) {
      stop("The date ", format(as_index_value(v, idx)), " is not in ", where,
           " (", format(idx[1]), " - ", format(idx[length(idx)]), ").", call. = FALSE)
    }
    k
  }
  n <- length(idx)
  cols <- c(lapply(imp, function(v) as.numeric(seq_len(n) == pos(v))),
            lapply(stp, function(v) as.numeric(seq_len(n) >= pos(v))))
  nm <- c(vapply(imp, function(v) paste0("I(", format(idx[pos(v)]), ")"), character(1)),
          vapply(stp, function(v) paste0("S(", format(idx[pos(v)]), ")"), character(1)))
  if (anyDuplicated(nm)) stop("A dummy date is repeated.", call. = FALSE)
  names(cols) <- nm
  as.data.frame(cols, check.names = FALSE)
}
