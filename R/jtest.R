# J-test for non-nested models (ADL3) -----------------------------------------
#
# Davidson and MacKinnon (1981): add the fitted values of the rival model to
# each model and test their coefficient, in both directions, on the sample
# common to the two models. The four outcomes are read as in ADL3.

#' J-test of two non-nested ADL models
#'
#' The J-test of Davidson and MacKinnon (1981), as in ADL3. Model A is
#' augmented with the fitted values of model B,
#' \deqn{y_t = X_{A,t}'\gamma_A + \theta_B\, \hat y_{B,t} + w_t,}
#' and \eqn{H_0: \theta_B = 0} is tested (A encompasses B); then the roles
#' are reversed. Both directions are run, and the pair of decisions gives one
#' of four outcomes:
#'
#' | Rejected | Reading |
#' |:--|:--|
#' | neither | both models are compatible with the data; the test cannot discriminate between them |
#' | only "A encompasses B" | evidence against A, in favour of B |
#' | only "B encompasses A" | evidence against B, in favour of A |
#' | both | neither model accounts for the other: both are mis-specified; consider a more general model |
#'
#' Both models are refitted on the observations they have in common (the
#' longer lag structure loses initial observations), and must have the same
#' dependent variable, so that the fitted values are comparable. The test is
#' for **non-nested** models: if one model's regressors span the other's, its
#' fitted values are a combination of the larger model's regressors and
#' \eqn{\theta} is not identified; then the restrictions are tested directly,
#' with [adl_joint_test()] on the larger model.
#'
#' The J statistic is asymptotically normal under the null. In small
#' samples it over-rejects, because \eqn{\hat y_B} is correlated with the
#' errors of A through the estimated \eqn{\hat\gamma_B}: treat marginal
#' rejections with care. `dist = "t"` gives the usual t reference.
#'
#' @param fit_a,fit_b Two fits of the same dependent variable on the same
#'   data: `"adl"` objects, or `dynlm` fits (which can also be pure
#'   autoregressions, such as the AR(1) of ADL3). Their observations are
#'   matched by date.
#' @inheritParams adl_vcov
#' @param alpha Significance level for the decisions.
#' @param dist `"normal"` (default) or `"t"`.
#' @param names Labels for the two models in the output; by default `"A"` and
#'   `"B"`.
#' @return An object of class `"j_test"`: a list with `tests` (one row per
#'   direction: `null_hypothesis`, `model`, `added`, `theta`, `std.error`,
#'   `statistic`, `p.value`, `reject_null`), `outcome` (`"neither"`,
#'   `"against_a"`, `"against_b"` or `"both"`), `interpretation`, `names`,
#'   `formulas`, `y`, `nobs`, `sample` (first and last date), and the
#'   settings `cov_type`, `maxlags`, `alpha` and `dist`.
#' @references Davidson, R. and MacKinnon, J. G. (1981). Several tests for
#'   model specification in the presence of alternative hypotheses.
#'   *Econometrica* 49(3), 781-793.
#' @examples
#' # Partial adjustment against DL(1): neither nests the other
#' pa <- adl(y ~ x, data = adl_sim, order = c(1, 0))
#' dl1 <- adl(y ~ x, data = adl_sim, order = c(0, 1))
#' j_test(pa, dl1, names = c("PA", "DL(1)"))
#' j_test(pa, dl1, cov_type = "HAC")
#'
#' # AR(1) against DL(1), as in ADL3: the AR(1) has no x, so fit it with dynlm
#' if (requireNamespace("dynlm", quietly = TRUE)) {
#'   library(dynlm)  # dynlm finds L() only when attached
#'   ar1 <- dynlm(y ~ L(y, 1), data = adl_sim)
#'   j_test(ar1, dl1, names = c("AR(1)", "DL(1)"))
#' }
#' @export
j_test <- function(fit_a, fit_b, cov_type = c("nonrobust", "HAC"),
                   maxlags = NULL, alpha = 0.05, dist = c("normal", "t"),
                   names = c("A", "B")) {
  pa <- j_parts(fit_a, "fit_a")
  pb <- j_parts(fit_b, "fit_b")
  cov_type <- match.arg(cov_type)
  dist <- match.arg(dist)
  check_alpha(alpha)
  if (!is.character(names) || length(names) != 2 || anyNA(names) ||
      names[1] == names[2]) {
    stop("`names` must be two different labels.", call. = FALSE)
  }

  # common sample
  ia <- pa$index
  ib <- pb$index
  num <- function(i) class(i)[1] %in% c("integer", "numeric")
  if (!identical(class(ia), class(ib)) && !(num(ia) && num(ib))) {
    stop("The two models have different kinds of time index; fit them on the ",
         "same data.", call. = FALSE)
  }
  common <- ia[ia %in% ib]
  ra <- match(common, ia)
  rb <- match(common, ib)
  Xa <- pa$X[ra, , drop = FALSE]
  Xb <- pb$X[rb, , drop = FALSE]
  ya <- pa$y[ra]
  yb <- pb$y[rb]
  if (!length(common) || !isTRUE(all.equal(unname(ya), unname(yb)))) {
    stop("The two models do not have the same dependent variable on a common ",
         "sample: fit them on the same data.", call. = FALSE)
  }
  n <- length(common)
  if (n <= ncol(Xa) + 1 || n <= ncol(Xb) + 1) {
    stop("Too few common observations (", n, ") for the J-test.", call. = FALSE)
  }

  # nesting: the test needs regressors of each model outside the other's span
  rk <- function(X) qr(X)$rank
  r_ab <- rk(cbind(Xa, Xb))
  a_in_b <- r_ab == rk(Xb)
  b_in_a <- r_ab == rk(Xa)
  if (a_in_b || b_in_a) {
    what <- if (a_in_b && b_in_a) "The two models span the same regressors" else
      if (b_in_a) paste0("Model ", names[2], " is nested in model ", names[1]) else
        paste0("Model ", names[1], " is nested in model ", names[2])
    stop(what, ": the J-test is for non-nested models. Test the restrictions ",
         "directly on the larger model with adl_joint_test().", call. = FALSE)
  }

  y <- unname(ya)
  fitted_ols <- function(X) drop(X %*% qr.coef(qr(X), y))
  yhat_a <- fitted_ols(Xa)
  yhat_b <- fitted_ols(Xb)
  lags <- hac_lags(maxlags, cov_type, n)

  one <- function(X, added) {
    Z <- cbind(X, J_fitted = added)
    m <- stats::lm(y ~ 0 + Z)
    k <- ncol(Z)
    th <- unname(stats::coef(m)[k])
    V <- if (cov_type == "HAC") {
      sandwich::NeweyWest(m, lag = lags, prewhite = FALSE, adjust = FALSE)
    } else stats::vcov(m)
    se <- sqrt(V[k, k])
    inf <- single_coef_inference(th, se, stats::df.residual(m), dist)
    c(theta = th, std.error = se, statistic = inf$statistic, p.value = inf$p.value)
  }
  ta <- one(Xa, yhat_b)
  tb <- one(Xb, yhat_a)
  tests <- data.frame(
    null_hypothesis = c(paste(names[1], "encompasses", names[2]),
                        paste(names[2], "encompasses", names[1])),
    model = names,
    added = paste0("fitted ", rev(names)),
    theta = c(ta[["theta"]], tb[["theta"]]),
    std.error = c(ta[["std.error"]], tb[["std.error"]]),
    statistic = c(ta[["statistic"]], tb[["statistic"]]),
    p.value = c(ta[["p.value"]], tb[["p.value"]]),
    stringsAsFactors = FALSE)
  tests$reject_null <- tests$p.value < alpha

  rej <- tests$reject_null
  outcome <- if (!rej[1] && !rej[2]) "neither" else if (rej[1] && !rej[2])
    "against_a" else if (!rej[1] && rej[2]) "against_b" else "both"
  interpretation <- switch(outcome,
    neither = paste0("Neither rejected: both models are compatible with the ",
                     "data; the test cannot discriminate between them."),
    against_a = paste0("Evidence against ", names[1], ", in favour of ", names[2],
                       ": ", names[1], " does not encompass ", names[2], "."),
    against_b = paste0("Evidence against ", names[2], ", in favour of ", names[1],
                       ": ", names[2], " does not encompass ", names[1], "."),
    both = paste0("Both rejected: neither model accounts for the other; both ",
                  "look mis-specified, so consider a more general model."))

  structure(list(
    tests = tests,
    outcome = outcome,
    interpretation = interpretation,
    names = names,
    formulas = c(pa$formula, pb$formula),
    y = pa$yname,
    nobs = n,
    sample = format(range(common)),
    cov_type = cov_type,
    maxlags = lags,
    alpha = alpha,
    dist = dist
  ), class = "j_test")
}

#' Print a J-test
#'
#' @param x A `"j_test"` object.
#' @param digits Significant digits.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.j_test <- function(x, digits = 4, ...) {
  title <- "J-test of non-nested models (Davidson-MacKinnon)"
  cat(title, "\n", strrep("-", nchar(title)), "\n", sep = "")
  for (i in 1:2) {
    cat("Model ", x$names[i], ": ", x$formulas[i], "\n", sep = "")
  }
  cat("Sample:      ", x$sample[1], " - ", x$sample[2], " (", x$nobs,
      " observations, common to both models)\n", sep = "")
  cov <- if (x$cov_type == "HAC") paste0("HAC (Bartlett, ", x$maxlags, " lags)")
         else "conventional OLS"
  cat("Covariance:  ", cov, "\n", sep = "")
  cat("Decisions:   ", if (x$dist == "t") "t" else "z", " tests at ",
      format(100 * x$alpha), "%\n\n", sep = "")
  t <- x$tests
  show <- data.frame(
    H0 = t$null_hypothesis,
    Regression = paste0(t$model, " + ", t$added),
    theta = signif(t$theta, digits),
    `Std. error` = signif(t$std.error, digits),
    stat = signif(t$statistic, digits),
    `p-value` = signif(t$p.value, digits),
    Decision = ifelse(t$reject_null, "reject", "do not reject"),
    check.names = FALSE)
  names(show)[5] <- if (x$dist == "t") "t" else "z"
  print(show, row.names = FALSE, right = FALSE)
  cat("\n", x$interpretation, "\n", sep = "")
  invisible(x)
}

#' Regressors, response and dates of a fit for the J-test
#' @noRd
j_parts <- function(fit, arg) {
  if (inherits(fit, "adl")) {
    m <- fit$lm
    return(list(X = stats::model.matrix(m),
                y = stats::model.response(stats::model.frame(m)),
                index = fit$index, formula = fit$dynlm, yname = fit$y))
  }
  if (inherits(fit, "dynlm")) {
    f <- stats::formula(fit)
    return(list(X = stats::model.matrix(fit),
                y = as.numeric(stats::model.response(stats::model.frame(fit))),
                index = fit$index,
                formula = paste(deparse(f, width.cutoff = 500L), collapse = " "),
                yname = deparse(f[[2]])))
  }
  stop("`", arg, "` must be an object from adl() or a dynlm fit.", call. = FALSE)
}
