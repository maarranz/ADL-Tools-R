# Covariance of the OLS estimates ---------------------------------------------
#
# Same conventions as forecastingtools: conventional OLS by default (as in
# EViews); HAC is Newey-West with Bartlett weights 1 - j/(L+1), no
# prewhitening and no small-sample adjustment, with automatic bandwidth
# floor(4 (n/100)^(2/9)) capped at n - 1.

#' Covariance matrix of an ADL fit
#'
#' @param fit An `"adl"` object.
#' @param cov_type `"nonrobust"` (conventional OLS, the default) or `"HAC"`
#'   (Newey-West, Bartlett kernel, no prewhitening, no small-sample
#'   adjustment).
#' @param maxlags HAC bandwidth. `NULL` uses the automatic rule
#'   \eqn{\lfloor 4 (n/100)^{2/9} \rfloor}{floor(4 (n/100)^(2/9))}, capped at
#'   \eqn{n - 1}. Must be `NULL` with `cov_type = "nonrobust"`.
#' @return The covariance matrix, with the coefficient names, and the
#'   bandwidth used in attribute `"maxlags"` (`NA` for nonrobust).
#' @examples
#' set.seed(1)
#' x <- cumsum(rnorm(100))
#' y <- stats::filter(0.5 * x + rnorm(100), 0.6, method = "recursive")
#' d <- ts(cbind(y = as.numeric(y), x = x))
#' fit <- adl(y ~ x, data = d, order = c(1, 1))
#' adl_vcov(fit)
#' adl_vcov(fit, cov_type = "HAC")
#' @export
adl_vcov <- function(fit, cov_type = c("nonrobust", "HAC"), maxlags = NULL) {
  check_adl(fit)
  cov_type <- match.arg(cov_type)
  lags <- hac_lags(maxlags, cov_type, stats::nobs(fit))
  V <- if (cov_type == "HAC") {
    sandwich::NeweyWest(fit$lm, lag = lags, prewhite = FALSE, adjust = FALSE)
  } else {
    stats::vcov(fit$lm)
  }
  V <- as.matrix(V)
  attr(V, "maxlags") <- lags
  V
}

#' Validate maxlags and return the bandwidth (NA for nonrobust)
#' @noRd
hac_lags <- function(maxlags, cov_type, n) {
  if (cov_type == "nonrobust") {
    if (!is.null(maxlags)) {
      stop("`maxlags` is only used with HAC inference; leave it NULL with ",
           "cov_type = \"nonrobust\".", call. = FALSE)
    }
    return(NA_integer_)
  }
  if (is.null(maxlags)) return(as.integer(min(n - 1, floor(4 * (n / 100)^(2 / 9)))))
  if (!is.numeric(maxlags) || length(maxlags) != 1 || is.na(maxlags) ||
      maxlags < 0 || maxlags != round(maxlags)) {
    stop("`maxlags` must be NULL or a non-negative integer.", call. = FALSE)
  }
  if (maxlags >= n) {
    stop("The HAC bandwidth `maxlags` must be smaller than the number of ",
         "observations (", n, ").", call. = FALSE)
  }
  as.integer(maxlags)
}

check_adl <- function(fit) {
  if (!inherits(fit, "adl")) stop("`fit` must be an object from adl().", call. = FALSE)
  invisible(fit)
}
