# Granger non-causality: Sargent's direct procedure (ADL3) -------------------
#
# Four regressions on a common sample:
#   1. y on lags of y and lags 1..r of x      tests x -> y   (lags of x = 0)
#   2. x on lags of x and lags 1..r of y      tests y -> x   (lags of y = 0)
#   3. model 1 plus x_t                       instantaneous  (x_t = 0)
#   4. model 2 plus y_t                       instantaneous  (y_t = 0)
# and the eight-way classification of ADL3.

#' Granger non-causality tests (Sargent's direct procedure)
#'
#' Tests Granger non-causality in both directions and instantaneous
#' causality between two stationary series, following the direct procedure
#' of Sargent (1976), and classifies the relation into one of the eight
#' types of Granger causal relations.
#'
#' Four regressions are estimated on a common sample:
#'
#' 1. \eqn{y_t} on \eqn{y_{t-1},\dots,y_{t-s}} and
#'    \eqn{x_{t-1},\dots,x_{t-r}} (no \eqn{x_t}); \eqn{H_0}: the lags of
#'    \eqn{x} are jointly zero (x does not cause y);
#' 2. the same with the roles of \eqn{x} and \eqn{y} exchanged (y does not
#'    cause x);
#' 3. model 1 plus \eqn{x_t}; \eqn{H_0}: the coefficient of \eqn{x_t} is zero
#'    (no instantaneous causality);
#' 4. model 2 plus \eqn{y_t}; \eqn{H_0}: the coefficient of \eqn{y_t} is zero.
#'
#' `granger_test()` uses the same order `p` for every lag polynomial.
#' [granger_test_general()] allows different orders in each model.
#'
#' Each test is reported with both the chi-square (Wald) and the F
#' statistic; `dist` chooses which one decides the classification. The
#' default is the chi-square: the regressions are dynamic, so even with the
#' conventional covariance the exact F distribution does not apply and both
#' are asymptotic approximations. With the
#' conventional covariance (the default, as in EViews and
#' `lmtest::grangertest()`) and equal orders, the two instantaneous tests
#' are numerically identical; with HAC or different orders they can differ,
#' and the classification is then flagged as ambiguous when they disagree.
#'
#' Both series are assumed stationary. With integrated series the
#' chi-square and F references do not apply in general (the Wald statistic
#' can have a nonstandard distribution); difference the series, or use a
#' procedure designed for integrated variables. "Instantaneous causality" is
#' a conditional contemporaneous association between the two series: its
#' direction cannot be identified without further (structural) assumptions,
#' and the whole analysis is conditional on the information set used (the
#' two series only).
#'
#' @param y,x Column names in `data`, or, when `data` is `NULL`, the two
#'   series themselves (numeric vectors, `ts` or `zoo`). Alternatively, `y`
#'   can be a formula `y ~ x`, or `y ~ x | z1 + z2` to condition on other
#'   variables (see **Conditioning variables**); then `x` is not given.
#' @param p Lag order used in all four models.
#' @param cond_lags Number of lags \eqn{1, \dots, k} of each conditioning
#'   variable in every model; by default `p` (in [granger_test_general()],
#'   the largest of the orders).
#' @param data Optional data (`ts`, `zoo`/`xts`, data frame, tibble or
#'   tsibble) containing `y` and `x` (and the conditioning variables);
#'   required with a formula.
#' @param alpha Significance level for the decisions.
#' @param dist `"chisq"` (default) or `"F"`: the statistic used for the
#'   decisions and the classification.
#' @inheritParams adl_vcov
#' @inheritParams adl
#' @section Conditioning variables:
#' With `y ~ x | z1 + z2`, lags \eqn{1, \dots,} `cond_lags` of every
#' conditioning variable enter all four regressions (never \eqn{z_t}: the
#' information set of Granger's definition is the past). The tests then ask
#' whether \eqn{x} helps to predict \eqn{y} given the past of \eqn{y}
#' **and** of the \eqn{z}'s, and the classification is conditional on them.
#' A third variable that drives both series can create apparent causality in
#' the bivariate test (ADL3), and conditioning on it can remove it; it can
#' also reveal causality that the bivariate test misses. With equal orders
#' and the conventional covariance the two instantaneous tests are still
#' numerically identical, because both regressions share the conditioning
#' set. Automatic order selection (`select`) is not available with
#' conditioning variables.
#'
#' @return An object of class `"granger_test"`: a list with `tests` (a data
#'   frame with one row per model), `classification` (type number 1-8,
#'   symbol, description and the `ambiguous` flag), `models` (the four
#'   `"adl"` fits), `conditioning` (the conditioning variables, if any) and
#'   `cond_lags`, and the settings used.
#' @references
#' Granger, C. W. J. (1969). Investigating causal relations by econometric
#' models and cross-spectral methods. *Econometrica*, 37(3), 424-438.
#'
#' Sargent, T. J. (1976). A classical macroeconometric model for the United
#' States. *Journal of Political Economy*, 84(2), 207-237.
#' @examples
#' set.seed(1)
#' n <- 200
#' x <- as.numeric(arima.sim(list(ar = 0.5), n))
#' y <- numeric(n)
#' for (t in 2:n) y[t] <- 0.5 * y[t - 1] + 0.4 * x[t - 1] + rnorm(1)
#' d <- ts(cbind(y = y, x = x))
#'
#' granger_test("y", "x", p = 2, data = d)
#' granger_test(y, x, p = 2)                      # the series directly
#' granger_test("y", "x", p = 2, data = d, cov_type = "HAC")
#' granger_test(y ~ x, p = 2, data = d)          # formula interface
#'
#' # a third variable z drives both: x seems to cause y until z is held fixed
#' z <- as.numeric(arima.sim(list(ar = 0.5), n))
#' x2 <- c(0, z[-n]) + rnorm(n)
#' y2 <- c(0, 0, z[-c(n - 1, n)]) + rnorm(n)
#' d2 <- ts(cbind(y = y2, x = x2, z = z))
#' granger_test(y ~ x, p = 2, data = d2)
#' granger_test(y ~ x | z, p = 2, data = d2)
#' @export
granger_test <- function(y, x, p, data = NULL, alpha = 0.05,
                         dist = c("chisq", "F"),
                         cov_type = c("nonrobust", "HAC"), maxlags = NULL,
                         time = NULL, start = NULL, end = NULL,
                         cond_lags = p) {
  ynm <- deparse(substitute(y))
  xnm <- deparse(substitute(x))
  cond <- character(0)
  if (inherits(y, "formula")) {
    if (!missing(x)) {
      if (missing(p) && is.numeric(x)) p <- x else
        stop("With a formula, give the order as `p = `.", call. = FALSE)
    }
    gf <- parse_granger_formula(y, data)
    y <- gf$y; x <- gf$x; cond <- gf$cond
  }
  check_whole(p, "p")
  if (length(p) != 1 || p < 1) stop("`p` must be a single integer >= 1.", call. = FALSE)
  orders <- list(y = c(p, p), x = c(p, p))
  granger_engine(y, x, data, orders, orders, alpha, match.arg(dist),
                 match.arg(cov_type), maxlags, time, start, end,
                 ynm = ynm, xnm = xnm, cond = cond,
                 cond_lags = check_cond_lags(cond_lags, cond))
}

#' Granger non-causality tests with model-specific lag orders
#'
#' As [granger_test()], but each of the four models has its own lag orders
#' \eqn{(s, r)}: \eqn{s} lags of the dependent variable and \eqn{r} lags of
#' the other variable.
#'
#' @inheritParams granger_test
#' @param orders Named list with the orders of the two Granger regressions:
#'   `y = c(s, r)` for the equation of `y` (model 1) and `x = c(s, r)` for
#'   the equation of `x` (model 2).
#' @param inst_orders Orders of the instantaneous-causality regressions
#'   (models 3 and 4), in the same form; by default those of `orders`. The
#'   contemporaneous term is added to these.
#' @param select Automatic order selection: `NULL` (default, use `orders`),
#'   `"bic"`, `"aic"`, `"hq"` or `"gets"`. Each Granger regression (models 1
#'   and 2) gets its own \eqn{(s, r)}, with \eqn{0 \le s \le} `max_order`
#'   and \eqn{1 \le r \le} `max_order`, chosen on a common sample: by the
#'   information criterion over the full grid (with [adl_ic()]), or by the
#'   general-to-specific search of ADL3 (with [adl_gets()]: from
#'   \eqn{(}`max_order`, `max_order`\eqn{)} the highest lags of both variables
#'   are dropped jointly while insignificant). The instantaneous models use
#'   the selected orders unless `inst_orders` is given.
#'
#'   Selection and testing use the same data. With the information criteria
#'   the lags of the other variable are judged by their contribution to the
#'   fit, together with the own lags; with `"gets"` they are dropped by their
#'   significance and then tested, a pre-test that favours rejection, so the
#'   information criteria are preferable. ADL3 notes the trade-off: more lags
#'   of the other variable capture its influence more completely, but cost
#'   power.
#'
#'   The orders are chosen on the sample of the largest model (`max_order`
#'   lags), but the four tests are then run on the common sample of the
#'   selected orders, which can be longer: the search and the reported test
#'   use different observations.
#' @param max_order Largest lag considered when `select` is used. By default
#'   4 for quarterly, 12 for monthly and 2 for other data.
#' @return An object of class `"granger_test"`; see [granger_test()].
#' @examples
#' set.seed(1)
#' n <- 200
#' x <- as.numeric(arima.sim(list(ar = 0.5), n))
#' y <- numeric(n)
#' for (t in 2:n) y[t] <- 0.5 * y[t - 1] + 0.4 * x[t - 1] + rnorm(1)
#' d <- ts(cbind(y = y, x = x))
#' granger_test_general("y", "x", data = d,
#'                      orders = list(y = c(1, 2), x = c(2, 1)))
#' # orders chosen by BIC for each equation
#' granger_test_general("y", "x", data = d, select = "bic", max_order = 4)
#' @export
granger_test_general <- function(y, x, orders, inst_orders = orders,
                                 data = NULL, alpha = 0.05,
                                 dist = c("chisq", "F"),
                                 cov_type = c("nonrobust", "HAC"),
                                 maxlags = NULL, select = NULL, time = NULL,
                                 start = NULL, end = NULL, max_order = NULL,
                                 cond_lags = NULL) {
  ynm <- deparse(substitute(y))
  xnm <- deparse(substitute(x))
  cond <- character(0)
  if (inherits(y, "formula")) {
    if (!missing(x)) {
      if (missing(orders) && is.list(x)) orders <- x else
        stop("With a formula, give the orders as `orders = `.", call. = FALSE)
    }
    gf <- parse_granger_formula(y, data)
    y <- gf$y; x <- gf$x; cond <- gf$cond
    if (missing(inst_orders) && !missing(orders)) inst_orders <- orders
  }
  if (length(cond) && !is.null(select)) {
    stop("Automatic order selection (`select`) is not available with ",
         "conditioning variables; give `orders`.", call. = FALSE)
  }
  selection <- NULL
  if (!is.null(select)) {
    select <- match.arg(tolower(select), c("bic", "aic", "hq", "gets"))
    if (!missing(orders)) {
      stop("Give either `orders` or `select` (with `max_order`), not both.",
           call. = FALSE)
    }
    pair <- granger_data(y, x, data, ynm, xnm)
    if (is.null(max_order)) {
      z <- to_zoo(pair$data, time = time)
      freq <- tryCatch(stats::frequency(z), error = function(e) 1)
      max_order <- if (freq == 4) 4L else if (freq == 12) 12L else 2L
    }
    check_whole(max_order, "max_order")
    if (length(max_order) != 1 || max_order < 1) {
      stop("`max_order` must be a single integer >= 1.", call. = FALSE)
    }
    pick <- function(dep, oth) {
      f <- stats::as.formula(paste(dep, "~", oth))
      if (select == "gets") {
        g <- suppressWarnings(adl_gets(f, data = pair$data, max_order = max_order,
                                       contemporaneous = "exclude", time = time,
                                       start = start, end = end))
        c(length(g$fit$lags[[1]]), max(g$fit$lags[[2]]))
      } else {
        tab <- suppressMessages(adl_ic(f, data = pair$data,
                                       max_order = c(max_order, max_order),
                                       contemporaneous = "exclude", time = time,
                                       start = start, end = end))
        best <- tab[which.min(tab[[select]]), ]
        c(best$p, best[[paste0("q_", oth)]])
      }
    }
    orders <- list(y = as.integer(pick(pair$y, pair$x)),
                   x = as.integer(pick(pair$x, pair$y)))
    if (missing(inst_orders)) inst_orders <- orders
    selection <- list(method = select, max_order = as.integer(max_order),
                      orders = orders)
  } else if (missing(orders)) {
    stop("Give `orders`, or choose them with `select`.", call. = FALSE)
  }
  orders <- check_orders(orders, "orders")
  inst_orders <- check_orders(inst_orders, "inst_orders")
  cond_lags <- cond_lags %||% max(unlist(orders), unlist(inst_orders))
  res <- granger_engine(y, x, data, orders, inst_orders, alpha, match.arg(dist),
                        match.arg(cov_type), maxlags, time, start, end,
                        ynm = ynm, xnm = xnm, cond = cond,
                        cond_lags = check_cond_lags(cond_lags, cond))
  res$selection <- selection
  res
}

# Engine --------------------------------------------------------------------

granger_engine <- function(y, x, data, orders, inst_orders, alpha, dist,
                           cov_type, maxlags, time, start, end, ynm, xnm,
                           cond = character(0), cond_lags = 0L) {
  check_alpha(alpha)
  pair <- granger_data(y, x, data, ynm, xnm, cond)
  d <- pair$data
  yv <- pair$y
  xv <- pair$x

  # Common sample: drop the largest lag used in any of the four models
  m <- max(unlist(orders), unlist(inst_orders), if (length(cond)) cond_lags)
  specs <- list(
    list(dep = yv, other = orders$y, inst = FALSE),
    list(dep = xv, other = orders$x, inst = FALSE),
    list(dep = yv, other = inst_orders$y, inst = TRUE),
    list(dep = xv, other = inst_orders$x, inst = TRUE)
  )
  models <- lapply(specs, function(s) {
    oth <- if (s$dep == yv) xv else yv
    lags <- stats::setNames(c(list(seq_len(s$other[1]),
                                   (if (s$inst) 0L else 1L):s$other[2]),
                              rep(list(seq_len(cond_lags)), length(cond))),
                            c(s$dep, oth, cond))
    rhs <- paste(c(oth, cond), collapse = " + ")
    suppressMessages(adl(stats::as.formula(paste(s$dep, "~", rhs)), data = d,
                         lags = lags, max_lag = m, time = time,
                         start = start, end = end))
  })
  names(models) <- c("granger_y", "granger_x", "inst_y", "inst_x")

  tests <- list(
    one_granger_test(models[[1]], terms = paste0(xv, "_lag", seq_len(orders$y[2])),
                     cov_type, maxlags),
    one_granger_test(models[[2]], terms = paste0(yv, "_lag", seq_len(orders$x[2])),
                     cov_type, maxlags),
    one_granger_test(models[[3]], terms = xv, cov_type, maxlags),
    one_granger_test(models[[4]], terms = yv, cov_type, maxlags)
  )
  tab <- do.call(rbind, tests)
  tab <- cbind(
    data.frame(model = 1:4,
               equation = c(yv, xv, yv, xv),
               null_hypothesis = c(
                 paste(xv, "does not Granger-cause", yv),
                 paste(yv, "does not Granger-cause", xv),
                 paste("no instantaneous causality (", xv, "_t in ", yv, " equation)",
                       sep = ""),
                 paste("no instantaneous causality (", yv, "_t in ", xv, " equation)",
                       sep = "")),
               own_lags = c(orders$y[1], orders$x[1], inst_orders$y[1], inst_orders$x[1]),
               other_lags = c(orders$y[2], orders$x[2], inst_orders$y[2], inst_orders$x[2]),
               stringsAsFactors = FALSE),
    tab)
  pcol <- if (dist == "F") "p.value_F" else "p.value_chisq"
  tab$reject_null <- tab[[pcol]] < alpha

  cls <- classify_granger(tab$reject_null, yv, xv)
  idx <- models[[1]]$index
  structure(list(
    tests = tab,
    classification = cls,
    models = models,
    y = yv, x = xv,
    conditioning = cond,
    cond_lags = if (length(cond)) as.integer(cond_lags) else NA_integer_,
    alpha = alpha, dist = dist, cov_type = cov_type,
    maxlags = attr(adl_vcov(models[[1]], cov_type, maxlags), "maxlags"),
    nobs = length(idx),
    sample = c(format(idx[1]), format(idx[length(idx)]))
  ), class = "granger_test")
}

#' Wald test of zero restrictions, reported in F and chi-square form
#' @noRd
one_granger_test <- function(fit, terms, cov_type, maxlags) {
  b <- stats::coef(fit)
  V <- adl_vcov(fit, cov_type, maxlags)
  R <- matrix(0, length(terms), length(b))
  R[cbind(seq_along(terms), match(terms, names(b)))] <- 1
  q <- length(terms)
  w <- wald_stat(b, V, R, rep(0, q))
  dfr <- df_resid(fit)
  data.frame(F = w / q, df_num = q, df_denom = dfr,
             p.value_F = stats::pf(w / q, q, dfr, lower.tail = FALSE),
             chisq = w, df = q,
             p.value_chisq = stats::pchisq(w, q, lower.tail = FALSE))
}

#' The eight causal types of ADL3
#' @noRd
classify_granger <- function(reject, y, x) {
  xy <- reject[1]
  yx <- reject[2]
  ambiguous <- reject[3] != reject[4]
  one <- function(inst) {
    type <- if (!xy && !yx) { if (inst) 2L else 1L }
            else if (xy && !yx) { if (inst) 5L else 3L }
            else if (!xy && yx) { if (inst) 6L else 4L }
            else { if (inst) 8L else 7L }
    sep <- c(", ", " - ", " -> ", " <- ", " => ", " <= ", " <-> ", " <=> ")[type]
    sym <- paste0("(", x, sep, y, ")")
    desc <- switch(type,
      paste(x, "and", y, "are independent"),
      "Only instantaneous causality",
      paste(x, "causes", y, "; no instantaneous causality"),
      paste(y, "causes", x, "; no instantaneous causality"),
      paste(x, "causes", y, ", with instantaneous causality"),
      paste(y, "causes", x, ", with instantaneous causality"),
      "Feedback; no instantaneous causality",
      "Feedback, with instantaneous causality")
    desc <- gsub(" ;", ";", gsub(" ,", ",", desc))
    list(type = type, symbol = sym, description = desc)
  }
  if (!ambiguous) {
    c(one(reject[3]), list(ambiguous = FALSE))
  } else {
    a <- one(FALSE)
    b <- one(TRUE)
    list(type = c(a$type, b$type), symbol = c(a$symbol, b$symbol),
         description = c(a$description, b$description), ambiguous = TRUE)
  }
}

#' Two series from names + data, or from the series themselves
#' @noRd
granger_data <- function(y, x, data, ynm, xnm, cond = character(0)) {
  if (length(cond) && is.null(data)) {
    stop("Conditioning variables need `data`.", call. = FALSE)
  }
  if (!is.null(data)) {
    if (!is.character(y) || !is.character(x) || length(y) != 1 || length(x) != 1) {
      stop("With `data`, `y` and `x` must be column names.", call. = FALSE)
    }
    yv <- y; xv <- x
  } else {
    if (is.character(y) || is.character(x)) {
      stop("`y` and `x` are names but no `data` was given.", call. = FALSE)
    }
    yv <- make.names(ynm); xv <- make.names(xnm)
    if (yv == xv) { yv <- "y"; xv <- "x" }
    if (inherits(y, "zoo") || inherits(x, "zoo")) {
      data <- zoo::merge.zoo(y, x, all = FALSE)
      colnames(data) <- c(yv, xv)
    } else if (stats::is.ts(y) && stats::is.ts(x)) {
      data <- stats::ts.intersect(y, x)
      colnames(data) <- c(yv, xv)
    } else {
      if (length(y) != length(x)) stop("`y` and `x` have different lengths.", call. = FALSE)
      data <- stats::ts(cbind(as.numeric(y), as.numeric(x)))
      colnames(data) <- c(yv, xv)
    }
  }
  if (yv == xv) stop("`y` and `x` must be different variables.", call. = FALSE)
  if (any(cond %in% c(yv, xv)) || anyDuplicated(cond)) {
    stop("The conditioning variables must differ from `y`, `x` and each ",
         "other.", call. = FALSE)
  }
  if (!identical(make.names(c(yv, xv)), c(yv, xv))) {
    stop("Variable names must be syntactic R names (got '", yv, "', '", xv,
         "'); rename the columns first.", call. = FALSE)
  }
  list(data = data, y = yv, x = xv)
}

#' Parse y ~ x | z1 + z2
#' @noRd
parse_granger_formula <- function(f, data) {
  if (is.null(data)) stop("A formula needs `data`.", call. = FALSE)
  if (length(f) != 3) stop("The formula must be of the form y ~ x | z.", call. = FALSE)
  y <- all.vars(f[[2]])
  rhs <- f[[3]]
  cond <- character(0)
  if (is.call(rhs) && identical(rhs[[1]], as.name("|"))) {
    cond <- vapply(split_plus(rhs[[3]]), function(e) {
      if (!is.name(e)) stop("Conditioning variables must be plain names.", call. = FALSE)
      as.character(e)
    }, character(1))
    rhs <- rhs[[2]]
  }
  if (length(y) != 1 || !is.name(f[[2]]) || !is.name(rhs)) {
    stop("The formula must be y ~ x | z1 + z2, with single variables y and x ",
         "(no lags or transformations).", call. = FALSE)
  }
  list(y = y, x = as.character(rhs), cond = unname(cond))
}

check_cond_lags <- function(k, cond) {
  if (!length(cond)) return(0L)
  check_whole(k, "cond_lags")
  if (length(k) != 1 || k < 1) {
    stop("`cond_lags` must be a single integer >= 1.", call. = FALSE)
  }
  as.integer(k)
}

check_orders <- function(o, what) {
  if (!is.list(o) || !setequal(names(o), c("y", "x"))) {
    stop("`", what, "` must be a list with elements `y` and `x`, each c(s, r).",
         call. = FALSE)
  }
  for (k in c("y", "x")) {
    v <- o[[k]]
    check_whole(v, paste0(what, "$", k))
    if (length(v) != 2 || v[1] < 0 || v[2] < 1) {
      stop("`", what, "$", k, "` must be c(s, r) with s >= 0 and r >= 1.",
           call. = FALSE)
    }
  }
  list(y = as.integer(o$y), x = as.integer(o$x))
}

# Printing ------------------------------------------------------------------

#' Print Granger non-causality tests
#'
#' @param x A `"granger_test"` object.
#' @param digits Significant digits.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.granger_test <- function(x, digits = 4, ...) {
  title <- "Granger non-causality tests (Sargent's direct procedure)"
  cat(title, "\n", strrep("-", nchar(title)), "\n", sep = "")
  cat("Sample:      ", x$sample[1], " - ", x$sample[2], " (", x$nobs,
      " observations, common to the four models)\n", sep = "")
  cov <- if (x$cov_type == "HAC") paste0("HAC (Bartlett, ", x$maxlags, " lags)")
         else "conventional OLS"
  cat("Covariance:  ", cov, "\n", sep = "")
  if (length(x$conditioning)) {
    cat("Conditioning: lags 1-", x$cond_lags, " of ",
        paste(x$conditioning, collapse = ", "), " in every model\n", sep = "")
  }
  cat("Decisions:   ", if (x$dist == "F") "F" else "Chi-square", " tests at ",
      format(100 * x$alpha), "%\n", sep = "")
  if (!is.null(x$selection)) {
    cat("Orders:      (s, r) selected by ",
        if (x$selection$method == "gets") "general-to-specific" else
          toupper(x$selection$method),
        " for each equation, up to ", x$selection$max_order, " lags\n", sep = "")
  }
  cat("\n")

  t <- x$tests
  show <- data.frame(
    H0 = t$null_hypothesis,
    `s,r` = paste0(t$own_lags, ",", t$other_lags),
    `Chi-sq` = signif(t$chisq, digits),
    df = t$df,
    `p(Chi-sq)` = signif(t$p.value_chisq, digits),
    F = signif(t$F, digits),
    `df(F)` = paste0(t$df_num, ",", t$df_denom),
    `p(F)` = signif(t$p.value_F, digits),
    Decision = ifelse(t$reject_null, "reject", "do not reject"),
    check.names = FALSE)
  print(show, row.names = FALSE, right = FALSE)

  cl <- x$classification
  cat("\nCausal relation: ")
  if (!cl$ambiguous) {
    cat(cl$symbol, " ", cl$description, " [type ", cl$type, " of 8]\n", sep = "")
  } else {
    cat("ambiguous - the two instantaneous-causality tests disagree:\n")
    cat("  ", cl$symbol[1], " ", cl$description[1], " [type ", cl$type[1],
        "], or\n", sep = "")
    cat("  ", cl$symbol[2], " ", cl$description[2], " [type ", cl$type[2],
        "]\n", sep = "")
  }
  if (length(x$conditioning)) {
    cat("The relation is conditional on the past of ",
        paste(x$conditioning, collapse = ", "), ".\n", sep = "")
  }
  cat(if (length(x$conditioning)) "All series" else "Both series",
      " are assumed stationary.\n", sep = "")
  invisible(x)
}
