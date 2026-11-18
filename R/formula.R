# Formula interface ---------------------------------------------------------
#
# Two forms (DESIGN.md section 3.3):
#   y ~ x1 + x2 | d1 + d2                 variables only; lags from order/lags
#   y ~ L(y, 1:2) + L(x1, 0:2) + x2 | d1  explicit lags, dynlm style
# A bare regressor in the second form means lag 0; L(x) means lag 1, as in
# dynlm. Terms after `|` enter without lags (dummies, other deterministics).

#' Parse an ADL formula
#'
#' @param formula A two-sided formula.
#' @return A list with `y`, `x`, `unlagged` and `lags` (named list of lag
#'   sets, or `NULL` when the formula has no `L()` terms).
#' @noRd
parse_adl_formula <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 3) {
    stop("`formula` must be two-sided, e.g. y ~ x or y ~ L(y, 1) + L(x, 0:1).",
         call. = FALSE)
  }
  env <- environment(formula) %||% parent.frame()
  lhs <- formula[[2]]
  if (!is.symbol(lhs)) {
    stop("The left-hand side must be a single variable name.", call. = FALSE)
  }
  y <- as.character(lhs)

  rhs <- formula[[3]]
  unlagged_expr <- NULL
  if (is.call(rhs) && identical(rhs[[1]], as.name("|"))) {
    unlagged_expr <- rhs[[3]]
    rhs <- rhs[[2]]
  }

  terms <- split_plus(rhs)
  lag_terms <- list()
  has_L <- FALSE
  for (tm in terms) {
    if (is.symbol(tm)) {
      v <- as.character(tm)
      lag_terms[[length(lag_terms) + 1]] <- list(var = v, lags = 0L, L = FALSE)
    } else if (is.call(tm) && identical(tm[[1]], as.name("L"))) {
      has_L <- TRUE
      if (!is.symbol(tm[[2]])) {
        stop("L() takes a variable name as first argument: ", deparse(tm),
             call. = FALSE)
      }
      k <- if (length(tm) >= 3) eval(tm[[3]], env) else 1L
      check_whole(k, paste0("lags in ", deparse(tm)))
      lag_terms[[length(lag_terms) + 1]] <- list(var = as.character(tm[[2]]),
                                                 lags = as.integer(k), L = TRUE)
    } else if (is.numeric(tm) && tm %in% c(0, 1)) {
      stop("Use the `intercept` argument instead of 0/1 in the formula.",
           call. = FALSE)
    } else {
      stop("Term not understood: ", deparse(tm),
           ". Use variable names or L(variable, lags).", call. = FALSE)
    }
  }

  vars <- unique(vapply(lag_terms, `[[`, "", "var"))
  x <- setdiff(vars, y)
  if (length(x) == 0) stop("The formula has no regressors.", call. = FALSE)

  if (y %in% vars && !has_L) {
    stop("Lags of '", y, "' must be written as L(", y, ", ...).", call. = FALSE)
  }

  lags <- NULL
  if (has_L) {
    lags <- stats::setNames(vector("list", length(vars)), vars)
    for (tm in lag_terms) lags[[tm$var]] <- c(lags[[tm$var]], tm$lags)
    if (y %in% names(lags) && any(lags[[y]] == 0)) {
      stop("The dependent variable cannot appear at lag 0 on the right-hand side.",
           call. = FALSE)
    }
  }

  unlagged <- character(0)
  if (!is.null(unlagged_expr)) {
    u <- split_plus(unlagged_expr)
    ok <- vapply(u, is.symbol, logical(1))
    if (!all(ok)) {
      stop("Only variable names are allowed after `|`.", call. = FALSE)
    }
    unlagged <- vapply(u, as.character, "")
    clash <- intersect(unlagged, c(y, x))
    if (length(clash)) {
      stop("Variable(s) both lagged and after `|`: ", paste(clash, collapse = ", "),
           call. = FALSE)
    }
  }

  list(y = y, x = x, unlagged = unlagged, lags = lags)
}

#' Split an expression on `+`
#' @noRd
split_plus <- function(e) {
  if (is.call(e) && identical(e[[1]], as.name("+")) && length(e) == 3) {
    c(split_plus(e[[2]]), split_plus(e[[3]]))
  } else if (is.call(e) && identical(e[[1]], as.name("("))) {
    split_plus(e[[2]])
  } else {
    list(e)
  }
}

#' Write a lag set compactly: 1:3, 2, c(1, 4)
#' @noRd
format_lags <- function(l) {
  if (length(l) == 1) return(as.character(l))
  if (all(diff(l) == 1)) return(paste0(l[1], ":", l[length(l)]))
  paste0("c(", paste(l, collapse = ", "), ")")
}

#' The dynlm formula equivalent to a lag specification
#' @noRd
dynlm_formula <- function(lags, unlagged = character(0), trend = FALSE,
                          intercept = TRUE) {
  y <- names(lags)[1]
  rhs <- character(0)
  for (v in names(lags)) {
    l <- lags[[v]]
    if (length(l) == 0) next
    rhs <- c(rhs, if (identical(l, 0L)) v else paste0("L(", v, ", ", format_lags(l), ")"))
  }
  rhs <- c(rhs, unlagged, if (trend) "trend(...)")
  if (!intercept) rhs <- c(rhs, "0")
  paste(y, "~", paste(rhs, collapse = " + "))
}

#' Short model label: ADL(p, q1, q2) when lags are contiguous from 1 / first_lag
#' @noRd
adl_label <- function(lags) {
  y <- lags[[1]]
  x <- lags[-1]
  simple_y <- length(y) == 0 || identical(y, seq_len(max(y)))
  simple_x <- all(vapply(x, function(l) all(diff(l) == 1) && l[1] %in% c(0, 1),
                         logical(1)))
  if (simple_y && simple_x) {
    p <- if (length(y)) max(y) else 0
    q <- vapply(x, max, numeric(1))
    lab <- paste0("ADL(", paste(c(p, q), collapse = ", "), ")")
    excl <- names(x)[vapply(x, function(l) l[1] == 1, logical(1))]
    if (length(excl)) lab <- paste0(lab, ", no contemporaneous ",
                                    paste(excl, collapse = ", "))
    return(lab)
  }
  "ADL (general lag sets)"
}

#' Build `y ~ x1 + x2 | u1 + u2` from variable names, as language objects
#' (symbols), so that any name -- spaces, backticks, backslashes -- keeps
#' its identity (no text is parsed)
#' @noRd
adl_formula <- function(y, x, unlagged = character(0)) {
  plus <- function(v) Reduce(function(a, b) call("+", a, b), lapply(v, as.name))
  rhs <- if (length(x)) plus(x) else 1
  if (length(unlagged)) rhs <- call("|", rhs, plus(unlagged))
  stats::as.formula(call("~", as.name(y), rhs), env = globalenv())
}
