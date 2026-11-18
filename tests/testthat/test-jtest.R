j_data <- function(n = 300, seed = 31) {
  set.seed(seed)
  x <- as.numeric(stats::arima.sim(list(ar = 0.6), n))
  # DGP: DL(1), no lag of y
  y <- 1 + 0.8 * x + 0.6 * c(0, x[-n]) + rnorm(n, sd = 0.5)
  stats::ts(cbind(y = y, x = x), start = c(1980, 1), frequency = 4)
}

# The J regressions written out on the common sample
j_by_hand <- function(Xa, Xb, y, cov_type = "nonrobust", lag = NULL) {
  one <- function(X, added) {
    m <- lm(y ~ 0 + X + added)
    V <- if (cov_type == "HAC") sandwich::NeweyWest(m, lag = lag, prewhite = FALSE,
                                                     adjust = FALSE) else vcov(m)
    k <- length(coef(m))
    c(coef(m)[[k]], sqrt(V[k, k]))
  }
  rbind(one(Xa, fitted(lm(y ~ 0 + Xb))), one(Xb, fitted(lm(y ~ 0 + Xa))))
}

test_that("the J statistics equal the regressions written out", {
  d <- j_data()
  fa <- adl(y ~ x, data = d, order = c(2, 0))   # loses two observations
  fb <- adl(y ~ x, data = d, order = c(0, 1))   # loses one
  for (ct in c("nonrobust", "HAC")) {
    j <- j_test(fa, fb, cov_type = ct)
    expect_equal(j$nobs, nobs(fa))
    Xa <- model.matrix(fa$lm)
    Xb <- model.matrix(fb$lm)[-1, ]
    y <- as.numeric(d[-(1:2), "y"])
    ref <- j_by_hand(Xa, Xb, y, ct, j$maxlags)
    expect_equal(j$tests$theta, ref[, 1], ignore_attr = TRUE)
    expect_equal(j$tests$std.error, ref[, 2], ignore_attr = TRUE)
    expect_equal(j$tests$statistic, ref[, 1] / ref[, 2], ignore_attr = TRUE)
  }
  jt <- j_test(fa, fb, dist = "t")
  expect_equal(jt$tests$p.value[1],
               2 * pt(abs(jt$tests$statistic[1]), nobs(fa) - 5, lower.tail = FALSE))
})

test_that("the four outcomes", {
  d <- j_data(n = 400)
  pa <- adl(y ~ x, data = d, order = c(1, 0))
  dl <- adl(y ~ x, data = d, order = c(0, 1))
  j <- j_test(pa, dl, names = c("PA", "DL"))
  expect_equal(j$outcome, "against_a")
  expect_match(j$interpretation, "against PA, in favour of DL")
  j2 <- j_test(dl, pa, names = c("DL", "PA"))
  expect_equal(j2$outcome, "against_b")
  # alpha so small that nothing is rejected, or so large that both are
  expect_equal(j_test(pa, dl, alpha = 1e-300)$outcome, "neither")
  expect_equal(j_test(pa, dl, alpha = 0.999999)$outcome, "both")
  expect_match(j_test(pa, dl, alpha = 0.999999)$interpretation, "more general")
})

test_that("dynlm fits: the AR(1) against the DL(1) of ADL3", {
  skip_if_not_installed("dynlm")
  suppressPackageStartupMessages(library(dynlm))
  # not `d`: dynlm evaluates `data` where its own d() (difference) is found
  dat <- j_data()
  ar1 <- dynlm::dynlm(y ~ L(y, 1), data = dat)
  dl <- adl(y ~ x, data = dat, order = c(0, 1))
  j <- j_test(ar1, dl, names = c("AR(1)", "DL(1)"))
  expect_equal(j$outcome, "against_a")
  expect_equal(j$nobs, nobs(ar1))
  expect_equal(j$y, "y")
})

test_that("nested, identical and incompatible models are refused", {
  d <- j_data()
  small <- adl(y ~ x, data = d, order = c(1, 0))
  big <- adl(y ~ x, data = d, order = c(1, 1))
  expect_error(j_test(small, big), "nested in model B")
  expect_error(j_test(big, small), "nested in model A")
  expect_error(j_test(big, big), "same regressors")
  d2 <- d
  d2[, "y"] <- d2[, "y"] + 1
  other <- adl(y ~ x, data = d2, order = c(0, 1))
  expect_error(j_test(small, other), "same dependent variable")
  expect_error(j_test(small, lm(y ~ x, data = as.data.frame(d))), "dynlm fit")
  expect_error(j_test(small, other, names = c("A", "A")), "different labels")
  d3 <- data.frame(y = as.numeric(d[, "y"]), x = as.numeric(d[, "x"]))
  dl3 <- adl(y ~ x, data = d3, order = c(0, 1))
  expect_error(j_test(small, dl3), "kinds of time index|same dependent")
})

test_that("printing and tidiers", {
  d <- j_data()
  j <- j_test(adl(y ~ x, data = d, order = c(1, 0)),
              adl(y ~ x, data = d, order = c(0, 1)), cov_type = "HAC")
  out <- capture.output(print(j))
  expect_true(any(grepl("Davidson-MacKinnon", out)))
  expect_true(any(grepl("A encompasses B", out)))
  expect_true(any(grepl("HAC", out)))
  expect_equal(nrow(tidy(j)), 2)
  expect_equal(glance(j)$outcome, j$outcome)
})
