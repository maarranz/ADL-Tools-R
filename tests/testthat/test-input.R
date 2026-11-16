make_ts <- function(n = 40) {
  set.seed(1)
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + rnorm(n)
  stats::ts(cbind(y = y, x = x), start = c(2000, 1), frequency = 4)
}

test_that("ts, zoo, data frame and tibble give the same data", {
  m <- make_ts()
  ref <- as_adl_data(m, c("y", "x"))
  expect_equal(ref$n, 40)

  z <- zoo::as.zoo(m)
  expect_equal(zoo::coredata(as_adl_data(z, c("y", "x"))$z),
               zoo::coredata(ref$z))

  df <- data.frame(date = zoo::as.Date(zoo::as.yearqtr(stats::time(m))),
                   y = m[, "y"], x = m[, "x"])
  expect_equal(unname(zoo::coredata(as_adl_data(df, c("y", "x"))$z)),
               unname(zoo::coredata(ref$z)))

  expect_message(r <- as_adl_data(df[, -1], c("y", "x")), "row order")
  expect_equal(r$n, 40)

  skip_if_not_installed("tibble")
  tb <- tibble::as_tibble(df)
  expect_equal(unname(zoo::coredata(as_adl_data(tb, c("y", "x"))$z)),
               unname(zoo::coredata(ref$z)))
})

test_that("data frames are sorted by the time column", {
  m <- make_ts()
  df <- data.frame(date = zoo::as.Date(zoo::as.yearqtr(stats::time(m))),
                   y = m[, "y"], x = m[, "x"])
  shuffled <- df[sample(nrow(df)), ]
  expect_equal(as_adl_data(shuffled, c("y", "x"))$z,
               as_adl_data(df, c("y", "x"))$z)
})

test_that("leading/trailing NA are trimmed, internal NA stop", {
  m <- make_ts()
  m[1:2, "x"] <- NA
  m[40, "y"] <- NA
  expect_message(r <- as_adl_data(m, c("y", "x")), "Trimmed 3")
  expect_equal(r$n, 37)

  m2 <- make_ts()
  m2[10, "x"] <- NA
  expect_error(as_adl_data(m2, c("y", "x")), "inside the sample")
})

test_that("start and end select the sample", {
  m <- make_ts()
  r <- as_adl_data(m, c("y", "x"), start = c(2001, 1), end = c(2004, 4))
  expect_equal(r$n, 16)
  expect_equal(format(r$index[1]), "2001 Q1")
})

test_that("bad inputs give clear errors", {
  m <- make_ts()
  expect_error(as_adl_data(m, c("y", "w")), "not found")
  expect_error(as_adl_data(list(1, 2), "y"), "Unsupported")
  expect_error(as_adl_data(m[, "y"], "y"), "single unnamed series")
})
