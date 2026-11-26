# Port of cointegrationtools 0.1.1, audit item H04 (26 Nov 2026)

test_that("H04: text dates must be ISO dates, and gaps are then detected", {
  set.seed(3)
  d <- data.frame(x = rnorm(120), y = rnorm(120))
  d$date <- format(seq(as.Date("2000-01-01"), by = "month", length.out = 120))
  f_text <- adl(y ~ x, data = d, order = c(1, 1), time = "date")
  f_date <- adl(y ~ x, data = transform(d, date = as.Date(date)), order = c(1, 1),
                time = "date")
  expect_equal(coef(f_text), coef(f_date))
  expect_s3_class(f_text$index, "Date")
  expect_error(adl(y ~ x, data = d[-35, ], order = c(1, 1), time = "date"), "gaps")
  # the same as a factor
  df <- transform(d, date = factor(date))
  expect_equal(coef(adl(y ~ x, data = df, order = c(1, 1), time = "date")), coef(f_text))
  # text that is not a set of ISO dates
  for (bad in list(sprintf("t%03d", 1:120),
                   format(seq(as.Date("2000-01-01"), by = "month", length.out = 120),
                          "%d/%m/%Y"))) {
    d$date <- bad
    expect_error(adl(y ~ x, data = d, order = c(1, 1), time = "date"), "ISO dates")
  }
  # values that parse but do not round-trip, or are impossible dates
  d$date <- format(seq(as.Date("2000-01-01"), by = "month", length.out = 120))
  d$date[10] <- "2000-10-1"
  expect_error(adl(y ~ x, data = d, order = c(1, 1), time = "date"), "ISO dates")
  d$date[10] <- "2000-02-30"
  expect_error(adl(y ~ x, data = d, order = c(1, 1), time = "date"), "ISO dates")
})

test_that("H04: an index of an unsupported class is refused", {
  expect_error(check_regular(factor(c("a", "b", "c"))), "Unsupported time index")
  expect_error(check_regular(as.difftime(1:5, units = "days")), "Unsupported time index")
})
