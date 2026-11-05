test_that("the package namespace loads", {
  expect_true(requireNamespace("adltools", quietly = TRUE))
})
