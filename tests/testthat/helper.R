#helpers

expect_equal_est <- function(object, expected, ..., ignore_attr = TRUE, tolerance = 1e-5) {
  expect_equal(object, expected, ignore_attr = ignore_attr, tolerance = tolerance, ...)
}

weighted_mean_diff <- function(y, t, w, subset = NULL) {
  if (is.null(subset)) subset <- seq_along(y)
  weighted.mean(y[subset][t[subset]==1], w[subset][t[subset]==1]) -
    weighted.mean(y[subset][t[subset]==0], w[subset][t[subset]==0])
}
