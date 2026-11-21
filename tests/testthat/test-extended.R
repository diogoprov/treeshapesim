# Extended tests (G5.10): run with TREESHAPESIM_EXTENDED_TESTS="true". They
# take a few minutes; see CONTRIBUTING.md.

skip_extended <- function() {
  testthat::skip_if_not(
    identical(Sys.getenv("TREESHAPESIM_EXTENDED_TESTS"), "true"),
    "extended tests not requested (TREESHAPESIM_EXTENDED_TESTS != 'true')"
  )
}

test_that("parameter recovery over many seeds (extended)", {
  # G5.6b
  skip_extended()
  grid <- expand.grid(beta = c(-1.5, -1, 0), c = c(-1, -0.5, 0, 0.5))
  for (j in seq_len(nrow(grid))) {
    est <- vapply(1:30, function(s) {
      set.seed(1000 * j + s)
      f <- fit_shape(sim_tree_shape(400, beta = grid$beta[j], c = grid$c[j]))
      c(f$beta, f$c)
    }, numeric(2))
    expect_lt(abs(mean(est[1, ]) - grid$beta[j]), 0.2)
    expect_lt(abs(mean(est[2, ]) - grid$c[j]), 0.1)
  }
})

test_that("the error of the estimates decreases with the number of tips (extended)", {
  # G5.7
  skip_extended()
  rmse <- vapply(c(50, 200, 800), function(n) {
    est <- vapply(1:40, function(s) {
      set.seed(n + s)
      f <- fit_shape(sim_tree_shape(n, beta = -1, c = -0.5))
      c(f$beta, f$c)
    }, numeric(2))
    sqrt(rowMeans((est - c(-1, -0.5))^2))
  }, numeric(2))
  expect_false(is.unsorted(rev(rmse[1, ])))
  expect_false(is.unsorted(rev(rmse[2, ])))
})
