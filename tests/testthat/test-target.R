cal <- local({
  set.seed(21)
  calibrate_shape(40, beta = c(-1.8, -1, 0, 2), c = c(-1.5, -0.5, 0, 0.5), reps = 60)
})

test_that("calibration curves are monotone in the expected direction", {
  expect_s3_class(cal, "shape_calibration")
  expect_true(all(diff(cal$gamma$gamma) > 0))
  expect_true(all(diff(cal$colless$colless) < 0))
  expect_output(print(cal), "n = 40 tips (60 draws per value)", fixed = TRUE)
})

test_that("shape_params inverts the curves and refuses targets outside them", {
  p <- shape_params(cal, gamma = cal$gamma$gamma[2], colless = cal$colless$colless[3])
  expect_equal(p$c, cal$gamma$c[2])
  expect_equal(p$beta, cal$colless$beta[3])
  expect_error(shape_params(cal, gamma = 100, colless = 0), "outside")
  expect_error(shape_params(list(), gamma = 0, colless = 0), "calibrate_shape")
  expect_error(shape_params(cal, gamma = NA_real_, colless = 0), "finite")
})

test_that("sim_tree_target returns trees within the tolerance", {
  set.seed(22)
  trees <- sim_tree_target(40, gamma = -1, colless = 0.3, nsim = 8, tol = c(0.3, 0.4),
                           calibration = cal, root_age = 10)
  expect_s3_class(trees, "multiPhylo")
  expect_length(trees, 8)
  shape <- t(vapply(trees, tree_shape, numeric(2)))
  expect_equal(unname(shape[, "gamma"]), attr(trees, "shape")$gamma, tolerance = 1e-10)
  expect_equal(unname(shape[, "colless"]), attr(trees, "shape")$colless, tolerance = 1e-10)
  expect_true(all(abs(shape[, "gamma"] + 1) <= 0.3))
  expect_true(all(abs(shape[, "colless"] - 0.3) <= 0.4))
  expect_true(all(vapply(trees, ape::is.ultrametric, TRUE)))
  expect_equal(max(ape::node.depth.edgelength(trees[[1]])), 10, tolerance = 1e-10)
  one <- sim_tree_target(40, gamma = 0, colless = 0, calibration = cal)
  expect_s3_class(one, "phylo")
  expect_error(sim_tree_target(30, gamma = 0, colless = 0, calibration = cal), "same `n`")
  # nothing within an impossible tolerance: a warning, then an error
  expect_warning(
    expect_error(sim_tree_target(40, gamma = 0, colless = 0, nsim = 3, tol = 1e-6,
                                 calibration = cal, max_tries = 5), "no tree"),
    "only")
})

test_that("noisy curves are made monotone and small n is refused", {
  set.seed(23)
  small <- calibrate_shape(6, reps = 5)
  expect_true(all(diff(small$gamma$gamma) >= 0))
  expect_true(all(diff(small$colless$colless) <= 0))
  expect_error(calibrate_shape(3), "at least 4")
  expect_error(calibrate_shape(10, beta = 0), "two distinct")
})

test_that("the class follows nsim even when fewer trees are found", {
  set.seed(24)
  expect_warning(trees <- sim_tree_target(40, gamma = -1, colless = 0.3, nsim = 3,
                                          tol = 0.05, calibration = cal, max_tries = 40),
                 "only")
  expect_s3_class(trees, "multiPhylo")
})

test_that("shape_grid crosses the levels", {
  grid <- shape_grid(40, gamma = cal$gamma$gamma[2:3], colless = cal$colless$colless[2:3],
                     calibration = cal)
  expect_equal(nrow(grid), 4)
  expect_setequal(grid$c, cal$gamma$c[2:3])
  expect_setequal(grid$beta, cal$colless$beta[2:3])
  expect_s3_class(attr(grid, "calibration"), "shape_calibration")
})
