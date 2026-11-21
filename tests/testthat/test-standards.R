#' Tests for the rOpenSci statistical software standards
#'
#' @srrstats {G5.0} Correctness is checked on a stored tree
#'   (`fixtures/tree80.nwk`, test-fit.R) with values from apTreeshape and
#'   laser; there are no standard reference data sets for tree shape.
#' @srrstats {G5.1} The test data are simulated with the exported functions,
#'   and the stored tree is in `tests/testthat/fixtures`.
#' @srrstats {G5.2, G5.2a, G5.2b} Every error and warning message is unique and
#'   is triggered in the tests.
#' @srrstats {G5.3} Outputs are tested for missing and undefined values.
#' @srrstats {G5.4, G5.4a} Probabilities and likelihoods are compared with
#'   explicit formulas.
#' @srrstats {G5.4b} Estimates match apTreeshape and laser (test-fit.R).
#' @srrstats {G5.5} All tests with random draws set a seed.
#' @srrstats {G5.6, G5.6a, G5.6b} Parameters are recovered within a tolerance
#'   over several seeds (more in test-extended.R).
#' @srrstats {G5.7} test-extended.R shows that the error of the estimates
#'   decreases with the number of tips.
#' @srrstats {G5.8, G5.8a, G5.8b, G5.8c, G5.8d} Edge conditions are tested.
#' @srrstats {G5.9, G5.9a, G5.9b} Noise susceptibility is tested.
#' @srrstats {G5.10, G5.12} Extended tests run with
#'   `TREESHAPESIM_EXTENDED_TESTS="true"`; their requirements are in
#'   CONTRIBUTING.md.
#' @srrstats {PD4.0, PD4.1} Numeric outputs are compared with explicit
#'   derivations.
#' @srrstats {PD4.3} Estimates are compared over alternative search intervals
#'   and optimizer tolerances.
#' @srrstats {PD4.4} Estimates are compared with an alternative optimizer.
#' @noRd
NULL

test_that("split probabilities match the beta-splitting formula", {
  # PD4.0, PD4.1: probabilities derived explicitly from Aldous (1996),
  # q(m, i) proportional to Gamma(beta + i + 1) Gamma(beta + m - i + 1) /
  # (Gamma(i + 1) Gamma(m - i + 1)), i = 1, ..., m - 1
  q <- function(beta, m) {
    i <- seq_len(m - 1)
    w <- gamma(beta + i + 1) * gamma(beta + m - i + 1) /
      (gamma(i + 1) * gamma(m - i + 1))
    w / sum(w)
  }
  # Yule (beta = 0): every split of a clade is equally likely
  expect_equal(exp(log_split_prob(0, 4)), rep(1 / 3, 3), tolerance = 1e-12)
  for (beta in c(-1.5, -1, 0.5, 3)) {
    for (m in c(3, 6, 15)) {
      expect_equal(exp(log_split_prob(beta, m)), q(beta, m), tolerance = 1e-10)
    }
  }
  # large clades do not overflow
  p <- exp(log_split_prob(-1.9, 2000))
  expect_true(all(is.finite(p)))
  expect_equal(sum(p), 1, tolerance = 1e-10)
})

test_that("the log-likelihood of beta sums the split log-probabilities", {
  # PD4.1: the caterpillar of 5 tips splits 5 -> (1, 4), 4 -> (1, 3),
  # 3 -> (1, 2); the 2-tip split has probability 1
  cat5 <- ape::stree(5, "left")
  sp <- split_sizes(cat5)
  expected <- log_split_prob(-1, 5)[sp[sp[, "m"] == 5, "i"]] +
    log_split_prob(-1, 4)[sp[sp[, "m"] == 4, "i"]] +
    log_split_prob(-1, 3)[sp[sp[, "m"] == 3, "i"]]
  expect_equal(loglik_beta(-1, sp), expected, tolerance = 1e-12)
})

test_that("the profile log-likelihood of c follows the DDX formula", {
  # PD4.1: log L = sum_k log(k * lambda_k) - sum_k k lambda_k g_k over the
  # events k = 2, ..., n - 1, with lambda_k = r k^c, minus the censored last
  # interval, and r at its maximum likelihood value
  g <- c(0.3, 0.25, 0.2, 0.15, 0.1)
  n <- length(g) + 1
  cc <- -0.7
  k <- 2:n
  ev <- 2:(n - 1)
  r <- (n - 2) / sum(k^(1 + cc) * g)
  expected <- sum(log(ev * r * ev^cc)) - sum(k * r * k^cc * g)
  expect_equal(loglik_c(cc, g), expected, tolerance = 1e-12)
})

test_that("estimates do not depend on the search interval or the tolerance", {
  # PD4.3: alternative values for every parameter of the optimizer
  set.seed(40)
  tr <- sim_tree_shape(150, beta = -0.5, c = -0.6)
  f <- fit_shape(tr)
  wide <- fit_shape(tr, beta_range = c(-1.9999, 50), c_range = c(-20, 20))
  narrow <- fit_shape(tr, beta_range = c(-1.5, 2), c_range = c(-2, 1))
  expect_equal(c(wide$beta, wide$c), c(f$beta, f$c), tolerance = 1e-4)
  expect_equal(c(narrow$beta, narrow$c), c(f$beta, f$c), tolerance = 1e-4)
  sp <- split_sizes(tr)
  g <- tree_intervals(tr)
  fine_b <- stats::optimize(loglik_beta, c(-1.999, 10), splits = sp,
                            maximum = TRUE, tol = 1e-10)$maximum
  fine_c <- stats::optimize(loglik_c, c(-5, 5), g = g, maximum = TRUE,
                            tol = 1e-10)$maximum
  expect_equal(c(fine_b, fine_c), c(f$beta, f$c), tolerance = 1e-4)
  # the optimizer is reported (PD3.3)
  expect_named(f$optimizer, c("algorithm", "tol", "evaluations"))
  expect_identical(f$optimizer$tol, 1e-6)
  expect_true(all(f$optimizer$evaluations > 0))
})

test_that("an alternative optimizer gives the same estimates", {
  # PD4.4: quasi-Newton search with box constraints instead of
  # golden-section search
  set.seed(41)
  tr <- sim_tree_shape(120, beta = -1, c = 0.3)
  f <- fit_shape(tr)
  sp <- split_sizes(tr)
  g <- tree_intervals(tr)
  b <- stats::optim(0, function(x) -loglik_beta(x, sp), method = "L-BFGS-B",
                    lower = -1.999, upper = 10)$par
  cc <- stats::optim(0, function(x) -loglik_c(x, g), method = "L-BFGS-B",
                     lower = -5, upper = 5)$par
  expect_equal(c(b, cc), c(f$beta, f$c), tolerance = 1e-3)
})

test_that("outputs hold no missing or undefined values", {
  # G5.3
  set.seed(42)
  tr <- sim_tree_shape(40, beta = -1, c = -0.5)
  expect_false(anyNA(tr$edge.length))
  expect_true(all(is.finite(tr$edge.length)))
  expect_true(all(is.finite(tree_shape(tr))))
  f <- fit_shape(tr, boot = 5)
  expect_true(all(is.finite(unlist(f[c("beta", "c", "r", "loglik_beta",
                                       "loglik_c")]))))
  expect_true(all(is.finite(as.matrix(f$boot))))
  expect_true(all(is.finite(intervals_lineage_rate(30, c = 1))))
  cal <- calibrate_shape(20, reps = 10)
  expect_false(anyNA(cal$gamma))
  expect_false(anyNA(cal$colless))
  expect_false(anyNA(extinct_tips(tr, 0.3)$edge.length))
})

test_that("edge conditions give clear errors", {
  # G5.8a zero-length input
  expect_error(sim_tree_shape(numeric(0)), "`n`")
  expect_error(fit_shape(tree = NULL), "`tree`")
  expect_error(shape_params(calibrate_shape(10, reps = 5), gamma = numeric(0),
                            colless = 0), "finite")
  # G5.8b unsupported types
  expect_error(sim_tree_shape("10"), "`n`")
  expect_error(sim_tree_shape(10, beta = "0"), "`beta`")
  expect_error(sim_tree_shape(10, c = 1i), "`c`")
  expect_error(tree_shape(matrix(1, 3, 3)), "`tree`")
  expect_error(assign_ages(list(), 1), "`topo`")
  # G5.8c all-NA input
  set.seed(43)
  tr <- sim_tree_shape(10)
  na_tree <- tr
  na_tree$edge.length[] <- NA_real_
  expect_error(fit_shape(na_tree), "branch lengths")
  expect_error(extinct_tips(na_tree, 0.3), "branch lengths")
  expect_error(sim_tree_shape(NA_real_), "`n`")
  expect_error(fit_shape(tr, boot = NA), "`boot`")
  # G5.8d data outside the scope of the algorithm
  expect_error(fit_shape(sim_tree_shape(3)), "at least 4")
  expect_error(fit_shape(ape::rtree(10)), "ultrametric")
  expect_error(extinct_tips(sim_tree_shape(2), 0.5), "at least 3")
})

test_that("every error and warning message is triggered", {
  # G5.2, G5.2b: messages not covered by the other test files
  set.seed(44)
  tr <- sim_tree_shape(20)
  cal <- calibrate_shape(20, reps = 10)
  expect_error(fit_shape(tr, tol = -1), "`tol`")
  expect_error(fit_shape(tr, tol = c(1, 2)), "`tol`")
  neg <- tr
  neg$edge.length[1] <- -0.1
  expect_error(fit_shape(neg), "non-negative, finite branch lengths")
  expect_error(sim_tree_target(20, gamma = NA_real_, colless = 0,
                               calibration = cal), "single finite numbers")
  expect_error(sim_tree_target(20, gamma = 0, colless = 0, tol = c(1, 1, 1),
                               calibration = cal), "`tol`")
  expect_error(sim_tree_target(20, gamma = 0, colless = 0, max_tries = 0,
                               calibration = cal), "`max_tries`")
  expect_error(calibrate_shape(20, reps = 0), "`reps`")
  expect_error(shape_grid(30, gamma = 0, colless = 0, calibration = cal),
               "same `n`")
  expect_error(sim_tree_shape(10, root_age = Inf), "`root_age`")
  # the calibration is computed when it is not given
  grid <- shape_grid(20, gamma = 0, colless = 0, reps = 20)
  expect_s3_class(attr(grid, "calibration"), "shape_calibration")
  # the warning printed when the corrected c is at the limit of the curve
  f <- fit_shape(sim_tree_shape(20, sampling_fraction = 0.5),
                 sampling_fraction = 0.5, correction_reps = 2)
  f$c_at_bound <- TRUE
  expect_output(print(f), "limit of the correction")
})

test_that("trivial noise in branch lengths does not change the estimates", {
  # G5.9a
  set.seed(45)
  tr <- sim_tree_shape(100, beta = -1, c = -0.5)
  noisy <- tr
  terminal <- noisy$edge[, 2] <= ape::Ntip(noisy)
  noisy$edge.length[terminal] <- noisy$edge.length[terminal] *
    (1 + stats::runif(sum(terminal), 0, 10 * .Machine$double.eps))
  f <- fit_shape(tr)
  fn <- fit_shape(noisy)
  expect_equal(c(fn$beta, fn$c), c(f$beta, f$c), tolerance = 1e-6)
})

test_that("the correction for sampling is stable across seeds", {
  # G5.9b: the correction is stochastic, but different seeds give close values
  set.seed(46)
  tr <- sim_tree_shape(120, c = -0.5, sampling_fraction = 0.5)
  est <- vapply(1:3, function(s) {
    set.seed(s)
    fit_shape(tr, sampling_fraction = 0.5, correction_reps = 30)$c
  }, 0)
  expect_lt(diff(range(est)), 0.2)
})

test_that("parameters are recovered across seeds", {
  # G5.6, G5.6a, G5.6b: recovery within a tolerance, over several seeds; the
  # extended version in test-extended.R uses more seeds and tree sizes
  est <- vapply(1:5, function(s) {
    set.seed(100 + s)
    f <- fit_shape(sim_tree_shape(300, beta = -1, c = -0.5))
    c(f$beta, f$c)
  }, numeric(2))
  expect_lt(abs(mean(est[1, ]) + 1), 0.25)
  expect_lt(abs(mean(est[2, ]) + 0.5), 0.15)
})
