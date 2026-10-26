test_that("fit_shape recovers the parameters of large trees", {
  set.seed(30)
  tr <- sim_tree_shape(400, beta = -1, c = -0.8)
  f <- fit_shape(tr)
  expect_s3_class(f, "shape_fit")
  expect_lt(abs(f$beta + 1), 0.3)
  expect_lt(abs(f$c + 0.8), 0.2)
  expect_false(f$beta_at_bound)
  expect_equal(f$shape, tree_shape(tr))
})

test_that("estimates are stable (values checked against apTreeshape and laser)", {
  # beta matched apTreeshape::maxlik.betasplit() and c matched laser::DDX()
  # to within 3e-5 on 30 random trees (both packages archived on CRAN)
  set.seed(31)
  f <- fit_shape(sim_tree_shape(80, -1.2, -0.6))
  expect_equal(c(f$beta, f$c, f$r), c(-1.2207825020, -0.4522262078, 13.0757228117),
               tolerance = 1e-6)
})

test_that("the profiled r maximizes the DDX likelihood", {
  set.seed(32)
  g <- intervals_lineage_rate(60, -0.5)
  n <- length(g) + 1
  full <- function(p) {
    r <- exp(p[1]); x <- p[2]; ev <- 2:(n - 1)
    sum(log(ev)) + (n - 2) * log(r) - x * sum(log(ev)) - r * sum((2:n)^(1 - x) * g)
  }
  opt <- stats::optim(c(0, 0), full, control = list(fnscale = -1, reltol = 1e-12))
  fc <- fit_c(g, c(-5, 5))
  expect_equal(fc$c, -opt$par[2], tolerance = 1e-4)
  expect_equal(fc$loglik, opt$value, tolerance = 1e-6)
  expect_equal(log(fc$r), opt$par[1], tolerance = 1e-4)
})

test_that("a caterpillar tree puts beta at the lower bound", {
  tr <- ape::compute.brlen(ape::stree(30, "left"), method = "Grafen")
  f <- fit_shape(tr)
  expect_true(f$beta_at_bound)
  expect_lt(f$beta, -1.99)
  expect_output(print(f), "imbalanced")
})

test_that("the bootstrap gives intervals", {
  set.seed(33)
  f <- fit_shape(sim_tree_shape(60, -1, -0.5), boot = 20, level = 0.9)
  expect_equal(dim(f$boot), c(20L, 2L))
  expect_equal(rownames(f$ci), c("beta", "c"))
  expect_true(all(f$ci[, 1] <= f$ci[, 2]))
  expect_output(print(f), "estimate")
})

test_that("fit_shape checks its arguments", {
  set.seed(34)
  tr <- sim_tree_shape(20)
  expect_error(fit_shape(ape::unroot(tr)), "rooted")
  expect_error(fit_shape(ape::di2multi(ape::collapse.singles(tr), tol = 0.05)), "binary")
  expect_error(fit_shape(sim_tree_shape(3)), "at least 4")
  no_len <- tr
  no_len$edge.length <- NULL
  expect_error(fit_shape(no_len), "ultrametric")
  expect_error(fit_shape(ape::rtree(20)), "ultrametric")
  expect_error(fit_shape(tr, boot = -1), "`boot`")
  expect_error(fit_shape(tr, level = 1), "`level`")
  expect_error(fit_shape(tr, beta_range = c(-2, 1)), "above -2")
  expect_error(fit_shape(tr, c_range = c(1, 0)), "`c_range`")
})

test_that("small trees, bounds and rounded branch lengths", {
  set.seed(35)
  small <- fit_shape(sim_tree_shape(4))
  expect_true(is.finite(small$beta) && is.finite(small$c))
  # a narrow c_range is hit and reported
  narrow <- fit_shape(sim_tree_shape(50, c = -1.5), c_range = c(-0.5, 0.5))
  expect_true(narrow$c_at_bound)
  expect_output(print(narrow), "c is at a bound")
  # the lower bound of beta_range is reported as such, whatever its sign
  cat_tree <- ape::compute.brlen(ape::stree(10, "left"), method = "Grafen")
  expect_output(print(fit_shape(cat_tree, beta_range = c(0.5, 10))), "lower bound")
  # rounded branch lengths need a larger tolerance
  tr <- sim_tree_shape(100, root_age = 100)
  tr$edge.length <- round(tr$edge.length, 3)
  expect_error(fit_shape(tr), "ultrametric")
  expect_s3_class(fit_shape(tr, tol = 1e-3), "shape_fit")
})
