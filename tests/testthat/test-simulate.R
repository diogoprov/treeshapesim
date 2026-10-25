tree_shape <- function(tree) {
  c(gamma = treestats::gamma_statistic(tree),
    colless = treestats::colless(tree, normalization = "yule"))
}

test_that("sim_tree_shape gives binary ultrametric trees with root depth 1", {
  set.seed(2)
  tr <- sim_tree_shape(60, beta = -1, c = -1)
  expect_equal(ape::Ntip(tr), 60)
  expect_true(ape::is.binary(tr))
  expect_true(ape::is.ultrametric(tr))
  expect_equal(max(ape::node.depth.edgelength(tr)), 1, tolerance = 1e-10)
  expect_true(all(tr$edge.length > 0))
})

test_that("gamma depends only on the intervals, not on the topology", {
  set.seed(3)
  g <- intervals_lineage_rate(40, c = -0.7)
  t1 <- assign_ages(rtopology_beta(40, beta = -1.8), g)
  t2 <- assign_ages(rtopology_beta(40, beta = 10), g)
  expect_equal(ape::gammaStat(t1), ape::gammaStat(t2), tolerance = 1e-10)
})

test_that("beta = 0 and c = 0 reproduce Yule expectations", {
  skip_if_not_installed("treestats")
  set.seed(4)
  s <- t(replicate(300, tree_shape(sim_tree_shape(50))))
  expect_lt(abs(mean(s[, "gamma"])), 0.2)   # gamma ~ N(0, 1) under Yule
  expect_lt(abs(mean(s[, "colless"])), 0.2) # Yule-normalized Colless ~ 0
})

test_that("Colless decreases with beta and gamma increases with c", {
  skip_if_not_installed("treestats")
  set.seed(5)
  m <- function(b, c) colMeans(t(replicate(100, tree_shape(sim_tree_shape(50, b, c)))))
  expect_gt(m(-1.2, 0)[["colless"]], m(1, 0)[["colless"]])
  expect_lt(m(0, -1)[["gamma"]], m(0, 0.5)[["gamma"]])
})

test_that("invalid arguments are rejected", {
  expect_error(rtopology_beta(1))
  expect_error(rtopology_beta(10, beta = -2))
  expect_error(assign_ages(rtopology_beta(10), intervals_lineage_rate(9)))
})
