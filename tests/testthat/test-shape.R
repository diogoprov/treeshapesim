test_that("tree_shape matches ape and treestats", {
  set.seed(11)
  for (r in 1:20) {
    n <- sample(3:150, 1)
    g <- intervals_lineage_rate(n, stats::runif(1, -1.5, 1))
    tp <- rtopology_beta(n, stats::runif(1, -1.95, 5))
    tr <- assign_ages(tp, g)
    s <- tree_shape(tr)
    expect_equal(s[["gamma"]], ape::gammaStat(tr), tolerance = 1e-10)
    expect_equal(gamma_intervals(g), s[["gamma"]], tolerance = 1e-10)
    expect_identical(colless_yule(tp), s[["colless"]])
    if (requireNamespace("treestats", quietly = TRUE)) {
      expect_equal(s[["colless"]], treestats::colless(tr, normalization = "yule"),
                   tolerance = 1e-10)
    }
  }
  expect_error(tree_shape(ape::rtree(2)), "at least 3")
})
