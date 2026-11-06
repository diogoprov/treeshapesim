# mean phylogenetic distance among the lost tips
mpd_lost <- function(tree, pruned) {
  d <- ape::cophenetic.phylo(tree)[attr(pruned, "lost"), attr(pruned, "lost")]
  mean(d[upper.tri(d)])
}

test_that("prune_tips removes the requested fraction and keeps the labels", {
  set.seed(40)
  tr <- sim_tree_shape(50, beta = -1, c = -0.5, root_age = 10)
  p <- prune_tips(tr, 0.3, signal = 0.5)
  expect_equal(ape::Ntip(p), 35)
  expect_length(attr(p, "lost"), 15)
  expect_setequal(c(p$tip.label, attr(p, "lost")), tr$tip.label)
  expect_true(ape::is.binary(p))
  expect_true(ape::is.ultrametric(p))
  # the tree is not rescaled: distances between survivors are unchanged
  sp <- p$tip.label
  expect_equal(ape::cophenetic.phylo(p)[sp, sp], ape::cophenetic.phylo(tr)[sp, sp])
})

test_that("signal clusters the losses on the tree", {
  set.seed(41)
  tr <- sim_tree_shape(80)
  m <- function(s) mean(replicate(100, mpd_lost(tr, prune_tips(tr, 0.3, signal = s))))
  expect_lt(m(1), m(0))
})

test_that("the liability has the covariance of Pagel's lambda", {
  set.seed(42)
  tr <- sim_tree_shape(6, root_age = 1)
  x <- replicate(20000, tip_liability(tr, 0.5))
  target <- ape::vcv.phylo(tr) * 0.5
  diag(target) <- 1
  expect_equal(unname(stats::cov(t(x))), unname(target[rownames(x), rownames(x)]),
               tolerance = 0.05)
  # lambda = 0 gives independent tips
  x0 <- replicate(5000, tip_liability(tr, 0))
  off <- stats::cor(t(x0))[upper.tri(diag(6))]
  expect_lt(max(abs(off)), 0.1)
})

test_that("risk links extinction to a trait", {
  set.seed(43)
  tr <- sim_tree_shape(60)
  trait <- stats::setNames(stats::rnorm(60), tr$tip.label)
  lost <- replicate(100, attr(prune_tips(tr, 0.25, risk = exp(2 * trait)), "lost"))
  expect_gt(mean(trait[lost]), mean(trait) + 0.5)
  # the order of the names does not matter, and zero risk protects a tip
  safe <- tr$tip.label[1:10]
  r <- stats::setNames(rep(1, 60), tr$tip.label)
  r[safe] <- 0
  p <- prune_tips(tr, 0.5, risk = rev(r))
  expect_true(all(safe %in% p$tip.label))
})

test_that("prune_tips checks its arguments", {
  set.seed(44)
  tr <- sim_tree_shape(10)
  expect_error(prune_tips(ape::unroot(tr), 0.5), "rooted")
  no_len <- tr
  no_len$edge.length <- NULL
  expect_error(prune_tips(no_len, 0.5), "branch lengths")
  expect_error(prune_tips(tr, 0), "`fraction`")
  expect_error(prune_tips(tr, 1), "`fraction`")
  expect_error(prune_tips(tr, 0.01), "at least one")
  expect_error(prune_tips(tr, 0.95), "leave at least two")
  expect_error(prune_tips(tr, 0.5, signal = 2), "`signal`")
  expect_error(prune_tips(tr, 0.5, risk = rep(1, 10)), "`risk`")
  r <- stats::setNames(rep(1, 10), tr$tip.label)
  expect_error(prune_tips(tr, 0.5, risk = r[1:9]), "`risk`")
  r[1] <- -1
  expect_error(prune_tips(tr, 0.5, risk = r), "`risk`")
  r[] <- 0
  r[1] <- 1
  expect_error(prune_tips(tr, 0.5, risk = r), "positive")
})
