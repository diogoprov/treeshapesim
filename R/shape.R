#' Gamma statistic and Colless index of a tree
#'
#' Computes Pybus and Harvey's (2000) gamma statistic, with [ape::gammaStat()],
#' and the Colless index normalized against its expectation under the Yule
#' model, \eqn{(I_C - n \log n - n(\gamma_E - 1 - \log 2)) / n}{(Ic - n log n -
#' n (gammaE - 1 - log 2)) / n}, where \eqn{n} is the number of tips and
#' \eqn{\gamma_E}{gammaE} is Euler's constant (the `"yule"` normalization of
#' `treestats::colless()`). Both are near zero for Yule trees: gamma is
#' negative when branching times are concentrated near the root, and the
#' normalized Colless index is positive for trees more imbalanced than Yule
#' trees.
#'
#' @param tree A binary, rooted `phylo` object with at least three tips. The
#'   gamma statistic assumes an ultrametric tree; this is not checked.
#' @return A named vector with `gamma` and `colless`.
#' @references Pybus, O. G., & Harvey, P. H. (2000). Testing macro-evolutionary
#'   models using incomplete molecular phylogenies. *Proceedings of the Royal
#'   Society B*, 267, 2267--2272. \doi{10.1098/rspb.2000.1278}
#' @examples
#' set.seed(1)
#' tree_shape(sim_tree_shape(100, beta = -1, c = -1))
#' @export
tree_shape <- function(tree) {
  if (!inherits(tree, "phylo") || !ape::is.binary(tree) || !ape::is.rooted(tree) ||
      ape::Ntip(tree) < 3) {
    stop("`tree` must be a binary, rooted `phylo` object with at least 3 tips.",
         call. = FALSE)
  }
  c(gamma = ape::gammaStat(tree), colless = colless_yule(tree))
}

# Colless index of a binary topology, normalized against its Yule expectation
colless_yule <- function(topo) {
  n <- ape::Ntip(topo)
  tips_below <- ape::node.depth(topo, method = 1)
  internal <- n + seq_len(topo$Nnode)
  kids <- split(topo$edge[, 2], factor(topo$edge[, 1], levels = internal))
  ic <- sum(vapply(kids, function(k) abs(tips_below[k[1]] - tips_below[k[2]]), 0))
  euler <- -digamma(1)
  (ic - n * log(n) - n * (euler - 1 - log(2))) / n
}

# Gamma statistic from the internode intervals g_2, ..., g_n of
# intervals_lineage_rate() (the same value as ape::gammaStat() on any tree with
# these intervals)
gamma_intervals <- function(g) {
  n <- length(g) + 1
  k <- 2:n
  weighted <- k * g
  total <- sum(weighted)
  inner <- sum(cumsum(weighted)[seq_len(n - 2)]) / (n - 2)
  (inner - total / 2) / (total * sqrt(1 / (12 * (n - 2))))
}
