#' Random topology under Aldous' beta-splitting model
#'
#' A clade of m tips splits into i and m - i tips with probability proportional
#' to Gamma(beta + i + 1) Gamma(beta + m - i + 1) / (Gamma(i + 1) Gamma(m - i + 1)).
#' `beta = 0` is Yule, `beta = -1.5` is PDA, `beta` close to -2 is
#' caterpillar-like and large `beta` is balanced.
#'
#' @param n Number of tips (at least 2).
#' @param beta Aldous' beta, greater than -2. Lower values give more imbalanced
#'   trees.
#' @return A `phylo` object without branch lengths, with tips `t1`, ..., `tn`.
#' @references Aldous, D. (1996). Probability distributions on cladograms. In
#'   D. Aldous & R. Pemantle (Eds.), *Random discrete structures* (pp. 1--18).
#'   Springer. \doi{10.1007/978-1-4612-0719-1_1}
#' @examples
#' set.seed(1)
#' rtopology_beta(10, beta = -1)
#' @export
rtopology_beta <- function(n, beta = 0) {
  stopifnot(n >= 2, beta > -2)
  split_clade <- function(tips) {
    m <- length(tips)
    if (m == 1) return(tips)
    i <- seq_len(m - 1)
    lw <- lgamma(beta + i + 1) + lgamma(beta + m - i + 1) - lgamma(i + 1) - lgamma(m - i + 1)
    k <- sample(i, 1, prob = exp(lw - max(lw)))
    paste0("(", split_clade(tips[seq_len(k)]), ",", split_clade(tips[-seq_len(k)]), ")")
  }
  ape::read.tree(text = paste0(split_clade(paste0("t", seq_len(n))), ";"))
}

#' Internode intervals of a pure-birth process with lineage-dependent rates
#'
#' With k lineages each splits at rate k^c, so the interval with k lineages is
#' Exp(k^(1 + c)). `c = 0` is Yule; `c < 0` slows diversification as lineages
#' accumulate (gamma < 0); `c > 0` accelerates it (gamma > 0). The final
#' interval (n lineages, up to the present) is drawn from the same law.
#'
#' @param n Number of tips.
#' @param c Lineage dependence of the split rate.
#' @return Intervals g_2, ..., g_n (length n - 1), scaled to sum 1.
#' @references Pybus, O. G., & Harvey, P. H. (2000). Testing macro-evolutionary
#'   models using incomplete molecular phylogenies. *Proceedings of the Royal
#'   Society B*, 267, 2267--2272. \doi{10.1098/rspb.2000.1278}
#' @examples
#' set.seed(1)
#' intervals_lineage_rate(10, c = -1)
#' @export
intervals_lineage_rate <- function(n, c = 0) {
  k <- 2:n
  g <- stats::rexp(n - 1, rate = k^(1 + c))
  g / sum(g)
}

#' Assign node ages to a topology from internode intervals
#'
#' Internal nodes are ranked by a uniform random linear extension of the tree
#' order (the ranking a Yule process implies given the topology): the next node
#' to split is drawn among the available ones with probability proportional to
#' the number of internal nodes in its subtree.
#'
#' @param topo Binary, rooted `phylo` topology.
#' @param g Internode intervals, as from [intervals_lineage_rate()]
#'   (length `Ntip(topo) - 1`).
#' @return `topo` with branch lengths; ultrametric, with root depth `sum(g)`.
#' @examples
#' set.seed(1)
#' assign_ages(rtopology_beta(10), intervals_lineage_rate(10))
#' @export
assign_ages <- function(topo, g) {
  n <- ape::Ntip(topo)
  stopifnot(length(g) == n - 1, ape::is.binary(topo))
  ages_sorted <- rev(cumsum(rev(g)))      # ages of the n - 1 splits, oldest first
  parent <- topo$edge[, 1]
  child <- topo$edge[, 2]
  n_int <- vapply(ape::prop.part(topo), length, 0) - 1   # internal nodes per subtree
  age <- numeric(n + topo$Nnode)
  available <- n + 1
  for (r in seq_len(n - 1)) {
    v <- if (length(available) == 1) available else
      sample(available, 1, prob = n_int[available - n])
    age[v] <- ages_sorted[r]
    kids <- child[parent == v]
    available <- c(setdiff(available, v), kids[kids > n])
  }
  topo$edge.length <- age[parent] - age[child]
  topo
}

#' Ultrametric tree with given imbalance and branching-time parameters
#'
#' Draws the topology with [rtopology_beta()] and the internode intervals with
#' [intervals_lineage_rate()], and combines them with [assign_ages()]. Because
#' the gamma statistic depends only on the intervals and the Colless index only
#' on the topology, `beta` mainly controls imbalance and `c` mainly controls
#' gamma.
#'
#' @inheritParams rtopology_beta
#' @inheritParams intervals_lineage_rate
#' @return A binary, ultrametric `phylo` object with root depth 1.
#' @examples
#' set.seed(1)
#' tr <- sim_tree_shape(50, beta = -1, c = -1)
#' ape::gammaStat(tr)
#' @export
sim_tree_shape <- function(n, beta = 0, c = 0) {
  assign_ages(rtopology_beta(n, beta), intervals_lineage_rate(n, c))
}
