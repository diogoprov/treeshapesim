#' Random topology under Aldous' beta-splitting model
#'
#' A clade of m tips splits into i and m - i tips with probability proportional
#' to Gamma(beta + i + 1) Gamma(beta + m - i + 1) / (Gamma(i + 1) Gamma(m - i + 1)).
#' `beta = 0` is Yule, `beta = -1.5` is PDA, `beta` close to -2 is
#' caterpillar-like and large `beta` is balanced.
#'
#' @param n Number of tips, a single integer of at least 2.
#' @param beta Aldous' beta, a single finite number greater than -2. Lower
#'   values give more imbalanced trees.
#' @return A `phylo` object without branch lengths, with tips `t1`, ..., `tn`.
#' @references Aldous, D. (1996). Probability distributions on cladograms. In
#'   D. Aldous & R. Pemantle (Eds.), *Random discrete structures* (pp. 1--18).
#'   Springer. \doi{10.1007/978-1-4612-0719-1_1}
#' @examples
#' set.seed(1)
#' rtopology_beta(10, beta = -1)
#' @export
rtopology_beta <- function(n, beta = 0) {
  check_n(n)
  check_beta(beta)
  n <- as.integer(n)
  # Clades are split depth first, left before right, in the order of a
  # recursive construction, with an explicit stack instead of recursion (which
  # fails for large, very imbalanced trees). Clade = tips start..start + m - 1.
  # Internal nodes are numbered and edges listed in preorder, as ape::read.tree()
  # does for the equivalent Newick string, so the same seed gives the same tree.
  edge <- matrix(0L, 2 * n - 2, 2)
  n_edge <- 0L
  next_node <- n + 1L
  stack <- list(c(start = 1L, m = n, parent = 0L))
  while (length(stack) > 0) {
    clade <- stack[[length(stack)]]
    stack[[length(stack)]] <- NULL
    m <- clade[["m"]]
    if (m == 1) {
      id <- clade[["start"]]
    } else {
      id <- next_node
      next_node <- next_node + 1L
    }
    if (clade[["parent"]] > 0) {
      n_edge <- n_edge + 1L
      edge[n_edge, ] <- c(clade[["parent"]], id)
    }
    if (m > 1) {
      i <- seq_len(m - 1)
      lw <- lgamma(beta + i + 1) + lgamma(beta + m - i + 1) - lgamma(i + 1) - lgamma(m - i + 1)
      k <- sample(i, 1, prob = exp(lw - max(lw)))
      start <- clade[["start"]]
      stack[[length(stack) + 1]] <- c(start = start + k, m = m - k, parent = id)
      stack[[length(stack) + 1]] <- c(start = start, m = k, parent = id)
    }
  }
  tree <- list(edge = edge, Nnode = n - 1L, tip.label = paste0("t", seq_len(n)))
  class(tree) <- "phylo"
  attr(tree, "order") <- "cladewise"
  tree
}

#' Internode intervals of a pure-birth process with lineage-dependent rates
#'
#' With k lineages each splits at rate k^c, so the interval with k lineages is
#' Exp(k^(1 + c)). `c = 0` is Yule; `c < 0` slows diversification as lineages
#' accumulate (gamma < 0); `c > 0` accelerates it (gamma > 0). The final
#' interval (n lineages, up to the present) is drawn from the same law.
#'
#' Any other dependence on the number of lineages can be given with `rate`, a
#' function of k that returns the per-lineage split rate; the interval with k
#' lineages is then Exp(k * rate(k)), and `c` is the special case
#' `rate = function(k) k^c`. Only the shape of `rate` matters: the intervals are
#' rescaled, so multiplying it by a constant does not change the tree. As with
#' `c`, the shape of `rate` describes the lineage-through-time curve and does
#' not identify the mechanism behind it (Pannetier et al. 2021).
#' [calibrate_shape()], [shape_params()] and [sim_tree_target()] work with `c`
#' only; with `rate`, the calibration is left to the user.
#'
#' @param n Number of tips, a single integer of at least 2.
#' @param c Lineage dependence of the split rate, a single finite number.
#' @param rate `NULL` (the default), to use `c`, or a vectorized function of the
#'   number of lineages k that returns the per-lineage split rate. It must
#'   return n - 1 positive, finite values for k = 2, ..., n, and cannot be
#'   combined with a `c` other than 0.
#' @return Intervals g_2, ..., g_n (length n - 1), scaled to sum 1.
#' @references Pannetier, T., Martinez, C., Bunnefeld, L., & Etienne, R. S.
#'   (2021). Branching patterns in phylogenies cannot distinguish
#'   diversity-dependent diversification from time-dependent diversification.
#'   *Evolution*, 75, 25--38. \doi{10.1111/evo.14124}
#'
#'   Pybus, O. G., & Harvey, P. H. (2000). Testing macro-evolutionary
#'   models using incomplete molecular phylogenies. *Proceedings of the Royal
#'   Society B*, 267, 2267--2272. \doi{10.1098/rspb.2000.1278}
#' @examples
#' set.seed(1)
#' intervals_lineage_rate(10, c = -1)
#' # logistic slowdown towards a carrying capacity K > n
#' intervals_lineage_rate(10, rate = function(k) 1 - k / 12)
#' # the rate drops once there are more than 5 lineages
#' intervals_lineage_rate(10, rate = function(k) ifelse(k <= 5, 1, 0.1))
#' @export
intervals_lineage_rate <- function(n, c = 0, rate = NULL) {
  check_n(n)
  check_c(c)
  g <- stats::rexp(n - 1, rate = total_rates(n, c, rate))
  g / sum(g)
}

# total split rates with k = 2, ..., n lineages, after checking `rate`
total_rates <- function(n, c, rate) {
  k <- 2:n
  if (is.null(rate)) {
    total <- k^(1 + c)
    if (any(!is.finite(total)) || any(total <= 0)) {
      stop(sprintf("`c` = %g is too extreme for %d tips: the split rates overflow or underflow.",
                   c, n), call. = FALSE)
    }
    return(total)
  }
  if (c != 0) stop("Use either `c` or `rate`, not both.", call. = FALSE)
  msg <- paste("`rate` must be a vectorized function that returns n - 1",
               "positive, finite values for k = 2, ..., n.")
  if (!is.function(rate)) stop(msg, call. = FALSE)
  r <- rate(k)
  if (!is.numeric(r) || length(r) != n - 1) stop(msg, call. = FALSE)
  total <- k * as.vector(r)
  if (any(!is.finite(total)) || any(total <= 0)) stop(msg, call. = FALSE)
  total
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
#'   (length `Ntip(topo) - 1`, positive and finite).
#' @return `topo` with branch lengths; ultrametric, with root depth `sum(g)`.
#' @examples
#' set.seed(1)
#' assign_ages(rtopology_beta(10), intervals_lineage_rate(10))
#' @export
assign_ages <- function(topo, g) {
  if (!inherits(topo, "phylo") || !ape::is.binary(topo) || !ape::is.rooted(topo)) {
    stop("`topo` must be a binary, rooted `phylo` object.", call. = FALSE)
  }
  n <- ape::Ntip(topo)
  if (!is.numeric(g) || length(g) != n - 1 || any(!is.finite(g)) || any(g <= 0)) {
    stop("`g` must hold Ntip(topo) - 1 positive, finite intervals.", call. = FALSE)
  }
  ages_sorted <- rev(cumsum(rev(g)))      # ages of the n - 1 splits, oldest first
  parent <- topo$edge[, 1]
  child <- topo$edge[, 2]
  internal <- n + seq_len(topo$Nnode)
  # internal nodes per subtree, and the children of each node in edge order
  n_int <- ape::node.depth(topo, method = 1)[internal] - 1
  kids_of <- split(child, factor(parent, levels = internal))
  age <- numeric(n + topo$Nnode)
  available <- n + 1
  for (r in seq_len(n - 1)) {
    v <- if (length(available) == 1) available else
      sample(available, 1, prob = n_int[available - n])
    age[v] <- ages_sorted[r]
    kids <- kids_of[[v - n]]
    available <- c(available[available != v], kids[kids > n])
  }
  topo$edge.length <- age[parent] - age[child]
  topo
}

#' Ultrametric trees with given imbalance and branching-time parameters
#'
#' Draws the topology with [rtopology_beta()] and the internode intervals with
#' [intervals_lineage_rate()], and combines them with [assign_ages()]. The
#' topology and the intervals are drawn independently, the gamma statistic
#' depends only on the intervals and the Colless index only on the topology, so
#' for a given `n` the realized gamma depends only on `c` (or `rate`) and the
#' realized Colless index only on `beta`. To draw trees near target values of
#' both statistics, see [sim_tree_target()].
#'
#' With `sampling_fraction` below 1, a complete tree of `round(n /
#' sampling_fraction)` tips is drawn, `beta` and `c` (or `rate`) apply to that
#' tree, and `n` tips are kept at random, so the result still has `n` tips and
#' root depth `root_age`. Incomplete sampling removes recent splits and so lowers
#' gamma (Pybus & Harvey 2000), and it changes the Colless index too, so the
#' realized statistics no longer follow the calibration of [calibrate_shape()].
#' A `rate` function must be valid for the number of tips of the complete tree.
#'
#' Very extreme values of `c` (or `rate`) make some intervals negligible relative
#' to the others, and the corresponding branches round to length zero; a warning
#' is given when this happens.
#'
#' @inheritParams rtopology_beta
#' @inheritParams intervals_lineage_rate
#' @param nsim Number of trees.
#' @param root_age Age of the root (positive).
#' @param sampling_fraction Fraction of the tips of the complete tree that are
#'   kept, a single number greater than 0 and at most 1.
#' @return With `nsim = 1`, a binary, ultrametric `phylo` object with root
#'   depth `root_age`; otherwise a `multiPhylo` object with `nsim` such trees.
#'   Attribute `params` holds `beta`, `c` (`NA` with `rate`), `root_age` and
#'   `sampling_fraction`, and attribute `rate` holds the `rate` function, if any.
#' @references Pybus, O. G., & Harvey, P. H. (2000). Testing macro-evolutionary
#'   models using incomplete molecular phylogenies. *Proceedings of the Royal
#'   Society B*, 267, 2267--2272. \doi{10.1098/rspb.2000.1278}
#' @section Reproducibility: The same seed gives the same trees on the same
#'   platform. Between platforms (for instance, arm64 and x86-64), differences
#'   in the last bit of floating-point functions can change a few draws of the
#'   topology, so trees may differ.
#' @examples
#' set.seed(1)
#' tr <- sim_tree_shape(50, beta = -1, c = -1)
#' ape::gammaStat(tr)
#' trees <- sim_tree_shape(50, beta = -1, c = -1, nsim = 10, root_age = 30)
#' # logistic slowdown towards a carrying capacity of 60 lineages
#' sim_tree_shape(50, rate = function(k) 1 - k / 60)
#' # 50 tips sampled from a complete tree of 100
#' sim_tree_shape(50, sampling_fraction = 0.5)
#' @export
sim_tree_shape <- function(n, beta = 0, c = 0, nsim = 1, root_age = 1, rate = NULL,
                           sampling_fraction = 1) {
  check_n(n)
  check_beta(beta)
  check_c(c)
  check_count(nsim, "nsim")
  check_positive(root_age, "root_age")
  check_fraction(sampling_fraction)
  n_full <- round(n / sampling_fraction)
  check_n(n_full)
  total_rates(n_full, c, rate)
  one <- function() {
    tr <- assign_ages(rtopology_beta(n_full, beta), intervals_lineage_rate(n_full, c, rate))
    if (n_full > n) tr <- sample_tips(tr, n)
    tr$edge.length <- tr$edge.length * root_age
    tr
  }
  trees <- replicate(nsim, one(), simplify = FALSE)
  n_zero <- sum(vapply(trees, function(tr) any(tr$edge.length <= 0), logical(1)))
  if (n_zero > 0) {
    warning(sprintf("%d of %d trees have branches of length zero: `c` or `rate` is too extreme.",
                    n_zero, nsim), call. = FALSE)
  }
  out <- if (nsim == 1) trees[[1]] else structure(trees, class = "multiPhylo")
  attr(out, "params") <- c(beta = beta, c = if (is.null(rate)) c else NA_real_,
                           root_age = root_age, sampling_fraction = sampling_fraction)
  attr(out, "rate") <- rate
  out
}

# keep n random tips, with root depth 1 and tips relabelled t1, ..., tn
sample_tips <- function(tree, n) {
  tree <- ape::drop.tip(tree, sample(tree$tip.label, ape::Ntip(tree) - n))
  tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
  tree$tip.label <- paste0("t", seq_len(n))
  tree
}

# argument checks shared by the exported functions
check_n <- function(n) {
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n < 2 || n != round(n) ||
      n > .Machine$integer.max) {
    stop("`n` must be a single integer of at least 2.", call. = FALSE)
  }
}
check_beta <- function(beta) {
  if (!is.numeric(beta) || length(beta) != 1 || !is.finite(beta) || beta <= -2) {
    stop("`beta` must be a single finite number greater than -2.", call. = FALSE)
  }
}
check_c <- function(c) {
  if (!is.numeric(c) || length(c) != 1 || !is.finite(c)) {
    stop("`c` must be a single finite number.", call. = FALSE)
  }
}
check_fraction <- function(x) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x <= 0 || x > 1) {
    stop("`sampling_fraction` must be a single number greater than 0 and at most 1.",
         call. = FALSE)
  }
}
check_count <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < 1 || x != round(x)) {
    stop(sprintf("`%s` must be a single positive integer.", name), call. = FALSE)
  }
}
check_positive <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x <= 0) {
    stop(sprintf("`%s` must be a single positive, finite number.", name), call. = FALSE)
  }
}
