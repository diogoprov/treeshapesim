#' Remove tips by non-random extinction
#'
#' Applies one extinction event to a tree and removes a fixed fraction of its
#' tips. Tips are lost according to `risk`, which links extinction to a trait or
#' to any other property of the species, or, without `risk`, according to a
#' latent liability that evolves on the tree, so that `signal` sets how
#' clustered the losses are: `signal = 0` is random extinction, and with
#' `signal = 1` close relatives tend to be lost together, up to whole clades.
#'
#' Under a homogeneous birth-death process the shape of the reconstructed tree is
#' the same as under pure birth (Gernhard 2008), so extinction changes the
#' topology only if it differs among lineages. The event acts on the tree as it
#' stands, at the present; it does not model extinction that differs among
#' lineages throughout the history of the clade, for which speciation and
#' extinction have to be simulated together (Maddison et al. 2007; Verboom et al.
#' 2020). Risk of extinction clustered on the tree does not always lead to more
#' loss of phylogenetic diversity than random extinction (Veron et al. 2017), so
#' the same `fraction` with `signal = 0` is the natural control.
#'
#' Unlike `sampling_fraction` in [sim_tree_shape()], the tree is not rescaled and
#' the tip labels are kept, so that traits stay aligned with the species that
#' remain. If a whole daughter clade of the root is lost, the root of the
#' remaining tree is younger. The separation between gamma and the Colless index
#' holds before the event, not after it, and [calibrate_shape()] does not apply
#' to the remaining tree.
#'
#' @param tree A binary, rooted `phylo` object with branch lengths.
#' @param fraction Fraction of the tips that are lost, a single number greater
#'   than 0 and less than 1. `round(fraction * Ntip(tree))` tips are removed, and
#'   at least two must remain.
#' @param risk `NULL`, or a numeric vector of non-negative weights named with the
#'   tip labels of `tree`. Tips are drawn without replacement with probability
#'   proportional to `risk`. To link extinction to a trait, use for instance
#'   `exp(b * trait)`.
#' @param signal Pagel's lambda of the latent liability, a single number from 0
#'   to 1, used when `risk` is `NULL`. The tips with the highest liability are
#'   lost.
#' @return The remaining tree, with the labels of the lost tips in attribute
#'   `lost`.
#' @references
#' Gernhard, T. (2008). The conditioned reconstructed process. *Journal of
#' Theoretical Biology*, 253, 769--778. \doi{10.1016/j.jtbi.2008.04.005}
#'
#' Maddison, W. P., Midford, P. E., & Otto, S. P. (2007). Estimating a binary
#' character's effect on speciation and extinction. *Systematic Biology*, 56,
#' 701--710. \doi{10.1080/10635150701607033}
#'
#' Verboom, G. A., Boucher, F. C., Ackerly, D. D., Wootton, L. M., & Freyman,
#' W. A. (2020). Species selection regime and phylogenetic tree shape.
#' *Systematic Biology*, 69, 774--794. \doi{10.1093/sysbio/syz076}
#'
#' Veron, S., Davies, T. J., Cadotte, M. W., Clergeau, P., & Pavoine, S. (2017).
#' Predicting loss of evolutionary history: Where are we? *Biological Reviews*,
#' 92, 271--291. \doi{10.1111/brv.12228}
#' @examples
#' set.seed(1)
#' tr <- sim_tree_shape(50, beta = -1)
#' clustered <- prune_tips(tr, 0.4, signal = 1)
#' random <- prune_tips(tr, 0.4, signal = 0)
#' attr(clustered, "lost")
#' # extinction linked to a trait
#' body_size <- setNames(rnorm(50), tr$tip.label)
#' prune_tips(tr, 0.4, risk = exp(2 * body_size))
#' @export
prune_tips <- function(tree, fraction, risk = NULL, signal = 1) {
  if (!inherits(tree, "phylo") || !ape::is.binary(tree) || !ape::is.rooted(tree) ||
      is.null(tree$edge.length)) {
    stop("`tree` must be a binary, rooted `phylo` object with branch lengths.",
         call. = FALSE)
  }
  if (!is.numeric(fraction) || length(fraction) != 1 || !is.finite(fraction) ||
      fraction <= 0 || fraction >= 1) {
    stop("`fraction` must be a single number between 0 and 1.", call. = FALSE)
  }
  tips <- tree$tip.label
  n_lost <- round(fraction * length(tips))
  if (n_lost < 1 || length(tips) - n_lost < 2) {
    stop("`fraction` must remove at least one tip and leave at least two.", call. = FALSE)
  }
  if (is.null(risk)) {
    if (!is.numeric(signal) || length(signal) != 1 || !is.finite(signal) ||
        signal < 0 || signal > 1) {
      stop("`signal` must be a single number from 0 to 1.", call. = FALSE)
    }
    liability <- tip_liability(tree, signal)
    lost <- names(sort(liability, decreasing = TRUE))[seq_len(n_lost)]
  } else {
    if (!is.numeric(risk) || is.null(names(risk)) || !setequal(names(risk), tips) ||
        anyDuplicated(names(risk)) > 0 || any(!is.finite(risk)) || any(risk < 0)) {
      stop("`risk` must hold one non-negative, finite value per tip, named with the tip labels.",
           call. = FALSE)
    }
    if (sum(risk > 0) < n_lost) {
      stop("`risk` must be positive for at least as many tips as are lost.", call. = FALSE)
    }
    risk <- risk[tips]
    lost <- sample(tips, n_lost, prob = risk)
  }
  out <- ape::drop.tip(tree, lost)
  attr(out, "lost") <- lost
  out
}

# Brownian liability at the tips under Pagel's lambda: internal branches are
# multiplied by lambda, and terminal branches extended to keep root-to-tip depths
tip_liability <- function(tree, lambda) {
  n <- ape::Ntip(tree)
  depth <- ape::node.depth.edgelength(tree)
  terminal <- tree$edge[, 2] <= n
  tree$edge.length <- tree$edge.length * lambda
  tree$edge.length[terminal] <- tree$edge.length[terminal] +
    (1 - lambda) * depth[tree$edge[terminal, 2]]
  x <- ape::rTraitCont(tree, model = "BM", sigma = 1)
  x[tree$tip.label]
}
