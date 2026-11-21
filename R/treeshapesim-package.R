#' treeshapesim: Simulate Trees with Separately Controlled Imbalance and
#' Branching Times
#'
#' Simulates rooted, binary, ultrametric trees in which topological imbalance
#' and the distribution of branching times are set by separate parameters, and
#' estimates those parameters from empirical trees.
#'
#' @section Terms:
#' \describe{
#'   \item{Topology}{The branching pattern of a tree, without branch lengths.}
#'   \item{Imbalance}{How unequally the tips are divided between the two
#'     daughter clades of the internal nodes, measured here by the Colless
#'     index normalized against its expectation under the Yule model (see
#'     [tree_shape()]).}
#'   \item{beta}{The parameter of Aldous' (1996) beta-splitting model of
#'     topologies: `beta = 0` gives the distribution of Yule topologies, lower
#'     values more imbalanced and higher values more balanced topologies.}
#'   \item{Internode intervals}{The times \eqn{g_k} during which the tree has
#'     \eqn{k} lineages, from the root (\eqn{k = 2}) to the present
#'     (\eqn{k = n}).}
#'   \item{c}{The exponent of the dependence of the per-lineage split rate on
#'     the number of lineages \eqn{k}, \eqn{k^c}; `c = 0` is a Yule process.}
#'   \item{gamma}{Pybus and Harvey's (2000) statistic of the internode
#'     intervals: negative when splits are concentrated near the root.}
#'   \item{Sampling fraction}{The proportion of the species of a clade that are
#'     in the tree.}
#' }
#'
#' @section Algorithms:
#' The combination of a beta-splitting topology with intervals from a
#' lineage-dependent pure-birth process, joined by a uniform random ranking of
#' the internal nodes, is, as far as we know, a new algorithm; this is its
#' first implementation. The estimators of `beta` (maximum likelihood under
#' beta-splitting) and of `c` (the density-dependent exponential model) were
#' implemented in the R packages apTreeshape and laser, both archived on CRAN;
#' [fit_shape()] reimplements them and adds a correction of `c` for incomplete
#' sampling, which is new.
#'
#' @section Life cycle:
#' The package is in active development and its interface is
#' \emph{experimental} (version 0.x): function names and arguments may still
#' change, with changes listed in `NEWS.md`. The core simulators
#' ([sim_tree_shape()], [rtopology_beta()], [intervals_lineage_rate()],
#' [assign_ages()]) are expected to be stable; planned work, including a
#' sampler with extinction and a CRAN release, is listed in `ROADMAP.md` in the
#' source repository.
#'
#' @references
#' Aldous, D. (1996). Probability distributions on cladograms. In D. Aldous &
#' R. Pemantle (Eds.), *Random discrete structures* (pp. 1--18). Springer.
#' \doi{10.1007/978-1-4612-0719-1_1}
#'
#' Pybus, O. G., & Harvey, P. H. (2000). Testing macro-evolutionary models
#' using incomplete molecular phylogenies. *Proceedings of the Royal Society
#' B*, 267, 2267--2272. \doi{10.1098/rspb.2000.1278}
#'
#' @srrstats {G1.0} Primary references are listed here and in each function.
#' @srrstats {G1.1} The algorithms are described as new or as
#'   reimplementations in the section Algorithms.
#' @srrstats {G1.2} The life cycle statement is in the section Life cycle and
#'   in the README.
#' @srrstats {G1.3} Statistical terms are defined in the section Terms.
#' @keywords internal
"_PACKAGE"
