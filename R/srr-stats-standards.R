#' Standards documented elsewhere
#'
#' rOpenSci statistical software standards (general and probability
#' distributions) that apply to the package as a whole. Standards that apply to
#' a function are tagged in its documentation, and testing standards in the
#' test files.
#'
#' @srrstatsVerbose TRUE
#'
#' @srrstats {G1.4} All functions are documented with roxygen2.
#' @srrstats {G1.4a} Internal functions are documented with roxygen2 and
#'   `@noRd`.
#' @srrstats {G1.5} The validation of the estimators and of the correction for
#'   sampling, and the benchmarks behind the performance claims, are scripts in
#'   `data-raw/` of the source repository (`validate_estimators.R`,
#'   `validate_sampling.R`, `benchmark_topology.R`), with their results in
#'   `data-raw/validation/`.
#' @srrstats {G1.6} `data-raw/validate_estimators.R` compares [fit_shape()]
#'   with `apTreeshape::maxlik.betasplit()` and `laser::DDX()`, and the README
#'   compares the package with TreeSim, TreeSimGM and apTreeshape.
#' @srrstats {G2.4} Inputs are numeric scalars, vectors, functions and `phylo`
#'   objects; the only conversion needed is that of whole numbers to integer
#'   (G2.4a, in [rtopology_beta()]).
#' @srrstats {G2.8} Trees pass through `check_tree()`, and numeric arguments
#'   through the shared checks in `R/simulate.R`, before any computation.
#' @srrstats {PD1.0} References for the beta-splitting model and the
#'   lineage-dependent pure-birth process are given in each function.
#' @noRd
NULL

#' NA_standards
#'
#' @srrstatsNA {G2.3, G2.3a, G2.3b} No argument takes character input.
#' @srrstatsNA {G2.4b, G2.4c, G2.4d, G2.4e} No argument needs conversion to or
#'   from continuous, character or factor types.
#' @srrstatsNA {G2.5} No argument takes factors.
#' @srrstatsNA {G2.7, G2.9, G2.10, G2.11, G2.12} No function takes tabular
#'   input; the data are trees (`phylo` objects).
#' @srrstatsNA {G2.14b, G2.14c} Missing values are always an error: a tree with
#'   a missing branch length has no defined shape, so values cannot be ignored
#'   or imputed.
#' @srrstatsNA {G3.1, G3.1a} No computation uses covariances.
#' @srrstatsNA {G4.0} No function writes files.
#' @srrstatsNA {G5.4c} Correctness is tested against previous implementations
#'   (G5.4b) and explicit formulas (G5.4a), not against published values.
#' @srrstatsNA {G5.11, G5.11a} The extended tests use simulated trees and need
#'   no downloads.
#' @srrstatsNA {PD2.0} The distributions are over trees (topologies and
#'   internode intervals), which general packages for probability
#'   distributions (such as distributional) do not represent.
#' @srrstatsNA {PD3.1, PD4.2} No function takes the name of a distribution as
#'   an argument: each model is fixed.
#' @srrstatsNA {PD3.4, PD3.5, PD3.5a} No function integrates a distribution.
#' @noRd
NULL
