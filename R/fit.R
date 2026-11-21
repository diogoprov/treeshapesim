#' Estimate the shape parameters of a tree
#'
#' Estimates, by maximum likelihood, the `beta` of Aldous' beta-splitting model
#' from the topology and the `c` of the lineage-dependent pure-birth process
#' from the branching times, which places an empirical tree on the same scale
#' as the trees of [sim_tree_shape()]. Optionally, a parametric bootstrap with
#' the package's own generators gives percentile intervals.
#'
#' @details
#' \strong{beta.} Each internal node with a clade of \eqn{m} tips split into
#' \eqn{i} and \eqn{m - i} tips contributes the probability of that split under
#' the beta-splitting model (see [rtopology_beta()]), and the log-likelihood is
#' the sum over nodes (Aldous 1996; this is the estimator applied by Blum &
#' Francois 2006 and implemented in `maxlik.betasplit()` of the archived package
#' apTreeshape). It is maximized over `beta_range`. In the package's validation
#' runs (`data-raw/validate_estimators.R` in the source repository), the
#' estimate was nearly unbiased for `beta` of -1 or less, but had a positive
#' bias and a large variance for balanced trees with few tips (`beta` of 0 or
#' more with 50 to 100 tips). An estimate at a bound of `beta_range` only says
#' that the topology is extreme for the model.
#'
#' \strong{c.} With \eqn{k} lineages, each lineage splits at rate
#' \eqn{r k^c}{r k^c}, which is the density-dependent exponential (DDX) model of
#' the R package laser (Rabosky 2006), whose per-lineage rate is
#' \eqn{r k^{-x}}{r k^(-x)}, so that \eqn{c = -x}. As in laser, the splits at
#' 2 to \eqn{n - 1} lineages are events and the last interval, from the most
#' recent split to the present, is not; \eqn{r} is profiled out (its maximum
#' likelihood value is \eqn{(n - 2) / \sum_k k^{1 + c} g_k}{(n - 2) / sum_k
#' k^(1 + c) g_k}, where \eqn{g_k} is the time with \eqn{k} lineages). The
#' estimate describes the shape of the lineage-through-time curve, not its
#' cause: branching times alone do not separate diversity-dependent from
#' time-dependent diversification (Pannetier et al. 2021).
#'
#' \strong{Incomplete sampling.} Trees with incomplete sampling give estimates
#' of `c` that are too low, because recent splits are missing; the estimate of
#' `beta` is little affected and is not corrected. With `sampling_fraction`
#' below 1, `c` is corrected by simulation: trees are drawn with
#' [sim_tree_shape()] at the estimated `beta`, the same `sampling_fraction` and
#' a grid of `c` from 1 below to 3 above the uncorrected estimate, and the mean
#' uncorrected estimate at each value (`correction_reps` trees) gives a curve,
#' made monotone by isotonic regression and inverted at the estimate from `tree`.
#' This assumes that tips are missing at random and that the fraction is known;
#' sampling biased towards old or distinct lineages is not corrected. The result
#' varies slightly with the seed.
#'
#' \strong{Bootstrap.} With `boot > 0`, topologies are drawn with
#' [rtopology_beta()] at the estimated `beta` and intervals with
#' [intervals_lineage_rate()] at the estimated `c` (with `sampling_fraction`
#' below 1, whole trees with [sim_tree_shape()], and `c` is corrected with the
#' same curve), and both parameters are re-estimated, following the recommendation of a parametric bootstrap for
#' diversification models by Etienne et al. (2016). The spread of the bootstrap
#' estimates around the point estimate includes the small-sample bias of the
#' estimators, and, for `c`, a small additional bias because
#' [intervals_lineage_rate()] draws the last interval from the same law,
#' whereas the likelihood treats it as censored. Percentile intervals need a
#' few hundred replicates.
#'
#' \strong{Optimization.} Each parameter is estimated separately with
#' [stats::optimize()], a golden-section search with successive parabolic
#' interpolation over the bounded interval `beta_range` or `c_range`, with a
#' tolerance of `1e-6` on the parameter. The search needs no starting value
#' and assumes a single maximum within the interval; an estimate at a bound is
#' flagged (`beta_at_bound`, `c_at_bound`). The name of the algorithm, the
#' tolerance and the number of evaluations of each log-likelihood are returned
#' in `optimizer`.
#'
#' @param tree A binary, rooted, ultrametric `phylo` object with at least four
#'   tips.
#' @param boot Number of parametric bootstrap replicates, a single
#'   non-negative integer (0, the default, for none).
#' @param level Coverage of the bootstrap percentile intervals, a single number
#'   between 0 and 1.
#' @param beta_range,c_range Search intervals for `beta` (lower bound above -2)
#'   and `c`: numeric vectors of two increasing finite values.
#' @param tol Tolerance of the ultrametricity check, a single non-negative
#'   number, as in
#'   [ape::is.ultrametric()] (scaled range of root-to-tip distances). Trees with
#'   rounded branch lengths may need a larger value, or can be made ultrametric
#'   first (for instance with `phytools::force.ultrametric()`).
#' @param sampling_fraction Fraction of the species of the clade that are in
#'   `tree`, a single number greater than 0 and at most 1. Below 1, `c` is
#'   corrected for incomplete sampling (see Details).
#' @param correction_reps Number of trees at each value of the grid used in the
#'   correction for incomplete sampling, a single positive integer.
#' @return An object of class `shape_fit`: a list with `n`, `beta`, `c`, `r`
#'   (speciation rate per lineage at one lineage, in the time units of the
#'   tree), `loglik_beta` (log-likelihood of `beta` up to an additive constant,
#'   so not comparable with other software), `loglik_c` (profile log-likelihood
#'   of `c`, equal to the log-likelihood of laser's DDX model at its maximum),
#'   `beta_at_bound` and `c_at_bound` (estimate at a bound of the search
#'   interval), `beta_range`, `c_range`, `shape` (output of [tree_shape()]),
#'   `optimizer` (list with the `algorithm`, its `tol` and the number of
#'   `evaluations` of each log-likelihood) and,
#'   with `boot > 0`, `boot` (data frame of bootstrap estimates) and `ci`
#'   (matrix of percentile intervals). With `sampling_fraction` below 1, `c` is
#'   the corrected estimate, `r` is `NA`, `loglik_c` refers to the uncorrected
#'   estimate, and the list also holds `c_uncorrected`, `sampling_fraction` and
#'   `correction` (data frame with the curve).
#' @references
#' Aldous, D. (1996). Probability distributions on cladograms. In D. Aldous &
#' R. Pemantle (Eds.), *Random discrete structures* (pp. 1--18). Springer.
#' \doi{10.1007/978-1-4612-0719-1_1}
#'
#' Blum, M. G. B., & Francois, O. (2006). Which random processes describe the
#' tree of life? A large-scale study of phylogenetic tree imbalance.
#' *Systematic Biology*, 55, 685--691. \doi{10.1080/10635150600889625}
#'
#' Etienne, R. S., Pigot, A. L., & Phillimore, A. B. (2016). How reliably can
#' we infer diversity-dependent diversification from phylogenies? *Methods in
#' Ecology and Evolution*, 7, 1092--1099. \doi{10.1111/2041-210X.12565}
#'
#' Pannetier, T., Martinez, C., Bunnefeld, L., & Etienne, R. S. (2021).
#' Branching patterns in phylogenies cannot distinguish diversity-dependent
#' diversification from time-dependent diversification. *Evolution*, 75,
#' 25--38. \doi{10.1111/evo.14124}
#'
#' Rabosky, D. L. (2006). LASER: a maximum likelihood toolkit for detecting
#' temporal shifts in diversification rates from molecular phylogenies.
#' *Evolutionary Bioinformatics*, 2. \doi{10.1177/117693430600200024}
#' @srrstats {G1.0} Primary references for both estimators are listed.
#' @srrstats {G1.3} The likelihoods, the correction for sampling and the
#'   bootstrap are defined in Details.
#' @srrstats {G2.0a, G2.1a} Expected lengths and types are given for every
#'   argument.
#' @srrstats {PD1.0} The beta-splitting and DDX models are justified with
#'   references (Aldous 1996; Rabosky 2006; Etienne et al. 2016).
#' @srrstats {PD3.2, PD3.3} The optimizer, its tolerance and the number of
#'   evaluations are documented and returned in `optimizer`.
#' @examples
#' set.seed(1)
#' tr <- sim_tree_shape(200, beta = -1, c = -0.8)
#' fit_shape(tr)
#' fit <- fit_shape(tr, boot = 50)
#' fit$ci
#' # a tree with half of the species of the clade
#' sampled <- sim_tree_shape(100, c = -0.5, sampling_fraction = 0.5)
#' fit_shape(sampled, sampling_fraction = 0.5, correction_reps = 20)
#' @export
fit_shape <- function(tree, boot = 0, level = 0.95, beta_range = c(-1.999, 10),
                      c_range = c(-5, 5), tol = .Machine$double.eps^0.5,
                      sampling_fraction = 1, correction_reps = 100) {
  check_fit_args(tree, boot, level, beta_range, c_range, tol, sampling_fraction,
                 correction_reps)
  n <- ape::Ntip(tree)
  splits <- split_sizes(tree)
  g <- tree_intervals(tree)
  fb <- fit_beta(splits, beta_range)
  fc <- fit_c(g, c_range)
  out <- list(n = n, beta = fb$beta, c = fc$c, r = fc$r, loglik_beta = fb$loglik,
              loglik_c = fc$loglik, beta_at_bound = fb$at_bound,
              c_at_bound = fc$at_bound, beta_range = beta_range, c_range = c_range,
              shape = tree_shape(tree),
              optimizer = list(algorithm = optimizer_name, tol = optimizer_tol,
                               evaluations = c(beta = fb$evaluations,
                                               c = fc$evaluations)))
  sampled <- sampling_fraction < 1
  if (sampled) {
    curve <- sampling_curve(n, fb$beta, fc$c, sampling_fraction, c_range,
                            correction_reps)
    out$c_uncorrected <- fc$c
    out$c <- invert_curve(curve, fc$c)
    out$r <- NA_real_
    out$c_at_bound <- at_bound(out$c, range(curve$c))
    out$sampling_fraction <- sampling_fraction
    out$correction <- curve
  }
  if (boot > 0) {
    draws <- boot_draws(boot, n, fb$beta, fc$c, out$c, beta_range, c_range,
                        sampling_fraction, out$correction)
    probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
    out$boot <- as.data.frame(draws)
    out$ci <- rbind(beta = stats::quantile(draws[, "beta"], probs),
                    c = stats::quantile(draws[, "c"], probs))
  }
  structure(out, class = "shape_fit")
}

#' @export
print.shape_fit <- function(x, digits = 3, ...) {
  cat("Shape estimates for a tree with", x$n, "tips\n")
  sampled <- !is.null(x$sampling_fraction)
  if (sampled) {
    cat("c corrected for a sampling fraction of ", x$sampling_fraction,
        " (uncorrected: ", format(x$c_uncorrected, digits = digits), ")\n", sep = "")
  }
  est <- c(beta = x$beta, c = x$c)
  if (!is.null(x$ci)) {
    tab <- cbind(estimate = est, x$ci)
    print(signif(tab, digits))
  } else {
    print(signif(est, digits))
  }
  if (isTRUE(x$beta_at_bound)) {
    upper <- abs(x$beta - x$beta_range[2]) < abs(x$beta - x$beta_range[1])
    cat("beta is at the", if (upper) "upper" else "lower",
        "bound of `beta_range`: the topology is extremely",
        if (upper) "balanced" else "imbalanced", "for this model.\n")
  }
  if (isTRUE(x$c_at_bound) && sampled) {
    cat("c is at the limit of the correction for sampling: treat it as",
        "unreliable.\n")
  } else if (isTRUE(x$c_at_bound)) {
    cat("c is at a bound of `c_range`: widen it, or the branching times carry",
        "little information (few tips).\n")
  }
  cat("gamma =", format(x$shape[["gamma"]], digits = digits),
      "| normalized Colless =", format(x$shape[["colless"]], digits = digits), "\n")
  invisible(x)
}

#' Optimizer used by fit_shape()
#'
#' Both parameters are estimated with [stats::optimize()], a golden-section
#' search with successive parabolic interpolation (Brent 1973) over a bounded
#' interval, with the tolerance below; it needs no starting value.
#' @srrstats {PD3.2} The optimizer, its tolerance and the search intervals are
#'   fixed here and documented in [fit_shape()].
#' @noRd
optimizer_name <- "stats::optimize (golden-section search, parabolic interpolation)"
optimizer_tol <- 1e-6

#' Check the arguments of fit_shape()
#'
#' @srrstats {G2.0, G2.1, G2.2} Lengths and types of all arguments are asserted.
#' @noRd
check_fit_args <- function(tree, boot, level, beta_range, c_range, tol,
                           sampling_fraction, correction_reps) {
  check_tree(tree, min_tips = 4, lengths = TRUE)
  check_number(tol, "`tol` must be a single non-negative number.", lower = 0)
  check_ultrametric(tree, tol)
  check_number(boot, "`boot` must be a single non-negative integer.",
               lower = 0, whole = TRUE)
  check_number(level, "`level` must be a single number between 0 and 1.",
               lower = 0, upper = 1, open = c(TRUE, TRUE))
  check_range(beta_range, "beta_range")
  check_range(c_range, "c_range")
  if (beta_range[1] <= -2) stop("`beta_range` must lie above -2.", call. = FALSE)
  check_fraction(sampling_fraction)
  check_count(correction_reps, "correction_reps")
}

#' Check that a tree is ultrametric within a tolerance
#'
#' @noRd
check_ultrametric <- function(tree, tol) {
  if (!ape::is.ultrametric(tree, tol = tol)) {
    stop("`tree` must be ultrametric (increase `tol` for rounded branch ",
         "lengths, or use e.g. phytools::force.ultrametric()).", call. = FALSE)
  }
}

#' Parametric bootstrap estimates of beta and c
#'
#' @param boot Number of replicates.
#' @param n Number of tips.
#' @param beta,c_hat Estimates used to simulate (uncorrected `c`).
#' @param c_sim Value of `c` used to simulate whole sampled trees (corrected).
#' @param curve Correction curve, or `NULL` without sampling.
#' @return A matrix with columns `beta` and `c`, one row per replicate.
#' @noRd
boot_draws <- function(boot, n, beta, c_hat, c_sim, beta_range, c_range,
                       sampling_fraction, curve) {
  one <- function(b) {
    if (is.null(curve)) {
      return(c(beta = fit_beta(split_sizes(rtopology_beta(n, beta)),
                               beta_range)$beta,
               c = fit_c(intervals_lineage_rate(n, c_hat), c_range)$c))
    }
    tr <- suppressWarnings(sim_tree_shape(n, beta, c_sim,
                                          sampling_fraction = sampling_fraction))
    c(beta = fit_beta(split_sizes(tr), beta_range)$beta,
      c = invert_curve(curve, fit_c(tree_intervals(tr), c_range)$c))
  }
  t(vapply(seq_len(boot), one, numeric(2)))
}

#' Clade sizes at the internal nodes
#'
#' @param topo A binary `phylo` object, with or without branch lengths.
#' @return A two-column matrix with the clade size `m` and the size `i` of the
#'   first child at every internal node.
#' @noRd
split_sizes <- function(topo) {
  n <- ape::Ntip(topo)
  tips_below <- ape::node.depth(topo, method = 1)
  internal <- n + seq_len(topo$Nnode)
  kids <- split(topo$edge[, 2], factor(topo$edge[, 1], levels = internal))
  i <- vapply(kids, function(k) tips_below[k[1]], 0)
  m <- tips_below[internal]
  cbind(m = m, i = i)
}

#' Log-probabilities of the splits of a clade under beta-splitting
#'
#' Analytic probabilities (Aldous 1996), normalized on the log scale with the
#' log-sum-exp trick so that large clades do not overflow.
#'
#' @param beta Value of Aldous' beta at which the probabilities are computed.
#' @param m Clade size (at least 2).
#' @return Log-probabilities of splits into `1, ..., m - 1` and `m - i` tips.
#' @srrstats {PD3.0} Probabilities are computed analytically.
#' @noRd
log_split_prob <- function(beta, m) {
  j <- seq_len(m - 1)
  lw <- lgamma(beta + j + 1) + lgamma(beta + m - j + 1) - lgamma(j + 1) -
    lgamma(m - j + 1)
  mx <- max(lw)
  lw - (mx + log(sum(exp(lw - mx))))
}

#' Log-likelihood of beta given the splits of a tree
#'
#' @param beta Value of Aldous' beta at which the likelihood is computed.
#' @param splits Output of `split_sizes()`.
#' @return The log-likelihood; splits of two tips have probability 1 whatever
#'   `beta` and are left out.
#' @noRd
loglik_beta <- function(beta, splits) {
  keep <- splits[, "m"] > 2
  sum(mapply(function(m, i) log_split_prob(beta, m)[i],
             splits[keep, "m"], splits[keep, "i"]))
}

#' Maximum likelihood estimate of beta
#'
#' @return A list with `beta`, `loglik`, `at_bound` and `evaluations` (number
#'   of evaluations of the log-likelihood).
#' @srrstats {PD3.3} The number of evaluations is returned with the estimate.
#' @noRd
fit_beta <- function(splits, beta_range) {
  counter <- new.env(parent = emptyenv())
  counter$calls <- 0L
  f <- function(beta) {
    counter$calls <- counter$calls + 1L
    loglik_beta(beta, splits)
  }
  opt <- stats::optimize(f, beta_range, maximum = TRUE, tol = optimizer_tol)
  list(beta = opt$maximum, loglik = opt$objective,
       at_bound = at_bound(opt$maximum, beta_range),
       evaluations = counter$calls)
}

#' Is an estimate at a bound of its search interval?
#'
#' @return `TRUE` if `x` is within 0.1% of the width of `range` from a bound.
#' @noRd
at_bound <- function(x, range) {
  eps <- 1e-3 * diff(range)
  x - range[1] < eps || range[2] - x < eps
}

#' Internode intervals of an ultrametric tree
#'
#' @return Intervals `g_2, ..., g_n`, where `g[k - 1]` is the time with `k`
#'   lineages; the last one runs from the most recent split to the present.
#' @noRd
tree_intervals <- function(tree) {
  bt <- sort(ape::branching.times(tree), decreasing = TRUE)
  c(-diff(bt), bt[length(bt)])
}

#' Profile log-likelihood of c
#'
#' The density-dependent exponential (DDX) log-likelihood of laser::DDX() with
#' `x = -c`, with the rate `r` at its maximum likelihood value.
#'
#' @param c Value of `c` at which the likelihood is computed.
#' @param g Internode intervals, as from `tree_intervals()`.
#' @return The profile log-likelihood.
#' @srrstats {PD3.0} The likelihood and the profiled rate are analytic.
#' @noRd
loglik_c <- function(c, g) {
  n <- length(g) + 1
  k <- 2:n
  ev <- 2:(n - 1)
  s <- sum(k^(1 + c) * g)
  r <- (n - 2) / s
  sum(log(ev)) + (n - 2) * log(r) + c * sum(log(ev)) - r * s
}

#' Maximum likelihood estimate of c
#'
#' @return A list with `c`, `loglik`, `r`, `at_bound` and `evaluations`.
#' @noRd
fit_c <- function(g, c_range) {
  counter <- new.env(parent = emptyenv())
  counter$calls <- 0L
  f <- function(c) {
    counter$calls <- counter$calls + 1L
    loglik_c(c, g)
  }
  opt <- stats::optimize(f, c_range, maximum = TRUE, tol = optimizer_tol)
  n <- length(g) + 1
  list(c = opt$maximum, loglik = opt$objective,
       r = (n - 2) / sum((2:n)^(1 + opt$maximum) * g),
       at_bound = at_bound(opt$maximum, c_range),
       evaluations = counter$calls)
}

#' Curve of the uncorrected estimate of c under incomplete sampling
#'
#' @param n Number of tips of the sampled tree.
#' @param beta,c_hat Estimates from the sampled tree.
#' @param reps Number of trees at each value of the grid.
#' @return A data frame with the true `c` on the grid and the mean uncorrected
#'   estimate, made monotone by isotonic regression.
#' @noRd
sampling_curve <- function(n, beta, c_hat, sampling_fraction, c_range, reps) {
  grid <- unique(pmin(pmax(c_hat + seq(-1, 3, by = 0.5), c_range[1]),
                      c_range[2]))
  estimate <- function(v) {
    tr <- suppressWarnings(sim_tree_shape(n, beta, v,
                                          sampling_fraction = sampling_fraction))
    fit_c(tree_intervals(tr), c_range)$c
  }
  mean_c <- vapply(grid, function(v) mean(replicate(reps, estimate(v))), 0)
  data.frame(c = grid, uncorrected = stats::isoreg(grid, mean_c)$yf)
}

#' Invert the correction curve at an uncorrected estimate
#'
#' @return The corrected `c`, clamped to the range of the grid.
#' @noRd
invert_curve <- function(curve, x) {
  stats::approx(curve$uncorrected, curve$c, xout = x, rule = 2, ties = mean)$y
}
