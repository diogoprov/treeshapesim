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
#' \strong{Bootstrap.} With `boot > 0`, topologies are drawn with
#' [rtopology_beta()] at the estimated `beta` and intervals with
#' [intervals_lineage_rate()] at the estimated `c`, and both parameters are
#' re-estimated, following the recommendation of a parametric bootstrap for
#' diversification models by Etienne et al. (2016). The spread of the bootstrap
#' estimates around the point estimate includes the small-sample bias of the
#' estimators, and, for `c`, a small additional bias because
#' [intervals_lineage_rate()] draws the last interval from the same law,
#' whereas the likelihood treats it as censored. Percentile intervals need a
#' few hundred replicates.
#'
#' @param tree A binary, rooted, ultrametric `phylo` object with at least four
#'   tips.
#' @param boot Number of parametric bootstrap replicates (0 for none).
#' @param level Coverage of the bootstrap percentile intervals.
#' @param beta_range,c_range Search intervals for `beta` (lower bound above -2)
#'   and `c`.
#' @param tol Tolerance of the ultrametricity check, as in
#'   [ape::is.ultrametric()] (scaled range of root-to-tip distances). Trees with
#'   rounded branch lengths may need a larger value, or can be made ultrametric
#'   first (for instance with `phytools::force.ultrametric()`).
#' @return An object of class `shape_fit`: a list with `n`, `beta`, `c`, `r`
#'   (speciation rate per lineage at one lineage, in the time units of the
#'   tree), `loglik_beta` (log-likelihood of `beta` up to an additive constant,
#'   so not comparable with other software), `loglik_c` (profile log-likelihood
#'   of `c`, equal to the log-likelihood of laser's DDX model at its maximum),
#'   `beta_at_bound` and `c_at_bound` (estimate at a bound of the search
#'   interval), `beta_range`, `c_range`, `shape` (output of [tree_shape()]) and,
#'   with `boot > 0`, `boot` (data frame of bootstrap estimates) and `ci`
#'   (matrix of percentile intervals).
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
#' @examples
#' set.seed(1)
#' tr <- sim_tree_shape(200, beta = -1, c = -0.8)
#' fit_shape(tr)
#' fit <- fit_shape(tr, boot = 50)
#' fit$ci
#' @export
fit_shape <- function(tree, boot = 0, level = 0.95, beta_range = c(-1.999, 10),
                      c_range = c(-5, 5), tol = .Machine$double.eps^0.5) {
  if (!inherits(tree, "phylo") || !ape::is.binary(tree) || !ape::is.rooted(tree) ||
      ape::Ntip(tree) < 4) {
    stop("`tree` must be a binary, rooted `phylo` object with at least 4 tips.",
         call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1 || !is.finite(tol) || tol < 0) {
    stop("`tol` must be a single non-negative number.", call. = FALSE)
  }
  if (is.null(tree$edge.length) || !ape::is.ultrametric(tree, tol = tol)) {
    stop("`tree` must be ultrametric, with branch lengths (increase `tol` for rounded branch lengths, or use e.g. phytools::force.ultrametric()).",
         call. = FALSE)
  }
  if (!is.numeric(boot) || length(boot) != 1 || !is.finite(boot) || boot < 0 ||
      boot != round(boot)) {
    stop("`boot` must be a single non-negative integer.", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1 || !(level > 0 && level < 1)) {
    stop("`level` must be a single number between 0 and 1.", call. = FALSE)
  }
  check_range(beta_range, "beta_range")
  check_range(c_range, "c_range")
  if (beta_range[1] <= -2) stop("`beta_range` must lie above -2.", call. = FALSE)

  n <- ape::Ntip(tree)
  splits <- split_sizes(tree)
  g <- tree_intervals(tree)
  fb <- fit_beta(splits, beta_range)
  fc <- fit_c(g, c_range)
  out <- list(n = n, beta = fb$beta, c = fc$c, r = fc$r, loglik_beta = fb$loglik,
              loglik_c = fc$loglik, beta_at_bound = fb$at_bound,
              c_at_bound = fc$at_bound, beta_range = beta_range, c_range = c_range,
              shape = tree_shape(tree))
  if (boot > 0) {
    draws <- t(vapply(seq_len(boot), function(b) {
      c(beta = fit_beta(split_sizes(rtopology_beta(n, fb$beta)), beta_range)$beta,
        c = fit_c(intervals_lineage_rate(n, fc$c), c_range)$c)
    }, numeric(2)))
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
  if (isTRUE(x$c_at_bound)) {
    cat("c is at a bound of `c_range`: widen it, or the branching times carry",
        "little information (few tips).\n")
  }
  cat("gamma =", format(x$shape[["gamma"]], digits = digits),
      "| normalized Colless =", format(x$shape[["colless"]], digits = digits), "\n")
  invisible(x)
}

# clade sizes (m) and size of one child (i) at every internal node
split_sizes <- function(topo) {
  n <- ape::Ntip(topo)
  tips_below <- ape::node.depth(topo, method = 1)
  internal <- n + seq_len(topo$Nnode)
  kids <- split(topo$edge[, 2], factor(topo$edge[, 1], levels = internal))
  i <- vapply(kids, function(k) tips_below[k[1]], 0)
  m <- tips_below[internal]
  cbind(m = m, i = i)
}

# log-likelihood of beta (up to a constant) given the splits
loglik_beta <- function(beta, splits) {
  log_q <- function(m, i) {
    j <- seq_len(m - 1)
    lw <- lgamma(beta + j + 1) + lgamma(beta + m - j + 1) - lgamma(j + 1) - lgamma(m - j + 1)
    mx <- max(lw)
    lw[i] - (mx + log(sum(exp(lw - mx))))
  }
  # splits of 2 tips have probability 1 whatever beta
  keep <- splits[, "m"] > 2
  sum(mapply(log_q, splits[keep, "m"], splits[keep, "i"]))
}

fit_beta <- function(splits, beta_range) {
  opt <- stats::optimize(loglik_beta, beta_range, splits = splits, maximum = TRUE,
                         tol = 1e-6)
  list(beta = opt$maximum, loglik = opt$objective,
       at_bound = at_bound(opt$maximum, beta_range))
}

at_bound <- function(x, range) {
  eps <- 1e-3 * diff(range)
  x - range[1] < eps || range[2] - x < eps
}

# internode intervals g_2, ..., g_n of an ultrametric tree (g[k - 1] is the
# time with k lineages; the last one runs from the most recent split to the
# present)
tree_intervals <- function(tree) {
  bt <- sort(ape::branching.times(tree), decreasing = TRUE)
  c(-diff(bt), bt[length(bt)])
}

# DDX profile log-likelihood of c (r at its maximum likelihood value), as in
# laser::DDX() with x = -c
loglik_c <- function(c, g) {
  n <- length(g) + 1
  k <- 2:n
  ev <- 2:(n - 1)
  s <- sum(k^(1 + c) * g)
  r <- (n - 2) / s
  sum(log(ev)) + (n - 2) * log(r) + c * sum(log(ev)) - r * s
}

fit_c <- function(g, c_range) {
  opt <- stats::optimize(loglik_c, c_range, g = g, maximum = TRUE, tol = 1e-6)
  n <- length(g) + 1
  list(c = opt$maximum, loglik = opt$objective,
       r = (n - 2) / sum((2:n)^(1 + opt$maximum) * g),
       at_bound = at_bound(opt$maximum, c_range))
}

check_range <- function(x, name) {
  if (!is.numeric(x) || length(x) != 2 || any(!is.finite(x)) || x[1] >= x[2]) {
    stop(sprintf("`%s` must hold two increasing finite numbers.", name), call. = FALSE)
  }
}
