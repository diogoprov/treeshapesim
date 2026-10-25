#' Calibrate the shape parameters for a number of tips
#'
#' Estimates, for a given number of tips, the mean gamma statistic as a
#' function of `c` and the mean normalized Colless index (see [tree_shape()])
#' as a function of `beta`. Because the topology and the internode intervals
#' are drawn independently, each curve is estimated from its own process alone:
#' gamma from intervals drawn with [intervals_lineage_rate()] and the Colless
#' index from topologies drawn with [rtopology_beta()], without building trees.
#' The expected values are monotone (gamma increases with `c` and the Colless
#' index decreases with `beta`), but the estimated means need not be, especially
#' for few tips or few draws, so each curve is made monotone by isotonic
#' regression ([stats::isoreg()]) before it is inverted with [shape_params()].
#'
#' @param n Number of tips, a single integer of at least 4 (with 3 tips there
#'   is a single topology, so the Colless index cannot be calibrated).
#' @param beta Values of `beta` at which the Colless index is estimated (at
#'   least two).
#' @param c Values of `c` at which gamma is estimated (at least two).
#' @param reps Number of topologies or interval sets per value.
#' @return An object of class `shape_calibration`: a list with `n`, `reps` and
#'   the data frames `colless` (columns `beta`, `colless` and `colless_mean`)
#'   and `gamma` (columns `c`, `gamma` and `gamma_mean`), where the `_mean`
#'   columns hold the mean statistic at each value and `colless` and `gamma`
#'   their monotone fit, which [shape_params()] inverts.
#' @examples
#' set.seed(1)
#' cal <- calibrate_shape(50, reps = 50)
#' cal$gamma
#' @export
calibrate_shape <- function(n, beta = c(-1.95, -1.8, -1.6, -1.4, -1.2, -1, -0.75,
                                        -0.5, 0, 0.5, 1, 2, 4, 8),
                            c = seq(-2, 1, by = 0.25), reps = 200) {
  check_n(n)
  if (n < 4) stop("`n` must be at least 4 to calibrate the Colless index.", call. = FALSE)
  beta <- sort(unique(beta))
  c <- sort(unique(c))
  if (length(beta) < 2 || length(c) < 2) {
    stop("`beta` and `c` must hold at least two distinct values each.", call. = FALSE)
  }
  for (b in beta) check_beta(b)
  for (v in c) check_c(v)
  check_count(reps, "reps")
  colless_mean <- vapply(beta, function(b)
    mean(replicate(reps, colless_yule(rtopology_beta(n, b)))), 0)
  gamma_mean <- vapply(c, function(v)
    mean(replicate(reps, gamma_intervals(intervals_lineage_rate(n, v)))), 0)
  # monotone fits: Colless decreasing in beta, gamma increasing in c
  colless <- -stats::isoreg(beta, -colless_mean)$yf
  gamma <- stats::isoreg(c, gamma_mean)$yf
  structure(list(n = as.integer(n), reps = as.integer(reps),
                 colless = data.frame(beta = beta, colless = colless,
                                      colless_mean = colless_mean),
                 gamma = data.frame(c = c, gamma = gamma, gamma_mean = gamma_mean)),
            class = "shape_calibration")
}

#' @export
print.shape_calibration <- function(x, ...) {
  cat("Shape calibration for n = ", x$n, " tips (", x$reps, " draws per value)\n", sep = "")
  cat("  mean gamma from", format(min(x$gamma$gamma), digits = 3), "to",
      format(max(x$gamma$gamma), digits = 3), "for c in [",
      min(x$gamma$c), ",", max(x$gamma$c), "]\n")
  cat("  mean Colless from", format(min(x$colless$colless), digits = 3), "to",
      format(max(x$colless$colless), digits = 3), "for beta in [",
      min(x$colless$beta), ",", max(x$colless$beta), "]\n")
  invisible(x)
}

#' Shape parameters that give target mean values of gamma and Colless
#'
#' Inverts the monotone curves of [calibrate_shape()] by linear interpolation.
#' Where a curve is flat (tied values of its monotone fit), the parameter is the
#' mean of the tied values.
#'
#' @param calibration Output of [calibrate_shape()].
#' @param gamma,colless Target mean gamma statistic and normalized Colless
#'   index (see [tree_shape()]); vectors that are recycled to a common length,
#'   as in [data.frame()].
#' @return A data frame with the targets and the corresponding `beta` and `c`.
#' @examples
#' set.seed(1)
#' cal <- calibrate_shape(50, reps = 50)
#' shape_params(cal, gamma = c(0, -2), colless = 0.5)
#' @export
shape_params <- function(calibration, gamma, colless) {
  if (!inherits(calibration, "shape_calibration")) {
    stop("`calibration` must come from calibrate_shape().", call. = FALSE)
  }
  if (!is.numeric(gamma) || !is.numeric(colless) || length(gamma) == 0 ||
      length(colless) == 0 || any(!is.finite(gamma)) || any(!is.finite(colless))) {
    stop("`gamma` and `colless` must be finite numbers.", call. = FALSE)
  }
  out <- data.frame(gamma = gamma, colless = colless)
  invert <- function(stat, value, target, what) {
    o <- order(stat)
    res <- stats::approx(stat[o], value[o], xout = target, rule = 1, ties = mean)$y
    if (anyNA(res)) {
      stop(sprintf("target %s outside the calibrated range [%.3g, %.3g]; widen the values of %s in calibrate_shape().",
                   what, min(stat), max(stat), if (what == "gamma") "`c`" else "`beta`"),
           call. = FALSE)
    }
    res
  }
  out$beta <- invert(calibration$colless$colless, calibration$colless$beta,
                     out$colless, "colless")
  out$c <- invert(calibration$gamma$gamma, calibration$gamma$c, out$gamma, "gamma")
  out
}

#' Simulate trees near target values of gamma and Colless
#'
#' Draws trees whose gamma statistic and normalized Colless index (see
#' [tree_shape()]) fall within a tolerance of target values. The parameters
#' `beta` and `c` whose mean statistics equal the targets are found with
#' [calibrate_shape()] and [shape_params()]. Topologies are then drawn with
#' [rtopology_beta()] and kept if their Colless index is within `tol` of the
#' target, interval sets are drawn with [intervals_lineage_rate()] and kept if
#' their gamma is within `tol` of the target, and the kept topologies and
#' interval sets are combined with [assign_ages()]. Because the topology and
#' the intervals are independent, rejecting them separately gives the same
#' distribution of trees as rejecting whole trees on both statistics, with a
#' much higher acceptance rate. Random numbers are used in a different order
#' than in [sim_tree_shape()], so the same seed does not give the same trees as
#' that function.
#'
#' @inheritParams rtopology_beta
#' @inheritParams shape_params
#' @param gamma,colless Target gamma statistic and normalized Colless index
#'   (single values).
#' @param nsim Number of trees.
#' @param tol Tolerance around the targets: one value for both statistics, or
#'   two values for gamma and Colless, in this order.
#' @param calibration Output of [calibrate_shape()] for the same `n`; if `NULL`,
#'   it is computed with `reps` draws per value.
#' @param reps Passed to [calibrate_shape()] when `calibration` is `NULL`.
#' @param max_tries Maximum number of topologies, and of interval sets, drawn.
#' @param root_age Age of the root.
#' @return With `nsim = 1`, a `phylo` object; otherwise a `multiPhylo` object
#'   (also when fewer trees than `nsim` are found).
#'   Attributes: `shape` (realized gamma and Colless of each tree), `params`
#'   (`beta` and `c` used) and `acceptance` (proportion of topologies and of
#'   interval sets kept). If fewer than `nsim` trees are found within
#'   `max_tries`, the trees found are returned with a warning.
#' @examples
#' set.seed(1)
#' cal <- calibrate_shape(50, reps = 50)
#' trees <- sim_tree_target(50, gamma = -1, colless = 0.5, nsim = 5,
#'                          calibration = cal)
#' attr(trees, "shape")
#' @export
sim_tree_target <- function(n, gamma, colless, nsim = 1, tol = 0.5, calibration = NULL,
                            reps = 200, max_tries = 1000 * nsim, root_age = 1) {
  check_n(n)
  if (!is.numeric(gamma) || length(gamma) != 1 || !is.finite(gamma) ||
      !is.numeric(colless) || length(colless) != 1 || !is.finite(colless)) {
    stop("`gamma` and `colless` must be single finite numbers.", call. = FALSE)
  }
  check_count(nsim, "nsim")
  check_count(max_tries, "max_tries")
  check_positive(root_age, "root_age")
  if (!is.numeric(tol) || !length(tol) %in% 1:2 || any(!is.finite(tol)) || any(tol <= 0)) {
    stop("`tol` must hold one or two positive numbers.", call. = FALSE)
  }
  tol <- rep(tol, length.out = 2)
  if (is.null(calibration)) {
    calibration <- calibrate_shape(n, reps = reps)
  } else if (!inherits(calibration, "shape_calibration") || calibration$n != n) {
    stop("`calibration` must come from calibrate_shape() with the same `n`.", call. = FALSE)
  }
  par <- shape_params(calibration, gamma = gamma, colless = colless)

  topologies <- list()
  colless_kept <- numeric(0)
  tries_topology <- 0
  while (length(topologies) < nsim && tries_topology < max_tries) {
    tries_topology <- tries_topology + 1
    tp <- rtopology_beta(n, par$beta)
    s <- colless_yule(tp)
    if (abs(s - colless) <= tol[2]) {
      topologies[[length(topologies) + 1]] <- tp
      colless_kept <- c(colless_kept, s)
    }
  }
  intervals <- list()
  gamma_kept <- numeric(0)
  tries_intervals <- 0
  while (length(intervals) < nsim && tries_intervals < max_tries) {
    tries_intervals <- tries_intervals + 1
    g <- intervals_lineage_rate(n, par$c)
    s <- gamma_intervals(g)
    if (abs(s - gamma) <= tol[1]) {
      intervals[[length(intervals) + 1]] <- g
      gamma_kept <- c(gamma_kept, s)
    }
  }
  k <- min(length(topologies), length(intervals))
  if (k < nsim) {
    warning(sprintf("only %d of %d trees within the tolerance after %d tries; increase `max_tries` or `tol`.",
                    k, nsim, max_tries), call. = FALSE)
  }
  if (k == 0) stop("no tree within the tolerance.", call. = FALSE)
  trees <- lapply(seq_len(k), function(i) {
    tr <- assign_ages(topologies[[i]], intervals[[i]])
    tr$edge.length <- tr$edge.length * root_age
    tr
  })
  out <- if (nsim == 1) trees[[1]] else structure(trees, class = "multiPhylo")
  attr(out, "shape") <- data.frame(gamma = gamma_kept[seq_len(k)],
                                   colless = colless_kept[seq_len(k)])
  attr(out, "params") <- c(beta = par$beta, c = par$c)
  attr(out, "acceptance") <- c(topology = length(topologies) / tries_topology,
                               intervals = length(intervals) / tries_intervals)
  out
}

#' Grid of tree shapes
#'
#' Crosses levels of the gamma statistic with levels of the normalized Colless
#' index and finds the `beta` and `c` that give each combination on average,
#' as a design for simulation studies. The levels can be chosen, for instance,
#' from quantiles of [tree_shape()] over empirical trees of similar size.
#'
#' @inheritParams sim_tree_target
#' @param gamma,colless Levels of each statistic.
#' @return A data frame with one row per combination of levels (`gamma`,
#'   `colless`, `beta`, `c`), with the calibration in attribute `calibration`.
#'   Each row can be passed to [sim_tree_target()] (targets) or to
#'   [sim_tree_shape()] (parameters).
#' @examples
#' set.seed(1)
#' grid <- shape_grid(50, gamma = c(-2, -1, 0), colless = c(0, 0.5, 1), reps = 50)
#' grid
#' @export
shape_grid <- function(n, gamma, colless, calibration = NULL, reps = 200) {
  check_n(n)
  if (is.null(calibration)) {
    calibration <- calibrate_shape(n, reps = reps)
  } else if (!inherits(calibration, "shape_calibration") || calibration$n != n) {
    stop("`calibration` must come from calibrate_shape() with the same `n`.", call. = FALSE)
  }
  grid <- expand.grid(gamma = gamma, colless = colless)
  out <- shape_params(calibration, gamma = grid$gamma, colless = grid$colless)
  attr(out, "calibration") <- calibration
  out
}
