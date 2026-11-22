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
#' @param n Number of tips, a single positive integer of at least 4 (with 3
#'   tips there is a single topology, so the Colless index cannot be
#'   calibrated).
#' @param beta Values of `beta` at which the Colless index is estimated (at
#'   least two).
#' @param c Values of `c` at which gamma is estimated (at least two).
#' @param reps Number of topologies or interval sets per value, a single
#'   positive integer.
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
calibrate_shape <- function(n,
                            beta = c(-1.95, -1.8, -1.6, -1.4, -1.2, -1, -0.75,
                                     -0.5, 0, 0.5, 1, 2, 4, 8),
                            c = seq(-2, 1, by = 0.25), reps = 200) {
  check_n(n)
  if (n < 4) {
    stop("`n` must be at least 4 to calibrate the Colless index.",
         call. = FALSE)
  }
  beta <- sort(unique(beta))
  c <- sort(unique(c))
  if (length(beta) < 2 || length(c) < 2) {
    stop("`beta` and `c` must hold at least two distinct values each.",
         call. = FALSE)
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
                 gamma = data.frame(c = c, gamma = gamma,
                                    gamma_mean = gamma_mean)),
            class = "shape_calibration")
}

#' @export
print.shape_calibration <- function(x, ...) {
  cat("Shape calibration for n = ", x$n, " tips (", x$reps,
      " draws per value)\n", sep = "")
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
  if (!is_finite_numeric(gamma) || !is_finite_numeric(colless)) {
    stop("`gamma` and `colless` must be finite numbers.", call. = FALSE)
  }
  out <- data.frame(gamma = gamma, colless = colless)
  out$beta <- invert_stat(calibration$colless$colless,
                          calibration$colless$beta, out$colless, "colless")
  out$c <- invert_stat(calibration$gamma$gamma, calibration$gamma$c,
                       out$gamma, "gamma")
  out
}

#' Invert a calibration curve by linear interpolation
#'
#' @param stat,value Mean statistic and the parameter value that gives it.
#' @param target Target values of the statistic.
#' @param what `"gamma"` or `"colless"`, for the error message.
#' @return The parameter values for `target`; an error outside the curve.
#' @noRd
invert_stat <- function(stat, value, target, what) {
  o <- order(stat)
  res <- stats::approx(stat[o], value[o], xout = target, rule = 1,
                       ties = mean)$y
  if (anyNA(res)) {
    param <- c(gamma = "`c`", colless = "`beta`")[[what]]
    stop(sprintf(paste("target %s outside the calibrated range [%.3g, %.3g];",
                       "widen the values of %s in calibrate_shape()."),
                 what, min(stat), max(stat), param),
         call. = FALSE)
  }
  res
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
#' @inheritParams sim_tree_shape
#' @param gamma,colless Target gamma statistic and normalized Colless index
#'   (single values).
#' @param tol Tolerance around the targets: one value for both statistics, or
#'   two values for gamma and Colless, in this order.
#' @param calibration Output of [calibrate_shape()] for the same `n`; if `NULL`,
#'   it is computed with `reps` draws per value.
#' @param reps Passed to [calibrate_shape()] when `calibration` is `NULL`.
#' @param max_tries Maximum number of topologies, and of interval sets, drawn,
#'   a single positive integer.
#' @return With `nsim = 1`, a `phylo` object; otherwise a `multiPhylo` object
#'   (also when fewer trees than `nsim` are found).
#'   Attributes: `shape` (realized gamma and Colless of each tree), `params`
#'   (`beta` and `c` used) and `acceptance` (proportion of topologies and of
#'   interval sets kept). If fewer than `nsim` trees are found within
#'   `max_tries`, the trees found are returned with a warning.
#' @srrstats {G2.0a, G2.1a} Expected lengths and types are given for every
#'   argument.
#' @examples
#' set.seed(1)
#' cal <- calibrate_shape(50, reps = 50)
#' trees <- sim_tree_target(50, gamma = -1, colless = 0.5, nsim = 5,
#'                          calibration = cal)
#' attr(trees, "shape")
#' @export
sim_tree_target <- function(n, gamma, colless, nsim = 1, tol = 0.5,
                            calibration = NULL, reps = 200,
                            max_tries = 1000 * nsim, root_age = 1) {
  check_n(n)
  if (!is_number(gamma) || !is_number(colless)) {
    stop("`gamma` and `colless` must be single finite numbers.", call. = FALSE)
  }
  check_count(nsim, "nsim")
  check_count(max_tries, "max_tries")
  check_positive(root_age, "root_age")
  tol <- check_tol(tol)
  calibration <- get_calibration(calibration, n, reps)
  par <- shape_params(calibration, gamma = gamma, colless = colless)
  topo <- draw_kept(function() rtopology_beta(n, par$beta), colless_yule,
                    colless, tol[2], nsim, max_tries)
  ints <- draw_kept(function() intervals_lineage_rate(n, par$c),
                    gamma_intervals, gamma, tol[1], nsim, max_tries)
  k <- min(length(topo$kept), length(ints$kept))
  if (k < nsim) {
    warning(sprintf(paste("only %d of %d trees within the tolerance after %d",
                          "tries; increase `max_tries` or `tol`."),
                    k, nsim, max_tries), call. = FALSE)
  }
  if (k == 0) stop("no tree within the tolerance.", call. = FALSE)
  combine_kept(topo, ints, k, nsim, root_age, par)
}

#' Check the tolerance of sim_tree_target()
#'
#' @return `tol` recycled to length 2 (gamma, Colless).
#' @noRd
check_tol <- function(tol) {
  if (!is_finite_numeric(tol) || length(tol) > 2 || any(tol <= 0)) {
    stop("`tol` must hold one or two positive numbers.", call. = FALSE)
  }
  rep(tol, length.out = 2)
}

#' Combine kept topologies and interval sets into trees
#'
#' @param topo,ints Outputs of draw_kept() for topologies and intervals.
#' @param k Number of trees to build.
#' @param par Parameters from shape_params().
#' @return A `phylo` (`nsim = 1`) or `multiPhylo` object with attributes
#'   `shape`, `params` and `acceptance`.
#' @noRd
combine_kept <- function(topo, ints, k, nsim, root_age, par) {
  trees <- lapply(seq_len(k), function(i) {
    tr <- assign_ages(topo$kept[[i]], ints$kept[[i]])
    tr$edge.length <- tr$edge.length * root_age
    tr
  })
  out <- if (nsim == 1) trees[[1]] else structure(trees, class = "multiPhylo")
  attr(out, "shape") <- data.frame(gamma = ints$stat[seq_len(k)],
                                   colless = topo$stat[seq_len(k)])
  attr(out, "params") <- c(beta = par$beta, c = par$c)
  attr(out, "acceptance") <- c(topology = topo$acceptance,
                               intervals = ints$acceptance)
  out
}

#' Draw objects until enough are within a tolerance of a target statistic
#'
#' @param draw Function with no arguments that draws one object.
#' @param stat Function that computes the statistic of an object.
#' @param target,tol Target value and tolerance.
#' @param nsim,max_tries Number of objects wanted and maximum number of draws.
#' @return A list with the `kept` objects, their `stat` and the `acceptance`
#'   rate.
#' @noRd
draw_kept <- function(draw, stat, target, tol, nsim, max_tries) {
  kept <- list()
  values <- numeric(0)
  tries <- 0
  while (length(kept) < nsim && tries < max_tries) {
    tries <- tries + 1
    x <- draw()
    s <- stat(x)
    if (abs(s - target) <= tol) {
      kept[[length(kept) + 1]] <- x
      values <- c(values, s)
    }
  }
  list(kept = kept, stat = values, acceptance = length(kept) / tries)
}

#' Check a calibration, or compute one
#'
#' @return `calibration`, or the output of [calibrate_shape()] if it is `NULL`.
#' @noRd
get_calibration <- function(calibration, n, reps) {
  check_count(reps, "reps")
  if (is.null(calibration)) return(calibrate_shape(n, reps = reps))
  if (!inherits(calibration, "shape_calibration") || calibration$n != n) {
    stop("`calibration` must come from calibrate_shape() with the same `n`.",
         call. = FALSE)
  }
  calibration
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
#' grid <- shape_grid(50, gamma = c(-2, -1, 0), colless = c(0, 0.5, 1),
#'                    reps = 50)
#' grid
#' @export
shape_grid <- function(n, gamma, colless, calibration = NULL, reps = 200) {
  check_n(n)
  calibration <- get_calibration(calibration, n, reps)
  grid <- expand.grid(gamma = gamma, colless = colless)
  out <- shape_params(calibration, gamma = grid$gamma, colless = grid$colless)
  attr(out, "calibration") <- calibration
  out
}
