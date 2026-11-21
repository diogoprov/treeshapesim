#' Argument checks shared by the exported functions
#'
#' `is_number()` is `TRUE` for a single finite number; `is_whole()` compares
#' with the nearest integer within a tolerance, not with `==`;
#' `is_finite_numeric()` is `TRUE` for a non-empty numeric vector of finite
#' values, optionally of a given length. `check_number()` stops with `msg`
#' unless `x` is a single finite number within `[lower, upper]` (each bound
#' excluded if the matching element of `open` is `TRUE`), and a whole number
#' if `whole = TRUE`.
#'
#' @srrstats {G2.0, G2.1, G2.2} Single-valued numeric arguments are checked for
#'   length, type and finiteness, so vectors are rejected where one value is
#'   expected.
#' @srrstats {G2.13, G2.14, G2.14a, G2.15, G2.16} Missing (`NA`) and undefined
#'   (`NaN`, `Inf`) values are rejected with an error in every numeric argument
#'   and in branch lengths; no function imputes or drops them, because a tree
#'   with a missing branch length has no defined shape.
#' @srrstats {G3.0} Integer checks use a tolerance rather than equality of
#'   floating-point numbers.
#' @noRd
is_number <- function(x) {
  is.numeric(x) && length(x) == 1 && is.finite(x)
}
is_whole <- function(x) {
  abs(x - round(x)) < sqrt(.Machine$double.eps)
}
is_finite_numeric <- function(x, len = length(x)) {
  is.numeric(x) && length(x) > 0 && length(x) == len && all(is.finite(x))
}
in_range <- function(x, lower, upper, open) {
  above <- if (open[1]) x > lower else x >= lower
  below <- if (open[2]) x < upper else x <= upper
  above && below
}
check_number <- function(x, msg, lower = -Inf, upper = Inf,
                         open = c(FALSE, FALSE), whole = FALSE) {
  ok <- is_number(x) && in_range(x, lower, upper, open) &&
    (!whole || is_whole(x))
  if (!ok) stop(msg, call. = FALSE)
}

#' Check that a tree is a binary, rooted phylo object
#'
#' @param lengths If `TRUE`, also require non-negative, finite branch lengths.
#' @param name Name of the argument, for the error message.
#' @noRd
check_tree <- function(tree, min_tips = 2, lengths = FALSE, name = "tree") {
  if (!is_tree(tree, min_tips)) {
    stop(sprintf(paste("`%s` must be a binary, rooted `phylo` object with at",
                       "least %d tips."), name, min_tips), call. = FALSE)
  }
  if (lengths && !valid_lengths(tree$edge.length)) {
    stop(sprintf("`%s` must have non-negative, finite branch lengths.", name),
         call. = FALSE)
  }
}
is_tree <- function(tree, min_tips) {
  inherits(tree, "phylo") && ape::Ntip(tree) >= min_tips &&
    ape::is.binary(tree) && ape::is.rooted(tree)
}
valid_lengths <- function(x) {
  is_finite_numeric(x) && all(x >= 0)
}

#' Checks of single arguments
#'
#' @noRd
check_n <- function(n) {
  check_number(n, "`n` must be a single integer of at least 2.",
               lower = 2, upper = .Machine$integer.max, whole = TRUE)
}
check_beta <- function(beta) {
  check_number(beta, "`beta` must be a single finite number greater than -2.",
               lower = -2, open = c(TRUE, FALSE))
}
check_c <- function(c) {
  check_number(c, "`c` must be a single finite number.")
}
check_fraction <- function(x) {
  check_number(x, paste("`sampling_fraction` must be a single number greater",
                        "than 0 and at most 1."),
               lower = 0, upper = 1, open = c(TRUE, FALSE))
}
check_count <- function(x, name) {
  check_number(x, sprintf("`%s` must be a single positive integer.", name),
               lower = 1, whole = TRUE)
}
check_positive <- function(x, name) {
  msg <- sprintf("`%s` must be a single positive, finite number.", name)
  check_number(x, msg, lower = 0, open = c(TRUE, FALSE))
}

#' Check a search interval
#'
#' @noRd
check_range <- function(x, name) {
  if (!is_finite_numeric(x, 2) || x[1] >= x[2]) {
    stop(sprintf("`%s` must hold two increasing finite numbers.", name),
         call. = FALSE)
  }
}
