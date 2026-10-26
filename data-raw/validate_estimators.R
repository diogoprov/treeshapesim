# Validation of two estimators on trees simulated with treeshapesim (not part of
# the package build; data-raw/ is in .Rbuildignore).
#   1. beta: apTreeshape::maxlik.betasplit() (maximum likelihood of Aldous'
#      beta-splitting model; Blum & Francois 2006) on topologies from
#      rtopology_beta(). An independent implementation recovering beta checks
#      that the generator follows the beta-splitting model.
#   2. c: maximum likelihood of the DDX model of laser (Rabosky 2006), in which
#      the per-lineage speciation rate with k lineages is r * k^(-x), so that
#      c = -x. The likelihood below follows laser::DDX(): splits at 2, ..., n - 1
#      lineages and no event in the last interval (n lineages to the present).
#      intervals_lineage_rate() also draws that last interval from the same law,
#      so a small bias is possible and is what this run measures.
# Writes data-raw/validation/beta.csv, c.csv and summary.csv.
# Usage (from the package root): Rscript data-raw/validate_estimators.R [reps]

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[[1]]) else 200L
lib <- Sys.getenv("TSS_VALIDATION_LIB", unset = "")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
if (!requireNamespace("apTreeshape", quietly = TRUE)) {
  stop("apTreeshape is needed (install.packages('apTreeshape')).")
}
pkgload::load_all(".", quiet = TRUE)
out_dir <- file.path("data-raw", "validation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(20261026)

# DDX: negative log-likelihood of the intervals g_2, ..., g_n (g[k - 1] is the
# time with k lineages) for log(r) and x
ddx_nll <- function(p, intervals) {
  g <- intervals
  r <- exp(p[1])
  x <- p[2]
  n <- length(g) + 1
  k <- 2:n
  split <- 2:(n - 1)
  -(sum(log(split * r * split^(-x))) - sum(r * k^(1 - x) * g))
}
fit_ddx <- function(g) {
  start <- c(log((length(g) - 1) / sum((2:(length(g) + 1)) * g)), 0)
  # (`intervals =`, not `g =`, which would partially match optim()'s `gr`)
  fit <- stats::optim(start, ddx_nll, intervals = g, method = "BFGS")
  c(c_hat = -fit$par[2], converged = fit$convergence == 0)
}

grid_n <- c(50, 100, 200)
grid_beta <- c(-1.8, -1.5, -1, -0.5, 0, 1, 3)
grid_c <- c(-1.5, -1, -0.5, 0, 0.5)

t0 <- Sys.time()
beta_res <- do.call(rbind, lapply(grid_n, function(n) do.call(rbind, lapply(grid_beta, function(b) {
  est <- vapply(seq_len(reps), function(i) {
    apTreeshape::maxlik.betasplit(rtopology_beta(n, b), up = 10)$max_lik
  }, 0)
  data.frame(n = n, beta = b, beta_hat = est)
}))))
utils::write.csv(beta_res, file.path(out_dir, "beta.csv"), row.names = FALSE)
message("beta done: ", format(Sys.time() - t0))

c_res <- do.call(rbind, lapply(grid_n, function(n) do.call(rbind, lapply(grid_c, function(cc) {
  est <- t(vapply(seq_len(reps), function(i) fit_ddx(intervals_lineage_rate(n, cc)), numeric(2)))
  data.frame(n = n, c = cc, c_hat = est[, "c_hat"], converged = est[, "converged"] == 1)
}))))
utils::write.csv(c_res, file.path(out_dir, "c.csv"), row.names = FALSE)
message("c done: ", format(Sys.time() - t0))

summarise <- function(d, true, est) {
  key <- interaction(d$n, d[[true]], drop = TRUE)
  do.call(rbind, lapply(split(d, key), function(x) data.frame(
    parameter = true, n = x$n[1], true = x[[true]][1],
    mean = mean(x[[est]]), bias = mean(x[[est]]) - x[[true]][1],
    rmse = sqrt(mean((x[[est]] - x[[true]][1])^2)),
    q025 = stats::quantile(x[[est]], 0.025), q975 = stats::quantile(x[[est]], 0.975))))
}
summary <- rbind(summarise(beta_res, "beta", "beta_hat"),
                 summarise(c_res[c_res$converged, ], "c", "c_hat"))
utils::write.csv(summary, file.path(out_dir, "summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE, digits = 3)
message("total: ", format(Sys.time() - t0))
