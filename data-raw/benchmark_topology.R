# Benchmark of the topology generator (not part of the package build).
# Compares rtopology_beta() in R with the C++ prototype in
# data-raw/prototype/rtopology_beta.cpp (which reproduces the R version exactly:
# same trees and RNG state for the same seed), and measures the share of
# rtopology_beta() and assign_ages() in sim_tree_shape(). Run on an idle
# machine: timings under load are not comparable.
# Writes data-raw/benchmark/topology.csv and session.txt.
# Usage (from the package root): Rscript data-raw/benchmark_topology.R [reps]

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[[1]]) else 3L
pkgload::load_all(".", quiet = TRUE)
Rcpp::sourceCpp("data-raw/prototype/rtopology_beta.cpp")
out_dir <- file.path("data-raw", "benchmark")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rtopology_cpp <- function(n, beta) {
  tree <- list(edge = rtopology_edges_cpp(as.integer(n), beta), Nnode = as.integer(n) - 1L,
               tip.label = paste0("t", seq_len(n)))
  class(tree) <- "phylo"
  attr(tree, "order") <- "cladewise"
  tree
}
elapsed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  force(expr)
  proc.time()[["elapsed"]] - t0
}

grid <- expand.grid(n = c(100, 500, 1000, 2000, 5000, 10000, 20000),
                    beta = c(-1.9, -1, 0, 2))
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(row) {
  n <- grid$n[row]
  b <- grid$beta[row]
  times <- t(vapply(seq_len(reps), function(r) {
    set.seed(r)
    t_r <- elapsed(tp_r <- rtopology_beta(n, b))
    set.seed(r)
    t_cpp <- elapsed(tp_cpp <- rtopology_cpp(n, b))
    stopifnot(identical(tp_r, tp_cpp))
    g <- intervals_lineage_rate(n, 0)
    t_ages <- elapsed(assign_ages(tp_r, g))
    c(r = t_r, cpp = t_cpp, ages = t_ages)
  }, numeric(3)))
  message(sprintf("n = %5d, beta = %4.1f: R %.3f s, C++ %.3f s, ages %.3f s", n, b,
                  stats::median(times[, "r"]), stats::median(times[, "cpp"]),
                  stats::median(times[, "ages"])))
  data.frame(n = n, beta = b, topology_r = stats::median(times[, "r"]),
             topology_cpp = stats::median(times[, "cpp"]),
             assign_ages = stats::median(times[, "ages"]))
}))
res$speedup <- res$topology_r / res$topology_cpp
utils::write.csv(res, file.path(out_dir, "topology.csv"), row.names = FALSE)
writeLines(c(paste("load average at start:", paste(Sys.getenv("TSS_LOAD"), collapse = "")),
             utils::capture.output(utils::sessionInfo()),
             utils::capture.output(print(Sys.info()[c("sysname", "release", "machine")]))),
           file.path(out_dir, "session.txt"))
print(res, row.names = FALSE, digits = 3)
