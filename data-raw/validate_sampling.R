# Correction of c for incomplete sampling in fit_shape(): bias and coverage of
# the estimates before and after the correction. Run from the package root;
# results go to data-raw/validation/sampling_summary.csv.
devtools::load_all(quiet = TRUE)
RNGkind("L'Ecuyer-CMRG")
set.seed(1)
cond <- expand.grid(n = 100, c = c(0, -1), rho = c(0.75, 0.5, 0.25))
run_one <- function(p) {
  tr <- sim_tree_shape(p$n, beta = 0, c = p$c, sampling_fraction = p$rho)
  f <- fit_shape(tr, sampling_fraction = p$rho, correction_reps = 50)
  c(uncorrected = f$c_uncorrected, corrected = f$c, at_bound = f$c_at_bound)
}
reps <- 40
res <- do.call(rbind, lapply(seq_len(nrow(cond)), function(i) {
  p <- cond[i, ]
  r <- do.call(rbind, parallel::mclapply(seq_len(reps), function(j) run_one(p),
                                         mc.cores = 4, mc.set.seed = TRUE))
  data.frame(p, reps = reps,
             bias_uncorrected = mean(r[, "uncorrected"]) - p$c,
             bias_corrected = mean(r[, "corrected"]) - p$c,
             rmse_uncorrected = sqrt(mean((r[, "uncorrected"] - p$c)^2)),
             rmse_corrected = sqrt(mean((r[, "corrected"] - p$c)^2)),
             at_bound = sum(r[, "at_bound"]))
}))
write.csv(res, "data-raw/validation/sampling_summary.csv", row.names = FALSE)
print(res, digits = 2, row.names = FALSE)
