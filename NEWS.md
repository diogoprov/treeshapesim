# treeshapesim (development version)

- `sim_tree_shape()` and `intervals_lineage_rate()` gain `rate`, a function of the number of lineages that returns the per-lineage split rate, for any dependence other than k^c (for instance, logistic or stepwise). The intervals remain independent of the topology, and with the default `rate = NULL` the trees are unchanged. `calibrate_shape()`, `shape_params()` and `sim_tree_target()` still work with `c` only.
- `sim_tree_shape()` gains `sampling_fraction`: it draws a complete tree of `round(n / sampling_fraction)` tips and keeps `n` at random, with root depth `root_age`. The default keeps the same trees.
- `sim_tree_shape()` stores `beta`, `c`, `root_age` and `sampling_fraction` in attribute `params`, and the `rate` function, if any, in attribute `rate`.
- `sim_tree_shape()` warns when branches round to length zero, which happens with very negative `c` (for instance, `c = -20` with 100 tips); `c` so extreme that the split rates overflow or underflow now gives a clear error instead of `NA` draws.
- `fit_shape()` gains `sampling_fraction` and `correction_reps`: with incomplete sampling, `c` is corrected by simulation, and the bootstrap uses the same correction. In a validation with 100 tips, the bias fell from -0.08 to -0.78 to within 0.03 for sampling fractions of 0.75 to 0.25 (`data-raw/validate_sampling.R`).
- `extinct_tips()` applies one extinction event to a tree: tips are lost at random, in clusters on the tree (a latent liability with Pagel's lambda, `signal`) or according to a trait (`risk`). The tip labels are kept, so that traits stay aligned with the remaining species. It was first added as `prune_tips()` and renamed before release, because treesliceR exports a function with that name.
- A new article shows trees from treeshapesim, with non-random extinction, as regional pools in comdynamite.
- The package follows the rOpenSci packaging guide and the general and probability-distribution statistical software standards, documented with srr tags (`srr::srr_report()`): a contributing guide, a Code of Conduct, a codemeta.json file, ORCID, status and lifecycle badges, and a package help page with definitions of terms, the novelty of the algorithms and a life cycle statement.
- `fit_shape()` returns `optimizer`, with the name of the algorithm, its tolerance and the number of evaluations of each log-likelihood.
- Branch lengths that are missing, infinite or negative are now an error in `fit_shape()` and `extinct_tips()`, and whole-number arguments are checked with a tolerance.
- New tests of split probabilities and likelihoods against explicit formulas, of alternative search intervals and optimizers, of noise, edge conditions and parameter recovery; extended tests run with `TREESHAPESIM_EXTENDED_TESTS="true"`.
- Test coverage on Codecov.
- The README compares treeshapesim with TreeSimGM too, and lists `fit_shape()` as estimation of tree shape.
- The README explains the use of trees from treeshapesim as regional pools in comdynamite.

# treeshapesim 0.2.1

- A test of `fit_shape()` used a simulated tree, whose topology differs between arm64 and x86-64 for the same seed; it now reads a fixed tree. `?sim_tree_shape` documents that the same seed gives the same trees on the same platform.

# treeshapesim 0.2.0

- `fit_shape()` estimates `beta` (maximum likelihood of the beta-splitting model from the topology) and `c` (maximum likelihood of the DDX model of laser, with `c = -x`) from an empirical tree, with an optional parametric bootstrap; on 30 random trees, `beta` matched the archived apTreeshape (`maxlik.betasplit()`) and `c` matched laser (`DDX()`) to within 3e-5.
- `sim_tree_target()` draws trees whose gamma statistic and normalized Colless index fall within a tolerance of target values: `beta` and `c` come from calibration curves (`calibrate_shape()`, inverted by `shape_params()`), and topologies and internode intervals are rejected separately, which is equivalent to rejecting whole trees but much more efficient.
- `shape_grid()` crosses levels of gamma and Colless and gives the parameters of each combination.
- `tree_shape()` returns the gamma statistic and the Colless index normalized against its Yule expectation, in pure R (same values as `ape::gammaStat()` and `treestats::colless(normalization = "yule")`).
- `sim_tree_shape()` gains `nsim` (returns a `multiPhylo`) and `root_age`; with the defaults it returns the same trees as before.
- `assign_ages()` finds subtree sizes and children in one pass instead of with `ape::prop.part()`, and checks its arguments; trees are unchanged.
- A vignette and a CITATION file.
- `rtopology_beta()` builds the topology iteratively instead of recursively, so large and very imbalanced trees (e.g. 3,000 tips with `beta` near -2) no longer fail; with the same seed it returns exactly the same trees as before.
- `sim_tree_shape()`, `rtopology_beta()` and `intervals_lineage_rate()` check their arguments: `n` must be a single integer of at least 2, `beta` a single finite number greater than -2 and `c` a single finite number (non-finite `c` used to return trees with `NaN` branch lengths).
- Package checks (R CMD check) on GitHub Actions.

# treeshapesim 0.1.0

- First release: `sim_tree_shape()`, `rtopology_beta()`,
  `intervals_lineage_rate()` and `assign_ages()`, taken from the pcps-sim project.
