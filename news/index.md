# Changelog

## treeshapesim (development version)

- [`sim_tree_shape()`](https://provetelab.org/treeshapesim/reference/sim_tree_shape.md)
  and
  [`intervals_lineage_rate()`](https://provetelab.org/treeshapesim/reference/intervals_lineage_rate.md)
  gain `rate`, a function of the number of lineages that returns the
  per-lineage split rate, for any dependence other than k^c (for
  instance, logistic or stepwise). The intervals remain independent of
  the topology, and with the default `rate = NULL` the trees are
  unchanged.
  [`calibrate_shape()`](https://provetelab.org/treeshapesim/reference/calibrate_shape.md),
  [`shape_params()`](https://provetelab.org/treeshapesim/reference/shape_params.md)
  and
  [`sim_tree_target()`](https://provetelab.org/treeshapesim/reference/sim_tree_target.md)
  still work with `c` only.
- [`sim_tree_shape()`](https://provetelab.org/treeshapesim/reference/sim_tree_shape.md)
  gains `sampling_fraction`: it draws a complete tree of
  `round(n / sampling_fraction)` tips and keeps `n` at random, with root
  depth `root_age`. The default keeps the same trees.
- [`sim_tree_shape()`](https://provetelab.org/treeshapesim/reference/sim_tree_shape.md)
  stores `beta`, `c`, `root_age` and `sampling_fraction` in attribute
  `params`, and the `rate` function, if any, in attribute `rate`.
- [`sim_tree_shape()`](https://provetelab.org/treeshapesim/reference/sim_tree_shape.md)
  warns when branches round to length zero, which happens with very
  negative `c` (for instance, `c = -20` with 100 tips); `c` so extreme
  that the split rates overflow or underflow now gives a clear error
  instead of `NA` draws.
- [`fit_shape()`](https://provetelab.org/treeshapesim/reference/fit_shape.md)
  gains `sampling_fraction` and `correction_reps`: with incomplete
  sampling, `c` is corrected by simulation, and the bootstrap uses the
  same correction. In a validation with 100 tips, the bias fell from
  -0.08 to -0.78 to within 0.03 for sampling fractions of 0.75 to 0.25
  (`data-raw/validate_sampling.R`).
- [`extinct_tips()`](https://provetelab.org/treeshapesim/reference/extinct_tips.md)
  applies one extinction event to a tree: tips are lost at random, in
  clusters on the tree (a latent liability with Pagel’s lambda,
  `signal`) or according to a trait (`risk`). The tip labels are kept,
  so that traits stay aligned with the remaining species. It was first
  added as `prune_tips()` and renamed before release, because treesliceR
  exports a function with that name.
- A new article shows trees from treeshapesim, with non-random
  extinction, as regional pools in comdynamite.
- The package follows the rOpenSci packaging guide and the general and
  probability-distribution statistical software standards, documented
  with srr tags (`srr::srr_report()`): a contributing guide, a Code of
  Conduct, a codemeta.json file, ORCID, status and lifecycle badges, and
  a package help page with definitions of terms, the novelty of the
  algorithms and a life cycle statement.
- [`fit_shape()`](https://provetelab.org/treeshapesim/reference/fit_shape.md)
  returns `optimizer`, with the name of the algorithm, its tolerance and
  the number of evaluations of each log-likelihood.
- Branch lengths that are missing, infinite or negative are now an error
  in
  [`fit_shape()`](https://provetelab.org/treeshapesim/reference/fit_shape.md)
  and
  [`extinct_tips()`](https://provetelab.org/treeshapesim/reference/extinct_tips.md),
  and whole-number arguments are checked with a tolerance.
- New tests of split probabilities and likelihoods against explicit
  formulas, of alternative search intervals and optimizers, of noise,
  edge conditions and parameter recovery; extended tests run with
  `TREESHAPESIM_EXTENDED_TESTS="true"`.
- Argument checks are gathered in internal helpers, so that no function
  has cyclomatic complexity above 10; error messages and results are
  unchanged.
- Test coverage on Codecov.
- The README compares treeshapesim with TreeSimGM too, and lists
  [`fit_shape()`](https://provetelab.org/treeshapesim/reference/fit_shape.md)
  as estimation of tree shape.
- The README explains the use of trees from treeshapesim as regional
  pools in comdynamite.

## treeshapesim 0.2.1

- A test of
  [`fit_shape()`](https://provetelab.org/treeshapesim/reference/fit_shape.md)
  used a simulated tree, whose topology differs between arm64 and x86-64
  for the same seed; it now reads a fixed tree.
  [`?sim_tree_shape`](https://provetelab.org/treeshapesim/reference/sim_tree_shape.md)
  documents that the same seed gives the same trees on the same
  platform.

## treeshapesim 0.2.0

- [`fit_shape()`](https://provetelab.org/treeshapesim/reference/fit_shape.md)
  estimates `beta` (maximum likelihood of the beta-splitting model from
  the topology) and `c` (maximum likelihood of the DDX model of laser,
  with `c = -x`) from an empirical tree, with an optional parametric
  bootstrap; on 30 random trees, `beta` matched the archived apTreeshape
  (`maxlik.betasplit()`) and `c` matched laser (`DDX()`) to within 3e-5.
- [`sim_tree_target()`](https://provetelab.org/treeshapesim/reference/sim_tree_target.md)
  draws trees whose gamma statistic and normalized Colless index fall
  within a tolerance of target values: `beta` and `c` come from
  calibration curves
  ([`calibrate_shape()`](https://provetelab.org/treeshapesim/reference/calibrate_shape.md),
  inverted by
  [`shape_params()`](https://provetelab.org/treeshapesim/reference/shape_params.md)),
  and topologies and internode intervals are rejected separately, which
  is equivalent to rejecting whole trees but much more efficient.
- [`shape_grid()`](https://provetelab.org/treeshapesim/reference/shape_grid.md)
  crosses levels of gamma and Colless and gives the parameters of each
  combination.
- [`tree_shape()`](https://provetelab.org/treeshapesim/reference/tree_shape.md)
  returns the gamma statistic and the Colless index normalized against
  its Yule expectation, in pure R (same values as
  [`ape::gammaStat()`](https://rdrr.io/pkg/ape/man/gammaStat.html) and
  `treestats::colless(normalization = "yule")`).
- [`sim_tree_shape()`](https://provetelab.org/treeshapesim/reference/sim_tree_shape.md)
  gains `nsim` (returns a `multiPhylo`) and `root_age`; with the
  defaults it returns the same trees as before.
- [`assign_ages()`](https://provetelab.org/treeshapesim/reference/assign_ages.md)
  finds subtree sizes and children in one pass instead of with
  [`ape::prop.part()`](https://rdrr.io/pkg/ape/man/boot.phylo.html), and
  checks its arguments; trees are unchanged.
- A vignette and a CITATION file.
- [`rtopology_beta()`](https://provetelab.org/treeshapesim/reference/rtopology_beta.md)
  builds the topology iteratively instead of recursively, so large and
  very imbalanced trees (e.g. 3,000 tips with `beta` near -2) no longer
  fail; with the same seed it returns exactly the same trees as before.
- [`sim_tree_shape()`](https://provetelab.org/treeshapesim/reference/sim_tree_shape.md),
  [`rtopology_beta()`](https://provetelab.org/treeshapesim/reference/rtopology_beta.md)
  and
  [`intervals_lineage_rate()`](https://provetelab.org/treeshapesim/reference/intervals_lineage_rate.md)
  check their arguments: `n` must be a single integer of at least 2,
  `beta` a single finite number greater than -2 and `c` a single finite
  number (non-finite `c` used to return trees with `NaN` branch
  lengths).
- Package checks (R CMD check) on GitHub Actions.

## treeshapesim 0.1.0

- First release:
  [`sim_tree_shape()`](https://provetelab.org/treeshapesim/reference/sim_tree_shape.md),
  [`rtopology_beta()`](https://provetelab.org/treeshapesim/reference/rtopology_beta.md),
  [`intervals_lineage_rate()`](https://provetelab.org/treeshapesim/reference/intervals_lineage_rate.md)
  and
  [`assign_ages()`](https://provetelab.org/treeshapesim/reference/assign_ages.md),
  taken from the pcps-sim project.
