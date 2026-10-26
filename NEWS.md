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
