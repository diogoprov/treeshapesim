# treeshapesim (development version)

- `rtopology_beta()` builds the topology iteratively instead of recursively, so large and very imbalanced trees (e.g. 3,000 tips with `beta` near -2) no longer fail; with the same seed it returns exactly the same trees as before.
- `sim_tree_shape()`, `rtopology_beta()` and `intervals_lineage_rate()` check their arguments: `n` must be a single integer of at least 2, `beta` a single finite number greater than -2 and `c` a single finite number (non-finite `c` used to return trees with `NaN` branch lengths).
- Package checks (R CMD check) on GitHub Actions.

# treeshapesim 0.1.0

- First release: `sim_tree_shape()`, `rtopology_beta()`,
  `intervals_lineage_rate()` and `assign_ages()`, taken from the pcps-sim project.
