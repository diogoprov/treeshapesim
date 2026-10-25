# Changelog

## treeshapesim (development version)

- [`rtopology_beta()`](https://diogoprov.github.io/treeshapesim/reference/rtopology_beta.md)
  builds the topology iteratively instead of recursively, so large and
  very imbalanced trees (e.g. 3,000 tips with `beta` near -2) no longer
  fail; with the same seed it returns exactly the same trees as before.
- [`sim_tree_shape()`](https://diogoprov.github.io/treeshapesim/reference/sim_tree_shape.md),
  [`rtopology_beta()`](https://diogoprov.github.io/treeshapesim/reference/rtopology_beta.md)
  and
  [`intervals_lineage_rate()`](https://diogoprov.github.io/treeshapesim/reference/intervals_lineage_rate.md)
  check their arguments: `n` must be a single integer of at least 2,
  `beta` a single finite number greater than -2 and `c` a single finite
  number (non-finite `c` used to return trees with `NaN` branch
  lengths).
- Package checks (R CMD check) on GitHub Actions.

## treeshapesim 0.1.0

- First release:
  [`sim_tree_shape()`](https://diogoprov.github.io/treeshapesim/reference/sim_tree_shape.md),
  [`rtopology_beta()`](https://diogoprov.github.io/treeshapesim/reference/rtopology_beta.md),
  [`intervals_lineage_rate()`](https://diogoprov.github.io/treeshapesim/reference/intervals_lineage_rate.md)
  and
  [`assign_ages()`](https://diogoprov.github.io/treeshapesim/reference/assign_ages.md),
  taken from the pcps-sim project.
