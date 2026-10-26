# Simulate trees near target values of gamma and Colless

Draws trees whose gamma statistic and normalized Colless index (see
[`tree_shape()`](https://diogoprov.github.io/treeshapesim/reference/tree_shape.md))
fall within a tolerance of target values. The parameters `beta` and `c`
whose mean statistics equal the targets are found with
[`calibrate_shape()`](https://diogoprov.github.io/treeshapesim/reference/calibrate_shape.md)
and
[`shape_params()`](https://diogoprov.github.io/treeshapesim/reference/shape_params.md).
Topologies are then drawn with
[`rtopology_beta()`](https://diogoprov.github.io/treeshapesim/reference/rtopology_beta.md)
and kept if their Colless index is within `tol` of the target, interval
sets are drawn with
[`intervals_lineage_rate()`](https://diogoprov.github.io/treeshapesim/reference/intervals_lineage_rate.md)
and kept if their gamma is within `tol` of the target, and the kept
topologies and interval sets are combined with
[`assign_ages()`](https://diogoprov.github.io/treeshapesim/reference/assign_ages.md).
Because the topology and the intervals are independent, rejecting them
separately gives the same distribution of trees as rejecting whole trees
on both statistics, with a much higher acceptance rate. Random numbers
are used in a different order than in
[`sim_tree_shape()`](https://diogoprov.github.io/treeshapesim/reference/sim_tree_shape.md),
so the same seed does not give the same trees as that function.

## Usage

``` r
sim_tree_target(
  n,
  gamma,
  colless,
  nsim = 1,
  tol = 0.5,
  calibration = NULL,
  reps = 200,
  max_tries = 1000 * nsim,
  root_age = 1
)
```

## Arguments

- n:

  Number of tips, a single integer of at least 2.

- gamma, colless:

  Target gamma statistic and normalized Colless index (single values).

- nsim:

  Number of trees.

- tol:

  Tolerance around the targets: one value for both statistics, or two
  values for gamma and Colless, in this order.

- calibration:

  Output of
  [`calibrate_shape()`](https://diogoprov.github.io/treeshapesim/reference/calibrate_shape.md)
  for the same `n`; if `NULL`, it is computed with `reps` draws per
  value.

- reps:

  Passed to
  [`calibrate_shape()`](https://diogoprov.github.io/treeshapesim/reference/calibrate_shape.md)
  when `calibration` is `NULL`.

- max_tries:

  Maximum number of topologies, and of interval sets, drawn.

- root_age:

  Age of the root.

## Value

With `nsim = 1`, a `phylo` object; otherwise a `multiPhylo` object (also
when fewer trees than `nsim` are found). Attributes: `shape` (realized
gamma and Colless of each tree), `params` (`beta` and `c` used) and
`acceptance` (proportion of topologies and of interval sets kept). If
fewer than `nsim` trees are found within `max_tries`, the trees found
are returned with a warning.

## Examples

``` r
set.seed(1)
cal <- calibrate_shape(50, reps = 50)
trees <- sim_tree_target(50, gamma = -1, colless = 0.5, nsim = 5,
                         calibration = cal)
attr(trees, "shape")
#>        gamma    colless
#> 1 -0.5955608 0.08390851
#> 2 -1.0764385 0.38390851
#> 3 -0.9091632 0.02390851
#> 4 -0.9009388 0.12390851
#> 5 -1.2630195 0.44390851
```
