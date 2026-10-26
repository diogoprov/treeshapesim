# Ultrametric trees with given imbalance and branching-time parameters

Draws the topology with
[`rtopology_beta()`](https://diogoprov.github.io/treeshapesim/reference/rtopology_beta.md)
and the internode intervals with
[`intervals_lineage_rate()`](https://diogoprov.github.io/treeshapesim/reference/intervals_lineage_rate.md),
and combines them with
[`assign_ages()`](https://diogoprov.github.io/treeshapesim/reference/assign_ages.md).
The topology and the intervals are drawn independently, the gamma
statistic depends only on the intervals and the Colless index only on
the topology, so for a given `n` the realized gamma depends only on `c`
and the realized Colless index only on `beta`. To draw trees near target
values of both statistics, see
[`sim_tree_target()`](https://diogoprov.github.io/treeshapesim/reference/sim_tree_target.md).

## Usage

``` r
sim_tree_shape(n, beta = 0, c = 0, nsim = 1, root_age = 1)
```

## Arguments

- n:

  Number of tips, a single integer of at least 2.

- beta:

  Aldous' beta, a single finite number greater than -2. Lower values
  give more imbalanced trees.

- c:

  Lineage dependence of the split rate, a single finite number.

- nsim:

  Number of trees.

- root_age:

  Age of the root (positive).

## Value

With `nsim = 1`, a binary, ultrametric `phylo` object with root depth
`root_age`; otherwise a `multiPhylo` object with `nsim` such trees.

## Examples

``` r
set.seed(1)
tr <- sim_tree_shape(50, beta = -1, c = -1)
ape::gammaStat(tr)
#> [1] -4.320064
trees <- sim_tree_shape(50, beta = -1, c = -1, nsim = 10, root_age = 30)
```
