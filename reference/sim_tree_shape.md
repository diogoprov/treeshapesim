# Ultrametric tree with given imbalance and branching-time parameters

Draws the topology with
[`rtopology_beta()`](https://diogoprov.github.io/treeshapesim/reference/rtopology_beta.md)
and the internode intervals with
[`intervals_lineage_rate()`](https://diogoprov.github.io/treeshapesim/reference/intervals_lineage_rate.md),
and combines them with
[`assign_ages()`](https://diogoprov.github.io/treeshapesim/reference/assign_ages.md).
Because the gamma statistic depends only on the intervals and the
Colless index only on the topology, `beta` mainly controls imbalance and
`c` mainly controls gamma.

## Usage

``` r
sim_tree_shape(n, beta = 0, c = 0)
```

## Arguments

- n:

  Number of tips, a single integer of at least 2.

- beta:

  Aldous' beta, a single finite number greater than -2. Lower values
  give more imbalanced trees.

- c:

  Lineage dependence of the split rate, a single finite number.

## Value

A binary, ultrametric `phylo` object with root depth 1.

## Examples

``` r
set.seed(1)
tr <- sim_tree_shape(50, beta = -1, c = -1)
ape::gammaStat(tr)
#> [1] -4.320064
```
