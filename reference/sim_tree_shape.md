# Ultrametric trees with given imbalance and branching-time parameters

Draws the topology with
[`rtopology_beta()`](https://provetelab.org/treeshapesim/reference/rtopology_beta.md)
and the internode intervals with
[`intervals_lineage_rate()`](https://provetelab.org/treeshapesim/reference/intervals_lineage_rate.md),
and combines them with
[`assign_ages()`](https://provetelab.org/treeshapesim/reference/assign_ages.md).
The topology and the intervals are drawn independently, the gamma
statistic depends only on the intervals and the Colless index only on
the topology, so for a given `n` the realized gamma depends only on `c`
(or `rate`) and the realized Colless index only on `beta`. To draw trees
near target values of both statistics, see
[`sim_tree_target()`](https://provetelab.org/treeshapesim/reference/sim_tree_target.md).

## Usage

``` r
sim_tree_shape(
  n,
  beta = 0,
  c = 0,
  nsim = 1,
  root_age = 1,
  rate = NULL,
  sampling_fraction = 1
)
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

- rate:

  `NULL` (the default), to use `c`, or a vectorized function of the
  number of lineages k that returns the per-lineage split rate. It must
  return n - 1 positive, finite values for k = 2, ..., n, and cannot be
  combined with a `c` other than 0.

- sampling_fraction:

  Fraction of the tips of the complete tree that are kept, a single
  number greater than 0 and at most 1.

## Value

With `nsim = 1`, a binary, ultrametric `phylo` object with root depth
`root_age`; otherwise a `multiPhylo` object with `nsim` such trees.
Attribute `params` holds `beta`, `c` (`NA` with `rate`), `root_age` and
`sampling_fraction`, and attribute `rate` holds the `rate` function, if
any.

## Details

With `sampling_fraction` below 1, a complete tree of
`round(n / sampling_fraction)` tips is drawn, `beta` and `c` (or `rate`)
apply to that tree, and `n` tips are kept at random, so the result still
has `n` tips and root depth `root_age`. Incomplete sampling removes
recent splits and so lowers gamma (Pybus & Harvey 2000), so gamma no
longer follows the calibration of
[`calibrate_shape()`](https://provetelab.org/treeshapesim/reference/calibrate_shape.md).
After sampling, gamma also depends slightly on the topology, because the
splits that remain depend on it, so the separation between gamma and the
Colless index is exact only for complete trees. A `rate` function must
be valid for the number of tips of the complete tree.

Very extreme values of `c` (or `rate`) make some intervals negligible
relative to the others, and the corresponding branches round to length
zero; a warning is given when this happens.

## Reproducibility

The same seed gives the same trees on the same platform. Between
platforms (for instance, arm64 and x86-64), differences in the last bit
of floating-point functions can change a few draws of the topology, so
trees may differ.

## References

Pybus, O. G., & Harvey, P. H. (2000). Testing macro-evolutionary models
using incomplete molecular phylogenies. *Proceedings of the Royal
Society B*, 267, 2267–2272.
[doi:10.1098/rspb.2000.1278](https://doi.org/10.1098/rspb.2000.1278)

## Examples

``` r
set.seed(1)
tr <- sim_tree_shape(50, beta = -1, c = -1)
ape::gammaStat(tr)
#> [1] -4.320064
trees <- sim_tree_shape(50, beta = -1, c = -1, nsim = 10, root_age = 30)
# logistic slowdown towards a carrying capacity of 60 lineages
sim_tree_shape(50, rate = function(k) 1 - k / 60)
#> 
#> Phylogenetic tree with 50 tips and 49 internal nodes.
#> 
#> Tip labels:
#>   t1, t2, t3, t4, t5, t6, ...
#> 
#> Rooted; includes branch length(s).
# 50 tips sampled from a complete tree of 100
sim_tree_shape(50, sampling_fraction = 0.5)
#> 
#> Phylogenetic tree with 50 tips and 49 internal nodes.
#> 
#> Tip labels:
#>   t1, t2, t3, t4, t5, t6, ...
#> 
#> Rooted; includes branch length(s).
```
