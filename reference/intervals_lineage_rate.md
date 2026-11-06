# Internode intervals of a pure-birth process with lineage-dependent rates

With k lineages each splits at rate k^c, so the interval with k lineages
is Exp(k^(1 + c)). `c = 0` is Yule; `c < 0` slows diversification as
lineages accumulate (gamma \< 0); `c > 0` accelerates it (gamma \> 0).
The final interval (n lineages, up to the present) is drawn from the
same law.

## Usage

``` r
intervals_lineage_rate(n, c = 0, rate = NULL)
```

## Arguments

- n:

  Number of tips, a single integer of at least 2.

- c:

  Lineage dependence of the split rate, a single finite number.

- rate:

  `NULL` (the default), to use `c`, or a vectorized function of the
  number of lineages k that returns the per-lineage split rate. It must
  return n - 1 positive, finite values for k = 2, ..., n, and cannot be
  combined with a `c` other than 0.

## Value

Intervals g_2, ..., g_n (length n - 1), scaled to sum 1.

## Details

Any other dependence on the number of lineages can be given with `rate`,
a function of k that returns the per-lineage split rate; the interval
with k lineages is then Exp(k \* rate(k)), and `c` is the special case
`rate = function(k) k^c`. Only the shape of `rate` matters: the
intervals are rescaled, so multiplying it by a constant does not change
the tree. As with `c`, the shape of `rate` describes the
lineage-through-time curve and does not identify the mechanism behind it
(Pannetier et al. 2021).
[`calibrate_shape()`](https://provetelab.org/treeshapesim/reference/calibrate_shape.md),
[`shape_params()`](https://provetelab.org/treeshapesim/reference/shape_params.md)
and
[`sim_tree_target()`](https://provetelab.org/treeshapesim/reference/sim_tree_target.md)
work with `c` only; with `rate`, the calibration is left to the user.

## References

Pannetier, T., Martinez, C., Bunnefeld, L., & Etienne, R. S. (2021).
Branching patterns in phylogenies cannot distinguish diversity-dependent
diversification from time-dependent diversification. *Evolution*, 75,
25–38. [doi:10.1111/evo.14124](https://doi.org/10.1111/evo.14124)

Pybus, O. G., & Harvey, P. H. (2000). Testing macro-evolutionary models
using incomplete molecular phylogenies. *Proceedings of the Royal
Society B*, 267, 2267–2272.
[doi:10.1098/rspb.2000.1278](https://doi.org/10.1098/rspb.2000.1278)

## Examples

``` r
set.seed(1)
intervals_lineage_rate(10, c = -1)
#> [1] 0.09121461 0.14272468 0.01759918 0.01688517 0.05267053 0.34966867 0.14851261
#> [8] 0.06518557 0.11553897
# logistic slowdown towards a carrying capacity K > n
intervals_lineage_rate(10, rate = function(k) 1 - k / 12)
#> [1] 0.01812595 0.12698669 0.05870831 0.08717482 0.30295889 0.07428034 0.07975727
#> [8] 0.17129897 0.08070879
# the rate drops once there are more than 5 lineages
intervals_lineage_rate(10, rate = function(k) ifelse(k <= 5, 1, 0.1))
#> [1] 0.05331123 0.06207473 0.18706249 0.04062533 0.15512381 0.25581133 0.04195830
#> [8] 0.02089944 0.18313334
```
