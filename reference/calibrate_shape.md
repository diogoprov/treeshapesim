# Calibrate the shape parameters for a number of tips

Estimates, for a given number of tips, the mean gamma statistic as a
function of `c` and the mean normalized Colless index (see
[`tree_shape()`](https://diogoprov.github.io/treeshapesim/reference/tree_shape.md))
as a function of `beta`. Because the topology and the internode
intervals are drawn independently, each curve is estimated from its own
process alone: gamma from intervals drawn with
[`intervals_lineage_rate()`](https://diogoprov.github.io/treeshapesim/reference/intervals_lineage_rate.md)
and the Colless index from topologies drawn with
[`rtopology_beta()`](https://diogoprov.github.io/treeshapesim/reference/rtopology_beta.md),
without building trees. The expected values are monotone (gamma
increases with `c` and the Colless index decreases with `beta`), but the
estimated means need not be, especially for few tips or few draws, so
each curve is made monotone by isotonic regression
([`stats::isoreg()`](https://rdrr.io/r/stats/isoreg.html)) before it is
inverted with
[`shape_params()`](https://diogoprov.github.io/treeshapesim/reference/shape_params.md).

## Usage

``` r
calibrate_shape(
  n,
  beta = c(-1.95, -1.8, -1.6, -1.4, -1.2, -1, -0.75, -0.5, 0, 0.5, 1, 2, 4, 8),
  c = seq(-2, 1, by = 0.25),
  reps = 200
)
```

## Arguments

- n:

  Number of tips, a single integer of at least 4 (with 3 tips there is a
  single topology, so the Colless index cannot be calibrated).

- beta:

  Values of `beta` at which the Colless index is estimated (at least
  two).

- c:

  Values of `c` at which gamma is estimated (at least two).

- reps:

  Number of topologies or interval sets per value.

## Value

An object of class `shape_calibration`: a list with `n`, `reps` and the
data frames `colless` (columns `beta`, `colless` and `colless_mean`) and
`gamma` (columns `c`, `gamma` and `gamma_mean`), where the `_mean`
columns hold the mean statistic at each value and `colless` and `gamma`
their monotone fit, which
[`shape_params()`](https://diogoprov.github.io/treeshapesim/reference/shape_params.md)
inverts.

## Examples

``` r
set.seed(1)
cal <- calibrate_shape(50, reps = 50)
cal$gamma
#>        c      gamma gamma_mean
#> 1  -2.00 -5.8035976 -5.8035976
#> 2  -1.75 -5.5579902 -5.5579902
#> 3  -1.50 -4.6638505 -4.6638505
#> 4  -1.25 -4.4493981 -4.4493981
#> 5  -1.00 -3.7145604 -3.7145604
#> 6  -0.75 -2.8344574 -2.8344574
#> 7  -0.50 -1.9903851 -1.9903851
#> 8  -0.25 -1.0701663 -1.0701663
#> 9   0.00 -0.1931103 -0.1931103
#> 10  0.25  1.5199069  1.5199069
#> 11  0.50  2.9293050  2.9293050
#> 12  0.75  4.0540840  4.0540840
#> 13  1.00  5.7964429  5.7964429
```
