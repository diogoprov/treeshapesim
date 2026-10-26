# Grid of tree shapes

Crosses levels of the gamma statistic with levels of the normalized
Colless index and finds the `beta` and `c` that give each combination on
average, as a design for simulation studies. The levels can be chosen,
for instance, from quantiles of
[`tree_shape()`](https://diogoprov.github.io/treeshapesim/reference/tree_shape.md)
over empirical trees of similar size.

## Usage

``` r
shape_grid(n, gamma, colless, calibration = NULL, reps = 200)
```

## Arguments

- n:

  Number of tips, a single integer of at least 2.

- gamma, colless:

  Levels of each statistic.

- calibration:

  Output of
  [`calibrate_shape()`](https://diogoprov.github.io/treeshapesim/reference/calibrate_shape.md)
  for the same `n`; if `NULL`, it is computed with `reps` draws per
  value.

- reps:

  Passed to
  [`calibrate_shape()`](https://diogoprov.github.io/treeshapesim/reference/calibrate_shape.md)
  when `calibration` is `NULL`.

## Value

A data frame with one row per combination of levels (`gamma`, `colless`,
`beta`, `c`), with the calibration in attribute `calibration`. Each row
can be passed to
[`sim_tree_target()`](https://diogoprov.github.io/treeshapesim/reference/sim_tree_target.md)
(targets) or to
[`sim_tree_shape()`](https://diogoprov.github.io/treeshapesim/reference/sim_tree_shape.md)
(parameters).

## Examples

``` r
set.seed(1)
grid <- shape_grid(50, gamma = c(-2, -1, 0), colless = c(0, 0.5, 1), reps = 50)
grid
#>   gamma colless        beta           c
#> 1    -2     0.0  0.07990851 -0.50284776
#> 2    -1     0.0  0.07990851 -0.22999948
#> 3     0     0.0  0.07990851  0.02818277
#> 4    -2     0.5 -0.40207838 -0.50284776
#> 5    -1     0.5 -0.40207838 -0.22999948
#> 6     0     0.5 -0.40207838  0.02818277
#> 7    -2     1.0 -0.67972320 -0.50284776
#> 8    -1     1.0 -0.67972320 -0.22999948
#> 9     0     1.0 -0.67972320  0.02818277
```
