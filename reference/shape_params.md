# Shape parameters that give target mean values of gamma and Colless

Inverts the monotone curves of
[`calibrate_shape()`](https://diogoprov.github.io/treeshapesim/reference/calibrate_shape.md)
by linear interpolation. Where a curve is flat (tied values of its
monotone fit), the parameter is the mean of the tied values.

## Usage

``` r
shape_params(calibration, gamma, colless)
```

## Arguments

- calibration:

  Output of
  [`calibrate_shape()`](https://diogoprov.github.io/treeshapesim/reference/calibrate_shape.md).

- gamma, colless:

  Target mean gamma statistic and normalized Colless index (see
  [`tree_shape()`](https://diogoprov.github.io/treeshapesim/reference/tree_shape.md));
  vectors that are recycled to a common length, as in
  [`data.frame()`](https://rdrr.io/r/base/data.frame.html).

## Value

A data frame with the targets and the corresponding `beta` and `c`.

## Examples

``` r
set.seed(1)
cal <- calibrate_shape(50, reps = 50)
shape_params(cal, gamma = c(0, -2), colless = 0.5)
#>   gamma colless       beta           c
#> 1     0     0.5 -0.4020784  0.02818277
#> 2    -2     0.5 -0.4020784 -0.50284776
```
