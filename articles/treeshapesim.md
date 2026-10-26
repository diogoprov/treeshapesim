# Simulating trees with controlled shape

``` r

library(treeshapesim)
set.seed(1)
```

## Two shape statistics, two processes

The gamma statistic of Pybus and Harvey (2000) summarizes the
distribution of branching times: it is near zero for Yule trees and
negative when branching times are concentrated near the root. The
Colless index summarizes topological imbalance; normalized against its
Yule expectation, it is near zero for Yule trees and positive for more
imbalanced ones.
[`tree_shape()`](https://diogoprov.github.io/treeshapesim/reference/tree_shape.md)
returns both.

[`sim_tree_shape()`](https://diogoprov.github.io/treeshapesim/reference/sim_tree_shape.md)
draws the topology from Aldous’ beta-splitting model (`beta`) and the
internode intervals from a pure-birth process in which each of the k
lineages splits at rate k^c (`c`). The two are drawn independently,
gamma depends only on the intervals and the Colless index only on the
topology, so `c` moves gamma without changing the Colless index and
`beta` does the opposite.

``` r

yule <- sim_tree_shape(100, nsim = 50)
early <- sim_tree_shape(100, c = -1, nsim = 50)
imbalanced <- sim_tree_shape(100, beta = -1.5, nsim = 50)
shapes <- lapply(list(yule = yule, early = early, imbalanced = imbalanced),
                 function(trees) colMeans(t(vapply(trees, tree_shape, numeric(2)))))
do.call(rbind, shapes)
#>                  gamma     colless
#> yule        0.31793823 -0.10063867
#> early      -5.58314659 -0.09283867
#> imbalanced  0.09541742 10.51216133
```

## Reaching target values

The mean of each statistic depends on its parameter and on the number of
tips.
[`calibrate_shape()`](https://diogoprov.github.io/treeshapesim/reference/calibrate_shape.md)
estimates both curves for a given number of tips, and
[`shape_params()`](https://diogoprov.github.io/treeshapesim/reference/shape_params.md)
inverts them. Since each curve is estimated from its own process,
without building trees, the calibration is fast. Use more draws per
value (`reps`, 200 by default) for a study; a few are enough here.

``` r

cal <- calibrate_shape(100, reps = 50)
cal
#> Shape calibration for n = 100 tips (50 draws per value)
#>   mean gamma from -8.46 to 9.44 for c in [ -2 , 1 ]
#>   mean Colless from -1.79 to 34.4 for beta in [ -1.95 , 8 ]
shape_params(cal, gamma = -2, colless = 1)
#>   gamma colless       beta          c
#> 1    -2       1 -0.5449625 -0.2866907
```

The parameters give the target values on average, but single trees vary
around them.
[`sim_tree_target()`](https://diogoprov.github.io/treeshapesim/reference/sim_tree_target.md)
keeps only trees whose statistics fall within a tolerance of the
targets. Because topology and intervals are independent, it rejects them
separately, which gives the same distribution of trees as rejecting
whole trees, with a much higher acceptance rate.

``` r

trees <- sim_tree_target(100, gamma = -2, colless = 1, nsim = 20, tol = 0.25,
                         calibration = cal, root_age = 50)
summary(attr(trees, "shape"))
#>      gamma           colless      
#>  Min.   :-2.242   Min.   :0.8008  
#>  1st Qu.:-2.085   1st Qu.:1.0058  
#>  Median :-1.957   Median :1.0608  
#>  Mean   :-1.971   Mean   :1.0398  
#>  3rd Qu.:-1.827   3rd Qu.:1.1058  
#>  Max.   :-1.756   Max.   :1.2008
attr(trees, "acceptance")
#>  topology intervals 
#> 0.1538462 0.1739130
```

## A grid of tree shapes

For simulation studies,
[`shape_grid()`](https://diogoprov.github.io/treeshapesim/reference/shape_grid.md)
crosses levels of both statistics. The levels can come, for instance,
from quantiles of
[`tree_shape()`](https://diogoprov.github.io/treeshapesim/reference/tree_shape.md)
over empirical trees of a similar size; here they are set by hand.

``` r

grid <- shape_grid(100, gamma = c(-2, -1, 0), colless = c(0, 0.5, 1),
                   calibration = cal)
grid
#>   gamma colless        beta           c
#> 1    -2     0.0  0.03840196 -0.28669067
#> 2    -1     0.0  0.03840196 -0.14838659
#> 3     0     0.0  0.03840196 -0.01769295
#> 4    -2     0.5 -0.26393708 -0.28669067
#> 5    -1     0.5 -0.26393708 -0.14838659
#> 6     0     0.5 -0.26393708 -0.01769295
#> 7    -2     1.0 -0.54496255 -0.28669067
#> 8    -1     1.0 -0.54496255 -0.14838659
#> 9     0     1.0 -0.54496255 -0.01769295
cell <- grid[grid$gamma == -1 & grid$colless == 0.5, ]
cell_trees <- sim_tree_target(100, gamma = cell$gamma, colless = cell$colless,
                              nsim = 10, calibration = cal)
```

## Placing an empirical tree on the same scale

[`fit_shape()`](https://diogoprov.github.io/treeshapesim/reference/fit_shape.md)
estimates `beta` from the topology and `c` from the branching times by
maximum likelihood (the beta-splitting model and the DDX model of the R
package laser, with `c = -x`), and a parametric bootstrap with the
package’s own generators gives percentile intervals. The estimate of
`beta` is imprecise for balanced trees with few tips, and `c` describes
the shape of the lineage-through-time curve, not its cause.

``` r

empirical <- sim_tree_shape(150, beta = -1, c = -0.7)   # stands in for a real tree
fit_shape(empirical, boot = 30)
#> Shape estimates for a tree with 150 tips
#>      estimate   2.5%  97.5%
#> beta   -1.200 -1.360 -0.835
#> c      -0.678 -0.856 -0.496
#> gamma = -7.5 | normalized Colless = 7.13
```

## References

Aldous, D. (1996). Probability distributions on cladograms. In D. Aldous
& R. Pemantle (Eds.), *Random discrete structures* (pp. 1–18). Springer.
<https://doi.org/10.1007/978-1-4612-0719-1_1>

Pybus, O. G., & Harvey, P. H. (2000). Testing macro-evolutionary models
using incomplete molecular phylogenies. *Proceedings of the Royal
Society B*, 267, 2267–2272. <https://doi.org/10.1098/rspb.2000.1278>
