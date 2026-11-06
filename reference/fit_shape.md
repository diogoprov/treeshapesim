# Estimate the shape parameters of a tree

Estimates, by maximum likelihood, the `beta` of Aldous' beta-splitting
model from the topology and the `c` of the lineage-dependent pure-birth
process from the branching times, which places an empirical tree on the
same scale as the trees of
[`sim_tree_shape()`](https://provetelab.org/treeshapesim/reference/sim_tree_shape.md).
Optionally, a parametric bootstrap with the package's own generators
gives percentile intervals.

## Usage

``` r
fit_shape(
  tree,
  boot = 0,
  level = 0.95,
  beta_range = c(-1.999, 10),
  c_range = c(-5, 5),
  tol = .Machine$double.eps^0.5,
  sampling_fraction = 1,
  correction_reps = 100
)
```

## Arguments

- tree:

  A binary, rooted, ultrametric `phylo` object with at least four tips.

- boot:

  Number of parametric bootstrap replicates (0 for none).

- level:

  Coverage of the bootstrap percentile intervals.

- beta_range, c_range:

  Search intervals for `beta` (lower bound above -2) and `c`.

- tol:

  Tolerance of the ultrametricity check, as in
  [`ape::is.ultrametric()`](https://rdrr.io/pkg/ape/man/is.ultrametric.html)
  (scaled range of root-to-tip distances). Trees with rounded branch
  lengths may need a larger value, or can be made ultrametric first (for
  instance with
  [`phytools::force.ultrametric()`](https://rdrr.io/pkg/phytools/man/force.ultrametric.html)).

- sampling_fraction:

  Fraction of the species of the clade that are in `tree`, a single
  number greater than 0 and at most 1. Below 1, `c` is corrected for
  incomplete sampling (see Details).

- correction_reps:

  Number of trees at each value of the grid used in the correction for
  incomplete sampling.

## Value

An object of class `shape_fit`: a list with `n`, `beta`, `c`, `r`
(speciation rate per lineage at one lineage, in the time units of the
tree), `loglik_beta` (log-likelihood of `beta` up to an additive
constant, so not comparable with other software), `loglik_c` (profile
log-likelihood of `c`, equal to the log-likelihood of laser's DDX model
at its maximum), `beta_at_bound` and `c_at_bound` (estimate at a bound
of the search interval), `beta_range`, `c_range`, `shape` (output of
[`tree_shape()`](https://provetelab.org/treeshapesim/reference/tree_shape.md))
and, with `boot > 0`, `boot` (data frame of bootstrap estimates) and
`ci` (matrix of percentile intervals). With `sampling_fraction` below 1,
`c` is the corrected estimate, `r` is `NA`, `loglik_c` refers to the
uncorrected estimate, and the list also holds `c_uncorrected`,
`sampling_fraction` and `correction` (data frame with the curve).

## Details

**beta.** Each internal node with a clade of \\m\\ tips split into \\i\\
and \\m - i\\ tips contributes the probability of that split under the
beta-splitting model (see
[`rtopology_beta()`](https://provetelab.org/treeshapesim/reference/rtopology_beta.md)),
and the log-likelihood is the sum over nodes (Aldous 1996; this is the
estimator applied by Blum & Francois 2006 and implemented in
`maxlik.betasplit()` of the archived package apTreeshape). It is
maximized over `beta_range`. In the package's validation runs
(`data-raw/validate_estimators.R` in the source repository), the
estimate was nearly unbiased for `beta` of -1 or less, but had a
positive bias and a large variance for balanced trees with few tips
(`beta` of 0 or more with 50 to 100 tips). An estimate at a bound of
`beta_range` only says that the topology is extreme for the model.

**c.** With \\k\\ lineages, each lineage splits at rate \\r k^c\\, which
is the density-dependent exponential (DDX) model of the R package laser
(Rabosky 2006), whose per-lineage rate is \\r k^{-x}\\, so that \\c =
-x\\. As in laser, the splits at 2 to \\n - 1\\ lineages are events and
the last interval, from the most recent split to the present, is not;
\\r\\ is profiled out (its maximum likelihood value is \\(n - 2) /
\sum_k k^{1 + c} g_k\\, where \\g_k\\ is the time with \\k\\ lineages).
The estimate describes the shape of the lineage-through-time curve, not
its cause: branching times alone do not separate diversity-dependent
from time-dependent diversification (Pannetier et al. 2021).

**Incomplete sampling.** Trees with incomplete sampling give estimates
of `c` that are too low, because recent splits are missing; the estimate
of `beta` is little affected and is not corrected. With
`sampling_fraction` below 1, `c` is corrected by simulation: trees are
drawn with
[`sim_tree_shape()`](https://provetelab.org/treeshapesim/reference/sim_tree_shape.md)
at the estimated `beta`, the same `sampling_fraction` and a grid of `c`
from 1 below to 3 above the uncorrected estimate, and the mean
uncorrected estimate at each value (`correction_reps` trees) gives a
curve, made monotone by isotonic regression and inverted at the estimate
from `tree`. This assumes that tips are missing at random and that the
fraction is known; sampling biased towards old or distinct lineages is
not corrected. The result varies slightly with the seed.

**Bootstrap.** With `boot > 0`, topologies are drawn with
[`rtopology_beta()`](https://provetelab.org/treeshapesim/reference/rtopology_beta.md)
at the estimated `beta` and intervals with
[`intervals_lineage_rate()`](https://provetelab.org/treeshapesim/reference/intervals_lineage_rate.md)
at the estimated `c` (with `sampling_fraction` below 1, whole trees with
[`sim_tree_shape()`](https://provetelab.org/treeshapesim/reference/sim_tree_shape.md),
and `c` is corrected with the same curve), and both parameters are
re-estimated, following the recommendation of a parametric bootstrap for
diversification models by Etienne et al. (2016). The spread of the
bootstrap estimates around the point estimate includes the small-sample
bias of the estimators, and, for `c`, a small additional bias because
[`intervals_lineage_rate()`](https://provetelab.org/treeshapesim/reference/intervals_lineage_rate.md)
draws the last interval from the same law, whereas the likelihood treats
it as censored. Percentile intervals need a few hundred replicates.

## References

Aldous, D. (1996). Probability distributions on cladograms. In D. Aldous
& R. Pemantle (Eds.), *Random discrete structures* (pp. 1–18). Springer.
[doi:10.1007/978-1-4612-0719-1_1](https://doi.org/10.1007/978-1-4612-0719-1_1)

Blum, M. G. B., & Francois, O. (2006). Which random processes describe
the tree of life? A large-scale study of phylogenetic tree imbalance.
*Systematic Biology*, 55, 685–691.
[doi:10.1080/10635150600889625](https://doi.org/10.1080/10635150600889625)

Etienne, R. S., Pigot, A. L., & Phillimore, A. B. (2016). How reliably
can we infer diversity-dependent diversification from phylogenies?
*Methods in Ecology and Evolution*, 7, 1092–1099.
[doi:10.1111/2041-210X.12565](https://doi.org/10.1111/2041-210X.12565)

Pannetier, T., Martinez, C., Bunnefeld, L., & Etienne, R. S. (2021).
Branching patterns in phylogenies cannot distinguish diversity-dependent
diversification from time-dependent diversification. *Evolution*, 75,
25–38. [doi:10.1111/evo.14124](https://doi.org/10.1111/evo.14124)

Rabosky, D. L. (2006). LASER: a maximum likelihood toolkit for detecting
temporal shifts in diversification rates from molecular phylogenies.
*Evolutionary Bioinformatics*, 2.
[doi:10.1177/117693430600200024](https://doi.org/10.1177/117693430600200024)

## Examples

``` r
set.seed(1)
tr <- sim_tree_shape(200, beta = -1, c = -0.8)
fit_shape(tr)
#> Shape estimates for a tree with 200 tips
#>   beta      c 
#> -0.864 -0.895 
#> gamma = -7.59 | normalized Colless = 3.03 
fit <- fit_shape(tr, boot = 50)
fit$ci
#>           2.5%      97.5%
#> beta -1.102879 -0.5561590
#> c    -1.046498 -0.7468201
# a tree with half of the species of the clade
sampled <- sim_tree_shape(100, c = -0.5, sampling_fraction = 0.5)
fit_shape(sampled, sampling_fraction = 0.5, correction_reps = 20)
#> Shape estimates for a tree with 100 tips
#> c corrected for a sampling fraction of 0.5 (uncorrected: -0.876)
#>   beta      c 
#>  1.030 -0.536 
#> gamma = -6.18 | normalized Colless = -1.06 
```
