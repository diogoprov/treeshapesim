# treeshapesim

Simulate rooted, binary, ultrametric trees in which topological
imbalance and the distribution of branching times are set by separate
parameters.

- **Topology:** Aldous’ (1996) beta-splitting model. `beta = 0` gives
  Yule trees, lower values give more imbalanced trees and higher values
  more balanced ones.
- **Branching times:** a pure-birth process in which each of the k
  lineages splits at rate k^c. `c = 0` is Yule, `c < 0` slows
  diversification towards the present (negative gamma; Pybus &
  Harvey, 2000) and `c > 0` speeds it up.

The gamma statistic depends only on the internode intervals and the
Colless index only on the topology, so `beta` mainly moves the Colless
index and `c` mainly moves gamma.

Documentation: <https://diogoprov.github.io/treeshapesim/>

## Installation

``` r

install.packages("treeshapesim", repos = c("https://diogoprov.r-universe.dev",
                                           "https://cloud.r-project.org"))
# or
remotes::install_github("diogoprov/treeshapesim")
```

## Example

``` r

library(treeshapesim)
set.seed(1)
tr <- sim_tree_shape(100, beta = -1, c = -0.5)   # imbalanced, early splits
tree_shape(tr)                                    # gamma and normalized Colless

# trees near target values of gamma and Colless
cal <- calibrate_shape(100)                       # calibration for 100 tips
trees <- sim_tree_target(100, gamma = -2, colless = 1, nsim = 50, tol = 0.25,
                         calibration = cal)

# a 3 x 3 grid of tree shapes and the parameters of each cell
grid <- shape_grid(100, gamma = c(-2, -1, 0), colless = c(0, 0.5, 1),
                   calibration = cal)
```

See
[`vignette("treeshapesim")`](https://diogoprov.github.io/treeshapesim/articles/treeshapesim.md)
for a walk-through. Many other tree statistics can be computed with
[treestats](https://github.com/thijsjanzen/treestats).

## References

Aldous, D. (1996). Probability distributions on cladograms. In D. Aldous
& R. Pemantle (Eds.), *Random discrete structures* (pp. 1–18). Springer.
<https://doi.org/10.1007/978-1-4612-0719-1_1>

Pybus, O. G., & Harvey, P. H. (2000). Testing macro-evolutionary models
using incomplete molecular phylogenies. *Proceedings of the Royal
Society B*, 267, 2267–2272. <https://doi.org/10.1098/rspb.2000.1278>
