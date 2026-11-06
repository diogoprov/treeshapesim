# Tree shape and extinction in metacommunity simulations with comdynamite

treeshapesim and [comdynamite](https://provetelab.org/comdynamite/)
divide the work of a simulation study. treeshapesim produces the
phylogeny of the regional species pool, with the imbalance and the
branching times chosen by the user. comdynamite takes that tree as a
fixed input, places traits on it and runs the dynamics of the
metacommunity. This article connects the two through a question in which
the tree is not just a backdrop: what happens to communities when the
pool has lost species, and does it matter which species were lost?

Neither package depends on the other, so comdynamite has to be installed
separately.

``` r

remotes::install_github("diogoprov/comdynamite")
```

``` r

library(treeshapesim)
library(comdynamite)
```

## Why the pattern of extinction matters

If all lineages speciate and go extinct at the same constant rates, the
shape of the reconstructed tree is the same as under pure birth
(Gernhard 2008). Such extinction changes the branching times, but not
the topology. Extinction moves the tree away from that expectation when
it differs among lineages.

The same contrast organizes the work on the loss of evolutionary
history. For a randomly evolving clade, Nee and May (1997) showed that
about 80% of the tree can survive the loss of about 95% of the species,
and that choosing the survivors to preserve the most history is not much
better than choosing them at random. Current extinction risk, however,
is not phylogenetically random, and its clumped distribution puts more
higher taxa and more evolutionary history at risk (Purvis et al. 2000).
In simulations with phylogenetically patterned speciation rates and
extinction risks, Heard and Mooers (2000) found that risks uncorrelated
with speciation rates removed only slightly more history than in random
trees, whereas risks correlated with speciation rates could increase or
decrease the loss strongly. Clustered risk has empirical support, but
does not always cause more loss of phylogenetic diversity than random
extinction: the shape of the tree also matters (Veron et al. 2017). That
dependence is what
[`sim_tree_shape()`](https://provetelab.org/treeshapesim/reference/sim_tree_shape.md)
can vary directly, and comdynamite adds the next question, which is what
the losses do to communities.

[`prune_tips()`](https://provetelab.org/treeshapesim/reference/prune_tips.md)
applies a single round of extinction to a complete tree, as in Nee and
May (1997), but the losses need not be random. It removes a fixed
fraction of the tips, so that scenarios can be compared at the same
richness, in one of three ways.

``` r

set.seed(1)
tree <- sim_tree_shape(60, beta = -1.5)
traits <- species_traits(tree, lambda = 1)
risk <- setNames(exp(3 * as.numeric(scale(traits$optimum))), rownames(traits))

random <- prune_tips(tree, 0.4, signal = 0)     # the control
clustered <- prune_tips(tree, 0.4, signal = 1)  # close relatives lost together
linked <- prune_tips(tree, 0.4, risk = risk)    # species with high optima at risk
head(attr(linked, "lost"))
#> [1] "t54" "t1"  "t59" "t58" "t10" "t12"
```

The tip labels are kept, so the traits of the survivors are found by
name. Traits are drawn on the complete tree, before the event, because
extinction acts on species that already have them.

## One simulation

The remaining tree is the regional pool of the metacommunity.

``` r

landscape <- landscape_implicit(20)
sim <- simulate_metacommunity(
  landscape = landscape,
  traits    = traits[linked$tip.label, ],
  alpha     = competition_matrix(linked, type = "phylo"),
  dispersal = 0.01,
  n_steps   = 300
)
dim(sim$abundance)    # sites x species
#> [1] 20 36
```

## A factorial design

The function below runs the whole sequence once and returns the shape of
the pool and a few properties of the communities. Environmental optima
and the environment of the sites share the range from 0 to 100, so the
species at risk under `risk` are those of the sites with high values.

``` r

one_run <- function(beta, scenario, n = 60, fraction = 0.4) {
  tree <- sim_tree_shape(n, beta = beta)
  traits <- species_traits(tree, lambda = 1)
  risk <- setNames(exp(3 * as.numeric(scale(traits$optimum))), rownames(traits))
  pool <- switch(scenario,
    none      = tree,
    random    = prune_tips(tree, fraction, signal = 0),
    clustered = prune_tips(tree, fraction, signal = 1),
    linked    = prune_tips(tree, fraction, risk = risk))
  landscape <- landscape_implicit(20)
  sim <- simulate_metacommunity(landscape, traits[pool$tip.label, ],
                                competition_matrix(pool, type = "phylo"),
                                dispersal = 0.01, n_steps = 300)
  present <- sim$abundance > 0
  richness <- rowSums(present)
  high <- landscape$env > 50
  c(tree_shape(pool),
    pd = sum(pool$edge.length) / sum(tree$edge.length),
    richness_low = mean(richness[!high]),
    richness_high = mean(richness[high]))
}

set.seed(2)
cells <- expand.grid(scenario = c("none", "random", "clustered", "linked"),
                     beta = c(0, -1.5), stringsAsFactors = FALSE)
means <- t(mapply(function(beta, scenario) {
  rowMeans(replicate(100, one_run(beta, scenario)))
}, cells$beta, cells$scenario))
cbind(cells, round(means, 2))
#>    scenario beta gamma colless   pd richness_low richness_high
#> 1      none  0.0 -0.12   -0.03 1.00        16.03         16.36
#> 2    random  0.0 -0.93    0.10 0.76        11.63         11.27
#> 3 clustered  0.0  0.21    0.14 0.68        11.29         10.95
#> 4    linked  0.0 -0.30    0.22 0.71        13.35          9.27
#> 5      none -1.5 -0.01    6.98 1.00        16.28         16.48
#> 6    random -1.5 -0.61    3.79 0.71        11.32         11.33
#> 7 clustered -1.5 -0.14    3.69 0.66        11.24         10.54
#> 8    linked -1.5 -0.38    4.08 0.67        13.06          9.58
```

The comparison that matters is among the three scenarios of extinction,
which have the same richness, and not between them and the complete
pool, because the Colless index and gamma also depend on the number of
tips. In this example:

- **Phylogenetic diversity.** At the same richness, clustered and
  trait-linked extinction kept less of the total branch length than
  random extinction (about 0.70 against 0.77 with `beta = 0`, and 0.66
  against 0.71 with `beta = -1.5`), and each scenario kept less in the
  imbalanced pools.
- **Gamma.** Random losses lowered gamma, as incomplete sampling does,
  because they remove recent splits throughout the tree. Clustered
  losses left it near zero: they remove whole clades, and the recent
  splits of the other clades remain.
- **Local richness.** Random and clustered extinction reduced richness
  evenly along the environmental gradient. Only trait-linked extinction
  created a difference (about 13 species at sites with low values
  against 9 at sites with high values), because the species adapted to
  one end of the gradient were lost.

So the same loss of species can leave the same number of species per
site and different phylogenetic structure, or a different distribution
of richness along the environment, depending on which species were lost.

## Limits of this design

- The event happens once, at the present. It is not extinction that
  differs among lineages throughout the history of a clade, for which
  speciation and extinction have to be simulated together (Maddison et
  al. 2007; Rabosky 2009; Verboom et al. 2020), and which is outside the
  scope of treeshapesim.
- Gamma and the Colless index vary independently only before the event.
  After it, they are outcomes of the design, and
  [`calibrate_shape()`](https://provetelab.org/treeshapesim/reference/calibrate_shape.md)
  does not apply.
- `signal` and `risk` describe the pattern of the losses, not their
  cause.
- With a trait that carries phylogenetic signal, as here, extinction
  linked to the trait is also clustered on the tree. To separate the
  two, use a lower `lambda` in
  [`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.html).
- The numbers above come from 100 replicates of a small example, and
  only illustrate the workflow. A study needs its own choice of pool
  size, fraction and dynamics.

## References

Gernhard, T. (2008). The conditioned reconstructed process. *Journal of
Theoretical Biology*, 253, 769–778.
<https://doi.org/10.1016/j.jtbi.2008.04.005>

Heard, S. B., & Mooers, A. Ø. (2000). Phylogenetically patterned
speciation rates and extinction risks change the loss of evolutionary
history during extinctions. *Proceedings of the Royal Society B*, 267,
613–620. <https://doi.org/10.1098/rspb.2000.1046>

Maddison, W. P., Midford, P. E., & Otto, S. P. (2007). Estimating a
binary character’s effect on speciation and extinction. *Systematic
Biology*, 56, 701–710. <https://doi.org/10.1080/10635150701607033>

Nee, S., & May, R. M. (1997). Extinction and the loss of evolutionary
history. *Science*, 278, 692–694.
<https://doi.org/10.1126/science.278.5338.692>

Purvis, A., Agapow, P.-M., Gittleman, J. L., & Mace, G. M. (2000).
Nonrandom extinction and the loss of evolutionary history. *Science*,
288, 328–330. <https://doi.org/10.1126/science.288.5464.328>

Rabosky, D. L. (2009). Heritability of extinction rates links
diversification patterns in molecular phylogenies and fossils.
*Systematic Biology*, 58, 629–640.
<https://doi.org/10.1093/sysbio/syp069>

Verboom, G. A., Boucher, F. C., Ackerly, D. D., Wootton, L. M., &
Freyman, W. A. (2020). Species selection regime and phylogenetic tree
shape. *Systematic Biology*, 69, 774–794.
<https://doi.org/10.1093/sysbio/syz076>

Veron, S., Davies, T. J., Cadotte, M. W., Clergeau, P., & Pavoine, S.
(2017). Predicting loss of evolutionary history: Where are we?
*Biological Reviews*, 92, 271–291. <https://doi.org/10.1111/brv.12228>
