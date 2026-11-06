# Remove tips by non-random extinction

Applies one extinction event to a tree and removes a fixed fraction of
its tips. Tips are lost according to `risk`, which links extinction to a
trait or to any other property of the species, or, without `risk`,
according to a latent liability that evolves on the tree, so that
`signal` sets how clustered the losses are: `signal = 0` is random
extinction, and with `signal = 1` close relatives tend to be lost
together, up to whole clades.

## Usage

``` r
prune_tips(tree, fraction, risk = NULL, signal = 1)
```

## Arguments

- tree:

  A binary, rooted `phylo` object with branch lengths.

- fraction:

  Fraction of the tips that are lost, a single number greater than 0 and
  less than 1. `round(fraction * Ntip(tree))` tips are removed, and at
  least two must remain.

- risk:

  `NULL`, or a numeric vector of non-negative weights named with the tip
  labels of `tree`. Tips are drawn without replacement with probability
  proportional to `risk`. To link extinction to a trait, use for
  instance `exp(b * trait)`.

- signal:

  Pagel's lambda of the latent liability, a single number from 0 to 1,
  used when `risk` is `NULL`. The tips with the highest liability are
  lost.

## Value

The remaining tree, with the labels of the lost tips in attribute
`lost`.

## Details

Under a homogeneous birth-death process the shape of the reconstructed
tree is the same as under pure birth (Gernhard 2008), so extinction
changes the topology only if it differs among lineages. The event acts
on the tree as it stands, at the present; it does not model extinction
that differs among lineages throughout the history of the clade, for
which speciation and extinction have to be simulated together (Maddison
et al. 2007; Rabosky 2009; Verboom et al. 2020).

Random extinction removes relatively little evolutionary history (Nee &
May 1997), but extinction risk is clumped on the tree (Purvis et al.
2000), and its effect depends on how risk relates to speciation rates
(Heard & Mooers 2000) and on the shape of the tree (Veron et al. 2017).
The same `fraction` with `signal = 0` is therefore the natural control.

Unlike `sampling_fraction` in
[`sim_tree_shape()`](https://provetelab.org/treeshapesim/reference/sim_tree_shape.md),
the tree is not rescaled and the tip labels are kept, so that traits
stay aligned with the species that remain. If a whole daughter clade of
the root is lost, the root of the remaining tree is younger. The
separation between gamma and the Colless index holds before the event,
not after it, and
[`calibrate_shape()`](https://provetelab.org/treeshapesim/reference/calibrate_shape.md)
does not apply to the remaining tree.

## References

Gernhard, T. (2008). The conditioned reconstructed process. *Journal of
Theoretical Biology*, 253, 769–778.
[doi:10.1016/j.jtbi.2008.04.005](https://doi.org/10.1016/j.jtbi.2008.04.005)

Heard, S. B., & Mooers, A. O. (2000). Phylogenetically patterned
speciation rates and extinction risks change the loss of evolutionary
history during extinctions. *Proceedings of the Royal Society B*, 267,
613–620.
[doi:10.1098/rspb.2000.1046](https://doi.org/10.1098/rspb.2000.1046)

Maddison, W. P., Midford, P. E., & Otto, S. P. (2007). Estimating a
binary character's effect on speciation and extinction. *Systematic
Biology*, 56, 701–710.
[doi:10.1080/10635150701607033](https://doi.org/10.1080/10635150701607033)

Nee, S., & May, R. M. (1997). Extinction and the loss of evolutionary
history. *Science*, 278, 692–694.
[doi:10.1126/science.278.5338.692](https://doi.org/10.1126/science.278.5338.692)

Purvis, A., Agapow, P.-M., Gittleman, J. L., & Mace, G. M. (2000).
Nonrandom extinction and the loss of evolutionary history. *Science*,
288, 328–330.
[doi:10.1126/science.288.5464.328](https://doi.org/10.1126/science.288.5464.328)

Rabosky, D. L. (2009). Heritability of extinction rates links
diversification patterns in molecular phylogenies and fossils.
*Systematic Biology*, 58, 629–640.
[doi:10.1093/sysbio/syp069](https://doi.org/10.1093/sysbio/syp069)

Verboom, G. A., Boucher, F. C., Ackerly, D. D., Wootton, L. M., &
Freyman, W. A. (2020). Species selection regime and phylogenetic tree
shape. *Systematic Biology*, 69, 774–794.
[doi:10.1093/sysbio/syz076](https://doi.org/10.1093/sysbio/syz076)

Veron, S., Davies, T. J., Cadotte, M. W., Clergeau, P., & Pavoine, S.
(2017). Predicting loss of evolutionary history: Where are we?
*Biological Reviews*, 92, 271–291.
[doi:10.1111/brv.12228](https://doi.org/10.1111/brv.12228)

## Examples

``` r
set.seed(1)
tr <- sim_tree_shape(50, beta = -1)
clustered <- prune_tips(tr, 0.4, signal = 1)
random <- prune_tips(tr, 0.4, signal = 0)
attr(clustered, "lost")
#>  [1] "t35" "t28" "t6"  "t47" "t14" "t46" "t29" "t43" "t44" "t42" "t13" "t45"
#> [13] "t48" "t25" "t5"  "t27" "t49" "t24" "t26" "t16"
# extinction linked to a trait
body_size <- setNames(rnorm(50), tr$tip.label)
prune_tips(tr, 0.4, risk = exp(2 * body_size))
#> 
#> Phylogenetic tree with 30 tips and 29 internal nodes.
#> 
#> Tip labels:
#>   t2, t3, t4, t5, t6, t7, ...
#> 
#> Rooted; includes branch length(s).
```
