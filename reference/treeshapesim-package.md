# treeshapesim: Simulate Trees with Separately Controlled Imbalance and Branching Times

Simulates rooted, binary, ultrametric trees in which topological
imbalance and the distribution of branching times are set by separate
parameters, and estimates those parameters from empirical trees.

## Terms

- Topology:

  The branching pattern of a tree, without branch lengths.

- Imbalance:

  How unequally the tips are divided between the two daughter clades of
  the internal nodes, measured here by the Colless index normalized
  against its expectation under the Yule model (see
  [`tree_shape()`](https://provetelab.org/treeshapesim/reference/tree_shape.md)).

- beta:

  The parameter of Aldous' (1996) beta-splitting model of topologies:
  `beta = 0` gives the distribution of Yule topologies, lower values
  more imbalanced and higher values more balanced topologies.

- Internode intervals:

  The times \\g_k\\ during which the tree has \\k\\ lineages, from the
  root (\\k = 2\\) to the present (\\k = n\\).

- c:

  The exponent of the dependence of the per-lineage split rate on the
  number of lineages \\k\\, \\k^c\\; `c = 0` is a Yule process.

- gamma:

  Pybus and Harvey's (2000) statistic of the internode intervals:
  negative when splits are concentrated near the root.

- Sampling fraction:

  The proportion of the species of a clade that are in the tree.

## Algorithms

The combination of a beta-splitting topology with intervals from a
lineage-dependent pure-birth process, joined by a uniform random ranking
of the internal nodes, is, as far as we know, a new algorithm; this is
its first implementation. The estimators of `beta` (maximum likelihood
under beta-splitting) and of `c` (the density-dependent exponential
model) were implemented in the R packages apTreeshape and laser, both
archived on CRAN;
[`fit_shape()`](https://provetelab.org/treeshapesim/reference/fit_shape.md)
reimplements them and adds a correction of `c` for incomplete sampling,
which is new.

## Life cycle

The package is in active development and its interface is *experimental*
(version 0.x): function names and arguments may still change, with
changes listed in `NEWS.md`. The core simulators
([`sim_tree_shape()`](https://provetelab.org/treeshapesim/reference/sim_tree_shape.md),
[`rtopology_beta()`](https://provetelab.org/treeshapesim/reference/rtopology_beta.md),
[`intervals_lineage_rate()`](https://provetelab.org/treeshapesim/reference/intervals_lineage_rate.md),
[`assign_ages()`](https://provetelab.org/treeshapesim/reference/assign_ages.md))
are expected to be stable; planned work, including a sampler with
extinction and a CRAN release, is listed in `ROADMAP.md` in the source
repository.

## References

Aldous, D. (1996). Probability distributions on cladograms. In D. Aldous
& R. Pemantle (Eds.), *Random discrete structures* (pp. 1–18). Springer.
[doi:10.1007/978-1-4612-0719-1_1](https://doi.org/10.1007/978-1-4612-0719-1_1)

Pybus, O. G., & Harvey, P. H. (2000). Testing macro-evolutionary models
using incomplete molecular phylogenies. *Proceedings of the Royal
Society B*, 267, 2267–2272.
[doi:10.1098/rspb.2000.1278](https://doi.org/10.1098/rspb.2000.1278)

## See also

Useful links:

- <https://provetelab.org/treeshapesim/>

- <https://github.com/diogoprov/treeshapesim>

- Report bugs at <https://github.com/diogoprov/treeshapesim/issues>

## Author

**Maintainer**: Diogo B. Provete <dbprovete@gmail.com>
([ORCID](https://orcid.org/0000-0002-0097-0651))

Authors:

- Diogo B. Provete <dbprovete@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-0097-0651))
