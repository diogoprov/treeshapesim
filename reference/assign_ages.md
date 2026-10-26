# Assign node ages to a topology from internode intervals

Internal nodes are ranked by a uniform random linear extension of the
tree order (the ranking a Yule process implies given the topology): the
next node to split is drawn among the available ones with probability
proportional to the number of internal nodes in its subtree.

## Usage

``` r
assign_ages(topo, g)
```

## Arguments

- topo:

  Binary, rooted `phylo` topology.

- g:

  Internode intervals, as from
  [`intervals_lineage_rate()`](https://diogoprov.github.io/treeshapesim/reference/intervals_lineage_rate.md)
  (length `Ntip(topo) - 1`, positive and finite).

## Value

`topo` with branch lengths; ultrametric, with root depth `sum(g)`.

## Examples

``` r
set.seed(1)
assign_ages(rtopology_beta(10), intervals_lineage_rate(10))
#> 
#> Phylogenetic tree with 10 tips and 9 internal nodes.
#> 
#> Tip labels:
#>   t1, t2, t3, t4, t5, t6, ...
#> 
#> Rooted; includes branch length(s).
```
