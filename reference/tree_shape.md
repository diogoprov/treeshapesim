# Gamma statistic and Colless index of a tree

Computes Pybus and Harvey's (2000) gamma statistic, with
[`ape::gammaStat()`](https://rdrr.io/pkg/ape/man/gammaStat.html), and
the Colless index normalized against its expectation under the Yule
model, \\(I_C - n \log n - n(\gamma_E - 1 - \log 2)) / n\\, where \\n\\
is the number of tips and \\\gamma_E\\ is Euler's constant (the `"yule"`
normalization of
[`treestats::colless()`](https://thijsjanzen.github.io/treestats/reference/colless.html)).
Both are near zero for Yule trees: gamma is negative when branching
times are concentrated near the root, and the normalized Colless index
is positive for trees more imbalanced than Yule trees.

## Usage

``` r
tree_shape(tree)
```

## Arguments

- tree:

  A binary, rooted `phylo` object with at least three tips. The gamma
  statistic assumes an ultrametric tree; this is not checked.

## Value

A named vector with `gamma` and `colless`.

## References

Pybus, O. G., & Harvey, P. H. (2000). Testing macro-evolutionary models
using incomplete molecular phylogenies. *Proceedings of the Royal
Society B*, 267, 2267–2272.
[doi:10.1098/rspb.2000.1278](https://doi.org/10.1098/rspb.2000.1278)

## Examples

``` r
set.seed(1)
tree_shape(sim_tree_shape(100, beta = -1, c = -1))
#>     gamma   colless 
#> -5.956557  2.440761 
```
