# Random topology under Aldous' beta-splitting model

A clade of m tips splits into i and m - i tips with probability
proportional to Gamma(beta + i + 1) Gamma(beta + m - i + 1) /
(Gamma(i + 1) Gamma(m - i + 1)). `beta = 0` is Yule, `beta = -1.5` is
PDA, `beta` close to -2 is caterpillar-like and large `beta` is
balanced.

## Usage

``` r
rtopology_beta(n, beta = 0)
```

## Arguments

- n:

  Number of tips (at least 2).

- beta:

  Aldous' beta, greater than -2. Lower values give more imbalanced
  trees.

## Value

A `phylo` object without branch lengths, with tips `t1`, ..., `tn`.

## References

Aldous, D. (1996). Probability distributions on cladograms. In D. Aldous
& R. Pemantle (Eds.), *Random discrete structures* (pp. 1–18). Springer.
[doi:10.1007/978-1-4612-0719-1_1](https://doi.org/10.1007/978-1-4612-0719-1_1)

## Examples

``` r
set.seed(1)
rtopology_beta(10, beta = -1)
#> 
#> Phylogenetic tree with 10 tips and 9 internal nodes.
#> 
#> Tip labels:
#>   t1, t2, t3, t4, t5, t6, ...
#> 
#> Rooted; no branch length.
```
