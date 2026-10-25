# Internode intervals of a pure-birth process with lineage-dependent rates

With k lineages each splits at rate k^c, so the interval with k lineages
is Exp(k^(1 + c)). `c = 0` is Yule; `c < 0` slows diversification as
lineages accumulate (gamma \< 0); `c > 0` accelerates it (gamma \> 0).
The final interval (n lineages, up to the present) is drawn from the
same law.

## Usage

``` r
intervals_lineage_rate(n, c = 0)
```

## Arguments

- n:

  Number of tips.

- c:

  Lineage dependence of the split rate.

## Value

Intervals g_2, ..., g_n (length n - 1), scaled to sum 1.

## References

Pybus, O. G., & Harvey, P. H. (2000). Testing macro-evolutionary models
using incomplete molecular phylogenies. *Proceedings of the Royal
Society B*, 267, 2267–2272.
[doi:10.1098/rspb.2000.1278](https://doi.org/10.1098/rspb.2000.1278)

## Examples

``` r
set.seed(1)
intervals_lineage_rate(10, c = -1)
#> [1] 0.09121461 0.14272468 0.01759918 0.01688517 0.05267053 0.34966867 0.14851261
#> [8] 0.06518557 0.11553897
```
