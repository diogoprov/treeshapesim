# Contributing to treeshapesim

Thank you for considering a contribution. This guide covers bug reports,
changes to the code and the tests.

## Bug reports and questions

Open an issue at <https://github.com/diogoprov/treeshapesim/issues> with a
minimal reproducible example (the [reprex](https://reprex.tidyverse.org)
package helps), including `set.seed()` when the problem involves random draws,
and the output of `sessionInfo()`. The same seed gives the same trees on the
same platform, but not necessarily across platforms (see `?sim_tree_shape`).

## Changes to the code

- Discuss larger changes in an issue first.
- Work on a branch and open a pull request against `main`.
- Follow the existing style: tidyverse style, lines of at most 80 characters,
  roxygen2 documentation for every function (internal functions with `@noRd`).
- Add tests for new behaviour, and a test for every new error or warning
  message, which must be unique.
- Add a bullet to `NEWS.md` under the development version.
- Run `devtools::document()`, `devtools::check()` and
  `pkgcheck::pkgcheck()` before you open the pull request.

### Statistical software standards

The package follows the rOpenSci general and probability-distribution
[statistical software standards](https://stats-devguide.ropensci.org).
Compliance is documented with `@srrstats` tags from the
[srr](https://docs.ropensci.org/srr/) package, next to the code or test that
meets each standard; standards that do not apply are tagged `@srrstatsNA` in
`R/srr-stats-standards.R`. If a change affects a standard, update its tag, and
check with `srr::srr_stats_pre_submit()`.

## Tests

The regular tests run with `devtools::test()` in less than a minute.

### Extended tests

Extended tests of parameter recovery over many seeds and of the decrease of the
error of the estimates with the number of tips are in
`tests/testthat/test-extended.R`. They are skipped unless the environment
variable `TREESHAPESIM_EXTENDED_TESTS` is `"true"`:

```r
Sys.setenv(TREESHAPESIM_EXTENDED_TESTS = "true")
devtools::test(filter = "extended")
```

They need no downloads or special hardware, take about 20 seconds on a
recent laptop, and produce no files to inspect.

### Validation scripts

The longer validation runs behind the performance claims in the documentation
are scripts in `data-raw/` (`validate_estimators.R`, `validate_sampling.R`,
`benchmark_topology.R`), with their results in `data-raw/validation/`. Some
need archived packages (apTreeshape, laser); see the header of each script.

## Code of Conduct

Please note that this project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing, you agree
to abide by its terms.
