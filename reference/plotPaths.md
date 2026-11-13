# Profile Likelihood: Plot of the Parameter Paths

Profile Likelihood: Plot of the Parameter Paths

## Usage

``` r
plotPaths(
  profs,
  ...,
  whichPar = NULL,
  sort = FALSE,
  relative = TRUE,
  scales = "fixed"
)
```

## Arguments

- profs:

  A [parframe](https://dmodverse.github.io/dMod2/reference/parframe.md)
  of profiles as returned by
  [profile()](https://dmodverse.github.io/dMod2/reference/profile.md),
  or a list of those.

- ...:

  Logical expressions to subset the plotted data.

- whichPar:

  Character or index vector, the profiled parameters whose paths are
  drawn. `NULL` (default) takes all.

- sort:

  Logical. `TRUE` sorts each pair of parameters so that a pair appears
  once; `FALSE` (default) keeps both orders.

- relative:

  Logical. `TRUE` (default) shifts each path to start at the optimum.

- scales:

  Character, `"fixed"` (default) or `"free"`.

## Value

A `ggplot` object with the plotted data frame as attribute `"data"`.

## See also

[profile()](https://dmodverse.github.io/dMod2/reference/profile.md),
[`plotProfile()`](https://dmodverse.github.io/dMod2/reference/plotProfile.md),
[`plotPathsMulti()`](https://dmodverse.github.io/dMod2/reference/plotPathsMulti.md)

## Examples

``` r
pars <- c(a = 1, b = 0.5)
obj <- constraintL2(mu = pars, sigma = 0.1)
profs <- profile(obj, pars, whichPar = c("a", "b"), limits = c(-1, 1),
                 cores = 1)
plotPaths(profs)

plotPaths(profs, whichPar = "a")
```
