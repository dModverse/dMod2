# Profile Likelihood: Plot All Parameter Paths of One Profile Together

Draws, for each profiled parameter, the changes of all other parameters
along the profile of that parameter in one panel. The `nPars` parameters
with the largest change are colored and named, the others are gray.
Requires the package cowplot.

## Usage

``` r
plotPathsMulti(profs, whichPar, nPars = 5, normalizePaths = FALSE)
```

## Arguments

- profs:

  A [parframe](https://dmodverse.github.io/dMod2/reference/parframe.md)
  of profiles as returned by
  [profile()](https://dmodverse.github.io/dMod2/reference/profile.md).

- whichPar:

  Character vector, the profiled parameters to draw.

- nPars:

  Integer, number of colored and named paths. Defaults to 5.

- normalizePaths:

  Logical. `TRUE` scales each path to a maximum absolute value of 1.
  Defaults to `FALSE`.

## Value

A `ggplot` object; for several `whichPar` the grid of
[`cowplot::plot_grid()`](https://wilkelab.org/cowplot/reference/plot_grid.html).

## See also

[`plotPaths()`](https://dmodverse.github.io/dMod2/reference/plotPaths.md),
[`plotProfilesAndPaths()`](https://dmodverse.github.io/dMod2/reference/plotProfilesAndPaths.md)

## Examples

``` r
pars <- c(a = 1, b = 0.5, c = 2)
obj <- constraintL2(mu = pars, sigma = 0.1)
profs <- profile(obj, pars, whichPar = c("a", "b"), limits = c(-1, 1),
                 cores = 1)
plotPathsMulti(profs, c("a", "b"), nPars = 2)
```
