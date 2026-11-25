# Profile Likelihood: Plot Profiles Along with Their Parameter Paths

For each profiled parameter, draws the profile above the paths of
[`plotPathsMulti()`](https://dmodverse.github.io/dMod2/reference/plotPathsMulti.md).
Requires the package cowplot.

## Usage

``` r
plotProfilesAndPaths(
  profs,
  whichPar,
  nPars = 5,
  ncol = 3,
  normalizePaths = FALSE,
  modes = c("data", "prior"),
  maxvalue = 5,
  threshold = .profileLines(),
  ...
)
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

- ncol:

  Number of columns of the plot grid. Defaults to 3.

- normalizePaths:

  Logical. `TRUE` scales each path to a maximum absolute value of 1.
  Defaults to `FALSE`.

- modes:

  Character vector, the contributions to the objective drawn in the
  profile plot besides the total. Defaults to `c("data", "prior")`.

- maxvalue, threshold:

  Cut-off of the objective difference and the horizontal lines above the
  optimum, as in
  [`plotProfile()`](https://dmodverse.github.io/dMod2/reference/plotProfile.md).

- ...:

  Further arguments passed to
  [`cowplot::plot_grid()`](https://wilkelab.org/cowplot/reference/plot_grid.html).

## Value

A `ggplot` object, the grid of
[`cowplot::plot_grid()`](https://wilkelab.org/cowplot/reference/plot_grid.html).

## See also

[`plotProfile()`](https://dmodverse.github.io/dMod2/reference/plotProfile.md),
[`plotPathsMulti()`](https://dmodverse.github.io/dMod2/reference/plotPathsMulti.md)

## Examples

``` r
pars <- c(a = 1, b = 0.5, c = 2)
obj <- constraintL2(mu = pars, sigma = 0.1)
profs <- profile(obj, pars, whichPar = c("a", "b"), limits = c(-1, 1),
                 cores = 1)
plotProfilesAndPaths(profs, c("a", "b"), nPars = 2, ncol = 2)
```
