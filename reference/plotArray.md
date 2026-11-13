# Plot an Array of Trajectories Along the Profile of a Parameter

Predicts the model at parameter sets taken along one profile, starting
at the optimum, and colors the trajectories by the value of the profiled
parameter.

## Usage

``` r
plotArray(
  par,
  profs,
  prd,
  times,
  ...,
  direction = c("up", "down"),
  covtable = NULL,
  nSim = 4,
  fixed = NULL
)
```

## Arguments

- par:

  Character, the profiled parameter.

- profs:

  A [parframe](https://dmodverse.github.io/dMod2/reference/parframe.md)
  of profiles as returned by
  [profile()](https://dmodverse.github.io/dMod2/reference/profile.md).

- prd:

  A prediction function, e.g. `g*x*p`, called as
  `prd(times, pars, deriv = FALSE)`.

- times:

  Numeric vector of time points for the prediction.

- ...:

  Logical expressions in the columns `time`, `name`, `value`,
  `condition` and those of `covtable` to subset the plotted data.

- direction:

  `"up"` (default) or `"down"`, the side of the profile that is traced
  from the optimum.

- covtable:

  Condition table, e.g. from
  [`covariates()`](https://dmodverse.github.io/dMod2/reference/covariates.md),
  merged with the predictions by condition so that `...` can refer to
  the table columns. Default `NULL`.

- nSim:

  Number of trajectories, evenly spaced along the profile side and
  starting at the optimum; at most the number of profile points on that
  side. Defaults to 4.

- fixed:

  Named numeric vector of the parameters held fixed in the profiles (the
  `fixed` of
  [profile()](https://dmodverse.github.io/dMod2/reference/profile.md)),
  added to every parameter set. Default `NULL`.

## Value

A `ggplot` object.

## See also

[`plotProfile()`](https://dmodverse.github.io/dMod2/reference/plotProfile.md),
[`plotPathsMulti()`](https://dmodverse.github.io/dMod2/reference/plotPathsMulti.md)

## Examples

``` r
times <- 0:5
grid <- data.frame(name = "A", time = times, row.names = paste0("A", times))
x <- Xd(grid)
pars <- structure(exp(-times / 2), names = getParameters(x))
obj <- constraintL2(mu = pars, sigma = 0.1)
profs <- profile(obj, pars, whichPar = "A2", limits = c(-1, 1), cores = 1)
plotArray("A2", profs, x, seq(0, 5, 0.1))
```
