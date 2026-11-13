# Profile Likelihood Plot

Profile Likelihood Plot

## Usage

``` r
# S3 method for class 'parframe'
plotProfile(
  profs,
  ...,
  maxvalue = 5,
  parlist = NULL,
  ncol = NULL,
  threshold = .profileLines()
)

# S3 method for class 'list'
plotProfile(
  profs,
  ...,
  maxvalue = 5,
  parlist = NULL,
  ncol = NULL,
  threshold = .profileLines()
)

plotProfile(profs, ...)
```

## Arguments

- profs:

  A [parframe](https://dmodverse.github.io/dMod2/reference/parframe.md)
  of profiles as returned by
  [profile()](https://dmodverse.github.io/dMod2/reference/profile.md),
  or a list of those.

- ...:

  Logical expressions to subset the plotted data.

- maxvalue:

  Numeric, the value of the objective difference at which profiles are
  cut off. Defaults to 5.

- parlist:

  Matrix or data frame with parameter columns, drawn as points. With a
  `value` column, the differences of `value` to the lowest profile value
  are used.

- ncol:

  Number of columns of the plot grid. `NULL` (default) lets
  [`ggplot2::facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html)
  choose.

- threshold:

  Numeric, the horizontal lines and y-axis breaks above the optimum,
  which is always drawn. Defaults to the chi-square thresholds for 68%,
  90% and 95%. Names are used as axis labels. Pass the value of
  [`profileThreshold()`](https://dmodverse.github.io/dMod2/reference/profileThreshold.md)
  to match
  [`confint.parframe()`](https://dmodverse.github.io/dMod2/reference/confint.parframe.md).

## Value

A `ggplot` object with the plotted data frame as attribute `"data"`.

## See also

[profile()](https://dmodverse.github.io/dMod2/reference/profile.md),
[`plotPaths()`](https://dmodverse.github.io/dMod2/reference/plotPaths.md)

## Examples

``` r
pars <- c(a = 1, b = 0.5)
obj <- constraintL2(mu = pars, sigma = 0.1)
profs <- profile(obj, pars, whichPar = c("a", "b"), limits = c(-1, 1),
                 cores = 1)
plotProfile(profs)
```
