# Plot Parameter Values of a Collection of Fits

Draws one box per parameter, colored by the step of the objective value
the fits belong to.

## Usage

``` r
# S3 method for class 'parlist'
plotPars(x, ..., tol = 1)

# S3 method for class 'parframe'
plotPars(x, ..., tol = 1)

plotPars(x, ...)

# S3 method for class 'em'
plotPars(x, ..., level = 0.95, which = NULL)

# S3 method for class 'list'
plotPars(x, ..., level = 0.95, which = NULL)
```

## Arguments

- x:

  A [parframe](https://dmodverse.github.io/dMod2/reference/parframe.md)
  or the result of
  [`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md).

- ...:

  Logical expressions to subset `x`.

- tol:

  A value more than `tol` above the lowest value of the current step
  starts a new step. Defaults to 1.

- level:

  With an `em` fit or a list of them: the level of the Wald intervals,
  see
  [confint()](https://dmodverse.github.io/dMod2/reference/confint.em.md).
  Default 0.95.

- which:

  With an `em` fit or a list of them: names of the parameters to show,
  by default all estimated ones.

## Value

A `ggplot` object with the plotted data frame as attribute `"data"`.

## Details

An `em` fit, see
[`em()`](https://dmodverse.github.io/dMod2/reference/em.md), or a list
of them, as fits of one model by several methods: the estimates of fixed
effects and hyperparameters with Wald intervals, where the fit has
standard errors, one row per parameter and one color per fit. The fits
of a list are labeled by the names of the list, else by their method.

## See also

[`plotValues()`](https://dmodverse.github.io/dMod2/reference/plotValues.md),
[`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md)

## Examples

``` r
fits <- parframe(data.frame(value = c(10, 10.2, 15, 40), converged = TRUE,
                            iterations = 20, a = c(1, 1.1, 2, 3),
                            b = c(0, 0.1, 1, 2)),
                 parameters = c("a", "b"))
plotPars(fits)
```
