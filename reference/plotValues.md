# Plot Objective Values of a Collection of Fits

Draws the waterfall plot of a fit collection: the objective values in
ascending order against their rank. A converged fit is colored by the
number of iterations, an unconverged one is red.

## Usage

``` r
# S3 method for class 'alignment'
plotValues(
  x,
  ...,
  tol = 1,
  showSteps = FALSE,
  pseudoLog = FALSE,
  relative = FALSE
)

# S3 method for class 'parframe'
plotValues(
  x,
  ...,
  tol = 1,
  showSteps = FALSE,
  pseudoLog = FALSE,
  relative = FALSE
)

# S3 method for class 'data.frame'
plotValues(
  x,
  ...,
  tol = 1,
  showSteps = FALSE,
  pseudoLog = FALSE,
  relative = FALSE
)

# S3 method for class 'parlist'
plotValues(
  x,
  ...,
  tol = 1,
  showSteps = FALSE,
  pseudoLog = FALSE,
  relative = FALSE
)

# S3 method for class 'list'
plotValues(
  x,
  ...,
  tol = 1,
  showSteps = FALSE,
  pseudoLog = FALSE,
  relative = FALSE
)

plotValues(x, ...)

# S3 method for class 'em'
plotValues(
  x,
  ...,
  tol = 1,
  showSteps = FALSE,
  pseudoLog = FALSE,
  relative = FALSE
)
```

## Arguments

- x:

  A collection of fits: a
  [parframe](https://dmodverse.github.io/dMod2/reference/parframe.md),
  the result of
  [`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md),
  a data frame with column `value` and optionally `converged` and
  `iterations`, a fit of
  [`msem()`](https://dmodverse.github.io/dMod2/reference/msem.md) or an
  alignment from
  [`align()`](https://dmodverse.github.io/dMod2/reference/align.md),
  whose starts are plotted. Or a named list of collections, one panel
  each.

- ...:

  Logical expressions to subset `x`, every collection of a list.

- tol:

  Maximal difference between neighboring objective values that still
  counts as one step. Defaults to 1.

- showSteps:

  Logical. `TRUE` marks the detected steps by dashed vertical lines
  labeled by their index. Defaults to `FALSE`.

- pseudoLog:

  Logical. `TRUE` draws the values on a pseudo-log10 scale, linear below
  1 in absolute value and defined for negative values. Defaults to
  `FALSE`, a linear scale.

- relative:

  `TRUE` draws every value above the best value in the same collection,
  as for collections whose objectives differ by a constant; a number
  draws every value above that number, as the best of all collections.
  Defaults to `FALSE`.

## Value

A `ggplot` object with the plotted data frame as attribute `"data"`; the
column `delta` is the distance to the best fit in the same collection,
`set` the collection of a list.

## See also

[`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md),
[`as.parframe()`](https://dmodverse.github.io/dMod2/reference/as.parframe.md),
[`plotPars()`](https://dmodverse.github.io/dMod2/reference/plotPars.md)

## Examples

``` r
fits <- parframe(data.frame(value = c(10, 10.2, 15, 40), converged = TRUE,
                            iterations = 20, a = c(1, 1.1, 2, 3)),
                 parameters = "a")
plotValues(fits)
```
