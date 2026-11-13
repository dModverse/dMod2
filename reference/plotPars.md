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
```

## Arguments

- x:

  A [parframe](https://dmodverse.github.io/dMod2/reference/parframe.md)
  or the result of
  [`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md).

- ...:

  Logical expressions to subset `x`.

- tol:

  Maximal difference between neighboring objective values that still
  counts as one step. Defaults to 1.

## Value

A `ggplot` object with the plotted data frame as attribute `"data"`.

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
