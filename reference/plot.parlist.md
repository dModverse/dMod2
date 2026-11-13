# Plot a Parameter List

Plots, for the converged fits, the converged against the initial value
of each parameter, or with `path = TRUE` the parameter values over the
iterations.

## Usage

``` r
# S3 method for class 'parlist'
plot(x, path = FALSE, ...)
```

## Arguments

- x:

  object of class `parlist`, e.g. from
  [`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md).

- path:

  logical, plot the path of the parameters from the initial values to
  convergence. Default `FALSE`. `path = TRUE` needs fits run with
  `trust(..., blather = TRUE)` or `mstrust(..., blather = TRUE)`.

- ...:

  not used.

## Value

A `ggplot` object.

## See also

[`parlist()`](https://dmodverse.github.io/dMod2/reference/parlist.md),
[`plotValues()`](https://dmodverse.github.io/dMod2/reference/plotValues.md)
