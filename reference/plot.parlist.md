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

  object of class `parlist` from
  [`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md).
  Both plots need the starting point `parinit` that
  [`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md)
  adds to each fit; a list of
  [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md)
  results has none.

- path:

  logical, plot the path of the parameters from the initial values to
  convergence. Default `FALSE`. `path = TRUE` needs fits run with
  `mstrust(..., blather = TRUE)`.

- ...:

  not used.

## Value

A `ggplot` object.

## See also

[`parlist()`](https://dmodverse.github.io/dMod2/reference/parlist.md),
[`plotValues()`](https://dmodverse.github.io/dMod2/reference/plotValues.md)
