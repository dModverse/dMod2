# Plot Fluxes Given a List of Flux Equations

Evaluates flux expressions along a model prediction and draws them
stacked per condition.

## Usage

``` r
plotFluxes(pouter, x, times, fluxEquations, legendTitle = "Fluxes:", ...)
```

## Arguments

- pouter:

  Named numeric vector of outer parameters.

- x:

  A prediction function, called as
  `x(times, pouter, deriv = FALSE, ...)`.

- times:

  Numeric vector of time points for the prediction.

- fluxEquations:

  Character vector or list of flux expressions in the states and inner
  parameters, e.g. the `rates` of an
  [eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md).
  Names are shown in the legend; without names the expressions are.

- legendTitle:

  Character, the legend title. Defaults to `"Fluxes:"`.

- ...:

  Further arguments passed to `x`, such as `fixed` or `conditions`.

## Value

A `ggplot` object with the flux data frame as attribute `"out"`.

## See also

[`subset.eqnlist()`](https://dmodverse.github.io/dMod2/reference/subset.eqnlist.md)
to select reactions.

## Examples

``` r
times <- 0:5
grid <- data.frame(name = "A", time = times, row.names = paste0("A", times))
x <- Xd(grid)
pars <- structure(exp(-times / 2), names = getParameters(x))
plotFluxes(pars, x, seq(0, 5, 0.1),
           c(production = "0.2", degradation = "0.5*A"))
```
