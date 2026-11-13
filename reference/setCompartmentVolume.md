# Set Compartment Volumes

Sets the volume expression, and optionally the volume rule, of one or
more compartments and updates the per-state volumes.

## Usage

``` r
setCompartmentVolume(eqnlist, ..., rules = NULL)
```

## Arguments

- eqnlist:

  object of class
  [eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md)

- ...:

  named arguments `compartment = "volume expression"`, or a single named
  character vector of the same shape. Compartments that do not exist yet
  are created.

- rules:

  Optional named character vector `compartment = "dV/dt"`, adding a
  volume rule. `NA` removes an existing rule.

## Value

An object of class
[eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md).

## See also

[`assignCompartment()`](https://dmodverse.github.io/dMod2/reference/assignCompartment.md)

## Examples

``` r
f <- eqnlist() |>
  addReaction("A", "B", "k*A", compartment = "cyt") |>
  setCompartmentVolume(cyt = "V_cyt")
f$volumes
#>       A       B 
#> "V_cyt" "V_cyt" 
```
