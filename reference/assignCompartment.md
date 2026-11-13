# Assign States to Compartments

Declares the compartment of a state, before or after the state is
introduced by
[`addReaction()`](https://dmodverse.github.io/dMod2/reference/addReaction.md).

## Usage

``` r
assignCompartment(eqnlist, ..., volume = NULL, rule = NULL)
```

## Arguments

- eqnlist:

  object of class
  [eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md)

- ...:

  named arguments `state = "compartment"`, or a single named character
  vector of the same shape. Compartments that do not exist yet are
  created with volume `"1"`.

- volume:

  Optional volume expression for the target compartment. Only allowed
  when `...` names a single compartment; use
  [`setCompartmentVolume()`](https://dmodverse.github.io/dMod2/reference/setCompartmentVolume.md)
  for several.

- rule:

  optional volume rule (`dV/dt`) for the target compartment, see
  [`setCompartmentVolume()`](https://dmodverse.github.io/dMod2/reference/setCompartmentVolume.md).
  Default `NULL`.

## Value

An object of class
[eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md).

## See also

[`setCompartmentVolume()`](https://dmodverse.github.io/dMod2/reference/setCompartmentVolume.md),
[`addReaction()`](https://dmodverse.github.io/dMod2/reference/addReaction.md)

## Examples

``` r
f <- eqnlist() |>
  assignCompartment(TGFb = "extraCell", volume = "V_ext") |>
  addReaction("R + TGFb", "R_TGFb", "k_on*R*TGFb", compartment = "Cell")
f$compartmentOf
#>           R        TGFb      R_TGFb 
#>      "Cell" "extraCell"      "Cell" 
```
