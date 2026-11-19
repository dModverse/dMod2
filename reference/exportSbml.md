# Export an Equation List to SBML

Writes an
[eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md) with
parameter values and initial values to an SBML Level 3 Version 2 file.
Kinetic laws are multiplied by the volume of their compartment. Symbolic
volumes are written as an `<initialAssignment>` of the compartment.
Requires the same Python environment as
[`importSbml()`](https://dmodverse.github.io/dMod2/reference/importSbml.md).

## Usage

``` r
exportSbml(
  eqnlist,
  parameters = NULL,
  inits = NULL,
  filepath,
  modelID = "dMod_export",
  events = NULL
)
```

## Arguments

- eqnlist:

  Object of class
  [eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md) to
  export.

- parameters:

  Named numeric vector of parameter values (including any
  compartment-size parameters referenced in `eqnlist$compartments`). An
  entry named like a compartment whose volume is the compartment symbol
  is the size of that compartment; entries named like other compartments
  or like species are not written. Default `NULL`: no `<parameter>` is
  written.

- inits:

  Named numeric or character vector of initial values keyed by state
  name. Numbers, also as character, are written as
  `initialConcentration`, or `initialAmount` for the states in
  `eqnlist$amountStates`; other expressions as `<initialAssignment>`
  formulas. Missing states default to 0. Default `NULL`.

- filepath:

  Path of the SBML file to write.

- modelID:

  SBML model identifier. Defaults to `"dMod_export"`.

- events:

  Object of class
  [eventlist](https://dmodverse.github.io/dMod2/reference/eventlist.md)
  to emit as `<listOfEvents>`, or `NULL`. Rows sharing a trigger become
  one `<event>`. A compartment an event assigns to is written as
  non-constant. Default `NULL`.

## Value

`filepath`, invisibly.

## See also

[`importSbml()`](https://dmodverse.github.io/dMod2/reference/importSbml.md),
[`exportPEtab()`](https://dmodverse.github.io/dMod2/reference/exportPEtab.md)

## Examples

``` r
f <- addReaction(eqnlist(), from = "A", to = "B", rate = "k1*A")
f <- addReaction(f, from = "B", to = "", rate = "k2*B")
file <- file.path(tempdir(), "model.xml")
exportSbml(f, parameters = c(k1 = 0.5, k2 = 0.1), inits = c(A = 1, B = 0),
           filepath = file)
model <- importSbml(file)
model$reactions
#> 
#>   Check Educt -> Product       Rate Description
#> 1           A ->       B 1 * k1 * A           1
#> 2           B ->         1 * k2 * B           2
model$pars
#>          k1          k2 defaultComp 
#>         0.5         0.1         1.0 
```
