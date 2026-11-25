# Import an SBML Model

Reads an SBML Level 3 file. Requires rjson and reticulate, which
provides a Python with `python-libsbml`; set the environment variable
`DMOD_LIBSBML_PYTHON` to the path of a Python interpreter to use that
one instead. Kinetic laws are divided by the volume of their
compartment.

## Usage

``` r
importSbml(modelpath, keep = NULL)
```

## Arguments

- modelpath:

  Path to the SBML file.

- keep:

  Parameter ids that keep their own symbol although the SBML sets them
  by an `<initialAssignment>`; a constant assignment becomes their
  default value. Any other initial assignment is substituted into the
  rates, initial values and observables. Default `NULL`.

## Value

A list with components

- `reactions`:

  The model as an
  [eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md).

- `pars`:

  Named numeric, the parameter values.

- `inits`:

  Named character, the initial value of each state.

- `observables`:

  Named list of observable formulas, possibly empty.

- `assignmentRules`:

  Named list of the formulas of the assignment rules and of the
  parameter initial assignments that `keep` does not name, already
  substituted into rates, initial values and observables.

- `events`:

  An
  [eventlist](https://dmodverse.github.io/dMod2/reference/eventlist.md)
  or `NULL`, including the rescaling of states when an event changes a
  compartment volume.

- `eventsSource`:

  An
  [eventlist](https://dmodverse.github.io/dMod2/reference/eventlist.md)
  or `NULL`, the event assignments of the SBML file only.

- `compartmentInits`:

  Named character, the size at the start of each compartment an event
  assigns to; such a compartment keeps the id as volume.

- `renamed`:

  Named character, the new ids of SBML ids that are no valid R or C++
  names, named by the old ids.

## See also

[`exportSbml()`](https://dmodverse.github.io/dMod2/reference/exportSbml.md),
[`importPEtab()`](https://dmodverse.github.io/dMod2/reference/importPEtab.md)

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
