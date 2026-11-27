# 10. PEtab import and export

PEtab is the exchange format for ODE parameter estimation problems: an
SBML model plus tables of measurements, observables, conditions and
parameters, joined by a YAML file. `petab.R` writes the guide model at
the best fit as a PEtab problem and reads it back, as for problems of
the [PEtab benchmark
collection](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab).

## PEtab export

[`exportPEtab()`](https://dmodverse.github.io/dMod2/reference/exportPEtab.md)
writes reactions as SBML, observables, errors, data, and the [parameter
transformation](https://dmodverse.github.io/dMod2/articles/methods-steady-states.html#def-parameter-transformation)
`p` with fixed values and values per condition. PEtab 2, the default,
has no parameter scales: values of parameters entering as
[`exp()`](https://rdrr.io/r/base/Log.html) are written on linear scale,
with a warning. The
[gauge](https://dmodverse.github.io/dMod2/articles/methods-symmetries.html#def-gauge)
of
[Symmetries](https://dmodverse.github.io/dMod2/articles/guide-symmetries.md)
becomes fixed values.

- **State at rest.** `steadyState = pSS` gives a pre-equilibration, a
  [steady
  state](https://dmodverse.github.io/dMod2/articles/methods-steady-states.html#def-steady-state)
  before time 0, shared by both conditions. One species per [conserved
  total](https://dmodverse.github.io/dMod2/articles/methods-steady-states.html#conserved-quantities)
  starts at the total (`NfkIkb`, `Ikk`), all others at 0.
- **Event.** `events = events` writes the
  [event](https://dmodverse.github.io/dMod2/articles/methods-sensitivities.html#def-event)
  at time 0 as the `TNF` value of the condition.
- **Errors.** Imaging readouts get the [error
  model](https://dmodverse.github.io/dMod2/articles/methods-objective.html#def-error-model),
  blots the data standard deviations as noise parameters per
  measurement.

``` r

# the estimated parameters at the best fit, the gauge as fixed values on the linear scale of the
# model
best <- as.parvec(readRDS("results/fits.rds"))
yaml <- exportPEtab(data, reactions, observables = g, p = p, pouter = best, errors = e,
                    fixed = exp(gauge), steadyState = pSS, events = events,
                    modelID = "oppelt", dir = "results/petab", overwrite = TRUE)
list.files("results/petab")
#> [1] "conditions_oppelt.tsv"   "experiments_oppelt.tsv" 
#> [3] "measurements_oppelt.tsv" "observables_oppelt.tsv" 
#> [5] "oppelt.xml"              "oppelt.yaml"            
#> [7] "parameters_oppelt.tsv"
```

The export uses the best fit of `fit.R`, loaded from `results/fits.rds`.
Each experiment starts with the pre-equilibration at \\t = -\infty\\,
then the experiment condition at 0:

``` r

read.delim("results/petab/experiments_oppelt.tsv")
#>                experimentId time conditionId
#> 1     preequilibration__TNF -Inf            
#> 2     preequilibration__TNF    0         TNF
#> 3 preequilibration__control -Inf            
#> 4 preequilibration__control    0     control
read.delim("results/petab/conditions_oppelt.tsv")
#>   conditionId targetId targetValue
#> 1         TNF      TNF           1
#> 2     control      TNF           0
```

## PEtab import

[`importPEtab()`](https://dmodverse.github.io/dMod2/reference/importPEtab.md)
builds `prdfn`, the
[`obsfn`s](https://dmodverse.github.io/dMod2/articles/methods-objective.html#def-obsfn)
of observables and errors, the `parfn`, data and `objfn`; a
pre-equilibration becomes a
[`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md). With
`compile = FALSE`,
[`compile()`](https://dmodverse.github.io/dMod2/reference/compile.md)
builds the sources.

``` r

petab <- importPEtab(yaml, backend = "cppDE", options = list(atol = 1e-10, rtol = 1e-8),
                     modelname = "oppelt_petab", compile = FALSE)
```

``` r

compile(petab$prd, petab$e, output = "oppelt_petab", cores = 6)
```

``` r

# the data term of model.R and the objective of the problem read back
c(model.R = attr(obj(c(best, gauge), deriv = FALSE), "data"),
  PEtab   = petab$obj(petab$bestfit, deriv = FALSE)$value)
#>   model.R     PEtab 
#> -511.1689 -511.1690
```

The imported problem reproduces the [data
term](https://dmodverse.github.io/dMod2/articles/methods-objective.html#data-term)
of `model.R` at the best fit. `petab$bestfit` holds the nominal values
on the problem scales, `petab$parlower` and `petab$parupper` the bounds.

## SBML models

[`exportSbml()`](https://dmodverse.github.io/dMod2/reference/exportSbml.md)
writes species with compartments, reactions, parameters and initial
assignments; concentration rates are multiplied by the [compartment
volume](https://dmodverse.github.io/dMod2/articles/methods-steady-states.html#def-compartment-volume),
[`importSbml()`](https://dmodverse.github.io/dMod2/reference/importSbml.md)
divides again.
[`importSbml()`](https://dmodverse.github.io/dMod2/reference/importSbml.md)
returns the equation list in `reactions`, parameter and initial values,
observables, assignment rules and events:

``` r

sbml <- importSbml("results/petab/oppelt.xml")
names(sbml)
#> [1] "reactions"        "pars"             "inits"            "observables"     
#> [5] "assignmentRules"  "events"           "eventsSource"     "compartmentInits"
#> [9] "renamed"
c(states = length(sbml$reactions$states), reactions = length(sbml$reactions$rates))
#>    states reactions 
#>        23        38
identical(sort(sbml$reactions$states), sort(reactions$states))
#> [1] TRUE
```

`exportSbml(reactions, parameters, inits, filepath)` writes an equation
list without PEtab tables.

## Formats and limits

- SBML and PEtab run through reticulate with `python-libsbml`;
  `DMOD_LIBSBML_PYTHON` selects the interpreter.
- PEtab 1 and 2 are read. `exportPEtab(formatVersion = "1")` keeps
  scales but writes condition expressions that only dMod2 reads.
- [`exportPEtab()`](https://dmodverse.github.io/dMod2/reference/exportPEtab.md)
  writes events at time 0 that replace a state and refuses others.

Reference:
[`exportPEtab()`](https://dmodverse.github.io/dMod2/reference/exportPEtab.md),
[`importPEtab()`](https://dmodverse.github.io/dMod2/reference/importPEtab.md),
[`exportSbml()`](https://dmodverse.github.io/dMod2/reference/exportSbml.md),
[`importSbml()`](https://dmodverse.github.io/dMod2/reference/importSbml.md)

[← 9. Alignment of relative data
measurements](https://dmodverse.github.io/dMod2/articles/guide-alignment.md)
\| [11. Remote computing
→](https://dmodverse.github.io/dMod2/articles/guide-remote-computing.md)
