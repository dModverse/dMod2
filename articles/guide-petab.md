# 10. PEtab import and export

PEtab is the exchange format for ODE parameter estimation problems: an
SBML model plus tables of measurements, observables with noise models,
conditions, and parameters with scales and bounds, joined by a YAML
file. This chapter writes the guide model at the best fit as a PEtab
problem and reads it back (`petab.R`); problems of the [PEtab benchmark
collection](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab)
and other tools are read the same way.

## PEtab export

[`exportPEtab()`](https://dmodverse.github.io/dMod2/reference/exportPEtab.md)
writes reactions as SBML, observables, errors, data, and the
transformation `p` with scales, fixed values and values per condition. A
parameter entering only as [`exp()`](https://rdrr.io/r/base/Log.html) is
on log scale. PEtab 2, the default, has no parameter scales: values are
written on linear scale, with a warning listing the log-scale
parameters. The gauge of
[Symmetries](https://dmodverse.github.io/dMod2/articles/guide-symmetries.md)
becomes fixed values, also linear.

- **State at rest.** With `steadyState = pSS`, each condition starts
  with a pre-equilibration to steady state before time 0. The totals of
  the conserved moieties (the
  [`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md)
  parameters) become initial-value parameters: one species per total
  starts at the total, here complex `NfkIkb` for NFκB and
  unphosphorylated `Ikk` for IKK, all others at 0; the simulation
  preserves the totals. Both conditions differ only in the dose, given
  by the event, and share one pre-equilibration condition.
- **Event.** With `events = events`, the event at time 0 becomes the
  `TNF` value in the condition table, set after pre-equilibration.
- **Errors.** Imaging readouts get their error model; the blots, outside
  the error model, get the data standard deviations as noise parameters
  per measurement, as in
  [`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md).

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
builds `prdfn`, the `obsfn`s of observables and errors, the `parfn` with
the scales of the parameter table, data and `objfn`. A pre-equilibration
becomes a
[`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md) with
totals from the initial values. `options` go to the solver; with
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
#> -511.1623 -511.1623
```

The imported problem reproduces the data term of `model.R` at the best
fit. `petab$bestfit` has the nominal values on the problem scales,
`petab$parlower` and `petab$parupper` the bounds; it is fitted and
profiled like any other.

## SBML models

[`exportSbml()`](https://dmodverse.github.io/dMod2/reference/exportSbml.md)
writes species with compartments and volumes, reactions with kinetic
laws, parameters and initial assignments, here one species per total set
to that total. SBML kinetic laws are amounts per time:
[`exportSbml()`](https://dmodverse.github.io/dMod2/reference/exportSbml.md)
multiplies concentration rates by the compartment volume,
[`importSbml()`](https://dmodverse.github.io/dMod2/reference/importSbml.md)
divides again.
[`importSbml()`](https://dmodverse.github.io/dMod2/reference/importSbml.md)
returns an equation list with compartments, volumes, parameter values,
initial values, assignment rules and events:

``` r

sbml <- importSbml("results/petab/oppelt.xml")
names(sbml)
#> [1] "reactions"       "pars"            "inits"           "observables"    
#> [5] "assignmentRules" "events"          "eventsSource"    "renamed"
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
- PEtab 1 and 2 are read.
  [`exportPEtab()`](https://dmodverse.github.io/dMod2/reference/exportPEtab.md)
  writes 2 by default, 1 with `formatVersion = "1"`, which keeps scales
  but puts expressions into the condition table that only dMod2 reads.
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
