# Export a dMod Problem to PEtab

Writes a PEtab v1 or v2 problem to `dir`: the parameter, observable,
condition and measurement tables, for v2 also the experiment table, an
SBML model and a YAML file. Parameter scales, fixed values and
per-condition values are read off the parameter transformation `p`.

## Usage

``` r
exportPEtab(
  data,
  reactions,
  observables,
  p,
  pouter,
  errors = NULL,
  lower = NULL,
  upper = NULL,
  fixed = NULL,
  steadyState = NULL,
  events = NULL,
  parameterScale = NULL,
  observableTransformation = "lin",
  noiseDistribution = "normal",
  modelID = "dMod_export",
  formatVersion = "2.0.0",
  dir,
  overwrite = FALSE
)
```

## Arguments

- data:

  A [datalist](https://dmodverse.github.io/dMod2/reference/datalist.md)
  (or a list of data.frames keyed by condition, each with `name`,
  `time`, `value` columns; or a long-format data.frame that
  [`as.datalist()`](https://dmodverse.github.io/dMod2/reference/datalist.md)
  accepts). Only used for the measurements;
  `attr(data, "condition.grid")` is ignored.

- reactions:

  An [eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md)
  describing the ODE network.

- observables:

  Observable formulas keyed by observableId. Accepts a named character
  vector, an
  [eqnvec](https://dmodverse.github.io/dMod2/reference/eqnvec.md), or an
  observation function produced by
  [`Y()`](https://dmodverse.github.io/dMod2/reference/Y.md) (formulas
  read from `attr(observables, "equations")`).

- p:

  A `parfn` produced by
  [`P()`](https://dmodverse.github.io/dMod2/reference/P.md). Required:
  the symbolic trafo is the source of truth for parameters.tsv,
  conditions.tsv, and `<initialAssignment>` formulas.

- pouter:

  Named numeric vector of estimated outer parameters, on the chosen
  `parameterScale`. Names become `parameterId`s in parameters.tsv.

- errors:

  Noise formulas keyed by observableId. Accepts a named character
  vector,
  [eqnvec](https://dmodverse.github.io/dMod2/reference/eqnvec.md), or a
  Y-built error function. An observable it leaves out takes the column
  `sigma` of the data, or the noise `"1"` without data; a data row of
  such an observable without a sigma is an error. If `NULL`, every
  observable takes the column `sigma` of the data, or has the noise
  `"1"` when the data hold no sigma.

- lower, upper:

  Named numeric vectors of bounds, on the same scale as `pouter`.
  Default `NULL`: five decades around `pouter` for log and log10
  parameters, `-Inf` and `Inf` for linear ones.

- fixed:

  Named numeric vector of parameters that are not estimated, on the
  linear scale. Default `NULL`.

- steadyState:

  The `parfn` of
  [`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md)
  between `p` and the prediction, or `NULL` (default). Every condition
  then starts from the steady state of `reactions`, written as a
  pre-equilibration: one species per conserved total starts at the total
  divided by the coefficient of that species, the other states it sets
  at 0, and the conditions with equal parameters share one
  pre-equilibration condition.

- events:

  An `eventlist` of the model, or `NULL` (default). Events at time 0
  that replace a state set it in the condition table, after the
  pre-equilibration where there is one; their values are evaluated with
  the trafo of each condition.

- parameterScale:

  `NULL` (default) reads the scale off `p`: a parameter entering only as
  `exp(X)` is `"log"`, one entering only as `10^(X)` or `exp10(X)` is
  `"log10"`, anything else `"lin"`. Otherwise a scalar
  `"lin"`/`"log"`/`"log10"` (broadcast to all names in `pouter`) or a
  named character vector keyed by parameterId.

- observableTransformation:

  `"lin"` (default), `"log"` or `"log10"`, a single value or one per
  observableId.

- noiseDistribution:

  `"normal"` (default), `"laplace"` or `"log-normal"`, a single value or
  one per observableId. For v2, `"laplace"` with a log transformation is
  written as `"log-laplace"`.

- modelID:

  SBML model identifier. Default `"dMod_export"`.

- formatVersion:

  `"2.0.0"` (default) or `"1"`.

- dir:

  Output directory, created if missing.

- overwrite:

  Logical, whether existing files in `dir` are overwritten. Default
  `FALSE`.

## Value

Path to the written YAML manifest, invisibly.

## Details

An outer parameter named like a state (e.g. `"A <- exp(A)"`) cannot
share the species id in SBML. It is written as `init_<state>` and the
condition table sets the species to it. An outer parameter named like an
inner parameter keeps the name when the trafo is a pure scale wrap; for
any other mapping it is written as `<name>_outer`.

## Limitations

Only events at time 0 that replace a state are written; any other event
in `events` is an error. Fixed parameters are always written with
`parameterScale = "lin"`.

## See also

[`exportPEtabObject()`](https://dmodverse.github.io/dMod2/reference/exportPEtabObject.md)
for a problem returned by
[`importPEtab()`](https://dmodverse.github.io/dMod2/reference/importPEtab.md),
[`importPEtab()`](https://dmodverse.github.io/dMod2/reference/importPEtab.md);
the guide chapter
[PEtab](https://dmodverse.github.io/dMod2/articles/guide-petab.html).

## Examples

``` r
f <- addReaction(eqnlist(), from = "A", to = "", rate = "k*A")
p <- P(eqnvec(A = "exp(logA)", k = "exp(logk)"), condition = "C1",
       modelname = "petab_export_p", compile = FALSE, outdir = tempdir())
data <- datalist(C1 = data.frame(name = "y", time = 0:4, sigma = 0.1,
                                 value = c(2.1, 1.2, 0.7, 0.5, 0.3)))
yaml <- exportPEtab(data, f, observables = c(y = "A"), p = p,
                    pouter = c(logA = log(2), logk = log(0.5)),
                    dir = file.path(tempdir(), "petab_decay"),
                    formatVersion = "1")
readPEtabTables(yaml)$parameters
#>   parameterId parameterScale lowerBound upperBound nominalValue estimate
#> 1        logA            log      2e-05      2e+05          2.0        1
#> 2        logk            log      5e-06      5e+04          0.5        1
```
