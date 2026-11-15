# Migrating from dMod

What a dMod [\[1\]](#ref-kaschek2019dmod) script has to change for
dMod2.

## Unchanged interface

Function objects and algebra are the same:
[`Xs()`](https://dmodverse.github.io/dMod2/reference/Xs.md),
[`Xf()`](https://dmodverse.github.io/dMod2/reference/Xf.md),
[`Y()`](https://dmodverse.github.io/dMod2/reference/Y.md),
[`P()`](https://dmodverse.github.io/dMod2/reference/P.md),
[`Pexpl()`](https://dmodverse.github.io/dMod2/reference/Pexpl.md),
[`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md),
[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md),
[`constraintL2()`](https://dmodverse.github.io/dMod2/reference/constraintL2.md),
[`datapointL2()`](https://dmodverse.github.io/dMod2/reference/datapointL2.md),
with `*` for concatenation, `+` for conditions or objective terms, and
`normL2(data, g * x * p)`.
[`eqnlist()`](https://dmodverse.github.io/dMod2/reference/eqnlist.md),
[`addReaction()`](https://dmodverse.github.io/dMod2/reference/addReaction.md),
[`define()`](https://dmodverse.github.io/dMod2/reference/define.md),
[`insert()`](https://dmodverse.github.io/dMod2/reference/define.md),
[`branch()`](https://dmodverse.github.io/dMod2/reference/define.md),
[`repar()`](https://dmodverse.github.io/dMod2/reference/repar.md),
[`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md),
[`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md),
[`profile()`](https://rdrr.io/r/stats/profile.html) and `plot*()` keep
names and main arguments. dMod2 has
[`getSymbols()`](https://dmodverse.github.io/dMod2/reference/getSymbols.md)
and
[`replaceSymbols()`](https://dmodverse.github.io/dMod2/reference/replaceSymbols.md),
so [`library(cOde)`](https://rdrr.io/r/base/library.html) is not needed.

## Substitution syntax

[`define()`](https://dmodverse.github.io/dMod2/reference/define.md),
[`insert()`](https://dmodverse.github.io/dMod2/reference/define.md),
[`repar()`](https://dmodverse.github.io/dMod2/reference/repar.md) and
`branch(apply = )` take `"lhs <- rhs"`; `"lhs ~ rhs"` is an error naming
the new form.

| dMod | dMod2 |
|----|----|
| `define(trafo, "x ~ x", x = pars)` | `define(trafo, "x <- x", x = pars)` |
| `insert(trafo, "x ~ exp(x)", x = pars)` | `insert(trafo, "x <- exp(x)", x = pars)` |
| `insert(trafo, "x ~ x*x_OE", dose == 10, ...)` | `insert(trafo, "x <- x*x_OE", dose == 10, ...)` |

Placeholders in `...`, branch table columns in scope, unnamed logical
condition selectors and `conditionMatch` are unchanged.

`~` now writes effect formulas: on a branched transformation,
[`insert()`](https://dmodverse.github.io/dMod2/reference/define.md)
takes an unquoted formula and replaces the symbol by a linear predictor
over the branch table columns, fixed and random effects in lme4
notation:

``` r

trafo <- branch(trafo, table = covariates)        # one row per condition, columns ID and wt
trafo <- insert(trafo, "x <- exp(x)", x = c("ka", "cl"))
trafo <- insert(trafo, x ~ 1 + (1 | ID), x = c("ka", "cl"))
trafo <- insert(trafo, v ~ wt + (1 | gr(ID, dist = spikeSlab())))
randomEffects(trafo)
```

See [Sparse
deviations](https://dmodverse.github.io/dMod2/articles/guide-sparse-deviations.md)
and [Nonlinear mixed-effects
models](https://dmodverse.github.io/dMod2/articles/guide-nlme.md).

## Backend and compilation

The default `odemodel(backend = "cppDE")` generates C++ code through
[cppDE](https://dmodverse.github.io/cppDE/); `"Sundials"` runs CVODES
through
[`cppDE::cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.html),
`"deSolve"` keeps the dMod path with code compiled by cOde.

Every function object is compiled. With `compile = FALSE` only sources
are written, and one
[`compile()`](https://dmodverse.github.io/dMod2/reference/compile.md)
builds them into one shared object:

**dMod**

``` r

x <- Xs(odemodel(r, modelname = "x"))
g <- Y(obs, r, attach.input = FALSE)
p <- P(trafo, condition = "C1")
```

**dMod2**

``` r

x <- Xs(odemodel(r, modelname = "x", compile = FALSE))
g <- Y(obs, x, modelname = "g", compile = FALSE)
p <- P(trafo, condition = "C1", modelname = "p", compile = FALSE)
compile(x, g, p, output = "model")
```

Files go to `outdir`, default `getOption("dMod2.outdir")`, else the
working directory.

| dMod | dMod2 |
|----|----|
| `odemodel(solver = )` | `odemodel(backend = c("cppDE", "Sundials", "deSolve"))` |
| `odemodel(outputs = , estimate = , gridpoints = )` | `...`, deSolve only |
| `Xs(optionsOde = , optionsSens = )` | `Xs(options = , optionsSens = )`, see [`cppDE::solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.html) |
| `Xs(fcontrol = )` | deSolve only |
| none | `deriv2 = TRUE` in [`odemodel()`](https://dmodverse.github.io/dMod2/reference/odemodel.md), [`Y()`](https://dmodverse.github.io/dMod2/reference/Y.md), [`P()`](https://dmodverse.github.io/dMod2/reference/P.md): exact second derivatives |

## Renamed arguments

Dotted argument names are lowerCamelCase.

| dMod | dMod2 | in |
|----|----|----|
| `attach.input` | `attachInput` | [`Y()`](https://dmodverse.github.io/dMod2/reference/Y.md), [`Pexpl()`](https://dmodverse.github.io/dMod2/reference/Pexpl.md) |
| `keep.root` | `keepRoot` | [`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md) |
| `attr.name` | `attrName` | [`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md), [`constraintL2()`](https://dmodverse.github.io/dMod2/reference/constraintL2.md), [`datapointL2()`](https://dmodverse.github.io/dMod2/reference/datapointL2.md) |
| `split.by`, `keep.covariates` | `splitBy`, `keepCovariates` | [`as.datalist()`](https://dmodverse.github.io/dMod2/reference/datalist.md) |
| `sort.by` | `sortBy` | [`as.parframe()`](https://dmodverse.github.io/dMod2/reference/as.parframe.md) of a parlist |
| `obj.attributes` | `objAttributes` | [`parframe()`](https://dmodverse.github.io/dMod2/reference/parframe.md) |
| `studyname` | `name`, default `"mstrust"` | [`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md) |
| `value` | `parameter` | [`datapointL2()`](https://dmodverse.github.io/dMod2/reference/datapointL2.md) |
| `nameFlux` | `legendTitle` | [`plotFluxes()`](https://dmodverse.github.io/dMod2/reference/plotFluxes.md) |
| `atr` | `which` | [`attrs()`](https://dmodverse.github.io/dMod2/reference/attrs.md) |
| `S` (stoichiometric matrix) | `eqnlist` | [`conservedQuantities()`](https://dmodverse.github.io/dMod2/reference/conservedQuantities.md) |

## Optimization

[`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md) drops
`parscale`, `fterm`, `mterm`; thresholds go to `tolControl`, step
settings to `stepControl`.

| dMod | dMod2 |
|----|----|
| `trust(..., fterm = , mterm = )` | `trust(..., tolControl = list(ftol = , mtol = , gtol = ))` |
| `trust(parscale = )` | removed |
| Gauss-Newton Hessian only | `hessianMethod = c("gn", "bfgs", "sr1", "exact")` |
| `mstrust(..., studyname = "fits")` | `mstrust(..., name = "fits", resultPath = ".")` |
| `mstrust(cores = )` | `cores = n`, or `c(fits = , conditions = )` |

## Profile likelihood

`profile(obj, pars, whichPar)` keeps `method`, `stepControl`,
`algoControl` and defaults. Changes:

| dMod | dMod2 |
|----|----|
| `optControl = list(fterm = , mterm = )` | any [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md) argument, e.g. `tolControl` |
| `stepControl = list(stop = "data")` | `stop`: also any numeric attribute of the objective |
| threshold `qchisq(1 - alpha, 1)` | `delta =`, e.g. `profileThreshold(method = "F", n = , p = )` |
| both sides | `side = c("both", "left", "right")` |
| `algoControl = list(W = "hessian")`: Gauss-Newton | exact Hessian with `deriv2 = TRUE` |
| `cores =` over parameters | `cores = n`, or `c(pars = , conditions = )` |
| none | attribute `"ends"`: steps and stop reason per side |
| none | `sparse = "auto"`: sparse Hessian where available |

[`confint()`](https://rdrr.io/r/stats/confint.html) uses the data term
by default, `val.column = "value"` the whole objective; see [Profile
likelihood](https://dmodverse.github.io/dMod2/articles/methods-profiles.md)
for the method. Stop reasons in `"ends"` are threshold, limit, step
limit or error.

## Steady states

[`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md) solves
by pseudo-transient continuation in C++ instead of
`rootSolve::multiroot()`. On an eqnlist, `flow = TRUE` (default) returns
the stable steady state reached by the flow; settings such as `positive`
go to `controlsPTC`. `P(eqnlist)` without `method` is implicit,
[`steadyStates()`](https://dmodverse.github.io/dMod2/reference/steadyStates.md)
derives an analytic transformation.

## Removed functions

| dMod | dMod2 |
|----|----|
| `priorL2(mu, lambda)` | `constraintL2(mu, sigma)` |
| `nll()` | `normL2(optBLOQ = c("M3", "M1", "M4NM", "M4BEAL"))` |
| `mname()` | [`modelname()`](https://dmodverse.github.io/dMod2/reference/modelname.md) |
| `funC0()` | [`cppDE::cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.html) |
| `ginv()` | [`MASS::ginv()`](https://rdrr.io/pkg/MASS/man/ginv.html) |
| `wide2long(prediction)` | `as.data.frame(prediction)` |

## Changed defaults

These run without error but give results different from dMod.

- [`odemodel()`](https://dmodverse.github.io/dMod2/reference/odemodel.md):
  cppDE (BDF) instead of deSolve `lsoda`; trajectories differ within
  tolerances.
- [`Y()`](https://dmodverse.github.io/dMod2/reference/Y.md):
  `attachInput = FALSE`, `g * x` returns observables only.
- [`P()`](https://dmodverse.github.io/dMod2/reference/P.md): an eqnlist
  as `trafo` gives an implicit transformation.
- [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md):
  `iterlim = 1000` instead of `100`, `tolControl` with
  `ftol = mtol = gtol = 1e-6`;
  [`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md)
  passes both on.
- [`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md):
  `retry = TRUE` restarts a fit that errors, up to `nTries = 10`
  attempts.
- [`compile()`](https://dmodverse.github.io/dMod2/reference/compile.md):
  `cores = detectFreeCores()` instead of `1`.
- [`confint()`](https://rdrr.io/r/stats/confint.html):
  `val.column = "data"` instead of `"value"`, a prior no longer narrows
  the interval.
- [`fitErrorModel()`](https://dmodverse.github.io/dMod2/reference/fitErrorModel.md):
  `plotting = FALSE`.

## Further reading

- [Get started](https://dmodverse.github.io/dMod2/articles/dMod2.md):
  one small model end to end.
- [Guide](https://dmodverse.github.io/dMod2/articles/guide-model.md): a
  published model from setup to
  [optimization](https://dmodverse.github.io/dMod2/articles/guide-optimization.md)
  and [profile
  likelihood](https://dmodverse.github.io/dMod2/articles/guide-profiles.md).
- [Installation and
  troubleshooting](https://dmodverse.github.io/dMod2/articles/installation.md):
  compilers, Python, SBML.
- [cppDE](https://dmodverse.github.io/cppDE/): solver options and the
  expression language.

## References

\[1\]

D. Kaschek, W. Mader, M. Fehling-Kaschek, M. Rosenblatt, and J. Timmer,
“dMod: Efficient parameter estimation and uncertainty analysis for
dynamic models,” *Journal of Statistical Software*, vol. 88, no. 10, pp.
1–32, 2019, doi:
[10.18637/jss.v088.i10](https://doi.org/10.18637/jss.v088.i10).
