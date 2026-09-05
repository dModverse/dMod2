# dMod2 0.7.1

* `trust()` takes an interchangeable Hessian source through `hessianMethod`.
  `"gn"` (default) is the Gauss-Newton `J^T J` as before; `"bfgs"` and `"sr1"`
  maintain a dense quasi-Newton update seeded from it; `"hybrid"` runs `"gn"`
  until it stagnates, then switches to `"bfgs"` once. The quasi-Newton phase
  consumes only the gradient. Reflective boundary only.
* Objective functions take a call-time `hessian` argument (default `TRUE`).
  With `hessian = FALSE` they return value and gradient but skip the Hessian
  entirely -- the `J^T J` contraction never runs and the result carries a `NULL`
  hessian. `trust()` uses this in the quasi-Newton phase; it propagates through
  objective composition (`+`).
* `trust()` reports `neval` and, under `blather`, the `hessianSource` per
  iteration, so a multi-start can be scored on gradient evaluations. These reach
  `as.parframe()` as columns.
* `compile()` builds model shared objects correctly on Windows. The per-source
  link and the combined-output object compile ran `system()` with a `2>&1`
  token, but R's `system()` on Windows has no shell: the token reached
  `R CMD SHLIB` as the override `PKG_LIBS=2>&1` and the compiler as an input
  file, so the object never built and the `.dll` never linked against BLAS,
  LAPACK or the Sundials solver. Both now go through `system2()`.
* Adds `inst/examples/example_Boehm_JProteomeRes2014.R`, which builds the Boehm
  et al. (2014) STAT5 model and compares the Hessian sources over a multi-start.

# dMod2 0.7.0

* PEtab import and export. `importPEtab()` reads a v1 or v2 problem and returns
  the composed `g * x * p`, the error model, the objective with its priors and
  the published parameter vector. `exportPEtabObject()` writes one back, and
  `exportPEtab()` builds a problem from a hand-written dMod model.
* SBML import and export through libsbml, reached over reticulate, so a problem
  that needs neither stays Python-free. Identifiers that R or C++ reserve are
  renamed throughout, including in the PEtab tables, and the MathML constants
  `pi`, `exponentiale` and `avogadro` are folded to their value.
* `normL2()` returns gradient and Hessian in the order of the parameter vector
  it was called with. They used to follow the union of the per-condition
  sensitivity blocks, while `trust()` and `optim()` read them positionally, so
  every fit optimised a permuted model.
* `normL2()` reports its sum of squares as a `chi2` attribute. Adding
  objectives pools the terms sharing an `attr.name` and splits the rest into
  `chi2_<attr.name>`.
* `compile(output = )` points the cOde models inside a function object at the
  batched shared object. Without that the deSolve backend looked for an entry
  point in a library that was never built.
* Pre-equilibration takes its conserved moieties from the same reduced network
  `Pequil()` uses, freezes a time dependent input at the start time, and leaves
  symbols no outer parameter reaches out of the sensitivity system. A failure
  now reports the solver error and the non-finite parameters.
* `Y()` accepts one observation function per condition and takes `cores`.
  `importPEtab()` takes `cores`, `deriv`, `outdir`, `optionsOde` and
  `optionsSens`, and derives a model name that steps aside for what is already
  loaded.
* `as.data.frame()` on a prediction joins the error model by observable rather
  than by row order, which used to hand each observable the wrong sigma.
* `wide2long()` on an unnamed list numbers its conditions instead of dropping
  the `condition` column, and forwards `keep` and `na.rm`.
* `print()` on an objective result from a `deriv = FALSE` call shows the value
  alone, where it used to fail on the absent Hessian.
* New data set `bachmann` with `inst/examples/example_BachmannMSB2011.R`
  building the model: 25 states, 36 conditions, 113 estimated parameters.
* All 31 cases of the PEtab v2 test suite reproduce the published likelihood
  and survive the export round trip. Of the 35 problems in the Benchmark-Models
  collection, 33 import, evaluate with derivatives, export and reimport.

# dMod2 0.6.5

* `compile(output = )` no longer links into a shared object the process
  already holds. Overwriting one is not portable: Windows may keep the file
  handle and macOS may keep the image resident, so the reload would serve the
  old code. A loaded name is replaced by `<name>_2` with a warning, the rule
  cppDE's model constructors already apply to their own names.

# dMod2 0.6.4

* New soft constraints next to `constraintL2()`: `constraintL1()` (Laplace),
  `constraintCauchy()`, `constraintGamma()`, `constraintExponential()`,
  `constraintChisq()` and `constraintRayleigh()`. Each is the `-2 log` density
  including its normalisation, unlike `constraintL2()`, which is the penalty
  form, and each carries the chain rule so `constraint * P()` is exact.
* The multivariate-normal path of `constraintL2()` is gone, together with the
  `penaltySpec` plumbing and the `plotIndivs()` / `plotHistIndivs()` generics.
  All of them need an `omegaSpec` or a `penaltyspec`, which only the NLME layer
  builds, so here they were unreachable.
* `constraintExp2()`, an unused box prior inherited from dMod 1.x, is removed.

# dMod2 0.6.3

* `normL2()` takes `t0`, the time at which initial values take effect. The time
  grid starts there rather than at 0, so a prediction may begin after 0.
* The error model keeps its fixed parameters under `deriv = FALSE`. They were
  derived from the sensitivity rows, which are absent then.
* `Xt()` reports its fixed parameters, as `Xs()` does, so an error model can
  read them off the prediction.
* `eqnlist()` accepts a model without species, which is a parameter-only
  problem rather than an error.
* The steady-state heuristic evaluates SBML's `piecewise()` and declares a
  state structurally zero only when its rate vanishes at zero alone.

# dMod2 0.6.2

* Reloading a model's shared object no longer invalidates anything. cppDE
  (>= 0.9.2) resolves entry points by name, so the batch handle cached by
  `Xs()`, and every objective built on it, survives an unload and reload. This
  replaces the symbol-cache flush counter of 0.6.1, which only caught reloads
  that went through `loadDLL()` or `compile()`.
* `loadDLL()` also searches the directories the sources were generated in, so a
  model compiled into a temporary folder is found rather than silently skipped.
* `normL2()` carries `compileInfo`, as `*` and `+` already did, so `loadDLL()`
  and `compile()` reach the shared objects of a composed objective.

# dMod2 0.6.1

* `loadDLL()` skips shared objects that are already loaded in the current
  process. Unloading them nulled the native symbol pointers held by the
  prediction, observation and parameter functions built from them, leaving
  every later call without a way to resolve them again.
* `Xs()` caches a prepared batch handle that carries such a pointer of its
  own, which no symbol-cache flush reached. It is now rebuilt whenever the
  cache is flushed.

# dMod2 0.6.0

* The core is complete: equations, compiled prediction, observation and
  parameter transformations, objectives, trust region optimisation and
  profile likelihood. The PEtab, symmetry, mixed effects and Bayesian layers
  live on their own development branches.
* Test suite over equations, compilation, prediction, objectives, the trust
  region and profiles.

# dMod2 0.5.21

* `plotProfile()`, `plotPaths()` and `plotValues()` build on one long-format
  frame; the `parlist` overlay is measured against the profile optimum.
* dMod palettes, `theme_dMod()` and the colour scales.

# dMod2 0.5.20

* `runbg()` runs jobs in the background on a local or remote host.
* `distributedComputing()` spreads fits and profiles over a cluster.

# dMod2 0.5.19

* `reml()` estimates the error model from the restricted likelihood, charging
  every data point its own leverage.
* `remlLeverage()` reports the hat values, their rank and the effective
  degrees of freedom per observable.

# dMod2 0.5.18

* `mstrust(cores = c(fits = , conditions = ))` and
  `profile(cores = c(pars = , conditions = ))` split the two parallel axes.
* `profileThreshold()` computes the chi-square or finite-sample F threshold
  that `profile()`, `plotProfile()` and `confint()` share.

# dMod2 0.5.17

* `steadyStates()` solves for symbolic steady states through sympy.

# dMod2 0.5.16

* `Pexpl()`, `Pimpl()` and `Pequil()` build explicit, implicit and
  equilibrated parameter transformations.
* `define()`, `insert()`, `branch()` and `repar()` assemble transformations.
* Warm starts are kept per condition, so a condition restarts from its own
  previous root.

# dMod2 0.5.15

* `Xs()` hands every condition to `cppDE::solveODEBatch()` in one call.
* `Y()` and `P(method = "explicit")` evaluate all conditions in one call.
* `cores` is a call-time argument on every function object and defaults to
  `getOption("dMod.cores", 1L)`.

# dMod2 0.5.14

* Fixed: composing a condition-less function with a multi-condition one kept
  only the first condition and returned an unnamed result.

# dMod2 0.5.13

* The loop over experimental conditions moved from outside a composed chain
  into its leaves, so a leaf sees every condition at once and can batch them.
* Predictions and objectives are identical across thread counts.

# dMod2 0.5.12

* Constraint and parameter-vector kernels in C++.

# dMod2 0.5.11

* `trust()` and `trustL1()` run on a C++ trust region kernel.

# dMod2 0.5.10

* `normL2()` evaluates its data term in a C++ kernel.

# dMod2 0.5.9

* Residual kernel in C++, with error models and BLOQ handling (`M1`, `M3`,
  `M4NM`, `M4BEAL`).

# dMod2 0.5.8

* `compile()` links model sources into one shared object, with compile and
  link flags taken per file from the backend that produced it.
* Objects are reused when source and command are unchanged, recorded in a
  `.dMod_objects` index.

# dMod2 0.5.7

* `odemodel()` generates C++ with first and second order sensitivities.
  Backends `cppDE`, `Sundials` and `deSolve`.

# dMod2 0.5.6

* An `eventlist` is declared on `odemodel()`, so the sensitivity equations are
  extended consistently.

# dMod2 0.5.5

* `eqnlist` carries compartments and volumes; `assignCompartment()` and
  `setCompartmentVolume()` set them.
* `conservedQuantities()` and `getTotals()` report the conserved moieties.

# dMod2 0.5.4

* Package entry point renamed to `dMod2`.

# dMod2 0.5.3

* Symmetry detection and SBML import moved out of the core.

# dMod2 0.5.2

* Removed `dModFrame`, `dfoptim` and `normIndiv`.

# dMod2 0.5.1

* The per-author `tools*.R` files are merged into `utils.R` and `accessors.R`.

# dMod2 0.5.0

* Fork of dMod 0.5.
