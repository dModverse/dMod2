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
