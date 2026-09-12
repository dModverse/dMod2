# dMod2 (development version)

* A cotangent carries directions. Every backward node now hands on a matrix or
  an array with a trailing direction axis: slice one is the cotangent, the rest
  are its derivatives along the directions the value pass carried. First order is
  one direction and the same code path, so nothing about it changed.
* An objective that declines to build a Hessian gets an answer rather than a
  crash. `trust()` read the `hessian` element as a matrix at three points
  outside the handler that turns an objective's failure into a rejected step,
  so a `NULL` surfaced as an Rcpp conversion error naming neither the objective
  nor the argument. Each point now decides for itself: a quasi-Newton start
  seeds the identity and says so, a handover keeps the approximation it had, and
  a trial point counts as a failed evaluation. `hessianMethod = "gn"` and
  `boundary = "clip"` need one at every iterate and stop at `parinit` instead.

* A backward solve that returns no adjoint is an error. It used to read as a
  cotangent of zero, so a backend that dropped the answer produced a gradient
  that was quietly zero in every direction it dropped rather than a failure.
  Found by the CVODES batch fault fixed in cppDE.
* A reverse gradient integrates the states once. The value pass of a reverse
  evaluation now runs on the reverse object itself and keeps its checkpoints,
  and the backward pass replays them instead of integrating the trajectory a
  second time. Nothing to switch on: it happens wherever `obj(pars, sweep =
  "reverse")` walks a prediction built with
  `odemodel(..., derivMode = c("forward", "reverse"))`. A store is matched on
  the times and parameters it was taken at, so it can never answer for another
  point.
* The backward pass can weight its own step size. `Xs(..., optionsReverse =
  list(gradtol = ))` hands the adjoint of the previous evaluation to the
  controller, which refines the grid where a step carries objective error. The
  term enters under a maximum, so the grid only ever becomes finer than
  `abstol` and `reltol` ask: a weight from a parameter the optimiser has since
  left costs steps and never accuracy. Off by default, and worth turning on
  only where a measurement says the steps saved beat the steps spent.
* `importPEtab(backend = "Sundials")` imports onto the CVODE backend, which
  `odemodel()` has carried all along. `derivMode` is matched at the door there
  too, and `"reverse"` on `backend = "deSolve"` says so rather than failing
  later: the deSolve backend goes forward only.
* A PEtab prior term honours `hessian = FALSE`. It did not: the argument fell
  into `...` and was ignored, so an objective carrying a prior handed back a
  zero Hessian to a caller that asked for none. Under `sweep = "reverse"`, which
  cannot produce a Hessian at all, that made the objective look as though it
  had, which is the invariant a caller uses to check the direction arrived.
* `backend = "Sundials"` goes backwards as well.
  `odemodel(..., backend = "Sundials", derivMode = c("forward", "reverse"))`
  compiles CVODES adjoint sensitivity analysis beside the forward object, and
  `obj(pars, sweep = "reverse")` walks the same chain through it. It refuses
  events, which the cppDE backend carries. Having two independent adjoints
  under one interface is what makes a systematic error in either visible;
  `inst/examples/example_AdjointComparison.R` puts the three routes to a
  gradient side by side on Bachmann.
* One argument names the derivative direction everywhere it is chosen at build
  time. `odemodel(reverse = TRUE)` becomes
  `odemodel(derivMode = c("forward", "reverse"))`, and the `derivMode` of `Y()`,
  `Pexpl()` and `importPEtab()` takes the same vocabulary: `"forward"`,
  `"reverse"`, `"symbolic"`, one or more at a time. `derivMode = "dual"` is gone
  and is spelled `"forward"`. There is no deprecation path; both old spellings
  are errors. `normL2(sweep = )` keeps its name, being a choice per call rather
  than a property of a compiled object.
* Gradients can be taken backwards through the whole chain.
  `obj(pars, sweep = "reverse")` walks `normL2 -> Y -> Xs -> P`, seeding the
  objective's cotangent at the data and pushing it back to the outer parameters
  in one sweep instead of propagating one tangent per parameter forward. It needs
  `odemodel(..., derivMode = c("forward", "reverse"))` and compiled observation
  and transformation functions, since the reverse path has no interpreted
  fallback. Its cost does
  not grow with the number of parameters.
* The composition algebra gained a backward pass. A `*` node cannot be walked in
  one recursion, since p2 has to be evaluated before p1 runs and differentiated
  after p1 has been, so the protocol is two phases with an explicit tape between
  them, holding each node's forward values so nothing is recomputed. Conditions,
  the `+` of several branches, the batch entry and `dMod.batch.check` all carry
  through.
* A reverse objective returns no Hessian, because there is no Jacobian to
  contract. That is what the quasi-Newton Hessian sources want.
* The reverse gradient belongs to the trajectory a value-only prediction
  produces, so value and gradient are consistent with each other. Under forward
  sensitivities they are not: the error norm takes the maximum over every tangent
  column, so a sensitivity solve steps finer than the value solve reported
  beside it.
* `importPEtab(..., reverse = TRUE)` builds the reverse object, and
  `inst/benchmarks/bench_gradientCost.R` measures the adjoint next to the forward
  gradient it was written to be compared against.
* `inst/examples/example_ReverseAD.R` walks the chain piece by piece on toy
  models; the Boehm example gained a reverse section.
* Every objective wrapper now declares `sweep` rather than forwarding it through
  `...`. A caller decides which direction a term supports by reading its formals,
  and a wrapper that forwarded the argument without naming it read as a term with
  no reverse path, so a reverse call quietly returned a forward gradient.
* `cppDE::funCpp()` is now `cppDE::cppFUN()`.

# dMod2 0.7.4

* A prepared ODE batch handle no longer outlives the shared object it was
  resolved from. `Xs()` caches that handle for the derivative path and keyed the
  cache on shapes and labels only, so a workspace shipped to a cluster node went
  on calling an entry point of a shared object that was never built there. Every
  `mstrust()` start then failed with `parinit not feasible`, and the real error,
  `"solve_x_s_batch" not available for .Call()`, showed up only when the
  objective was evaluated by hand. The cache now also checks that the shared
  object is still loaded, and `modelname<-` drops it. The rename loop in the
  scripts `runbg()` and `distributedComputing()` generate covers objective
  functions too, which it had skipped.
* Needs cppDE 0.9.5, where a Windows install detects OpenMP. Without it a
  batched solve stays serial there.

# dMod2 0.7.3

* SBML and PEtab import work again in a fresh session. `.dmod_libsbml_python()`
  read `reticulate::py_exe()` before Python was initialised, where it falls back
  to discovery, misses the reticulate-managed environment and returns an empty
  string. The `system2()` that followed then failed with status 127, which the
  resolver reported as a missing `python-libsbml`, so a working environment
  looked like a broken one. It now initialises Python first and rejects an empty
  path.
* The libsbml integration tests did not catch that, because `.libsbml_works()`
  turns any error into `FALSE` and every test behind it into a skip. A test now
  covers the interpreter resolution itself, outside that guard.

# dMod2 0.7.2

* `trust()` reports the iterate the run actually reached. The best-iterate
  bookkeeping introduced with non-monotone acceptance took its snapshot at the
  head of an iteration, so a run ending on `fvalue`, `preddiff`, `step` or the
  iteration limit reported the second to last accepted point. `argument`,
  `value`, `gradient`, `hessian` and `atBound` were affected, and with them
  `vcov()`, `confint()` and profiles. The objective value moved by less than
  `ftol` on the three soft stops, but the gradient could be orders of magnitude
  too large, and an iteration-limit stop lost a full step. Runs ending on
  `gradient`, `stagnation`, `radius` or `objfun` were never affected.
* `trust()` takes `hessianFallback` and `fallbackLimit`. The fallback takes
  over at the first value, model, step or stagnation stop instead of ending the
  run; `fallbackLimit` above one alternates back to the primary source at the
  next such stop, re-seeding a Gauss-Newton phase from a fresh Hessian at the
  cost of one evaluation. `hessianMethod = "hybrid"` is gone, not deprecated: a
  run is a method and, optionally, a fallback, and no pair carries a name of its
  own. It was `"gn"` with `hessianFallback = "bfgs"`, and `trust()` says so when
  it is passed.
* `as.parframe()` carries `nSwitch`, the number of Hessian source handovers, so
  a multi-start can be scored on them.
* `trustL1()` is gone. It carried its own trust-region driver, shared no
  acceptance semantics with `trust()` after the changes above, and had no caller
  left here: the L1 penalty, its clustering and the EM layer that use it live on
  `devel-EM`, and it returns with them once that layer lands, aligned with the
  interface above. `constraintL1()` stays as the Laplace prior it always was.
* `createExample()`, `exmpextr()` and `extractExamples()` are gone. They served
  the pre-testthat `inst/tests` layout, and `createExample()` opened an editor
  on a template to write a unit test, which is maintainer tooling rather than
  package API.
* `distributedComputing()` and `runbg()` attach `dMod2` on the remote machine
  and nothing else. They used to replicate every package attached in the
  submitting session, plus a hard-coded `tidyverse`, which made a job depend on
  the session that sent it and filled the node logs with failed loads. A job
  that needs another package should attach it in its own expression.
* `trust()` groups its arguments into `tolControl`, `qnControl` and
  `stepControl`. The flat `ftol`, `mtol`, `gtol`, `xtol`, `rmin`, `boundary`,
  `theta.max`, `hessianInit`, `qnMemory` and `qnCautious` are gone, not
  deprecated; `trust()` and `mstrust()` name the control member a moved
  argument became. `trustL1` keeps its flat tolerances.
* `trust()` takes `qnControl$qnRejected` (SR1 also updates from a rejected
  trial point, Nocedal and Wright sec. 6.2; default `TRUE` wherever an SR1
  phase can occur) and `stepControl$nonmonotone` (Zhang-Hager acceptance,
  default `0`, the monotone rule).

# dMod2 0.7.1

* `trust()` ends a quasi-Newton run on the gradient. `fvalue`, `preddiff` and
  `step` read the prediction of the quadratic model, which during a
  quasi-Newton phase describes the approximation being built rather than the
  iterate, so they no longer terminate such a run; `gradient`, `stagnation` and
  the radius still do. This follows the usual recommendation for quasi-Newton
  methods. On the Boehm model it raises the starts that reach the published
  optimum from 10 to 13 for `"bfgs"` and from 12 to 18 for `"sr1"`, and lowers
  the evaluations spent per such start in both.
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
* `mstrust()` drops the deprecated `studyname` argument. Use `name`.
* `compile()` builds model shared objects correctly on Windows. The per-source
  link and the combined-output object compile ran `system()` with a `2>&1`
  token, but R's `system()` on Windows has no shell: the token reached
  `R CMD SHLIB` as the override `PKG_LIBS=2>&1` and the compiler as an input
  file, so the object never built and the `.dll` never linked against BLAS,
  LAPACK or the Sundials solver. Both now go through `system2()`.
* Adds `inst/examples/example_Boehm_JProteomeRes2014.R`, which builds the Boehm
  et al. (2014) STAT5 model and compares the Hessian sources over a multi-start.
* `define()` and `insert()` resolve the values passed through `...` in the frame
  they were called from. Resolution reached the global environment only, so the
  same call that worked at top level failed with an "object not found" inside a
  function or a knitr chunk. Condition columns and `.currentSymbols` keep
  precedence over the calling frame.  The internal names of `insert()` are no longer visible to those
  values either; `.currentTrafo` and `.currentSymbols` are the documented way to
  reach the branch being rewritten.
* `plotValues()` marks a converged fit with a circle and an unconverged one with
  a triangle, in every plot, and keeps both in the legend. The mapping followed
  the levels present in the data, so it could differ between two plots of the
  same kind.
* `subset()` on a `parframe` evaluates its condition against the columns of the
  frame and then in the frame it was written in. The second step reached the
  package namespace instead, so a condition naming a local variable failed from
  inside a function.
* Adds the `Optimisation` vignette on Hessian sources, seeding and the economics
  of a multi-start, shipped pre-rendered, with
  `inst/benchmarks/bench_hessianSource.R` behind its numbers.
* Attaching the package reports that BLAS is pinned to one thread inside the
  forked workers of `mstrust()` and `profile()`, so a serial BLAS there is not a
  surprise, and says so when no thread-control entry point was found and the
  deadlock is still reachable. `options(dMod.quiet = TRUE)` suppresses the line.
  The pin itself is cppDE's and needs no cooperation here; requires cppDE 0.9.4.

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
