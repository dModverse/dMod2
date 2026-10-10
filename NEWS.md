# dMod2 (devel-EM)

* `scanL1()` gains `lambda = "em"`: the strength of the L1 or Lq penalty is
  estimated per family (gates, reference parameters) by an EM on the marginal
  likelihood, as a variance component, instead of scanned over a grid. One
  multistart of the EM replaces the multistart per `lambda`. With
  `control$em$adaptive`, a reference parameter is penalised relative to its
  full estimate. Refits along the terms ordered by their penalised size make
  the final choice by likelihood ratio test.
* New `selectEB()`: selects which gated parameters and which reference
  parameters are present, with the penalty estimated instead of scanned. The
  slab is a g-prior on the effects, the share of present candidates is
  integrated under a beta-binomial prior, and a multistart over structures
  alternates single switches of candidates with the closed-form update of
  `g`. Per candidate it reports the gain in -2 log L, the threshold and the
  inclusion probability; `rule = "alpha"` thresholds at the chi-square
  quantile instead. Refits with a likelihood ratio test against the full model
  make the final choice, on the data term with `lrt = "data"`. Substitutes are
  reported as alternative structures and near-collinear groups. Vignette
  "Structure selection with an estimated penalty".
* New `scanL1()`: L1 selection over a grid of penalty strengths in the
  manner of Hauber, Rosenblatt and Timmer (2023). A penalised multistart fit
  per `lambda`, an unpenalised refit of each distinct structure, and the
  choice by likelihood ratio test against the full model or by BIC. Penalties
  pull gates to zero, fold changes to zero, or fuse parameters of a block
  pairwise; `q < 1` by reweighted L1. `gateL1()` inserts the gates, so a
  parameter in a log parametrisation can reach exactly zero.
* `scanL1()` gains `ssl`, the spike-and-slab lasso of Rockova and George
  (2018): every penalised term comes from a slab or a spike, an EM alternates
  closed-form inclusion probabilities with weighted `trustL1()` fits, and one
  warm-started chain over the spike strength replaces the multistart per
  `lambda`. `select = "plateau"` takes the structure the chain settles on.
  `pathFits` sets the starts per `lambda` and per refit apart from `fits`, the
  starts of the full model; a refit also starts from its penalised optimum.
  `plot()` shows the waterfall of the full model, the groups of every block
  over `lambda`, and the inclusion probabilities.
* `scanL1()` gains the control entries `hits`, `tolHits` and `maxFits`: a
  multistart of the full model, of a `lambda` or of a refit adds batches of
  starts until its best value is reached by `hits` runs. Starts and hits are
  reported per `lambda`, per refit and for the full model; `arguments` holds
  all free parameters of every penalised fit.
* The spike-and-slab lasso of `scanL1()` puts its block terms on the gaps
  between neighbours of the sorted values (Ke, Fan and Wu 2015), so that a
  split costs one slab term whatever the cluster sizes; spike terms within
  `control$snap` of their kink are put onto it.
* With `q < 1`, every start of a `scanL1()` multistart is reweighted and
  compared on the Lq objective; when the downward pass finds a better optimum,
  starts are added until the waterfall reaches it.
* `control$trust` of `scanL1()` takes `ftol`, `mtol` and `gtol` for both
  `trust()` and `trustL1()`.
* The spike-and-slab path of `scanL1()` runs the EM at each `lambda` from the
  previous mode and from the sparse point and keeps the smaller `-2 log`
  posterior, followed by the downward pass that `q < 1` uses.
* `scanL1()` returns `level` and `levelFits`: per `lambda`, the runs on the
  lowest level of the waterfall, within `tolHits` of the best value, and the
  structures they take.

# dMod2 0.11.3

* New `trustL1()`: the trust region of `trust()` with an L1 penalty on chosen
  parameters around reference values `mu`, one-sided with `one.sided`. `gate`
  keeps a parameter on or above its kink, so it can reach the reference
  exactly, and `fuse` penalises the pairwise differences within blocks, whose
  members move as one while equal. The penalised step runs in a C++ kernel.

# dMod2 0.11.2

* `symmetryDetection()` and `symmetryReduction()` install symident from PyPI on
  first use; under `RETICULATE_PYTHON`, `pip install symident`.

# dMod2 0.11.1

* `distributedComputing()` runs zstd as a pipe stage of its own and moves the
  uploaded files one by one. macOS reads `tar -I` as a file list and its `mv`
  has no `-t`, so the upload and the collection of results failed there.

# dMod2 0.11.0

* Breaking: `Xs()` on cppDE takes `refine` and `gradtol` in `optionsReverse`;
  `floor` is gone, and so is the weighting of the backward grid by the
  previous evaluation's adjoint. `refine = TRUE` holds each step of the
  backward sweep to the error test of the CVODES backward problem under
  `abstol` and `reltol`, through `cppDE::adjointControl(refine = TRUE)`, and
  `gradtol` is then the absolute tolerance on each step's share of the
  gradient. The reverse value pass takes the grid of a value solve.
* New `symmetryDetection()`, `reconstControl()` and `symmetryReduction()`:
  structural non-identifiabilities of an ODE model, as generators of its
  symmetries, and the parameter transformation that removes them. The
  computation runs in the Python package symident through reticulate, which
  installs it on first use from <https://github.com/dModverse/symident>; R
  arguments and fields are camelCase, `print()` and `summary()` are the
  reports of symident. symident's switches are read from `SYMIDENT_*` and from
  the options `dMod.sym.*`. msolve, for coupled steady states, comes from
  `symident.install_msolve()` or `SYMIDENT_MSOLVE`.
* reticulate no longer provisions symengine.
* `Xs()` on cppDE and Sundials stops on a failed solve (`onFailure = "stop"`),
  as `Xf()` does; a fit takes the error as a rejected step.
* `Xs()` on cppDE: `optionsSens = list(sensErrCon = FALSE)` takes the
  sensitivities out of the error test.
* `odemodel(backend = "Sundials", derivMode = "reverse")` takes events: the
  reverse objective takes the adjoint and the gradient through each jump.
* Reverse objectives run their solves through prepared batch handles,
  evaluate the error model batched and keep the data-to-prediction indices
  between evaluations. Under second order a transformation that no tangent
  reaches skips its second derivatives.
* `normL2()` with one condition evaluates the error model without derivatives
  when none are requested.
* `importPEtab()` gains `sparse`, which pins the sparse (KLU) or dense linear
  solver of the cppDE and Sundials backends.
* `benchmarks/`, outside the built package: value, forward and reverse
  gradients of PEtab benchmark problems by tier (`run-benchmarks.R --tier
  tiny|medium|full`), with the scripts of `inst/benchmarks/` under
  `benchmarks/scripts/`.
* A sum of objectives keeps the `sweep` attribute its terms agree on, so a
  reverse objective with a prior says it was evaluated in reverse.
* `constraintL1()` and `constraintL2()` are S3 generics; the default methods
  are the priors on named parameters.
* The option `dMod.outdir` sets where `odemodel()`, `P()`, `Pimpl()`, `Y()` and the
  PEtab import write generated sources and shared objects; unset, it is the
  working directory as before.
* `getSymbols()` and `replaceSymbols()` are reexported from cOde.
* Breaking: `Xs()` takes its solver options as `options`, which applies to
  every solve, and `optionsSens`, which now holds only the entries that
  override `options` for the solves with sensitivities. `Xf()` and
  `importPEtab()` take `options` as well, and the controls are named
  `options`, `optionsSens` and `optionsReverse`. `optionsOde` is still
  accepted as an argument and by `controls()`, with a deprecation warning.
  On deSolve the default methods stay `"lsoda"` and, for the solves with
  sensitivities, `"lsodes"`.
* `Xs()` on cppDE warns about an unknown entry of `optionsReverse`, given to
  the constructor or set through `controls()`.
* The solver option for the limit of consecutive rejected steps is
  `maxattempts`, as in cppDE.
* Fix: `constraintL2() * P()` takes `fixed`, and with it `trust(fixed = )` and
  `profile()` on such priors.
* Fix: `Xt()` returns its sensitivities by the outer parameters, so
  `g * Xt() * p` has derivatives and `normL2()` on it takes `fixed`.
* Fix: `Y()` differentiates by its own parameters when no transformation
  precedes the prediction: `g * x` lacked the derivatives by observation and
  error parameters that are not parameters of `x`.
* Fix: `normL2()` takes an error model whose parameters are all fixed.
* Fix: `as.datalist()` stores a `sigma` or `lloq` column of `NA` alone as
  numeric.
* Fix: `mstrust(cautiousMode = TRUE)` removes the `.Rda` fits of its own run
  rather than files of that name in the working directory.
* Breaking: arguments in dot case are camelCase. `Y()` and `Pexpl()` take
  `attachInput` (was `attach.input`), `Pimpl()` `keepRoot` (was
  `keep.root`); both are also the names of the controls.
* Breaking: `normL2()`, `datapointL2()` and the `constraint*()` priors take
  `attrName` (was `attr.name`), also as control of `datapointL2()`;
  `normL2()` and `evalConditionResidual()` take `optBLOQ` (was `opt.BLOQ`).
* Breaking: `datapointL2()` takes `parameter` (was `value`), the name of the
  parameter holding the data value.
* Breaking: `as.datalist()` takes `splitBy` and `keepCovariates` (were
  `split.by` and `keep.covariates`), `as.parframe()` `sortBy` (was `sort.by`),
  `parframe()` `objAttributes` (was `obj.attributes`; the attribute keeps its
  name).
* Breaking: `reml()` and `remlLeverage()` take `rankTol` (was `rank.tol`).
* Breaking: `trust(stepControl = )` takes `thetaMax` (was `theta.max`);
  `mstrust()` takes `startFromCenter` (was `start1stfromCenter`).
* Breaking: `distributedComputing()` takes `memPerCore`, `sshPasswd`,
  `varValues`, `nRep`, `purgeLocal` and `customFolders` (were `mem_per_core`,
  `ssh_passwd`, `var_values`, `no_rep`, `purge_local` and `custom_folders`);
  its `purge()` takes `purgeLocal`.
* Breaking: `plotPathsMulti()` and `plotProfilesAndPaths()` take `whichPar`
  and `nPars` (were `whichpars` and `npars`), `plotProfilesAndPaths()` `ncol`
  (was `ncols`), `plotArray()` `nSim` (was `nsimus`), `plotFluxes()`
  `legendTitle` (was `nameFlux`), `attrs()` `which` (was `atr`) and
  `profileParsPerNode()` `parsPerNode` (was `fits_per_node`).
* Breaking: `reconstControl()` takes `minSupportCandCap`, `perPrimeCap` and
  `perPrimeMinPrimes` (were `minsupportCandCap`, `perprimeCap` and
  `perprimeMinPrimes`).
* Breaking: `readPEtabYaml()` and `readPEtabTables()` replace
  `readPetabYaml()` and `readPetabTables()`. An imported problem prints as
  `<PEtab problem ...>`.
* Breaking: `repar(trafo, expr)` takes the transformation first, as
  `define()` and `insert()` do. The order `repar(expr, trafo)` is recognised
  by an expression in first place.
* Breaking: `profile()` and `vcov()` are methods of the generics of stats,
  which dMod2 no longer masks. `profile()` dispatches on objective functions
  and plain R functions, first argument `fitted` (was `objfun`). `trust()`
  returns a list of class `trustfit`, kept by `mstrust()`, on which `vcov()`
  dispatches, first argument `object` (was `fit`); a fit stored before needs
  `class(fit) <- c("trustfit", "list")`.
* Breaking: the operands of `+` and `*` on dMod functions, objectives,
  objlists and datalists are named `e1` and `e2`.
* The old names above are accepted with a deprecation warning; giving the old
  and the new name together is an error.
* Breaking: `distributedComputing()` submits by default (`recover = FALSE`).
  The `get()` of `distributedComputing()` and `runbg()` returns the results
  instead of assigning `cluster_result` or `.runbgOutput` in the global
  environment; `runbg(wait = TRUE)` returns them too. The former behaviour is
  `cluster_result <- job$get()`.
* Breaking: `steadyStates()` writes nothing unless `file` is given; the model
  for the solver goes to a temporary file. Its unused argument `rates` is
  deprecated and ignored, and `verbose` defaults to `FALSE`.
* Breaking: `fitErrorModel()` plots only with `plotting = TRUE`;
  `resetWarmStarts()` and `symmetryDetection()` print only with
  `verbose = TRUE`.
* Breaking: `symmetryDetection()` requires `f` and `g`.
* `ggopen()` opens the PDF with the viewer of the platform by default:
  `open` on macOS, `xdg-open` on Linux, `shell.exec()` on Windows.
* `trust(parscale = )` is deprecated and warns; it still scales the trust
  region. Fit on log scale instead.
* `msParframe()` leaves the global random number generator as it was.
  `mstrust()` takes `samplefun` as a function as well as its name.
* `rref()` returns a list named `rref` and `pivots`; its unused argument
  `fractions` is deprecated and ignored.
* Breaking: `Pimpl(controlsPTC)` takes the tolerances as `reltol` and `abstol`,
  as `cppDE::ptc()` does; `rtol` and `atol` are deprecated.
* `compile()` restores `PKG_*` flags set before the call without an error,
  and leaves no `symbols.rds` in the working directory under `R CMD check`.
* `reconstControl(homotopy = )` switches the last reconstruction route of
  symident, which fits directions still open under `equilibrate = TRUE` along
  lines in parameter space; on by default.
* Needs cppDE 0.12.0.

# dMod2 0.10.3

* The help pages of the constructors (`odemodel()`, `Xs()`, `Xf()`, `Xd()`,
  `Xt()`, `Y()`, `P()`, `Pexpl()`, `Pimpl()`, `normL2()`, `datapointL2()`,
  the `constraint*()` priors) and of the function classes are rewritten, with
  examples for `Xs()`, `P()`, `Pexpl()` and `Pimpl()`.
* `Xs()` and `Xf()` on cppDE and Sundials return the forcings after the
  states, so observables can contain them, and take forcing names as factor.
  A forcing holds its first and last value outside its points and may be a
  single point.
* `Xf()` on cppDE and Sundials starts states missing from `pars` at 0.
* `datapointL2()` returns an `objlist`.
* `constraintL2()` stops on a `sigma` that mixes numbers and parameter names.
* `normL2()` names the data conditions `x` does not have.
* `plotValues()` puts its breaks at 0 and at the decades from 10 on and no
  longer fixes the limits of the y axis.
* Needs cppDE 0.11.3.

# dMod2 0.10.2

* `runbg()` and `distributedComputing()` with `compile = TRUE` take SUNDIALS'
  include path, and MPI's if that SUNDIALS was built with it, from the cppDE
  installed on the remote machine instead of the submitting one. A model using
  SUNDIALS' LAPACK dense solver links it there too, and the build stops with a
  message if the remote SUNDIALS lacks it.

# dMod2 0.10.1

* `Pimpl()` draws its random starts uniformly over `startRange = c(-5, 5)` on
  `startScale`, log10 by default and linear if `positive = FALSE`. Each
  warm-start cache, i.e. each condition in each fit of `mstrust()`, draws from its
  own stream; the global RNG is read, not advanced. `maxit` defaults to
  `ceiling(70 * log(n + 1))` for `n` states (was 400).
* Needs cppDE 0.11.2: `ptc()` scales a row by its largest partial derivative, so
  a row whose own rate is tiny next to its partner's no longer stalls it, and it
  builds on macOS.

# dMod2 0.10.0

* `normL2()` loses `t0`. `times` is a vector for all conditions or a list named
  by condition, and a condition's prediction starts at the first time of its
  grid, its data times and its `times`: `times = list(C1 = -10, C2 = 5)` is what
  `t0 = c(C1 = -10, C2 = 5)` was. Without `times` the grid starts at the first
  data point, which matters only for a model built with
  `odemodel(..., includeTimeZero = FALSE)`; with the default the grid holds 0.
  PEtab start times go through `times`, in one objective for all conditions.
* **Multiple shooting.**
  `normL2(data, g * x * p, multipleShootingControl = TRUE)`, or a list of
  settings, makes `trust()` and `mstrust()` fit by Bock's generalised
  Gauss-Newton method. The time axis is cut into segments with node values of
  their own; continuity is a constraint, linearised and eliminated at every
  iterate (condensing), so the subproblem keeps the single-shooting size and
  box bounds work as before. At convergence the result is a stationary point
  of the single-shooting objective, which the objective still is when called
  as a function. The ODE model needs `odemodel(..., includeTimeZero = FALSE)`;
  `profile()` re-optimises by single shooting.
* The trust region bounds parameters and node values together. Steps are
  accepted by a filter or an l2 merit function (`stepControl$acceptance`),
  with a second-order correction (`stepControl$soc`) and a gap tolerance
  (`tolControl$ctol`).
* `stepControl$anneal` (default) first penalises the gaps, `f + w ||c||^2`
  with `w` rising tenfold per stage, then hands over to the exact method. A
  stage ends once a step gains less than a thousandth; the annealing ends
  early where the penalty dominates and a stage no longer closes the gaps.
* `nodes = "auto"` (default) starts from ten segments per condition, or one
  per half oscillation of the data where that is more, and cuts a segment
  where its propagation matrix grows by more than `growth`, where it misses its
  data by more than `misfit` and by twice the median of its condition, or where
  its solve fails. A misfit everywhere is left to the parameters, and each half
  of a cut keeps `minPoints` data points, by default one more than the states.
  Unobserved node values start from a run synchronised to the data. `charts`
  puts node values on a log10, linear or angle scale per state; the gap of an
  angle is taken modulo 2 pi. The scale of an unobserved linear state is the
  largest range it covers within one segment. `hessianMethod = "bfgs"` or
  `"sr1"` keeps one quasi-Newton block per segment.
* The cost of an evaluation is linear in the number of segments: segment
  parameters are built per condition, the annealing model is assembled in one
  pass, and the values at new nodes come from one batched solve.
  `fit$multipleShooting$cuts` counts the cuts by growth, misfit and failed
  solve.
* `c()` of one parameter vector, alone or beside `NULL` or an empty one, as
  every prediction builds its `parameters`, skips the list machinery.
* Methods of the Freiburg group (Bock; Horbelt, Timmer and Voss 2002; Peifer
  and Timmer 2007; Voss, Timmer and Kurths 2004), off by default:
  `stepControl$acceptance = "natural"` (damped Gauss-Newton on the natural
  level function, needs a fixed sigma), `stepControl$twoPhase` (parameters
  first, then nodes), `stepControl$regularise`, `stepControl$restore`
  (restoration phase of the filter), `init = "spline"`, `minPoints`,
  `nodes = "transitions"` (a node just before each fast change of the data)
  and `breaks` (nodes without continuity).
* Spiking data need tight integrator tolerances (1e-10); at 1e-8 the objective
  is too noisy far from the optimum and the trust region collapses.
* The objective returns the weighted residuals of each segment and their
  Jacobian on request (`residuals = TRUE`). A trial point after a rejected
  step is evaluated without sensitivities.
* The reverse walk takes several first-order seeds at once, which `Xs()`
  answers in one backward solve. Checkpoints serve any number of backward
  sweeps.
* New benches `bench/multipleShooting_preBotC.R` and
  `bench/multipleShooting_macrospin.R`.

# dMod2 0.9.2

* `Pimpl()` runs a multistart when the warm starts and the initial guess fail:
  `controlsPTC(nStarts = 20, startSd = 2, seed = 1)`, log-normal around the guess,
  from a fixed seed and without touching the global RNG. If every start fails it
  is an error.
* `plotValues()` draws the objective values minus the best one on a
  pseudo-log10 scale, linear below 1; `"data"` holds the distance in `delta`.

# dMod2 0.9.1

* `exportPEtab()` reads the parameter scale off `p`: a parameter entering only
  as `exp(X)` is written with `parameterScale = log`, one entering only as
  `10^(X)` or `exp10(X)` with `log10`, with linear nominal values and bounds
  five decades around them. `parameterScale = NULL` is the new default. v2
  writes the linear values.
* `exportPEtab()` no longer writes self-referencing files for a trafo like
  `A ~ exp(A)`. An outer parameter named like a state goes out as
  `init_<state>` and the condition table sets the species to it; one named
  like an inner parameter keeps its name under a pure scale wrap and becomes
  `<name>_outer` otherwise. Symbolic initial values go to the condition table
  instead of SBML `<initialAssignment>`. An export that would still refer an
  assignment to its own target stops with an error.
* `exportPEtabObject()` to v1 keeps preequilibration conditions in the
  condition table, leaves condition targets and compartment sizes out of the
  parameter table, writes observables in a v2 noise formula as their formula,
  and refuses condition switches during a simulation, which v1 cannot hold.
  A v2 experiment names no condition where the table has none.
* `exportPEtabObject()` to v1 moves an expression that sets a species in every
  condition into an SBML initial assignment and warns about per-condition
  expressions, which v1 does not allow; start times after 0 stop the export.
  Measurements name their conditions by id, as the condition table does.
* `exportPEtabObject()` to v1 keeps the SBML constants of an imported problem
  in the SBML, as v2 does. Listed in parameters.tsv they lost their role in a
  preequilibration, where a constant at 0 idles its reaction, and the
  reimport of Isensee_JCB2018 equilibrated to another state.
* `exportPEtabObject()` declares condition targets and compartment sizes that
  nothing else states, with the SBML default recorded at import (a surface
  compartment of Lang_PLOSComputBiol2024 went out sized by itself).
* `exportPEtabObject()` writes a single measured sigma per observable as its
  noise formula instead of per-row noise parameters.
* `exportPEtabObject()` translates priors between the versions: v1 spells the
  log families `logNormal`/`logLaplace`, and a v1 prior on the parameter scale
  goes to v2 as the matching prior of the linear parameter. `importPEtab()`
  reads the v1 spellings and truncates a prior on the parameter scale at the
  bounds on that scale; it used the linear bounds (Schwen_PONE2014).
* `importSbml()` substitutes the initial assignment of a constant parameter
  wherever the parameter appears, with species at their initial values. It
  kept the SBML default, so a condition that set the source parameter did not
  reach the rates (Laske_PLOSComputBiol2019: `ModelValue_80 := k_syn_P`;
  Bertozzi_PNAS2020: `beta_N := R0_*gamma_/N_`).
* `repar()`, `insert()` and `define()` keep a trailing underscore in an
  identifier. `gamma_` became `gamma`, so a PEtab condition that set `gamma_`,
  `N_`, `I0_` or `R0_` in Bertozzi_PNAS2020 changed nothing.
* `Y()` returns a NaN observable instead of stopping, and `normL2()` stops
  only on a NaN at a data point. A ratio of states that all start at 0 is
  undefined at t0 alone, where Laske_PLOSComputBiol2019 has no data, and the
  objective could not be evaluated.
* `importPEtab(derivMode = c("forward", "reverse"))` builds the observation,
  error model and parameter transformation for the reverse sweep as well;
  only the ODE model had it, so `obj(sweep = "reverse")` failed.
* `exportSbml()` declares species that enter a rate without being consumed or
  produced as modifiers, as SBML requires, and stops on a compartment sized
  by its own symbol without a value.

# dMod2 0.9.0

* `Pimpl()` solves steady states by pseudo-transient continuation: implicit
  Euler steps under local error control that turn into Newton steps near the
  root. It follows the flow to a stable steady state and checks the
  eigenvalues on the conservation manifold; `flow = FALSE` returns any regular
  root of plain equations. The iteration runs in C++, `cppDE::ptc()`.
* `Pimpl()` keeps every species and adds the conservation rows `C x = T`
  instead of eliminating pivot species. Iterates stay positive and on the
  manifold. States without influx at the given parameter values are 0.
* `Pimpl()` starts deterministically: the nearest kept root of the condition
  corrected by its sensitivities, then the initial guess. Conditions with
  identical inputs are solved once. A failed solve is an error.
* `Pimpl()`: `controlsPTC` replaces `controlsMS` and `controlsNleqslv`,
  `expressInTotals` is gone. nleqslv is no longer a dependency.
* `Pequil()` is removed. `P()` builds an `eqnlist` with `Pimpl()`;
  `method = "equilibrate"` is gone.
* `importPEtab()` computes preequilibration with `Pimpl()`. Events do not act
  during preequilibration (PEtab v2 case 0023).
* `importPEtab()`: species a preequilibration cannot move keep their initial
  values. They take part only in reactions idled by SBML constants at 0 or by
  species that start and stay at 0, and the steady-state equations leave them
  undetermined (Isensee_JCB2018: AC/pAC and PDE/pPDE).

# dMod2 0.8.4

* `steadyStates()` 1.3 and 1.4 detect structurally zero clusters that span
  compartments, such as a ligand in the medium bound by a receptor in the
  cell. The sink-cluster test runs in amounts, over reactions rather than the
  rows per volume ratio the backend reads. As in 0.8.3, a conserved moiety in
  the support of the test stays nonzero.

# dMod2 0.8.3

* Structurally zero states: a conserved moiety listed before a leaking cluster
  was declared zero together with it, since the sink-cluster LP returns their
  union. Only species that drain into a leaking reaction are zero now.

# dMod2 0.8.2

* `steadyStates(version = "1.4")` prints the resolved expressions through sympy
  again, as 1.3 does: no redundant brackets from the substitution, nothing
  expanded.
* `steadyStates(version = "1.4", solveQuadratic = TRUE)`: a state whose own
  balance, denominators cleared, is quadratic with positive production and
  quadratic consumption takes the unique positive root, tried before any rate
  constant. A root that would take a balance another unknown needs is dropped
  for the next attempt.
* `steadyStates(version = "1.4", givenCQs = )`: each conserved quantity keeps
  one of its states free, the first one whose balance can be spent elsewhere.
  `customTotals()` reach it as before.
* `steadyStates(version = "1.4")` never solves for a rate constant whose fluxes
  all contain a `neglect`ed symbol, which would divide by that symbol.
* `steadyStates(version = "1.4")` prints one line per attempt instead of the
  intermediate expressions, notes after its summary, and on failure the
  balances left.
* `steadyStates()` handles species named like sympy objects (`Ci`, `E`, `S`,
  `Q`, `gamma`).

# dMod2 0.8.1

* `steadyStates(version = "1.4")`: a new core on the same interface. Balances
  are linear forms over the flux terms, and each unknown comes from a
  combination of balances a linear program finds, so every solution is a ratio
  of positive sums. Where no single unknown is left, one side of a balance
  shares its sum by new flux ratios `r_*`. Needs `positive = TRUE`; volumes are
  never unknowns, and `solveQuadratic`, `branches` and `givenCQs` do not apply.
* `steadyStates(verbose = )`: `TRUE` (default) prints a few progress lines and
  the result, `FALSE` the result only, `"full"` every step.
* `steadyStates()` 1.3 solves linear balances directly instead of by
  `solve()`, reads signs off the expression tree after a numeric screen, and
  searches a block balance by balance when no pivot combination is positive.
* `steadyStates()` solves a strongly connected block of states jointly where no
  single balance stays positive, e.g. receptors recycling through endosomes.
* `steadyStates()` solves every rate constant at most once. A symbol defined
  twice, by itself or in a cycle is an error.
* `steadyStates()` factors out common terms instead of cancelling the resolved
  expressions, which stalled large models.
* `steadyStates(testSteady = "fast")` checks mod p only. A residual it cannot
  evaluate is an error with its reason, not a fallback to the symbolic test.
* `normL2(t0 = )` takes one start per condition, named by condition.
* `Pexpl()` and `Y()` build `derivMode = "forward"` by default and accept
  `"forward-reverse"`, which the reverse sweep with `deriv2 = TRUE` needs.
* `controls()` reaches every option an object uses at run time, and
  `controls<-` changes it for the next call. `Pimpl()` (`keep.root`,
  `controlsMS`, `controlsNleqslv`) and `Pexpl()` (`attach.input`) have controls
  now; `Pequil()` reads `keep.root` and `attach.input` from them, and its memo
  no longer answers a call made with other controls. `Xs()` and `Xf()` merge
  the solver defaults at each call, so replacing `optionsOde` keeps the options
  it does not name. The accessor walks `g * x * p`, the summands of a sum
  (attribute `terms`) and the objective inside `%.*%` and `objfn * parfn`;
  `condition = NULL` sets every condition, and an unknown name or condition is
  an error.
* `mstrust()` routes its arguments by the formals of `optmethod` instead of
  those of `trust()`, and no longer passes `cautiousMode` to the objective.
* **Bug fixes.** `normL2()` no longer reads past a prediction the solver cut
  short, which gave a wrong value or a segfault; it is an error naming the time
  the solver reached, which `trust()` takes as a failed step.
  `Pexpl(attach.input = TRUE)` gives the inputs it passes through their forward
  derivatives, so the forward gradient and Hessian along them are no longer
  zero; inputs in `fixed` are handed on as constants. The name clash check of
  `mstrust()` measures the intersection.
* Requires cppDE 0.10.2.

# dMod2 0.8.0

* **Reverse mode.** `obj(pars, sweep = "reverse")` walks `normL2 -> Y -> Xs -> P`
  backwards in one sweep, at a cost that does not grow with the number of
  parameters; with `deriv2 = TRUE` it returns the exact Hessian at a cost linear
  in it. Needs `odemodel(..., derivMode = c("forward", "reverse"))`, or
  `"forward-reverse"` for second order; `backend = "Sundials"` uses the CVODES
  adjoint, and `importPEtab(backend = "Sundials")` imports onto it. The reverse
  gradient belongs to the trajectory of the value solve and reuses its
  checkpoints; `Xs(optionsReverse = list(gradtol = ))` weights the backward step
  size.
* **Derivative arguments.** `deriv`, `hessian`, `deriv2` and `sweep` each ask
  for one thing, and `curvature` is gone. `derivMode` names the direction at
  build time in `odemodel()`, `Y()`, `Pexpl()` and `importPEtab()`;
  `reverse = TRUE` and `"dual"` are errors. `"symbolic"` is gone, so
  transformation and observation functions evaluate only after `compile()`.
* **Optimisation.** `trust()` asks per evaluation for a value, a gradient, a
  Gauss-Newton or an exact Hessian. `qnControl$hessianInit = "exact"` seeds a
  quasi-Newton run, `hessianReseed = "stall"` fetches a fresh Hessian when it
  stalls, and `hessianMethod = "exact"` is a Newton run.
* Requires cppDE 0.10.0 and follows its names: `tangent`, `hessian`,
  `cotangent` and `curvature`, and `cppFUN()` for `funCpp()`.
* `distributedComputing()` and `runbg()` gain `libs`, library paths put in
  front on the remote side.
* `compile()` reports the toolchain and reused objects through `message()`
  rather than `cat()`, so `suppressMessages()` silences it.
* **Bug fixes.** A summed objective keeps the curvature of every term. An
  objective that declines a Hessian no longer crashes `trust()`. A backward
  solve without an answer is an error instead of a zero gradient. A PEtab prior
  honours `hessian = FALSE`. The batched backward path propagates more than one
  direction. `compile()` stays within the 8191-character command line of
  Windows.
* New examples `inst/examples/example_ReverseAD.R`,
  `example_AdjointComparison.R` and `example_BachmannReverse.R`; the tests
  compile the models of each file jointly.

# dMod2 0.7.5

* `runbg(compile = TRUE)` builds on Debian and Ubuntu machines again. The
  generated build script looked for R's headers in `$(R RHOME)/include`, which
  Debian does not have, so every source failed with `R.h: No such file or
  directory` and the job never started. The script now asks the remote R for
  `R.home("include")`, as `compile()` does locally.
* `runbg(compile = TRUE)` compiles in parallel again. It exports
  `OMP_NUM_THREADS=1` before the build, and `nproc` honours that, so the build
  counted one core and compiled every source in turn. The count now ignores the
  OpenMP variables. `runbg()` also gains `buildCores` and `buildBundle`, as
  `distributedComputing()` has.

# dMod2 0.7.4

* A prepared ODE batch handle no longer outlives the shared object it was
  resolved from. `Xs()` keyed its cache on shapes and labels only, so a
  workspace shipped to a cluster node called an entry point that was never built
  there and every `mstrust()` start failed with `parinit not feasible`. The
  cache now checks that the object is loaded, `modelname<-` drops it, and the
  rename loop in `runbg()` and `distributedComputing()` covers objectives.
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
  run is a method and, optionally, a fallback, and no pair has a name of its
  own. It was `"gn"` with `hessianFallback = "bfgs"`, and `trust()` says so when
  it is passed.
* `as.parframe()` includes `nSwitch`, the number of Hessian source handovers, so
  a multi-start can be scored on them.
* `trustL1()` is gone. It had its own trust-region driver, shared no
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
* Objective functions take a call-time `hessian` argument. With
  `hessian = FALSE` the `J^T J` contraction never runs and the result has a
  `NULL` hessian. `trust()` uses this in the quasi-Newton phase, and it
  propagates through objective composition.
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
  form, and each applies the chain rule so `constraint * P()` is exact.
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
* `normL2()` attaches `compileInfo`, as `*` and `+` already did, so `loadDLL()`
  and `compile()` reach the shared objects of a composed objective.

# dMod2 0.6.1

* `loadDLL()` skips shared objects that are already loaded in the current
  process. Unloading them nulled the native symbol pointers held by the
  prediction, observation and parameter functions built from them, leaving
  every later call without a way to resolve them again.
* `Xs()` caches a prepared batch handle that holds such a pointer of its
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

* `eqnlist` stores compartments and volumes; `assignCompartment()` and
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
