# dMod2 0.15.3

* `trust()` ends with `stopReason = "radius"` once a step no longer changes the parameters in
  floating point, also with the default `rmin = 0`. Before, an objective whose value at the iterate
  differs between evaluations, as the Laplace marginal of `em()` and `msem()` can, kept the run
  evaluating the iterate until `iterlim`. Under `qnControl = list(hessianReseed = "stall")` a
  quasi-Newton phase reseeds first, as at a stall.

* With sparse random effects, `em()` and `msem()` drop the BFGS correction of the model Hessian,
  in the search of the random effect modes and in the joint fit, when a step shows no positive
  curvature, so the search returns to the model Hessian instead of creeping on a stale correction.

# dMod2 0.15.2

* The refit in each step of `profile()` starts at the radius `rinit` of `optControl`. Before, the
  radius was at most three times the length of the step, which for a parameter on a small scale,
  such as the validation parameter of a prediction profile, kept the other parameters in place
  and returned the profile of the validation point alone.

* `subset()` of an equation list keeps a state that a kept reaction does not change but whose rate
  depends on it, as a catalyst.

* `plotCombined(errfn = )` and `plotPrediction(errfn = )` draw the sigma band only where the error
  model gives a sigma, so an observable without one no longer stops the plot.

* `plotArray(fixed = )` adds the parameters held fixed in the profiles to every parameter set, so
  it plots the predictions along profiles computed with `fixed`.

# dMod2 0.15.1

* American spelling throughout: `trust(stepControl = list(regularize = ))` replaces `regularise`,
  the attribute `run` of the result of `mstrust()` names `optimizer` and `optimizerArgs` instead
  of `optimiser` and `optimiserArgs`, and the sequential palette `"grey"` of [dMod_palettes()] is `"gray"`.

* `wide2long()`, `%bmm%`, `rref()`, `nullZ()` and `submatrix()` are no longer exported, and
  `ggopen()` is removed. `as.data.frame()` of a prediction returns the long format of
  `wide2long()`, with a column `sigma`.

# dMod2 0.15.0

* New `constraintStability(pSS, p)`, a term of the objective that keeps the steady state of a
  `Pimpl()` stable: for every eigenvalue of the Jacobian at the steady state on the manifold of the
  conserved quantities, the square of how far the real part passes an onset, scaled to `weight`
  where the real part reaches the `stabilityTol` of `pSS`. Zero below the onset, with exact
  gradient and Gauss-Newton Hessian; the onset defaults to a decay slower than the data last.
  `Pimpl(deriv2 = TRUE)` provides the second derivatives it needs. With
  `controlsPTC = list(stability = FALSE)` the steady state may be unstable, and the term goes on
  rising there instead of the evaluation failing.

* `plotValues()` takes a named list of fit collections, one panel each and each ranked separately,
  and with `relative` draws the values above the best of each collection or above a given value.
  `plotValues()`, `plotPars()` and `as.parvec()` take the result of `mstrust()` as it is, and
  `as.parframe()` returns a parframe unchanged.

* `align(reml = TRUE)` takes the standard errors of the values from the information of the mean
  parameters alone, the error parameters held at their estimates (Lindstrom and Bates 1990).

* The integral of `em()` over sparse random effects is formed in natural parameters, without the
  center of the quadratic expansion, in C++ group by group where no Dirichlet-process random effect
  joins groups. The spike-and-slab random effects of a group are integrated exactly, by a
  depth-first sum over their subsets that leaves out subtrees below a bound, within a budget of
  nodes instead of a limit on their number. Where the data couple factors, or the sum would exceed
  the budget, the expectation consistent (EC) approximation runs as a convergent double loop with
  Newton steps in both loops; EC converges where the data tie coordinates closely or do not see
  them. Mean field and `control$laplace$coupling` are gone, with the fallback to mean field.

* The native evaluation chain takes observables and error models that read `time`; such models
  ran on the R path before.

* `em()` with Laplace, spike-and-slab or Dirichlet-process random effects follows the parameter
  transformations, steady states, ODE, events, observables and error model from every such random
  effect to the measured names of the group. One that reaches none of them rests at 0 with the
  prior: the group searches the mode without it, the integral and the hyperparameter updates
  leave it out, and `$random` marks it in a column `data`. A hyperparameter all of
  whose random effects reach no data is an error. Every expression is read with the names that are 0
  in the condition set to 0, so that a state that starts at 0 and that nothing drives passes no
  random effect on, and a root of a steady state depends only on the inputs that the equations and
  conservation rows of the steady state reach.

* The second derivatives of the Laplace marginal and of the FOCEI volume term by forward over
  reverse take the gradients of their traces direction by direction, each from a forward pass
  along one random effect, where that is cheaper than a backward sweep along all of them per
  random effect.

* `msem(sd = )` takes a vector named by the parameters it perturbs; an unnamed vector of several
  entries is an error. The table of `msem()` has the largest absolute entry of the gradient of every
  start (`gradient`).

* `em()` with sparse random effects finds the modes of the random effects group by group, by
  `trust()` where a group has no Laplace coordinates and by `trustL1()` where it has, in place of
  one trust region over the random effects of all groups. The search differentiates the model in
  the random effects of the group only, adds a BFGS correction to the Gauss-Newton Hessian and ends
  at the `tolControl` of `control$laplace$trustControl`, where given. A group starts from the lowest
  of the predicted modes and the minima found at the last four distinct accepted points, as the
  inner problem of a group may have several local minima; at the start of the fit and every 25 accepted steps the
  search also runs from 0 and from two starts around the modes.

* The model Hessian of the outer trust region of `em()` with sparse random effects under forward
  over reverse is built from Gauss-Newton blocks throughout, which keeps it positive semidefinite;
  the exact rows of the data Hessian in the random effects enter the gradient and the shift of the
  modes.

* The fits of `em()` with sparse random effects count why evaluations of their trust region failed
  (`failures`, also a column of the table of `msem()`). Every start of `msem()` runs to
  convergence, all starts at once; `msem(budget = )` is gone. `msem(which = )`
  fits a part of the starts, so that tasks of `distributedComputing()` share one multistart, and
  `get()` joins their fits into one.

* `residuals(type = "NPDE")` of an `em()` fit draws a reading at or below the limit given the other
  readings of the same group, on the Gaussian copula of the simulations, which stay uncensored. The NPDE
  of a group are then independent under the model also with censored readings. The pde is the
  mid-rank of a reading among the simulations.

* `detectFreeCores()` inside a SLURM job takes the cores the job was given and bound to.

* Multiple shooting with `nodes = "auto"` counts a turning point of the data where the replicate
  means reverse by more than the range that their noise reaches over the series in one draw of a
  thousand. The noise is that of the replicate means, from `sigma` where the data have it, else
  from the replicates or the scatter of the series. Noise on a plateau adds no segments, and
  oscillations are also counted without `sigma`.

# dMod2 0.14.8

* `Pimpl()`: a state whose reactions all vanish at the steady state rests. It takes the value
  from the input, 0 where the input has none, and is not among the parameters.

* Alignment models compile into `align/<host>_<process>/` of the output directory, so processes
  that build the same model at once share no file.

* `plotVPC(pc = TRUE)` leaves out the readings whose population prediction is not positive.
  `confint(method = "profile")` of an `em()` fit refits with the fit's derivative mode. The plots
  of `em()` fits order their panels as the model and the data do; `plotTrace()` orders the stages as
  they ran.

* `msem(retry = TRUE)`, the default, draws a start again where the fit from it fails, for instance where the
  steady state of `Pimpl()` has no stable root at the start, up to `nTries` attempts; the table of
  starts has their `tries` and `stop`, the test that ended each trust region. `keepFits = TRUE`
  returns the fit of every start. The fits of `em()` have `stopReason`.

* `align()` takes the standard errors of a fit without random terms from the exact Hessian at the
  fit, and the fallback where the marginal Hessian is not positive definite from the expected
  information, positive semidefinite by construction.

* `compile()` reuses the objects of earlier calls in the same directory.

* The help pages of the mixed-effects functions have an example each; what takes longer than
  seconds is not run.

# dMod2 0.14.7

* The website <https://dmodverse.github.io/dMod2/> replaces the vignettes. The site guide has eleven
  chapters: setting up models, symmetries, sensitivities, optimization, profile likelihood,
  prediction uncertainty, nonlinear mixed-effects models, group-specific effects, alignment of
  relative data measurements, PEtab import and export, remote computing. The chapters on
  optimization and on mixed-effects models derive the methods of `trust()`, `mstrust()`,
  `trustL1()`, multiple shooting, FOCEI, quadrature, Laplace EM and SAEM and compare them. The
  site also has the reference by topic and the example projects as articles.

* The applications that were demos are example projects in `examples/` of the repository: each
  folder holds the scripts, the saved fits and a README of one project. New example projects `pharmacometrics`
  (theophylline, warfarin and mavoglurant, FOCEI and SAEM multistarts with the diagnostics of
  `em()` fits, nlmixr2 on the same model where installed) and `TMDD` (target-mediated drug
  disposition with the target at steady state, censored readings, FOCEI and quadrature).

* License GPL (>= 3).

# dMod2 0.14.6

* `align(reml = TRUE)`, the default: the error parameters, and the hyperparameters of normal
  random terms, by restricted likelihood. The parameters of the symbols of the prediction are
  integrated out under a flat prior as by Lindstrom and Bates: in rounds from the maximum
  likelihood fit, the log determinant of their Gauss-Newton information at the parameters and
  modes of the previous round enters the objective with an exact gradient, through FOCEI where
  the model has random terms. Fits with Laplace, spike-and-slab or Dirichlet-process terms and
  `method = "quadrature"` stay maximum likelihood; `print()` and the attribute `estimation` of
  `fit$parameters` say which, `fit$remlRounds` counts the rounds. Where the information is
  singular, the fit is by maximum likelihood, with a warning. `anova()` of alignments refuses
  restricted likelihoods of fits with different parameters of the prediction.

* `align()` goes on from the current fit after a variance is set to 0, and the results have those etas
  among `parameters`. Where the Hessian of a fit with random terms is not positive definite at
  the optimum, the standard errors come from the Gauss-Newton Schur complement over the random
  effects (`fit$seSource`).

* `alignmentModel(technical = )` names the columns that tell technical replicates apart, such as
  lanes. Readings that agree in all other columns but `value` and `lloq` are replaced by their
  mean before the alignment: geometric, on the log scale, for an error model in `prediction`,
  arithmetic otherwise; a sample with censored readings is censored at the mean of the sample readings.

* Alignment models compile into `getOption("dMod2.outdir")`, and processes that build the same
  model at once take turns.

# dMod2 0.14.5

* `em()` sets a variance at the boundary to 0: a normal random effect that forms a block alone and
  whose standard deviation ends below 0.1 is dropped, the etas of that effect fixed at 0, and the other parameters are
  fitted again; where that raises -2 log L by at most `control$boundary` (default 0.1), the fit
  without it is returned and names the hyperparameter in `boundary`, which `msem()` lists per
  start. The fit stores the random effects in `design` and the values given as `fixed`, and
  `em()` and `msem()` take a fit as `init` and go on from it.

* FOCEI ends where the trust-region model can gain less than `mtol` in the largest region the fit allows
  (`rmax`), and after five consecutive rejected trials that leave the objective flat within
  `ftol`. The inner search for the modes, of FOCEI and of the quadrature, stops where the Newton
  step is below `innerControl$steptol` (default 1e-8), so that -2 log L at given parameters does
  not depend on where the search started. FOCEI over several groupings solves the modes of all
  components in lock step, with one batched evaluation of the model per round; `summary()` prints
  the shrinkage per grouping.

* `em()` and `msem()` stop with an error where -2 log L of the data is not finite at the start,
  and FOCEI where a Hessian is not finite. A trial point where the prediction fails, or where the
  value or derivatives there are not finite, is rejected, in the search for the modes for that subject
  alone; SAEM keeps the parameters of an iteration where the objective fails at the draws.

* `em(method = "quadrature")` and Laplace EM factor the Hessian of the modes instead of inverting
  it, so that a random effect whose standard deviation collapses to 0 ends the fit. A stage of
  the quadrature ends where the Hessian is not positive definite, where a Newton step on rebuilt
  nodes would lower -2 log L by less than `tolDecrement` or than the nodes resolve, or where a
  step raised -2 log L by more than a Newton step would gain; the stage trace has the resolution.

* SAEM evaluates the data term of the draws, and the quadrature that of the nodes, in C++ where
  the model allows it, as FOCEI does.

* `msem(budget = )` caps the outer trust region of FOCEI, also the one before quadrature, at
  `budget` times the iterations of the first fit; quadrature, Laplace EM and SAEM keep their own
  limits.

* Spike-and-slab random effects have no upper bound per group: up to 12 of a group are
  integrated jointly and exactly, more each alone, coupled by the expectation consistent (EC)
  approximation or mean field (`coupling` in `control`). With sparse random effects one coupling
  holds for the whole fit, EC where it converges at the start, else mean field, reported as
  `"meanFieldFallback"`.

* An effect formula leaves a condition whose grouping is missing (`NA` in the branch table)
  without random effect, at the fixed effects.

* `residuals(fit, type = "NPDE")`: normalized prediction distribution errors from simulated data
  sets, decorrelated per group, censored readings by a uniform draw below their simulated
  fraction, with the tests of mean, variance and normality as attribute `tests`. New
  `plotVPC()`: visual predictive check of an `em()` fit, prediction corrected with `pc = TRUE`,
  stratified by a column of the condition grid, with the fraction of censored readings.
  `plotResiduals()` of an `em()` fit takes `type`, any of IWRES, CWRES and NPDE, and
  `plotRandomEffects(fit, covariates = )` plots the random effects against covariates of the
  groups. The plots of `em()` fits order the conditions naturally.

* New data sets `warfarin` (PK/PD after an oral dose, 32 subjects) and `mavoglurant` (PK after
  an intravenous infusion, 120 subjects, 78 on two occasions), from nlmixr2data.

# dMod2 0.14.4

* `profile()` returns the attribute `ends`, one row per side with the number of steps and why
  the side ended: the threshold, a limit, the step limit, or the error of the objective. It
  solves the path equation with the exact Hessian where the objective returns one with
  `deriv2 = TRUE`; `deriv2 = FALSE` keeps the Gauss-Newton matrix. `verbose = TRUE` shows the
  progress toward the threshold of the contribution named by `stop`.

* `mstrust()` warns with the first error when every fit failed.

* `plotCombined(errfn = )` draws one standard deviation of the error model as a band around the
  prediction. `plotData()`, `plotCombined()` and `plotPrediction()` order the conditions
  naturally (L2 before L10). A color legend of a single condition, and a legend of
  `plotProfile()` or `plotPaths()` with a single entry, are left out; `plotProfile()` and
  `plotPaths()` keep the order of a list of profiles.

* `symmetryReduction(allowChecked = )` reduces on directions of level `"checked"`, which hold at
  random points without a proof.

* `distributedComputing()` writes the result of a task under a temporary name and renames it
  once written, so that `check()` and `get()` read complete results only; the remote build uses
  `/tmp` where `TMPDIR` names a missing directory. `terminate()` cancels a SLURM job by job ID.

# dMod2 0.14.3

* `trust()` also ends on a rejected Newton step from inside the trust region that predicts a
  decrease below `mtol`.

* The sparse trust-region step, and the step of multiple shooting, count a numerically singular
  factor or a failed solve as not positive definite and take the norm of the step and the
  curvature along it without under- or overflow: a Hessian with flat directions gives a step on
  the boundary of the region.

# dMod2 0.14.2

* `normL2(data, g * x * p, errmodel = e)` evaluates in C++ where `x` is `Xs()` of the cppDE
  backend, `Xt()`, or `Xf()` for values alone, `g` and `e` come from `Y()`, and `p` is any
  product of `Pexpl()` and `Pimpl()`, sums per condition included, or absent: per condition the
  transformations, the steady states with their derivatives by the implicit function theorem,
  the batch solve of the model and the data term, in parallel over the conditions on
  `getOption("dMod2.cores")`. Value, gradient and Hessian, Gauss-Newton or exact, and the adjoint
  gradient of `sweep = "reverse"` equal those of the R evaluation, which every other model takes,
  as does every model under `options(dMod2.evaluation = "R")`.

* The parallel regions over conditions run on `getOption("dMod2.cores")` threads, at most the
  processors, also for fewer conditions than cores.

* The checkpoints `Xs()` keeps for the reverse sweep are limited to twice the widest batch of
  conditions, at least 4, and freed with `cppDE::releaseStore()` once they leave the cache or the
  backward solve of the native chain used them.

* A sum of objectives in which one term has a value that is not finite returns that value.

* The option `dMod2.batch.check` is removed.

# dMod2 0.14.1

* The steady states of `Pimpl()` are solved in C++ by the solver of `cppDE::ptc()`, in the
  native chain for all conditions of an evaluation at once and in parallel on
  `getOption("dMod2.cores")`, with the same result on any number of cores. A root depends on the
  parameters, the states in `pars` and the controls only, so that an objective is the same
  function of the parameters in every call; identical values are answered from memory, which
  `keepRoot = FALSE` turns off. The control `archive` and `resetWarmStarts()` are removed.

* `Pimpl(controlsPTC = list(stabilityTol = ))`: a root counts as stable while the real part of
  the leading eigenvalue on the manifold of the conserved quantities is at most `stabilityTol`,
  in the units of the rates; default 0.

* `conservedQuantities(eqnlist)` replaces `getTotals()` and `conservedQuantities(S)` on a
  stoichiometric matrix. A conserved quantity is an amount: the concentrations enter with the
  volumes of their compartments, also symbolic ones (`1*A+1*V*B`). A quantity of a single
  species is named after it. `customTotals()` checks the same, and `Pimpl()` reads the
  coefficients at the parameters of each call. A state that no reaction changes is a parameter
  of `Pimpl()`, not a conserved total.

* `Pimpl()` and the PEtab pre-equilibration keep a conserved moiety at the total where the
  moiety reactions rest, and totals keep the names of the full network.

* `exportPEtab()` takes the observables and the error model as observation functions of `Y()`.
  `exportPEtab(steadyState = )` takes the `Pimpl()` of the chain and writes the steady state as a
  pre-equilibration, conditions with equal parameters sharing one; `exportPEtab(events = )`
  writes the events at time 0 that replace a state into the condition table and refuses others.
  An observable the error model leaves out takes the standard deviations of the observable data, one
  without data the noise 1, and a fixed parameter that enters `p` as `exp()` or `10^()` has the
  value it has in the model.

* `log10Transform(states = )` transforms the states named, all of them by default.

# dMod2 0.14.0

* Breaking: `define()`, `insert()`, `repar()` and `branch(apply = )` write a substitution as
  `"lhs <- rhs"`, e.g. `insert(trafo, "x <- exp(x)", x = pars)`. A `"lhs ~ rhs"` string or a
  formula is an error that names the new form; `~` is left to model formulas.

* `evalConditionResidual()` with `df` gives the Student-t likelihood of `normT()`.

* Random effects as formulas: `insert()` on a branched transformation takes effect formulas
  such as `cl ~ wt + (1 + dose | p | ID)`. The symbol is replaced by its linear predictor on
  the scale the transformation gives it, so after `insert(trafo, "x <- exp(x)", x = "cl")`
  the random effects are normal on the log scale. Fixed effects are coded as by
  `model.matrix()` over the columns of the branch table, random terms follow lme4 with the
  IDs of brms for correlation across formulas, `gr(ID, dist = )` sets the distribution.
  `x ~ 1 + (1 | ID), x = pars` and `ka + v ~ ...` apply one right-hand side to many symbols.
  `randomEffects()` returns the structure from the transformation or from the `P()` built
  from it.

* New `normal()`, `laplace(lambda)`, `spikeSlab(prob, tau)` and
  `dirichletProcess(concentration, tau)`, the distributions of a random term. A spike-and-slab
  deviation is exactly 0 or, with probability `prob`, drawn from N(0, tau^2). Under
  `dirichletProcess()` the groups share the population value or a common deviation in
  clusters, with the Ewens prior over the partitions and N(0, tau^2) as base measure. Every
  hyperparameter left `NULL` is estimated by marginal maximum likelihood with the other
  parameters.

* New `em()` and `msem()`, mixed-effects fits as `trust()` and `mstrust()`:
  `em(obj, random = randomEffects(p), init, fixed, method)` takes the objective of the data and
  builds the prior itself. A group may own several conditions, conditions without a group
  enter the fit, random terms over several groupings, nested or crossed, are integrated per
  connected component, and normal and sparse random effects combine in one fit. Methods:
    * `"focei"`, the Laplace approximation of the marginal with exact gradient, censored data
      by M3 and random effects in the error model; Laplace, spike-and-slab and
      Dirichlet-process random effects are integrated on the Schur complement of the normal
      ones, the partitions of a Dirichlet process summed exactly (at most 14 groups), and
      sparse random effects of one group that the data couple by expectation consistent
      (`control$laplace$coupling`, or mean field as a bound). Their hyperparameters move with
      the other parameters in one trust region: the Schur complement of the fixed parameters
      takes the sparse random effects at the curvature of their posterior, the curvature of a
      hyperparameter puts its Newton step at its M-step, and a BFGS correction learns the
      missing information. Such fits have no Hessian of the marginal: with Laplace and
      spike-and-slab random effects `vcov()` inverts the empirical information of the exact
      scores of the groups, `fit$information`, which needs more groups than parameters; with
      Dirichlet-process ones, whose groups are not independent, `vcov()` is an error and
      `confint(method = "profile")` gives intervals;
    * `"quadrature"`, from the FOCEI fit, on nested sparse grids whose level per group follows
      `control$quadrature$tol`, with an importance-sampling certificate per group;
    * `"laplaceEM"`, the same nodes in an ECM;
    * `"saem"`, stochastic approximation with the mu-referenced fixed effects and Omega in
      closed form; for Dirichlet-process random effects of any number of groups by Gibbs
      sweeps over the clusters, -2 log L by sequential importance sampling with `valueSE`.

  `msem()` perturbs the starts on the scale on which each parameter is unbounded and returns
  the best fit with the table of its starts, which `plotValues()` draws as a waterfall;
  `cores` splits as in `mstrust()`. The fit from `init` runs first, and every other start ends
  after `budget` (default 5) times its outer iterations, with `converged = FALSE`.

* `em()` fits report `omega` (standard deviations and correlations with standard errors),
  `random` (the modes per group with conditional standard errors, for the sparse
  distributions the cluster and the probability of 0), `strengths`, `slab`, `logLik()` (so
  `AIC()` and `BIC()`), `residuals()` of type IRES, IWRES, CWRES, PRES or PWRES,
  `predict(type = )`, `confint(method = "profile")`, `anova()`, `simulate()`, `summary()`,
  `coef()`, `vcov()`, `plot()`, `plotResiduals()`, `plotTrace()`, `plotCombined()` (data,
  IPRED and PRED per condition) and the new generic `plotRandomEffects()`, which compares
  draws from the conditional distributions of the random effects with their estimated
  normal, so the shrinkage of the modes does not distort it (`type = "modes"` for the modes).

* New `normT()`: the Student-t likelihood with fixed degrees of freedom, with exact
  derivatives and censored data by M3, through `trust()` and every method of `em()`.

* `em()` and `msem()` take `sparse = "auto"` (default), `TRUE` or `FALSE`, as `trust()` does.

* `em()` by FOCEI over normal random effects runs its outer trust region on the sparse
  Hessian where `sparse` chooses it: the joint Gauss-Newton Hessian minus the Schur complement
  of every group plus the volume term, each local to the parameters of the group, for one
  grouping and for several. `summary()` and the standard errors of Omega solve the sparse
  Hessian for the columns they need.

* The normal prior of `em()` returns a sparse Hessian on `hessian = "sparse"`.

* FOCEI, with normal and with sparse random effects, gets the second derivatives of its volume
  term by forward over reverse (`control$derivMode = "forward-reverse"`, the default, for models
  built with `derivMode = c("forward", "reverse", "forward-reverse")`): one pass forward along
  the random effects and a few backward sweeps give the exact gradient, the mode shift and the
  derivatives of the curvature of every group, at a cost that does not grow with the square of
  the parameters. The outer trust region takes the Gauss-Newton Schur complement with the volume
  curvature of the variance parameters from these derivatives; the reported Hessian is exact
  forward over forward where the model has `deriv2`, else the exact Schur complement of the
  joint Hessian with that volume curvature, exact for models linear in the random effects, and
  `summary()` says so. `"forward-forward"` keeps the second-order sensitivities; a model with
  only one of the two takes it. In both modes, dense and sparse, the outer model of FOCEI gets
  an SR1 correction from the changes of the full outer gradient between the last three accepted
  points: it adds the curvature Gauss-Newton leaves out, negative curvature included, so fewer
  fits stop at saddles of the Gauss-Newton model.

* FOCEI in both derivative modes: the inner trust region of the random effects steps on the
  Gauss-Newton Hessian with an SR1 correction from its accepted steps where that stays positive
  definite; the outer gradient has the first-order correction for the eta gradient the inner
  solve leaves at the modes; a rejected step shrinks the radius to a quarter of its own length
  where that is shorter; the fit ends on the evaluation of its last accepted point. Forward over
  reverse answers all second-order seeds of a step in one backward sweep.

* `em()` with Laplace, spike-and-slab or Dirichlet-process random effects takes the terms of
  `obj` besides the data, such as priors and constraints; a
  term on a hyperparameter or a random effect is an error.

* Sparse random effects by forward over reverse: the marginal alone needs no second-order
  sensitivities, so models built without `deriv2` also return to expectation consistent after
  a fall back to mean field.

* `em()` says in a message when its models give second derivatives forward over forward only and
  names `derivMode = c("forward", "reverse", "forward-reverse")`; `msem()` says it once.

* Removed from the devel-EM layer: `EM()` and `msEM()` (now `em()` and `msem()`), `omega()`,
  `penaltyL1()`, `emInit()`, `updateOmegaChol()`, `makeSubjectNodes()`, `emObjfn()`,
  `sparsify()`, `scanL1()`, `gateL1()`, `selectEB()` and the vignette "Structure selection
  with an estimated penalty".

* New `alignmentModel()` and `align()`: alignment of replicate readings with free expressions
  for the prediction and its standard deviation. Every symbol is a linear predictor given by a
  formula in `model.matrix()` coding with random terms, a column with one level per value,
  numeric or not, and `I(dose)` a covariate. `trafo = c("sj <- 10^sj")` sets the scale of a
  symbol, `error = NULL` an error model without symbols, `prior` a `constraintL2()` on any
  parameter, the spread of a random term included; censored readings by M3, a Student-t option.
  As in blotIt, the experiments connected by shared values form blocks, and per block the mean
  of the predictors of the first biological symbol is held at 0. Absolute readings, marked in a
  logical column `absolute`, such as concentrations by targeted mass spectrometry, are predicted
  by the expression `absolute` in the biological symbols and data columns (the first biological
  symbol by default, units by the expression), depend on no scaling parameter and may be `NA` in
  the columns of the scaling formulas; their block is not normalized but comes out in their
  units. Every experiment within its connected component of the groups of the random terms is a
  condition of its own, so the Hessian is sparse over them; designs and coefficients are sparse
  matrices and the standard errors of the values come from a sparse solve. `align()` starts the
  coefficients of the first scaling symbol where its predictor is, per experiment, the median of
  the predictors at which the prediction meets the relative readings, runs the fits of a
  multistart in parallel on `cores`, one core each, and stops at tolerances of
  `sqrt(.Machine$double.eps)`, inner fits included. Methods `print()`, `as.datalist()`,
  `plotValues()` (the starts by rank), `plotResiduals()`, `anova()` and `simulate()`, plots
  `plotAligned()` and `plotScaled()`. The demo "alignment" aligns qPCR experiments with
  undetected wells, and Western blots with a random background per gel together with targeted
  mass spectrometry, which puts the values in nM.

* `alignmentModel()` takes random terms of several groupings, nested or crossed, such as a lane
  effect per gel and lane shared by the antibodies of a gel beside a background per gel and
  antibody. They are normal and integrated per connected component; a condition is an experiment
  within its component, and the normalization blocks follow the fixed parts of the scaling
  formulas. `fit$omega` and `fit$random` have a column `grouping`, and `plotRandomEffects()`
  draws a panel per grouping and random effect.

* The objective of an alignment answers all its conditions by one prediction and one error-model
  function: the compiled model is evaluated once at the readings of all conditions, and the
  derivatives of every condition in its parameters and the contractions of the backward sweep
  run in C++, in parallel over the conditions. The error model of an alignment is taken from the
  prediction, which evaluates it alongside. Its backward sweeps answer several second-order
  seeds at once, so FOCEI fits it by forward over reverse.

* An alignment holds the evaluation that the backward sweeps of prediction and error model share
  only within one sweep, so a fitted alignment keeps none of it; an error model that depends on
  no symbol is not evaluated backwards.

* The reverse mode: the Hessian of `normL2()` with `sweep = "reverse"` has the curvature of an
  estimated error model, an error model that is a sum of per-condition ones goes backwards, and
  so does an observation function without parameters of its own at second order. `Y()`, `P()`
  and `Xs()` answer several second-order seeds in one backward sweep.

* A prediction or observation function whose kernel has a `condition` argument gets the
  condition of every call, its batch entries get `conditions`; an observation function
  answers the backward sweep through its batch entry where it has one. `normL2()` with `times`
  as a list evaluates a function of several conditions at the times of each, and its backward
  sweep sums the cotangents of the conditions in one pass. An error model whose output the
  prediction holds is taken from there.

* A cppDE model keeps the batch handles of its last eight request sets.

* A sum of objectives holds its terms, not the frame it is built in.

* `distributedComputing(build = "compile")` builds a model of eight or more C++ sources that
  lie outside the working directory, such as in `getOption("dMod2.outdir")`.

# dMod2 0.13.4

* `evalConditionResidual()` has the parameters that enter only the error model in its
  gradient and Hessian.

* An objective called with `hessian = "sparse"` returns its Gauss-Newton Hessian as a sparse
  `Matrix` assembled from the blocks of the conditions; `trust(sparse = TRUE)` factorises it by
  sparse Cholesky and `vcov(which = )` solves it for the columns named. `datapointL2()` returns
  a sparse Hessian there too.

* `trust()` and `profile()` take `sparse = "auto"` (default), `TRUE` or `FALSE`. `"auto"` takes
  the sparse Gauss-Newton Hessian for a minimisation without phases, bounds or exact Hessian
  when there are at least 500 free parameters and the Hessian at the start fills at most 5% of
  its matrix (`options(dMod2.sparseAuto = )`), and starts from that evaluation.
  `trust(sparse = TRUE)` on an objective that returns a dense Hessian is an error that names the
  terms of its sum that do. `mstrust()` passes `sparse` to `trust()`.

* `profile()` on a sparse Hessian decides by the pivots of its LU whether the path equation is
  singular, by the rule of `vcov()`: a pivot at most 1e-14 times the largest. A singular step
  then takes the pseudo-inverse and reoptimises as on a dense Hessian, up to 5000 parameters;
  above, the error names the non-identifiable parameters. The pseudo-inverse step keeps all its
  components. `algoControl$W = "identity"` moves the profiled parameter alone and reoptimises
  the others at every step, without a matrix in the path equation and without a warning on
  singularity.

* Requires cppDE 0.13.4: a second-order solve whose parameters have no second
  derivatives crashed R with an older cppDE.

# dMod2 0.13.3

* `distributedComputing()` keeps the manifest of a job also in its folder on every remote
  node. `recover = TRUE` with `on` given takes it from there, so a job submitted from
  another machine is recovered, and a job of the same name that still runs on the nodes
  is an error also when it was submitted from another machine.

* `distributedComputing()` runs on Windows, with the local commands in the sh of Rtools.
  `detectFreeCores()` measures the processor load of the local Windows machine, without
  a warning.

* `compile()` compiles in parallel on Windows too, with or without a combined `output`,
  and `P()` and `Y()` pass their `cores` to it there.

# dMod2 0.13.2

* `distributedComputing()` reports the start of a job as a message, on a SLURM cluster
  with the job ID and the array tasks. A build writes its output to `build.log` in the job
  folder, and a failed build is an error that shows its end instead of starting the tasks.

* `plotProfile()` leaves more space between its panels, and a narrow panel leaves out the
  x labels that would overlap.

# dMod2 0.13.1

* Breaking: `distributedComputing()` takes where its tasks run as `on`: `localNode()`,
  `sshNodes()` or `slurmCluster()`, each with the `cores` of a task. A node is one worker,
  and a SLURM job has `array` workers, its array tasks; the tasks are dealt to the workers
  in turn. `slurmCluster()` takes `walltime`, `memPerCore`, the job script and
  `partition`, `account`, `maxRunning` and further options of `sbatch`. `build =
  c("auto", "compile", "link", "copy")` replaces `compile` and `link`, and `"auto"`
  compiles on remote nodes. The arguments `machine`, `slurm`, `jobscript`, `memPerCore`,
  `walltime`, `compile` and `link` are gone, and `slurmTemplate()` names its first
  argument `login`.

* `distributedComputing()` ships the model files of the objects in `input` instead of
  every file of `dMod2.outdir`. A job of the same name that still runs where the new one
  would run is an error. The manifest is written before the start, and every job is
  registered until `purge(local = TRUE)`; `dMod2Jobs()` lists them, and `recover` finds a
  job from any working directory. The job has `wait()` and `log(i)`, and `get()` joins
  parlists with `c()` and data frames with `rbind()`, leaving out tasks without a result
  with a warning.

* `rbind()` of parameter frames keeps their attributes and joins the parts of a profile,
  such as its two sides, into one. `confint()` sorts the right side of a profile as it
  does the left.

* Breaking: `trust()` takes `boundary` as an argument instead of `stepControl$boundary`.
  `trustL1()` takes its settings as `trust()` does: `ftol`, `mtol`, `gtol`, `xtol` and
  `rmin` in `tolControl`, `thetaMax` (was `theta.max`) in `stepControl`, and `parupper`,
  `parlower` and `boundary` after `iterlim`; a flat name is refused with the list it moved
  to, and `iterlim` defaults to 1000.

* The `check()` of a `distributedComputing()` job counts a task whose process remains as a
  zombie as finished, as in a container whose init does not reap.

* `distributedComputing(slurm = TRUE)` reads the job ID of `sbatch` also when the build on
  the login node prints before it.

* `plotCombined(aesthetics = )` replaces a default aesthetic also under its other spelling,
  such as `colour` for `color`.

* `wide2long()` of a matrix builds its data frame without a row-name warning.

# dMod2 0.13.0

* Breaking: `steadyStates()` has one solver, the former version `"1.4"`: every unknown is
  solved from a combination of balances, so every result is a ratio of positive sums. The
  versions `"1.0"` to `"1.3"` and the arguments `version`, `sparsifyLevel`, `outputFormat`,
  `positive` and `branches` are gone; `testSteady` is `"fast"` or `"skip"`.

* Breaking: the options are named `dMod2.*`: `dMod2.cores`, `dMod2.outdir`,
  `dMod2.batch.check` and `dMod2.compile.cmdlimit`; the `dMod.*` names are no longer read.
  `dMod.quiet` and `dMod.sym.<name>` are gone: `suppressPackageStartupMessages()` silences
  the startup message, and `SYMIDENT_<NAME>` sets a switch of symident. `?dMod2` lists the
  options and environment variables.

* The solver options of `Xs()` and `Xf()` on the `cppDE` and `Sundials` backends take
  `maxOrder`, an upper bound on the order of the multistep methods, and `stald`, stability
  limit detection for `"bdf"`.

* Breaking: `distributedComputing()` is rebuilt and takes over `runbg()`; `runbg()` and
  `profileParsPerNode()` are removed. The code runs once per task: `tasks` is a number, a
  list or a data frame, and the code reads its element or row as `task` and its position as
  `taskId`, and the cores it may use as `cores`. Without `slurm` the tasks run as detached
  background processes on `"localhost"`, without ssh, or on machines reached by ssh, one
  after the other per machine and each with the machine's cores: those `machine` names, as
  in `c(server1 = 20, server2 = 10)`, or its free cores, and never more than are free. With
  `slurm = TRUE` they form an array job with the cores `machine` names, submitted through
  the job script template of the login node, with `{{ }}` fields. `input` is taken from the
  calling environment, every task gets the seed `seed + taskId`, `recover` needs only the
  `jobname`, and the job gains `terminate()`; `check()` returns `TRUE` once no task runs any
  more and warns about failed tasks. `cores`, `nRep`, `varValues`, `partition`, `nodes`,
  `sshPasswd`, `customFolders`, `returnAll`, `resetSeeds` and `purgeLocal` are gone.
  `demo("BA_transport")` sends its multistarts this way, and
  `examples/distributedComputing.md` describes it.

* `slurmTemplate()` sets up the job script template of a SLURM cluster once per login node,
  kept across sessions, from a job script of the cluster or from the one dMod2 ships for the
  bwForCluster Helix.

* `mstrust()` draws every random start, those of retries included, before it spreads the
  fits over `cores`, so `set.seed()` reproduces a run whatever `cores`.

* `unique()` on a parframe returns the rows with distinct values sorted by value, also for a
  parframe not sorted by value; rows are selected with `subset()` first. dMod2 needs
  ggplot2 3.5.0.

* `rref()` no longer takes `verbose`, which it did not use.

* `detectFreeCores()` names its values by machine and reaches `"localhost"` without ssh.

* `symmetryDetection(unknownInputs = )` declares unknown functions of time; an exact
  direction lists the change of each input as `inputChange`. Needs symident 0.12.0.

* Loading dMod2 no longer requires symident in the Python environment; symmetryDetection()
  and symmetryReduction() require it on first use, so a symident that cannot be resolved
  leaves the code generation of models intact.

* Breaking: `getElements()`, `blockdiagSymb()`, `getLinVars()`, `long2wide()`, `lsdMod()`,
  `compare()`, `forcingsSymb()` and `dot()` are removed, as is the `character` method of
  `reduceReplicates()`; `reduceReplicates(read.csv(file))` replaces it, and openxlsx and purrr
  are no longer suggested. `progressBar()` is no longer exported.

* Breaking: a dMod function refuses a named argument it does not take, more unnamed
  arguments than it has inputs and a missing input, so `x(times, pars, condition = "C1")`
  is an error instead of a prediction for all conditions. `Xs()` and `Xf()` refuse
  arguments they do not take.

* Breaking: the subset expressions of `plotPars()`, `plotValues()` and `plotArray()` come
  right after the data, `tol`, `direction`, `covtable` and `nSim` after them; `covtable`
  defaults to `NULL`, and the subset applies with and without it. `plotArray()` draws
  `nSim` trajectories evenly spaced from the optimum, at most one per profile point,
  grouped by parameter set.

* `Xd()` applies the chain rule of a parameter transformation and differentiates by the
  free parameters only; a missing grid parameter is an error.

* `predict()` on a parframe labels each prediction with its own row, also for a parframe
  not sorted by value, and adds a condition grid with one covariate; `plotResiduals()`
  uses the row of each fit.

* `addReaction()` reads any numeric coefficient, such as `3/2*A`, `(2)*A` or `2*3*A`, and
  refuses a term that is no `coefficient*species`, a reaction without species and a state
  named `Rate` or `Description`. `getReactions()` writes the coefficients as they are.
  `subset()` of an eqnlist takes `%in%` on a vector as base R does.

* `rref()` sets entries below `tol` to zero and returns the pivots as an integer vector;
  `nullZ()` passes `tol` on and returns the identity for a matrix without pivots.
  Conserved quantities have no rounding residues.

* `lloq = NA` in data means no limit of quantification.

* `as.datalist()` takes data without a condition or covariate column, as one condition
  `"1"`, tibbles and data.tables. Renaming a datalist renames its `condition` column,
  `subset()` of a datalist keeps the condition grid, and `covariates()` of a data frame
  returns a data frame for any number of covariates.

* `reduceReplicates()` names conditions by the levels of factor columns.

* `profile()` checks `whichPar`, `stepControl` and `algoControl`, names a failed profile
  with its error, reports an objective error that ended a side and singular Hessians once
  per profile, and returns the optimum alone when no step is accepted. Each path step
  solves the linear system directly and takes the pseudoinverse only for an
  ill-conditioned one.

* `confint()` takes a profile without a `data` column through `val.column = "value"`.

* `mstrust()` takes an anonymous `samplefun` or one with a long name, passes the arguments
  of an optimiser other than `trust()`, refuses
  `fits < 1`, and writes one log entry per fit also from PSOCK workers; `load.parlist()`
  reads the fits in the order of their numbers. `mstrust()` and `profile()` share one
  parallel backend.

* A scaled objective `a %.*% obj` keeps an absent gradient or Hessian absent and needs a
  finite scalar `a`.

* On the `deSolve` backend, `Xs()` and `Xf()` return the requested times only, warn when
  the solver stops early, name missing parameters, and `Xs()` checks `names` and refuses
  an odemodel without sensitivities.

* `[` on a parframe keeps a single column as a parframe and selects columns by a single
  index, as for a data frame; `as.parvec()` of a parframe checks `index`; `parframe()`
  takes as parameters the columns that are neither meta names nor objective attributes.
  `summary()` and `plot()` of a parlist skip aborted fits, `summary()` also takes an empty
  parlist and fits without convergence flag, and `c()` of parlists concatenates.

* `eqnvec()`, `parvec()` and `prdlist()` name unnamed entries by their position; a parvec
  with `NA` prints; `wide2long()` refuses duplicated names.

* `plotPaths()` and `plotPathsMulti()` build their paths in one pass; `plotPathsMulti()`
  normalizes each path to the range -1 to 1 and colors the `nPars` largest paths, the
  others in "Others". `plotPars()` detects steps on a frame in any order and names the
  column `step`. A layer added to `plotCombined()` uses its data. `coordTransform()`
  propagates errors only where `sigma` is given. `plotFluxes()` takes any number of
  fluxes.

* `wide2long()`, `lbind()` and `combine()` bind through data.table, which speeds up the
  data frames of predictions with many conditions.

* `write.eqnlist()` and `as.eventlist()` handle empty lists; `is.eqnlist()` checks that
  the columns of the stoichiometric matrix are the states.

* `Pexpl()` marks an output as fixed only when it depends on no varying input, and
  `Pimpl()` keeps the rows of every solved state, so an output whose first derivative
  vanishes at the evaluation point keeps its second derivative. `Pimpl()` treats a species
  that no reaction changes, such as a catalyst, as a constant state, not as zero.

* In a sum of objectives, `datapointL2()` finds the predictions of every `normL2()` term,
  also behind a prior or a term that `conditions` skips. `obj + NULL` returns `obj`; adding
  anything else that is not an objfn is an error.

* `constraintL2()` and the other priors apply a named `sigma` to the parameter it names and
  refuse a `sigma` without an entry for a parameter of `mu`. `trustL1(minimize = FALSE)`
  subtracts the penalty from the objective it maximises. `evalConditionResidual()` passes
  `fixed` parameters to the error model also without sensitivities.

* `odemodel(deriv = FALSE)` builds no derivatives, whatever `derivMode`. `define()` and
  `insert()` refuse `conditionMatch` for a trafo without a condition tree.

* `eqnlist()` takes `compartments` without `compartmentOf`, with every state in
  `defaultComp`, and refuses `compartmentOf` without its compartments. `c()` of eqnlists
  keeps the layout of an eqnlist without reactions and the custom totals, and `subset()`
  returns an empty eqnlist when no reaction matches. `log10Transform()` keeps names that
  sympy reads as constants, such as `E` and `I`. `format.eqnvec()` keeps powers and
  numbers such as `1e-3` whole. `as.eqnvec()` refuses empty and duplicated names.

* `lbind()` and `as.data.frame()` of a datalist fill a column that some conditions lack
  with `NA`. `c()` of datalists returns class `list` as well and refuses a condition named
  twice. `as.datalist()` refuses covariate values whose condition names coincide.

* `subset()` of a parframe and the subsets of `plotPars()`, `plotValues()` and
  `plotProfile()` take several expressions and see the frame of the caller. `as.parframe()`
  of a parlist fills a parameter that some fits lack with `NA`, and `c()` appends a single
  fit to a parlist. `plotProfile(parlist = )` draws the profiled parameters only,
  `plotProfilesAndPaths()` draws the total profile, `plot()` of a prdframe takes a data
  frame, and `coordTransform()` keeps the order of the rows.

* `modelname<-` assigns a vector of names condition by condition and takes `conditions`;
  `getParameters(conditions = )` of a function for every condition returns its
  parameters.

* The code that dMod 1 contributors other than the authors wrote is reimplemented from its
  documented behaviour: the print methods of eqnvec, eqnlist, objlist and parvec,
  `summary()`, `plot()` and `as.parframe()` of a parlist, `fitErrorModel()`,
  `reduceReplicates()`, the argument routing and fit log of `mstrust()`, `Xd()`,
  `define()`, `insert()`, `conservedQuantities()`, the datalist methods and the prediction,
  data, flux, residual and array plots. With the old `steadyStates()` versions gone,
  `DESCRIPTION` lists the two authors only.

* `Xd()` takes a state with a single grid point as constant; `plotFluxes()` labels an
  unnamed flux by its expression also when others are named; `eqnlist()` takes a plain
  vector as a one-column stoichiometric matrix; `reduceReplicates()` groups by the exact
  times instead of times rounded to seven digits.

* In a nested sum of objectives, a term without conditions is evaluated for any requested
  conditions, as in a flat sum.

* `summary()` of a parlist heads its table "Best and worst fit"; the messages of
  `as.parvec()` on an unconverged fit, `as.eqnvec()` on unparsable equations and
  `steadyStates()` on named `givenCQs` are reworded.

# dMod2 0.12.7

* `loadDLL()` also loads the shared objects of a model's auxiliary sources, such as the
  influx function of `Pimpl()`, so `mstrust()` and `profile()` with a condition axis run
  on PSOCK workers also for models compiled source by source.
* `trust()` and `trustL1()` append the objective's error message to
  `parinit not feasible: objfun failed`.
* `define()`, `insert()` and `repar()` on `trafo = NULL` map each left-hand side to its
  right-hand side, not to itself.

# dMod2 0.12.6

* `controls()` on an objective also reaches the functions it evaluates: the prediction and
  error model of each `normL2()` term and the parameter transformation of an
  `objfn * parfn`. Solver tolerances and other controls can be changed through the
  objective alone. `condition` restricts the lookup as for an `fn`.
* An unknown entry of `options`, `optionsSens` or `optionsReverse` of `Xs()` and `Xf()` is
  an error, also when set later with `controls()`.
* `Y()` names a state the prediction does not return, e.g. after the control `names` of
  `Xs()` left it out.
* A failed batch solve names the conditions that did not complete.
* `plotValues()` shows unconverged fits in red with a legend entry of their own instead of
  a caption.

# dMod2 0.12.5

* `plotValues()` colors converged fits by their iterations on the dMod ramp and draws
  unconverged ones in red, with a caption only when there are any. With
  `pseudoLog = TRUE`, a range within one decade gets breaks at 1-2-5 or finer.
* The demo `KokMSB2020` shows the model scheme of the publication and the statistics of
  200 fits, which ship in `inst/extdata/kok/`.

# dMod2 0.12.4

* `symmetryDetection()` and `symmetryReduction()` need symident 0.11.0. A direction
  checked without an error bound has level `"checked"`, and `info$assumptions` lists the
  declared domain and the declared starts of the switches, with the error bound.
* `symmetryReduction()` names the claims of a section that are not proven under
  `blocks[[i]]$unproven`, with level `"open"`, and lists `assumptions`; the level `"assumed"`,
  `blocks[[i]]$assumptions` and the argument `verbose` are removed. A section on a face
  holds `faceWhy` and `faceLevel`; a resting state checked only at random points of the face
  makes the block `"checked"`.
* In the blocks of `symmetryReduction()`, `newParameters`, `solvedDomain` and `choices`
  replace `survivorMeaning`, `carrierDomain` and `admissible`.
* `trust()` takes `phases`, a list of phases, each with its own `hessianMethod`,
  `iterlim`, `sweep` and `qnControl`, run in turn and cycled through for up to
  `switches` handovers. A phase that uses up its own `iterlim` hands over like a
  soft stop. `fallbackControl` is removed:
  `fallbackControl = list(hessianMethod = "sr1", limit = 3, iterlim = 200)` is
  `phases = list(list(hessianMethod = "gn"), list(hessianMethod = "sr1", iterlim = 200)), switches = 3`.
  A run without `phases` is unchanged.
* The result of `trust()` has `phaseIterations`, the iterations of each phase run
  in run order, and with `phases` the sweeps of the phases as `phaseSweep` in place
  of `fallbackSweep`. `as.parframe()` adds a column `iterations.<phase>` per phase.
* Each run of a quasi-Newton phase has its own allowance of reseeds by
  `qnControl$hessianReseed = "stall"`; `nReseed` counts over the run.
* `mstrust(cores = c(fits = , conditions = ))` runs both axes as given, each
  capped at the cores available, and records that layout in `attr(, "run")$cores`.
* `mstrust()` and `profile()` stop the cluster of a conditions axis and restore
  `options(dMod.cores)` on exit, on every platform.
* The cores available are those of the R process, `parallelly::availableCores()`,
  which follows the CPU affinity and a SLURM or PBS allocation. dMod2 imports
  parallelly.
* A sum of objective functions keeps the `compileInfo` of its terms, so `loadDLL()`
  on it searches the directories its sources were compiled in.

# dMod2 0.12.3

* `trustL1()` drops `parscale`, as `trust()` did in 0.12.2. A scaling of the
  parameters belongs in `P()`.

# dMod2 0.12.2

* `trust()` gains `fallbackControl`, the settings of the fallback phase:
  `hessianMethod` (default `"none"`), `limit`, `iterlim` (per fallback
  phase), `sweep` and `qnControl`. With `sweep = NULL`, the default, the
  fallback phase evaluates the objective backwards when the chain can.
  `hessianFallback`, `fallbackLimit` and `parscale` are removed.
* `trust()` defaults to `iterlim = 1000`, counted over all phases.
* The deprecated argument names and aliases are removed, among them
  `attr.name`, `opt.BLOQ`, `obj.attributes`, `sort.by`, `split.by`,
  `keep.covariates`, `rank.tol`, `attach.input`, `keep.root`, `optionsOde`,
  `start1stfromCenter`, `fterm`, `mterm`, `stepControl$theta.max`, the
  `controlsPTC` entries `rtol` and `atol`, the `cores` of `normL2()` and
  `readPetabYaml()`, `readPetabTables()`. `repar(expr, trafo)` is an error.
* `mstrust()` drops `hessian` from every fit and records `elapsed`, the wall
  time of each fit, and the attribute `"run"` with the start time, total wall
  time, cores, optimiser, sampler, their arguments and the package versions.
  `as.parframe()` adds the column `elapsed` and keeps `"run"`.
* `vcov()` takes an evaluation of an objective function, `vcov(obj(pars))`,
  with `fixed` for parameters held at a bound. `vcov()` on a fit without
  Hessian is an error instead of a matrix of `NA`.
* The `.R` files of `mstrust(cautiousMode = TRUE)` restore the fits exactly.

# dMod2 0.12.1

* Needs symident 0.10.1.
* When symident needs msolve and finds none, `symmetryDetection()` and
  `symmetryReduction()` ask in an interactive session whether to build it into the
  symident cache, on Windows with Rtools, and then repeat the call.

# dMod2 0.12.0

* Needs cppDE 0.13.0 and symident 0.10.0.
* `getSymbols()` and `replaceSymbols()` are dMod2's own, no longer re-exported
  from cOde, which is now needed for the deSolve backend only.
  `replaceSymbols()` brackets a replacement that is not a single operand,
  replaces variables only, not function or argument names, and keeps the
  spacing.
* `pi` is the number in equations, as in cppDE, and `getSymbols()` does not
  list it.
* `symmetryDetection()` gains `regimes`, the regimes an experiment visits in a
  model with switches, overall or per condition.
* A direction of `symmetryDetection()` has a `level` (`"exact"`, `"verified"`
  or `"open"`) in place of `certified`, and `stateComponents`. `info` gains
  `rankStatus`, `proof`, `bound`, `gapCertified`, `limits`, `errors` and
  `branches`.
* `reconstControl()` gains `homotopyDegree`, `degreeCapMax`, `exactProof`,
  `exactWorkCap` and `seed`; `NULL`, the default, keeps the value of symident.
* `symmetryReduction()` stops at a direction of level `"open"` unless
  `allowOpen = TRUE`. It gains `workCap`, `timeout` defaults to `Inf`, and
  `check = TRUE` runs the detection again on the reduced model. The result
  gains `level`, `limits` and `checks`, and each block its `level`.
* `plotValues()` draws the objective values on a linear scale, not their distance to
  the best fit on a pseudo-log scale; `pseudoLog = TRUE` switches to a pseudo-log10 scale.
* The data set `kok` and the demo `KokMSB2020`. The Bachmann, Boehm, bile acid and
  Hodgkin-Huxley scripts are demos as well, `demo(package = "dMod2")` lists them; the
  tutorials left in `inst/examples` from dMod and the R2C templates in `inst/templates`
  are removed, and the gradient comparisons moved to `bench/gradient`.
* `bench/` and `benchmarks/` are one directory, `bench/`, by topic: `harness/` (the PEtab
  benchmark suite), `gradient/`, `hessian/`, `symmetry/`, `multipleShooting/` and
  `numSteady/`.

# dMod2 0.11.4

* symident and python-libsbml are found in `~/.virtualenvs/r-reticulate` too,
  installed on first use.
* `summary(x, detailed = TRUE)` on a `symmetryReduction()` result explains each
  section. Needs symident 0.9.3.

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
