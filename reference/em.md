# Mixed-effects fit by marginal maximum likelihood

Fits a mixed-effects model: the random effects that
[`insert()`](https://dmodverse.github.io/dMod2/reference/define.md) with
effect formulas put into the parameter transformation are integrated
out, and the fixed parameters and the hyperparameters of the random
effects are estimated by maximizing the marginal likelihood.

A group, a level of the grouping of a random term such as a patient or a
cell line, may own several conditions, which share the random effects of
the group. Normal, spike-and-slab and Dirichlet-process random effects
of one grouping may be combined. Random terms over several groupings,
nested or crossed, are normal and integrated per connected component of
the conditions.

Methods:

- `"focei"`:

  Laplace approximation at the conditional modes with the exact
  derivative of the volume term. The default, and the method for several
  groupings; also for spike-and-slab and Dirichlet-process random
  effects (at most 14 groups).

- `"quadrature"`:

  adaptive sparse-grid quadrature around the modes, started from the
  FOCEI fit; section Quadrature.

- `"laplaceEM"`:

  ECM on the Laplace approximation with the covariance-corrected update
  of `Omega`.

- `"saem"`:

  stochastic-approximation EM, also for spike-and-slab and
  Dirichlet-process random effects; section SAEM.

## Usage

``` r
em(
  obj,
  random,
  init = NULL,
  fixed = NULL,
  method = c("focei", "quadrature", "laplaceEM", "saem"),
  control = list(),
  cores = 1L,
  verbose = TRUE,
  sparse = "auto"
)
```

## Arguments

- obj:

  Objective of the data,
  [`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md) or
  [`normT()`](https://dmodverse.github.io/dMod2/reference/normT.md) of
  the model `g * x * p`, without terms for the random effects.

- random:

  The random effects,
  [`randomEffects()`](https://dmodverse.github.io/dMod2/reference/randomEffects.md)
  of the parameter transformation or of `p`.

- init:

  Named numeric start with the parameters of `obj` that are neither
  random effects nor created by the effect formulas (those start at 0).
  It may also set those and the hyperparameters: `omega_<a>_<a>`, the
  log of the diagonal entry of the Cholesky factor of `Omega` (the log
  standard deviation of an uncorrelated random effect), and
  `omega_<a>_<b>`, the entry in row `a` and column `b` below the
  diagonal, with `a` and `b` the stems of the random effects without
  `eta_`, as in `omega_cl_cl`; `prob_*` and `tau_*` (spike and slab),
  `concentration_*` and `tau_*` (Dirichlet process), named
  `<kind>_<grouping>` and numbered `_1`, `_2`, ... where a grouping has
  several blocks. An `em` fit as `init` gives the start: the estimates,
  `fixed` values and random effects of that fit, without those it set to
  0, and with sparse random effects the inner modes.

- fixed:

  Named numeric of parameters held fixed.

- method:

  `"focei"` (default), `"quadrature"`, `"laplaceEM"` or `"saem"`.

- control:

  List of options; other entries are an error.

  `sd`

  :   start of the standard deviations of the random effects, default
      0.3.

  `inner`

  :   `"nested"` (default) solves the modes at every evaluation of the
      outer search; `"sqp"` takes `sqpIterlim` iterations of the inner
      search from the modes predicted along d eta / d theta, with value
      and gradient corrected to first order for the eta gradient left
      over, until the outer search stops, and then goes on nested from
      there, so that the estimates are those of the nested search.
      Applies to FOCEI and to spike-and-slab and Dirichlet-process
      random effects.

  `innerControl`

  :   the trust region of the modes: `rinit` and `rmax`, the initial and
      largest radius (default 1 and 10); `iterlim`, default 30 and over
      several groupings 100; `steptol`, default 1e-8, the Newton step at
      which it stops; `mtol`, default 1e-12, the predicted change of -2
      log L, relative to -2 log L, at which a step on the boundary of
      the region stops it; `eigenFloorRelative`, default 1e-10, the
      floor of the eigenvalues of the inner Hessian relative to their
      mean (one grouping); `sqpIterlim`, default 2, the inner iterations
      per evaluation with `inner = "sqp"`. With spike-and-slab and
      Dirichlet-process random effects only `sqpIterlim` applies. Other
      entries are an error.

  `trustControl`

  :   the outer trust region; it stops where the predicted change of -2
      log L is below `mtol`, default 1e-5. Takes `ftol`, `mtol`, `gtol`
      and `xtol` flat or in `tolControl`. With random effects over
      several groupings the outer search is
      [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md),
      and every argument of
      [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md)
      passes, `hessianMethod`, `phases` and `qnControl` included; FOCEI
      of one grouping reads `rinit`, `rmax`, `iterlim`, `ftol` and
      `mtol`. Other entries are an error.

  `secondOrderCorrection`

  :   the derivative of the volume term, default `TRUE`.

  `derivMode`

  :   how FOCEI gets the second derivatives of the volume term:
      `"forward-reverse"` (default), sweeps for models built with
      `derivMode = c("reverse", "forward-reverse")`, or
      `"forward-forward"`, the second-order sensitivities of models
      built with `derivMode = "forward-forward"`. Without the mode asked
      for the fit takes the other, with a message.

  `boundary`

  :   the largest rise of -2 log L at which a variance is set to 0,
      default 0.1, `FALSE` for none. With one grouping and several
      blocks, a normal random effect that forms a block alone and whose
      standard deviation ends below 0.1 is dropped and the fit repeated
      from the estimates; where -2 log L rises by at most `boundary`,
      that fit is returned and names the hyperparameter in `boundary`.

  `quadrature`, `saem`, `sparseFit`

  :   lists of the sections Quadrature, SAEM and Sparse random effects;
      other entries are an error.

- cores:

  Number of cores over the groups and conditions inside the fit; with 1
  (default) no threads are started.

- verbose:

  Logical; print the progress of the iterative methods. Default `TRUE`.

- sparse:

  `"auto"` (default), `TRUE` or `FALSE`: the outer trust region of FOCEI
  over normal random effects on the sparse Gauss-Newton Hessian.
  `"auto"` decides as
  [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md)
  does.

## Value

An object of class `em` with, among others:

- `argument`, `value`, `gradient`, `hessian`:

  the estimates, -2 log L of the marginal likelihood with gradient and
  Hessian (a sparse `Matrix` on the sparse path, `NULL` with
  spike-and-slab or Dirichlet-process random effects); on the sparse
  path `lowrank`, the part of the Hessian kept apart (see
  [`objlist()`](https://dmodverse.github.io/dMod2/reference/objlist.md)).

- `information`:

  with spike-and-slab random effects the empirical information that
  gives the standard errors: `G / (G - 1)` times the sum of the outer
  products of the centered scores of the `G` groups, of rank at most
  `G - 1`, plus half the curvature of the terms without random effects
  (conditions outside every group, priors) by their expected information
  (`hessian = "expected"`, see
  [dModfn](https://dmodverse.github.io/dMod2/reference/dModfn.md)),
  positive semidefinite. It is not the Hessian of the marginal -2 log L
  and is singular where the scores and those terms span fewer directions
  than there are parameters.
  [vcov()](https://dmodverse.github.io/dMod2/reference/vcov.em.md)
  inverts it where it is positive definite.

- `converged`, `iterations`:

  the convergence of the outer search and the iterations; with several
  groupings also `phaseIterations` of
  [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md).

- `innerIterations`, `sqpIterations`:

  with FOCEI and sparse random effects the iterations of the inner
  searches summed over the groups and evaluations, and with
  `control$inner = "sqp"` the outer iterations on inexact modes.

- `Omega`, `omega`:

  the covariance of the normal random effects, their standard deviations
  and correlations with standard errors.

- `random`:

  per group the modes with conditional standard errors; for
  spike-and-slab and Dirichlet-process random effects also the cluster,
  0 for the population value, and the posterior probability `p0` of a
  deviation of 0; a column `data` where a sparse random effect reaches
  no data of the group; such an effect rests at 0 with the prior.

- `slab`:

  the hyperparameters of the spike-and-slab and Dirichlet-process
  blocks.

- `valueSE`:

  the standard error of `value` where it comes from importance sampling.

- `coupling`:

  `"expectationConsistent"` where a group needed it, else `NULL`.

- `boundary`:

  the hyperparameters set to 0, see `control`.

- `etaModes`, `etaSE`, `shrinkage`:

- `certificate`:

  with method `"quadrature"` per group the importance-sampling check at
  the estimate, section Quadrature.

Methods: `print`, `summary`, `coef`, `vcov`,
[confint](https://dmodverse.github.io/dMod2/reference/confint.em.md),
[logLik](https://dmodverse.github.io/dMod2/reference/logLik.em.md),
[residuals](https://dmodverse.github.io/dMod2/reference/residuals.em.md),
[predict](https://dmodverse.github.io/dMod2/reference/predict.em.md),
`plot`.

## Quadrature

Method `"quadrature"` alternates between placing the nodes around the
modes at the current parameters and one trust step over every parameter
on these frozen nodes. The nodes of the mu-referenced fixed effects are
frozen in `phi = X beta + eta`. Entries of `control$quadrature`:

- `tol`:

  allowed change of -2 log L by the next level, default 0.1, split over
  the groups; each group takes the lowest level from the current stage
  of the group up within `tol / G` and never a lower one than before.
  Where the next level exceeds `maxNodes`, the change from the level
  below is the estimate. `Inf` keeps the levels of the stages.

- `levels`:

  the stages, default `K + 1` with `K` the number of random effects of a
  group; level `K + s` integrates polynomials of total degree `2 s + 1`
  exactly.

- `maxLevel`, `maxNodes`:

  highest level, default `K + 8`, and the largest node count of a level
  per group, default 1000. A group whose error exceeds `tol / G` at the
  highest level within both is capped. Column `capped` of `stageTrace`
  counts the capped groups; attribute `"levels"` of `stageTrace` and the
  `certificate` list level, level error and cap per group.

- `tolDecrement`:

  a stage ends when a Newton step on the rebuilt nodes would lower -2
  log L by less than this, default 1e-4 (`stopReason` `"decrement"`), or
  when -2 log L rose after a step on nodes of unchanged levels, so the
  step lies below what the levels resolve (`"resolution"`). `converged`
  means stationary at the levels reached, capped ones included; the fit
  returns the iterate of lowest -2 log L on the highest levels.

- `maxEcmPerStage`, `cm1Control`, `maxCm1Iter`:

  rebuilds per stage (default 10, `stopReason` `"iterlim"` where they
  run out),
  [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md)
  settings and iterations of a step (default 30).

- `pruneTol`:

  nodes whose weight lies this much in log below the largest are
  dropped, default `Inf`.

- `certificate`:

  `FALSE`, or a list of `points` (default 256) and `seed` (default 1):
  importance sampling at the estimate from a t5 proposal at the mode of
  each group. A group is flagged where the two -2 log L differ by more
  than three standard errors plus the error the level of the group
  allows, or the Pareto shape of the weights exceeds 0.7.

- `epsEcm`, `maxEcm`:

  for `"laplaceEM"`: the extrapolated distance to the fixed point at
  which it stops, default 1e-5, and the iteration limit.

## SAEM

Metropolis draws of the random effects by random-walk and independence
proposals at the conditional modes, Robbins-Monro averages of their
sufficient statistics, the mu-referenced fixed effects and `Omega` in
closed form, the other parameters by Newton steps on averaged quadratic
models. Entries of `control$saem`:

- `nBurnin`, `nEM`:

  at most that many burn-in and averaging iterations, default 200 and
  300.

- `nMcmc`, `nLaplace`:

  random-walk and independence steps per iteration, default 3 each;
  `nLaplace` with normal random effects alone.

- `laplaceEvery`, `laplaceScale`:

  iterations between refreshes of the modes, default 5, with normal
  random effects alone, and the factor of the proposal covariance,
  default 1.5, also of the slab value of an indicator move.

- `window`, `mcTol`:

  the burn-in ends when no parameter has a significant slope over the
  last `window` iterations (default 25), the averaging when the expected
  rise of -2 log L from the Monte Carlo error is below `mcTol` (default
  0.05).

- `particles`:

  the particles of the sequential importance sampling of -2 log L with
  Dirichlet-process random effects, default 500.

- `cm1Control`:

  [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md)
  settings of the Newton steps of the other parameters, default one
  iteration per SAEM iteration.

- `nFlip`:

  with spike-and-slab random effects the sweeps over the indicators per
  iteration, default 1.

- `nodes`:

  Gauss-Hermite nodes of the slab mass of an indicator, 3 (default) or
  5.

- `integral`:

  `"importance"` takes -2 log L by importance sampling over the
  partitions also where Dirichlet-process random effects beside
  spike-and-slab ones allow the exact sum.

With normal random effects alone, -2 log L at the estimate is the
quadrature of level `K + 2`, and `maxLevel`, `maxNodes` and `pruneTol`
of `control$quadrature` apply to it. With spike-and-slab or
Dirichlet-process random effects the averaging runs `nEM` iterations. A
sweep of Metropolis-Hastings moves visits every indicator of a group
once, in random order and for all groups at once, with the slab value
integrated out: the other slab values of the group follow along the
regression of the local quadratic model of the data, the slab mass on
that line is a Gauss-Hermite sum with the exact model at the nodes, and
the proposal draws the indicator from the resulting conditional
inclusion probability and the slab value from the local Gaussian; the
independence Metropolis-Hastings step keeps the posterior exact. `nMcmc`
random-walk steps move the slab values of every group jointly. `prob`
and `tau` follow in closed form from the averaged conditional inclusion
probabilities and slab second moments of the sweeps, not from the drawn
indicators; an estimated `tau` starts at least at the largest group mode
under the slab alone. An `em` fit as `init` starts the chain at the
selected random effects of that fit instead, the way to check a FOCEI
optimum by SAEM. The acceptance rates of births, deaths and random-walk
steps per iteration are columns of `stageTrace`. A Gibbs sweep draws the
clusters of the Dirichlet-process random effects (Neal 2000, algorithm
8). The reported clusters are the most frequent partition of the late
draws, a spike-and-slab random effect is selected where the mean
conditional inclusion probability of the late sweeps is above 1/2, `p0`
one minus that mean. -2 log L is that of `method = "focei"` at the
estimate; with Dirichlet-process random effects alone, or beside
spike-and-slab ones on more than 14 groups, by importance sampling over
the partitions, with `valueSE`.

## Sparse random effects

With spike-and-slab or Dirichlet-process random effects the parameters
and hyperparameters move in one trust region. The integral of a group is
exact where the data couple none of the factors and the sum over the
subsets of the spike-and-slab random effects of the group stays within a
budget of nodes; else it is the expectation consistent (EC)
approximation, and a trial point at which EC does not converge is a
failed evaluation. Entries of `control$sparseFit`:

- `maxOuter`:

  iterations of the trust region, default 200.

- `trustControl`:

  the [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md)
  settings of the outer search. The modes are found as by FOCEI, by the
  group kernel, where the Gauss-Newton step is below 1e-8; of
  `control$innerControl` only `sqpIterlim` applies.

## See also

[`msem()`](https://dmodverse.github.io/dMod2/reference/msem.md),
[`randomEffects()`](https://dmodverse.github.io/dMod2/reference/randomEffects.md),
[`insert()`](https://dmodverse.github.io/dMod2/reference/define.md); the
guide chapters [Nonlinear mixed-effects
models](https://dmodverse.github.io/dMod2/articles/guide-nlme.html) and
[Sparse
deviations](https://dmodverse.github.io/dMod2/articles/guide-sparse-deviations.html),
the methods articles [Marginal likelihood of mixed-effects
models](https://dmodverse.github.io/dMod2/articles/methods-nlme.html)
and [Sparsity
priors](https://dmodverse.github.io/dMod2/articles/methods-sparsity-priors.html).

## Examples

``` r
# \donttest{
## Theophylline, 12 subjects: one-compartment model with first-order absorption in closed
## form, log parameters with a normal random intercept per subject, additive error
th <- datasets::Theoph
th$Subject <- as.character(th$Subject)
subjects <- unique(th$Subject)
data <- as.datalist(data.frame(name = "y", time = th$Time, value = th$conc, sigma = NA_real_,
                               condition = th$Subject))
## second derivatives forward over reverse, the default of em()
dm <- c("forward", "reverse", "forward-reverse")
x <- Xt()
g <- Y(c(y = "amt * Ka / (V * Ka - Cl) * (exp(-Cl / V * time) - exp(-Ka * time))"),
       states = "time", modelname = "theoph_obs", compile = FALSE, derivMode = dm,
       outdir = tempdir())
#> Overwriting: /tmp/RtmpnLJd0N/theoph_obs.cpp
err <- Y(c(y = "sigma_add"), g, attachInput = FALSE, modelname = "theoph_err", compile = FALSE,
         derivMode = dm, outdir = tempdir())
#> Overwriting: /tmp/RtmpnLJd0N/theoph_err.cpp
## one condition per subject with the subject dose; the effect formula replaces each log parameter
## by an intercept plus a random effect per subject
doses <- tapply(th$Dose * th$Wt, th$Subject, `[`, 1)[subjects]
trafo <- branch(define(NULL, "x <- x", x = c("Ka", "V", "Cl", "amt", "sigma_add")),
                table = data.frame(ID = subjects, dose = doses, row.names = subjects))
trafo <- insert(trafo, "amt <- dose", dose = dose)
trafo <- insert(trafo, "x <- exp(x)", x = c("Ka", "V", "Cl", "sigma_add"))
p <- P(insert(trafo, Ka + V + Cl ~ 1 + (1 | ID)), method = "explicit", compile = FALSE,
       modelname = "theoph_p", derivMode = dm, outdir = tempdir())
compile(g, err, p, output = tempfile("theoph"))
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC]
#> reusing 14 unchanged object(s)
obj <- normL2(data, g * x * p, errmodel = err)
init <- c(Ka = 0, V = 3, Cl = 0.7, sigma_add = log(0.5))

## FOCEI fit, the estimates and the random effects per subject
fit <- em(obj, randomEffects(p), init = init)
fit
#> EM (focei): -2 log L 359.403186, converged TRUE, 10 iterations, 12 groups
#>          Ka           V          Cl   sigma_add omega_Ka_Ka   omega_V_V 
#>   0.4638825   3.4597479   1.0121568  -0.3647014  -0.4552240  -1.9774269 
#> omega_Cl_Cl 
#>  -1.3361153 
#> Random effects:
#>  type eta1 eta2  estimate         se
#>    sd   Ka   Ka 0.6343058 0.14972302
#>    sd    V    V 0.1384250 0.04075699
#>    sd   Cl   Cl 0.2628648 0.06586477
summary(fit)
#> EM summary (focei): -2 log L 359.403186, converged TRUE, 10 iterations
#>   without the 2 pi constant of the data (NONMEM OBJ): 116.803413
#> 
#> Parameters:
#>             estimate      se rsePct
#> Ka            0.4639 0.19180 41.350
#> V             3.4600 0.04718  1.364
#> Cl            1.0120 0.08333  8.233
#> sigma_add    -0.3647 0.07268 19.930
#> omega_Ka_Ka  -0.4552 0.23600 51.850
#> omega_V_V    -1.9770 0.29440 14.890
#> omega_Cl_Cl  -1.3360 0.25060 18.750
#> Standard errors without the curvature of the prediction in the random effects within the volume term (forward over reverse; derivMode = "forward-forward" gives the full one)
#> 
#> Random effects:
#>  type eta1 eta2 estimate      se
#>    sd   Ka   Ka   0.6343 0.14970
#>    sd    V    V   0.1384 0.04076
#>    sd   Cl   Cl   0.2629 0.06586
#> 
#> Shrinkage: eta_Ka = 0.018, eta_V = 0.103, eta_Cl = 0.039 
fit$omega
#>   type eta1 eta2  estimate         se
#> 1   sd   Ka   Ka 0.6343058 0.14972302
#> 2   sd    V    V 0.1384250 0.04075699
#> 3   sd   Cl   Cl 0.2628648 0.06586477
head(fit$random)
#>   group eta parameter        mode        se
#> 1     1  Ka  eta_Ka_1  0.08621976 0.1515854
#> 2     2  Ka  eta_Ka_2  0.19818408 0.1725709
#> 3     3  Ka  eta_Ka_3  0.37186008 0.1810720
#> 4     4  Ka  eta_Ka_4 -0.28173371 0.1640986
#> 5     5  Ka  eta_Ka_5 -0.04271272 0.1461175
#> 6     6  Ka  eta_Ka_6 -0.38568935 0.2043190

## quadrature started from the FOCEI estimates, and SAEM
fitQ <- em(obj, randomEffects(p), init = fit, method = "quadrature", verbose = FALSE)
fitS <- em(obj, randomEffects(p), init = init, method = "saem", verbose = FALSE)
c(focei = fit$value, quadrature = fitQ$value, saem = fitS$value)
#>      focei quadrature       saem 
#>   359.4032   360.7370   360.7588 
# }
if (FALSE) { # \dontrun{
## Ka and Cl with correlated random effects, V with a spike-and-slab one: a second transformation
p2 <- P(insert(trafo, Ka + Cl ~ 1 + (1 | p | ID), V ~ 1 + (1 | gr(ID, dist = spikeSlab()))),
        method = "explicit", compile = FALSE, modelname = "theoph_p2", derivMode = dm,
        outdir = tempdir())
compile(p2, output = tempfile("theoph"))
fit2 <- em(normL2(data, g * x * p2, errmodel = err), randomEffects(p2), init = init)
fit2$omega
fit2$slab
AIC(fit, fit2)
} # }
```
