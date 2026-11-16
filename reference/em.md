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
  derivative of the volume term. The default, and the method for
  spike-and-slab and Dirichlet-process random effects (at most 14
  groups) and for several groupings.

- `"quadrature"`:

  adaptive sparse-grid quadrature around the modes, started from the
  FOCEI fit; section Quadrature.

- `"laplaceEM"`:

  ECM on the Laplace approximation with the covariance-corrected update
  of `Omega`.

- `"saem"`:

  stochastic-approximation EM, also for Dirichlet-process random
  effects; section SAEM.

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

  `innerControl`

  :   the trust region of the modes; it stops where the Newton step is
      below `steptol`, default 1e-8.

  `trustControl`

  :   the outer trust region; it stops where the predicted change of -2
      log L is below `mtol`, default 1e-5.

  `secondOrderCorrection`

  :   the derivative of the volume term, default `TRUE`.

  `derivMode`

  :   how FOCEI gets the second derivatives of the volume term:
      `"forward-reverse"` (default), sweeps for models built with
      `derivMode = c("forward", "reverse", "forward-reverse")`, or
      `"forward-forward"`, the second-order sensitivities of models
      built with `deriv2 = TRUE`. Without the mode asked for the fit
      takes the other, with a message.

  `boundary`

  :   the largest rise of -2 log L at which a variance is set to 0,
      default 0.1, `FALSE` for none. With one grouping and several
      blocks, a normal random effect that forms a block alone and whose
      standard deviation ends below 0.1 is dropped and the fit repeated
      from the estimates; where -2 log L rises by at most `boundary`,
      that fit is returned and names the hyperparameter in `boundary`.

  `quadrature`, `saem`, `sparseFit`

  :   lists of the sections Quadrature, SAEM and Sparse random effects.

- cores:

  Number of cores over the groups and conditions inside the fit; with 1
  (default) no threads are started.

- verbose:

  Logical; print the progress of the iterative methods.

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
  products of the centered scores of the `G` groups, plus half the
  Hessian of the conditions without random effects.

- `converged`, `iterations`:

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
  `Inf` keeps the levels of the stages.

- `levels`:

  the stages, default `K + 1` with `K` the number of random effects of a
  group; level `K + s` integrates polynomials of total degree `2 s + 1`
  exactly.

- `maxLevel`:

  highest level, default `K + 8`.

- `tolDecrement`:

  a stage ends when a Newton step on the rebuilt nodes would lower -2
  log L by less than this, default 1e-4, or than the nodes resolve.

- `maxEcmPerStage`, `cm1Control`, `maxCm1Iter`:

  rebuilds per stage (default 10),
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

  random-walk and independence steps per iteration, default 3 each.

- `laplaceEvery`, `laplaceScale`:

  iterations between refreshes of the modes, default 5, and the factor
  of the proposal covariance, default 1.5.

- `window`, `mcTol`:

  the burn-in ends when no parameter has a significant slope over the
  last `window` iterations (default 25), the averaging when the expected
  rise of -2 log L from the Monte Carlo error is below `mcTol` (default
  0.05).

- `particles`:

  with Dirichlet-process random effects the particles of the sequential
  importance sampling of -2 log L, default 500.

- `cm1Control`:

  [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md)
  settings of the Newton steps of the other parameters, default one
  iteration per SAEM iteration.

With Dirichlet-process random effects a Gibbs sweep draws the clusters
of the groups (Neal 2000, algorithm 8), the averaging runs `nEM`
iterations, and the reported clusters are the most frequent partition of
the late draws.

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
  settings; `ftol`, `mtol` and `gtol` of the `tolControl` there also end
  the search of the modes (else 1e-10).

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
#> Overwriting: /tmp/RtmpZujll1/theoph_obs.cpp
err <- Y(c(y = "sigma_add"), g, attachInput = FALSE, modelname = "theoph_err", compile = FALSE,
         derivMode = dm, outdir = tempdir())
#> Overwriting: /tmp/RtmpZujll1/theoph_err.cpp
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
#> Standard errors without the curvature of the prediction in the random effects within the volume term (forward over reverse; deriv2 = TRUE gives the full one)
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
#>   359.4032   360.7548   360.7588 
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
