# Multistart mixed-effects fit

Runs [`em()`](https://dmodverse.github.io/dMod2/reference/em.md) from
`fits` starts, the first at `init`, the others perturbed by normal noise
of standard deviation `sd` (on the log scale for concentrations and slab
standard deviations, on the logit scale for slab probabilities), and
returns the best fit.

## Usage

``` r
msem(
  obj,
  random,
  init = NULL,
  fixed = NULL,
  method = c("focei", "quadrature", "laplaceEM", "saem"),
  control = list(),
  fits = 20L,
  cores = 1L,
  sd = 0.5,
  verbose = FALSE,
  sparse = "auto",
  retry = TRUE,
  nTries = 10L,
  keepFits = FALSE,
  which = NULL
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

- fits:

  Number of starts. Default 20.

- cores:

  Number of fits run in parallel, or `c(fits = , conditions = )` as in
  [`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md):
  fits in parallel and cores inside each fit.

- sd:

  Standard deviation of the perturbation of the starts: one number for
  every parameter, or a vector named by the parameters it perturbs, the
  others keeping their start.

- verbose:

  Logical. If `TRUE`, print the table of the starts; the fits print no
  progress. Default `FALSE`.

- sparse:

  `"auto"` (default), `TRUE` or `FALSE`: the outer trust region of FOCEI
  over normal random effects on the sparse Gauss-Newton Hessian.
  `"auto"` decides as
  [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md)
  does.

- retry:

  Logical. If `TRUE` (default), a fit that fails with an error, for
  instance at a start where a steady state of
  [`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md) has
  no stable root, is repeated from a new perturbed start, up to `nTries`
  attempts.

- nTries:

  Maximum number of attempts per start, including the first. Default
  `10`.

- keepFits:

  Logical. If `TRUE`, the result has `fits`, the fit of every start, or
  the error message of a failed one, without the model. Default `FALSE`.

- which:

  Positions of the starts to fit, among the `fits` drawn. `NULL`
  (default) fits all. Tasks of
  [`distributedComputing()`](https://dmodverse.github.io/dMod2/reference/distributedComputing.md)
  that set the same seed and fit disjoint parts fit the starts of one
  call, and [`get()`](https://rdrr.io/r/base/get.html) joins their
  results into one.

## Value

The fit with the smallest `value`, of class `em`, with `multistart`, a
data frame with a row per start:

- `start`:

  the position of the start among the `fits` drawn.

- `value`, `converged`, `iterations`:

  of the fit.

- `gradient`:

  the largest absolute entry of the gradient at the end of the fit.

- `tries`:

  the attempts it took.

- `stop`:

  the test that ended the trust region, `"error"` where every attempt
  failed.

- `failures`:

  why evaluations in the trust region failed, and how often.

- `error`, `boundary`:

  the error of a failed start, the variances the fit set to 0.

[`plotValues()`](https://dmodverse.github.io/dMod2/reference/plotValues.md)
draws it as a waterfall.

## See also

[`em()`](https://dmodverse.github.io/dMod2/reference/em.md)

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
#> Overwriting: /tmp/RtmpriQ23f/theoph_obs.cpp
err <- Y(c(y = "sigma_add"), g, attachInput = FALSE, modelname = "theoph_err", compile = FALSE,
         derivMode = dm, outdir = tempdir())
#> Overwriting: /tmp/RtmpriQ23f/theoph_err.cpp
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

## four starts: the first at init, the others perturbed
best <- msem(obj, randomEffects(p), init = init, fits = 4)
best$multistart
#>   start    value converged iterations    gradient tries stop failures error
#> 1     1 359.4032      TRUE         10 0.002975526     1 <NA>           <NA>
#> 2     2 359.4032      TRUE         11 0.002859698     1 <NA>           <NA>
#> 3     3 359.4032      TRUE         14 0.000496391     1 <NA>           <NA>
#> 4     4 359.4032      TRUE         11 0.001513915     1 <NA>           <NA>
#>   boundary
#> 1         
#> 2         
#> 3         
#> 4         
plotValues(best)

# }
if (FALSE) { # \dontrun{
## 20 starts, four fits at a time
best <- msem(obj, randomEffects(p), init = init, fits = 20, cores = c(fits = 4, conditions = 1))
plotValues(best)
} # }
```
