# Fit an alignment model

Estimates the parameters of an
[`alignmentModel()`](https://dmodverse.github.io/dMod2/reference/alignmentModel.md)
by maximum likelihood, the marginal one with random terms
([`em()`](https://dmodverse.github.io/dMod2/reference/em.md)), with the
normalization of the values per block of experiments (see
[`alignmentModel()`](https://dmodverse.github.io/dMod2/reference/alignmentModel.md))
in the objective; the error parameters and the hyperparameters of normal
random terms by default by restricted likelihood. The Hessian is sparse
over the experiments, with the normalization as a low-rank part (see
[`objlist()`](https://dmodverse.github.io/dMod2/reference/objlist.md));
[`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md) and
[`em()`](https://dmodverse.github.io/dMod2/reference/em.md) take the
sparse path as `sparse = "auto"` decides.

## Usage

``` r
align(
  model,
  fits = 10L,
  cores = 1L,
  method = c("focei", "quadrature"),
  sd = 0.5,
  iterlim = 1000L,
  init = NULL,
  reml = TRUE
)
```

## Arguments

- model:

  An `alignmentmodel`.

- fits:

  Number of fits from perturbed starts. Default 10.

- cores:

  Number of fits run in parallel; every fit runs on one core. Default 1.

- method:

  With random terms: `"focei"` (default) or `"quadrature"`, see
  [`em()`](https://dmodverse.github.io/dMod2/reference/em.md).

- sd:

  Standard deviation of the start perturbation on the parameter scale.
  Default 0.5.

- iterlim:

  Iteration limit of the optimizer; with spike-and-slab or
  Dirichlet-process terms the outer iterations of
  [`em()`](https://dmodverse.github.io/dMod2/reference/em.md). Default
  1000.

- init:

  Optional named start of some parameters; the intercepts of the
  biological symbols start where the symbol is the median reading, the
  other intercepts where their symbol is 1, the other parameters at 0.
  The coefficients of the first scaling symbol then start where the
  predictor of that symbol is, per experiment, the median of the
  predictors at which the prediction meets the relative readings.

- reml:

  `TRUE` (default): the error parameters, and the hyperparameters of
  normal random terms, by restricted likelihood, the parameters of the
  symbols of the prediction (`model` and `absolute` of the
  [`alignmentModel()`](https://dmodverse.github.io/dMod2/reference/alignmentModel.md))
  integrated out under a flat prior (Lindstrom and Bates 1990), in at
  most 20 rounds from the best maximum likelihood fit until the error
  parameters and hyperparameters move by less than 1e-4. Without random
  terms or with normal ones and `method = "focei"`; else maximum
  likelihood. Where the information of the values and scales is
  singular, at the maximum likelihood fit or at a round, the restricted
  likelihood is not defined, and the fit is by maximum likelihood with a
  warning. `FALSE`: maximum likelihood.

## Value

Object of class `alignment`: `aligned` (one row per combination of the
columns of the biological formulas: the columns, `symbol`, `estimate`
and `se` of the predictor, `value` and `sigma` on the scale of the
symbol; with REML `se` from the information of the mean parameters
alone, the error parameters held at their estimates), `scaled` (the
readings with `prediction`, `sigma`, `residual` and `bloq`),
`parameters` (with attribute `estimation`, `"REML"` or `"ML"`), `value`
(-2 log L, of the restricted likelihood with REML, with the
normalization and the prior), `npar` (one less per block without
absolute readings), `AIC` (with REML comparable only between fits with
the same parameters of the prediction), `values` (-2 log L of all fits,
sorted), `multistart` (`value`, `converged` and `iterations` per fit,
see
[`plotValues()`](https://dmodverse.github.io/dMod2/reference/plotValues.md);
with REML the maximum likelihood fits the rounds start from),
`remlRounds` (the rounds of the restricted likelihood, 0 without),
`seSource` (`"Hessian"`, or `"Gauss-Newton"` where the reported Hessian
of a fit with random terms is not positive definite: then the
Gauss-Newton Schur complement of the expected information over the
random effects; `"none"` where that fails too, without standard errors),
`fit`, `objfn`, `omega`, `random` and `randomEffects` with random terms,
and the `model`. `scaled` also has `scaled`, the value of the biological
symbol at which the prediction meets the reading with the other symbols
at the fit (`NaN` where none does), and `scaledSigma`, the error model
over the slope of the prediction there; `aligned` has `lloq`, the
smallest such value of the censored readings of a combination.

## Details

The objective of an alignment is exact to rounding, so the fits stop at
changes of the objective, predicted decreases and gradients below
`sqrt(.Machine$double.eps)`, the inner fits of the random effects too.

## References

Kemmer, S., Bang, S., Rosenblatt, M., Timmer, J. and Kaschek, D. (2022).
BlotIt: Optimal alignment of Western blot and qPCR experiments. *PLOS
ONE* 17(8), e0264295.
[doi:10.1371/journal.pone.0264295](https://doi.org/10.1371/journal.pone.0264295)

Lindstrom, M. J. and Bates, D. M. (1990). Nonlinear mixed effects models
for repeated measures data. *Biometrics* 46(3), 673-687.
[doi:10.2307/2532087](https://doi.org/10.2307/2532087)

## See also

[`alignmentModel()`](https://dmodverse.github.io/dMod2/reference/alignmentModel.md),
[`plotAligned()`](https://dmodverse.github.io/dMod2/reference/plotAligned.md);
the guide chapter [Alignment of relative data
measurements](https://dmodverse.github.io/dMod2/articles/guide-alignment.html).

## Examples

``` r
## Five time points on six gels, each gel with a separate scale and an additive background
set.seed(1)
d <- expand.grid(time = c("0", "10", "30", "60", "120"), gel = paste0("G", 1:6),
                 stringsAsFactors = FALSE)
truth <- c(0.02, 1, 4, 2, 0.5)[match(d$time, c("0", "10", "30", "60", "120"))]
d$value <- (rlnorm(6, 7, 0.5)[match(d$gel, unique(d$gel))] * truth + 50) * rlnorm(nrow(d), 0, 0.1)

# \donttest{
old <- options(dMod2.outdir = tempdir())
## every symbol positive, on the log scale of the predictor; one value per time point,
## one scale per gel
logs <- c("yi <- exp(yi)", "sj <- exp(sj)", "bj <- exp(bj)", "sigmaR <- exp(sigmaR)")
m <- alignmentModel(d, model = "yi * sj + bj", errmodel = "sigmaR * prediction",
                    biological = yi ~ time, scaling = sj ~ gel, error = bj + sigmaR ~ 1,
                    trafo = logs)
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC -fopenmp]
## a background varying between gels: a list of formulas for one role
m2 <- alignmentModel(d, model = "yi * sj + bj", errmodel = "sigmaR * prediction",
                     biological = yi ~ time, scaling = list(sj ~ gel, bj ~ 1 + (1 | gel)),
                     error = sigmaR ~ 1, trafo = logs)
m
#> alignmentModel: 30 readings, model yi * sj + bj, error sigmaR * prediction
#>   yi         biological    5 parameters, yi = exp(yi)
#>   sj         scaling       6 parameters, sj = exp(sj)
#>   bj         error         1 parameter, bj = exp(bj)
#>   sigmaR     error         1 parameter, sigmaR = exp(sigmaR)
#>   mean of yi normalized in 1 block of experiments
m2
#> alignmentModel: 30 readings, model yi * sj + bj, error sigmaR * prediction
#>   yi         biological    5 parameters, yi = exp(yi)
#>   sj         scaling       6 parameters, sj = exp(sj)
#>   bj         scaling       1 parameter, bj = exp(bj), random over gel
#>   sigmaR     error         1 parameter, sigmaR = exp(sigmaR)
#>   mean of yi normalized in 1 block of experiments
fit <- align(m, fits = 4)
fit
#> alignment: 30 readings, 5 biological values, 12 parameters
#>   restricted likelihood (REML, 2 rounds): -2 log L 389.06, AIC 413.06
#>   maximum likelihood fits, REML from the best: -2 log L 329.79 329.79 329.79 329.79
fit$aligned
#>   time symbol   estimate         se      value       sigma
#> 1    0     yi -3.2127997 0.18302237 0.04024378 0.007365512
#> 2   10     yi  0.4499769 0.05497631 1.56827596 0.086218023
#> 3   30     yi  1.7783403 0.05924187 5.92002292 0.350713229
#> 4   60     yi  1.1659908 0.05670480 3.20910085 0.181971415
#> 5  120     yi -0.1815083 0.05135481 0.83401135 0.042830490
plotAligned(fit)

options(old)
# }
```
