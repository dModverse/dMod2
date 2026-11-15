# Estimate the Error Model by Restricted Maximum Likelihood

Alternates between fitting the mean parameters at a fixed error model
and updating the error parameters from the REML stationarity condition,
until the error parameters stop moving. Every data point enters with the
leverage of that point, see
[`remlLeverage()`](https://dmodverse.github.io/dMod2/reference/remlLeverage.md).
Every L2 term of a composed objective is used, so
`normL2(d1, ...) + normL2(d2, ...)` is handled as one dataset.

## Usage

``` r
reml(
  objfun,
  pars,
  errpars = NULL,
  fixed = NULL,
  iterlim = 25L,
  tol = 1e-06,
  rankTol = 1e-08,
  optimizer = trust,
  control = NULL,
  cores = getOption("dMod2.cores", 1L),
  ...
)
```

## Arguments

- objfun:

  Objective function built by
  [`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md),
  with an error model.

- pars:

  Named numeric starting vector of mean and error parameters.

- errpars:

  Character, the error-model parameters. Default `NULL`: the symbols of
  the error equations that appear in `pars`.

- fixed:

  Named numeric, held fixed throughout. Default `NULL`.

- iterlim:

  Maximum number of rounds. Default `25`.

- tol:

  Convergence threshold on the largest change of an error parameter
  between two rounds. Default `1e-6`.

- rankTol:

  Passed to
  [`remlLeverage()`](https://dmodverse.github.io/dMod2/reference/remlLeverage.md).
  Default `1e-8`.

- optimizer:

  Optimizer for the mean step, called as
  `optimizer(objfun, parinit, fixed = , cores = , ...)`. Default
  [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md).

- control:

  Named list of arguments for `optimizer`, merged into
  `list(rinit = 0.1, rmax = 10)`. Default `NULL`.

- cores:

  Passed to the prediction and the objective. Default
  `getOption("dMod2.cores", 1)`.

- ...:

  Further arguments for `objfun` in the mean step, passed through
  `optimizer`. `value.plain` is evaluated without them.

## Value

A list with components

- `argument`:

  The full parameter vector.

- `value`:

  The REML criterion, `value.plain + logdet`.

- `value.plain`:

  The objective at `argument`.

- `logdet`:

  Log determinant of \\J^\top W J\\, see
  [`remlLeverage()`](https://dmodverse.github.io/dMod2/reference/remlLeverage.md).

- `errpars`:

  The error-model parameters.

- `leverage`:

  The data frame returned by
  [`remlLeverage()`](https://dmodverse.github.io/dMod2/reference/remlLeverage.md).

- `dof`:

  Effective degrees of freedom per observable.

- `rank`:

  Effective number of mean parameters.

- `iterations`:

  Number of rounds.

- `converged`:

  Logical.

- `fit`:

  The last mean fit.

## Details

Only first-order sensitivities are used. A term acting on the error
parameters from outside the data term, such as a prior, is not part of
the REML condition; `reml()` warns when it finds one.

## See also

[`remlLeverage()`](https://dmodverse.github.io/dMod2/reference/remlLeverage.md),
[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md)

## Examples

``` r
# \donttest{
old <- options(dMod2.outdir = tempdir())
f <- addReaction(eqnlist(), from = "A", to = "", rate = "k*A")
x <- Xs(odemodel(f, modelname = "reml_x", compile = FALSE))
g <- Y(c(y = "A"), f = x, attachInput = FALSE, modelname = "reml_g",
       compile = FALSE)
e <- Y(c(y = "sigma_y"), f = g, attachInput = FALSE, condition = "C1",
       modelname = "reml_e", compile = FALSE)
p <- P(eqnvec(A = "exp(logA)", k = "exp(logk)", sigma_y = "sigma_y"),
       condition = "C1", modelname = "reml_p", compile = FALSE)
compile(x, g, e, p, output = "reml_example", cores = 1)
#> using C++ compiler: g++
#>   every source : -O2 -DNDEBUG -w -fPIC
#>   2 of 5 also  : -fopenmp

set.seed(1)
times <- seq(0, 5, by = 0.25)
data <- datalist(C1 = data.frame(
  name = "y", time = times, sigma = NA_real_,
  value = 2 * exp(-0.5 * times) + rnorm(length(times), sd = 0.1)))
obj <- normL2(data, g * x * p, errmodel = e)

fit <- reml(obj, c(logA = 0, logk = 0, sigma_y = 1))
fit$argument
#>        logA        logk     sigma_y 
#>  0.67832557 -0.75124935  0.09166024 
fit$dof
#>  y 
#> 19 
head(remlLeverage(obj, fit$argument))
#>   condition time name      sigma    residual   leverage
#> 1        C1 0.00    y 0.09166024  0.03322076 0.39880136
#> 2        C1 0.25    y 0.09166024 -0.03202037 0.24442887
#> 3        C1 0.50    y 0.09166024  0.08245573 0.15525247
#> 4        C1 0.75    y 0.09166024 -0.15077817 0.10691867
#> 5        C1 1.00    y 0.09166024 -0.01658392 0.08322519
#> 6        C1 1.25    y 0.09166024  0.10417243 0.07360126
options(old)
# }
```
