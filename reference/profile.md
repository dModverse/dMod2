# Profile Likelihood

Computes the profile likelihood of one or more parameters around an
optimum. `method = "integrate"` follows the profile path by Euler steps
of the path equation of the Lagrangian; `method = "optimize"`
re-optimizes the other parameters at every step. A method of
[`stats::profile()`](https://rdrr.io/r/stats/profile.html) for dMod
objective functions and plain R functions.

## Usage

``` r
# S3 method for class 'objfn'
profile(
  fitted,
  pars,
  whichPar,
  alpha = 0.05,
  limits = c(lower = -Inf, upper = Inf),
  method = c("integrate", "optimize"),
  stepControl = NULL,
  algoControl = NULL,
  optControl = NULL,
  verbose = FALSE,
  cores = 1,
  cautiousMode = FALSE,
  side = c("both", "left", "right"),
  delta = NULL,
  sparse = "auto",
  ...
)

# S3 method for class '`function`'
profile(fitted, ...)
```

## Arguments

- fitted:

  Objective function `fitted(pars, fixed, ...)` returning a list with
  `value`, `gradient` and `hessian`. Numeric attributes of the result,
  such as `"data"` and `"prior"`, are kept as columns of the profile,
  except `chi2` and `chi2_*`.

- pars:

  Named numeric parameter vector at the optimum.

- whichPar:

  Numeric or character vector, the parameters to profile.

- alpha:

  Significance level; the profile stops where the objective has risen by
  `qchisq(1 - alpha, 1)`. Default `0.05`.

- limits:

  Numeric of length 2, the lower and upper deviation from
  `pars[whichPar]` at which a profile stops. Default
  `c(lower = -Inf, upper = Inf)`.

- method:

  `"integrate"` (default) or `"optimize"`. Selects the defaults of
  `stepControl`, `algoControl` and `optControl`.

- stepControl:

  List of step-size settings, merged into the defaults. Defaults are
  given as `"integrate"` / `"optimize"`.

  `stepsize`

  :   Initial step in the profiled parameter. `1e-4` / `1e-2`.

  `min`, `max`

  :   Bounds on the step size. `1e-4` and `Inf` for both methods.

  `atol`

  :   The step shrinks by a factor 1.5 when actual and predicted change
      of the objective differ by more than `atol`. `1e-2` / `1e-1`.

  `rtol`

  :   The step doubles when the difference is below `0.3 * atol`, or
      below `0.3 * rtol` relative to the actual change. `1e-2` / `1e-1`.

  `limit`

  :   Maximum number of steps per side. `500` / `100`.

  `stop`

  :   Component of the objective whose rise is compared with `delta`:
      `"value"` or the name of a numeric attribute such as `"data"`.
      `"value"` for both methods.

- algoControl:

  List of path settings, merged into the defaults. Defaults are given as
  `"integrate"` / `"optimize"`.

  `gamma`

  :   Initial factor of the correction toward the optimum of the other
      parameters. `1` / `0`.

  `W`

  :   Matrix in the path equation, `"hessian"` or `"identity"`.
      `"hessian"` is the exact Hessian where `fitted` returns one with
      `deriv2 = TRUE`, as a model built with second derivatives does,
      otherwise the one `fitted` returns, for
      [`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md)
      the Gauss-Newton matrix; `deriv2 = FALSE` in `...` keeps that one.
      With `"identity"` a step moves the profiled parameter alone and
      the others are reoptimized at every step. `"hessian"` /
      `"identity"`.

  `reoptimize`

  :   Whether every step is followed by a fit of the other parameters
      with
      [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md).
      `TRUE` for both methods.

  `correction`

  :   Relative correction above which `gamma` is halved; below half of
      it `gamma` doubles up to the initial value. `1` for both methods.

  `reg`

  :   Where the Hessian bordered by the constraint is numerically
      singular, the path equation is solved by a pseudoinverse without
      the singular values below `reg` times the largest one and the
      matrix dimension. `.Machine$double.eps` / `0`.

- optControl:

  List of arguments for
  [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md) in
  the re-optimization, merged into
  `list(rinit = 0.1, rmax = 10, iterlim = 10)` for `"integrate"` and
  `list(rinit = 0.1, rmax = 10, iterlim = 100)` for `"optimize"`. Any
  other argument of
  [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md) can
  be added.

- verbose:

  Logical, print progress. Default `FALSE`.

- cores:

  Number of parallel workers over parameters, or
  `c(pars = , conditions = )` to also give each worker a condition axis,
  which selects a PSOCK backend. Keep the product with the threads of
  each objective evaluation below the number of cores. Default `1`.

- cautiousMode:

  Logical. If `TRUE`, every step is written to
  `profiles-interRes/<parameter>-<side>.R` in the working directory, and
  the files are kept. Default `FALSE`.

- side:

  `"both"` (default), `"left"` or `"right"`, the side of the optimum to
  profile.

- delta:

  Rise in the objective at which a profile stops. Default `NULL`, which
  uses `qchisq(1 - alpha, 1)`. Pass
  [`profileThreshold()`](https://dmodverse.github.io/dMod2/reference/profileThreshold.md)
  with the arguments later given to
  [`confint.parframe()`](https://dmodverse.github.io/dMod2/reference/confint.parframe.md),
  so the profile reaches that threshold.

- sparse:

  `"auto"` (default), `TRUE` or `FALSE`: the path equation and the
  reoptimizations on the sparse Hessian of `fitted`, chosen at `pars` as
  in [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md).
  `hessian = "sparse"` in `...` is `TRUE`.

- ...:

  Further arguments passed to `fitted`, for example `fixed`.

## Value

A
[`parframe()`](https://dmodverse.github.io/dMod2/reference/parframe.md)
with one row per profile point of all profiled parameters, with
metanames `value` (the objective), `constraint` (deviation of the
profiled parameter from the value in `pars`), `stepsize`, `gamma` and
`whichPar` (the profiled parameter), columns for the numeric attributes
of the objective other than `chi2` and `chi2_*`, and one column per
parameter. Profiles that fail are dropped with a message; when all fail,
`profile()` stops with an error. The attribute `"ends"` is a data frame
with one row per side and the columns `whichPar`, `side`, `steps` and
`reason`: `"threshold"`, `"limit"`, `"step limit"` or the error of the
objective that ended the side.

## See also

[`confint.parframe()`](https://dmodverse.github.io/dMod2/reference/confint.parframe.md)
for confidence intervals,
[`plotProfile()`](https://dmodverse.github.io/dMod2/reference/plotProfile.md)
and
[`plotPaths()`](https://dmodverse.github.io/dMod2/reference/plotPaths.md)
for plots,
[`profileThreshold()`](https://dmodverse.github.io/dMod2/reference/profileThreshold.md);
the guide chapters [Profile
likelihood](https://dmodverse.github.io/dMod2/articles/guide-profiles.html),
[Prediction
uncertainty](https://dmodverse.github.io/dMod2/articles/guide-prediction-uncertainty.html)
and [Remote
computing](https://dmodverse.github.io/dMod2/articles/guide-remote-computing.html),
which explains `cores`.

## Examples

``` r
# \donttest{
old <- options(dMod2.outdir = tempdir())

## Decay model, observation and parameter transformation, compiled together
f <- addReaction(eqnlist(), from = "A", to = "", rate = "k*A")
x <- Xs(odemodel(f, modelname = "profile_x", compile = FALSE))
g <- Y(c(y = "A"), f = x, attachInput = FALSE, modelname = "profile_g",
       compile = FALSE)
p <- P(eqnvec(A = "exp(logA)", k = "exp(logk)"), condition = "C1",
       modelname = "profile_p", compile = FALSE)
compile(x, g, p, output = "profile_example", cores = 1)
#> using C++ compiler: g++
#>   every source : -O2 -DNDEBUG -w -fPIC
#>   2 of 4 also  : -fopenmp

## Data, objective and fit
set.seed(1)
times <- seq(0, 5, by = 0.5)
data <- datalist(C1 = data.frame(
  name = "y", time = times, sigma = 0.05,
  value = 2 * exp(-0.5 * times) + rnorm(length(times), sd = 0.05)))
obj <- normL2(data, g * x * p)
myfit <- trust(obj, c(logA = 0, logk = 0), rinit = 1, rmax = 10)

## Profiles by integration and by repeated optimization
profiles.approx <- profile(obj, myfit$argument, whichPar = c("logA", "logk"))
profiles.exact <- profile(obj, myfit$argument, whichPar = c("logA", "logk"),
                          method = "optimize")

## Plots and confidence intervals
plotProfile(list(approx = profiles.approx, exact = profiles.exact))

plotPaths(profiles.approx, whichPar = "logk")

confint(profiles.approx, val.column = "value")
#>      name      value      lower      upper
#> logA logA  0.6774738  0.6369417  0.7167936
#> logk logk -0.7369431 -0.8051900 -0.6699037

options(old)
# }
```
