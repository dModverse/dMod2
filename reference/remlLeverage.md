# Leverage of the Data Points in the Mean Model

Hat values of the weighted mean model, \\h = \mathrm{diag}(W^{1/2} J
(J^\top W J)^{-1} J^\top W^{1/2})\\ with \\J =
\partial\mu/\partial\theta\\ the prediction sensitivities and \\W =
\mathrm{diag}(1/\sigma_i^2)\\: the leverage of each data point. Their
sum is the numerical rank of the weighted sensitivity matrix, so a
non-identifiable direction contributes nothing.

## Usage

``` r
remlLeverage(
  objfun,
  pars,
  meanpars = NULL,
  fixed = NULL,
  rankTol = 1e-08,
  cores = getOption("dMod2.cores", 1L)
)
```

## Arguments

- objfun:

  Objective function built by
  [`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md),
  with prediction function, data and error model.

- pars:

  Named numeric parameter vector, usually a fit.

- meanpars:

  Character, the mean parameters spanning \\J\\. Default `NULL`:
  `names(pars)` without the error-model parameters.

- fixed:

  Named numeric, passed to the prediction. Default `NULL`.

- rankTol:

  Relative threshold on the singular values of the weighted sensitivity
  matrix below which a direction counts as absent. Default `1e-8`.

- cores:

  Passed to the prediction. Default `getOption("dMod2.cores", 1)`.

## Value

A data frame with one row per data point and columns `condition`,
`time`, `name`, `sigma`, `residual` and `leverage`. Attributes:
`"rank"`, the effective number of mean parameters; `"logdet"`, the log
determinant of \\J^\top W J\\ over the directions present; `"dof"`, per
observable the number of data points minus the summed leverage.

## See also

[`reml()`](https://dmodverse.github.io/dMod2/reference/reml.md)

## Examples

``` r
# \donttest{
old <- options(dMod2.outdir = tempdir())
f <- addReaction(eqnlist(), from = "A", to = "", rate = "k*A")
x <- Xs(odemodel(f, modelname = "reml_x", compile = FALSE))
#> Overwriting existing file: /tmp/RtmptyY1gm/reml_x.cpp
#> Overwriting existing file: /tmp/RtmptyY1gm/reml_x_s.cpp
g <- Y(c(y = "A"), f = x, attachInput = FALSE, modelname = "reml_g",
       compile = FALSE)
#> Overwriting: /tmp/RtmptyY1gm/reml_g.cpp
e <- Y(c(y = "sigma_y"), f = g, attachInput = FALSE, condition = "C1",
       modelname = "reml_e", compile = FALSE)
#> Overwriting: /tmp/RtmptyY1gm/reml_e_C1.cpp
p <- P(eqnvec(A = "exp(logA)", k = "exp(logk)", sigma_y = "sigma_y"),
       condition = "C1", modelname = "reml_p", compile = FALSE)
#> Overwriting: /tmp/RtmptyY1gm/reml_p_C1.cpp
compile(x, g, e, p, output = "reml_example", cores = 1)
#> Warning: A shared library named 'reml_example' is already loaded; overwriting it is not portable. Using 'reml_example_2' instead.
#> using C++ compiler: g++
#>   every source : -O2 -DNDEBUG -w -fPIC
#>   2 of 5 also  : -fopenmp
#> reusing 5 unchanged object(s)

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
