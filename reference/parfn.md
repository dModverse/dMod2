# Parameter Transformation Function

Builds a function that transforms one parameter vector into another and
propagates the Jacobian.
[`P()`](https://dmodverse.github.io/dMod2/reference/P.md) calls it
internally; use it directly for a transformation of your own.

## Usage

``` r
parfn(p2p, parameters = NULL, condition = NULL)
```

## Arguments

- p2p:

  transformation for one condition, a function
  `p2p(pars, fixed, deriv, deriv2)` that maps the parameter vector
  `pars` and the fixed values `fixed` to a new parameter vector. With
  `deriv = TRUE` the result has the attribute `deriv`, the Jacobian of
  the transformation, with `deriv2 = TRUE` also `deriv2`, the second
  derivatives. A `p2p` with an argument `condition` receives the
  condition name.

- parameters:

  character, the parameters accepted by the function. Default `NULL`.

- condition:

  character, the condition for which the transformation is defined.
  Default `NULL`, every condition.

## Value

Object of class `parfn`, a function
`p(..., fixed, deriv, deriv2, hessian, conditions, env, cores, sweep)`
that takes the parameter vector `pars` through `...`. The attributes of
the function are `"mappings"`, the list of `p2p` functions,
`"parameters"`, the union of their parameters, and `"conditions"`.

## Calling a dMod function

The leading arguments depend on the class and are matched by position or
by name:

- prediction function: `x(times, pars, ...)`, returning a
  [prdlist](https://dmodverse.github.io/dMod2/reference/prdlist.md);

- observation function: `g(out, pars, ...)` with a
  [prdlist](https://dmodverse.github.io/dMod2/reference/prdlist.md)
  `out`, returning a
  [prdlist](https://dmodverse.github.io/dMod2/reference/prdlist.md);

- parameter transformation: `p(pars, ...)`, returning a list with one
  [parvec](https://dmodverse.github.io/dMod2/reference/parvec.md) per
  condition;

- objective function: `obj(pars, ...)`, returning an
  [objlist](https://dmodverse.github.io/dMod2/reference/objlist.md).

`pars` is a named numeric vector. If it has a `deriv` attribute, as the
output of a parameter transformation does, the derivatives are taken
with respect to the parameters behind it. A composition is called
according to the resulting class: `g * x * p` is a prediction function,
`(g * x * p)(times, pars)`.

The other arguments are named:

- `fixed`: named numeric vector of parameters that enter the evaluation
  without derivatives. Default `NULL`.

- `deriv`: `TRUE` (default) returns first derivatives with respect to
  `pars`.

- `deriv2`: `TRUE` returns second derivatives as well. Default `FALSE`.
  Every function of the chain has to be built for them, see the
  arguments `deriv2` and `derivMode` of
  [`odemodel()`](https://dmodverse.github.io/dMod2/reference/odemodel.md),
  [`Y()`](https://dmodverse.github.io/dMod2/reference/Y.md) and
  [`P()`](https://dmodverse.github.io/dMod2/reference/P.md).

- `hessian`: objective functions only, see the table below. Default
  `NULL`.

- `sweep`: objective functions only. `"forward"` (default) propagates
  the sensitivities through the chain; `"reverse"` computes the gradient
  in one backward pass, whose cost does not grow with the number of
  parameters. `"reverse"` needs every function of the chain built with
  `derivMode = "reverse"`, and with `"forward-reverse"` for
  `deriv2 = TRUE`.

- `conditions`: character vector, the conditions to evaluate. Default
  all conditions of the function.

- `env`: environment for intermediate results. Default `NULL`, a new
  one. An objective function stores the prediction there as `prediction`
  and returns the environment as attribute `env`.

- `cores`: number of threads for the ODE solves and the residual kernel.
  Default `getOption("dMod2.cores", 1L)`.

Prediction, observation and parameter functions accept `hessian` and
`sweep` and ignore them. With `deriv = TRUE` an objective function
returns this Hessian:

|             |                |                      |                      |
|-------------|----------------|----------------------|----------------------|
| `sweep`     | `hessian`      | `deriv2 = FALSE`     | `deriv2 = TRUE`      |
| `"forward"` | `NULL`, `TRUE` | Gauss-Newton         | exact                |
| `"forward"` | `FALSE`        | none                 | none, with a warning |
| `"reverse"` | `NULL`         | none                 | exact                |
| `"reverse"` | `TRUE`         | none, with a warning | exact                |
| `"reverse"` | `FALSE`        | none                 | none, with a warning |

`hessian = "sparse"` returns the Gauss-Newton Hessian of the forward
sweep as a sparse `Matrix`, assembled from the blocks of the conditions;
[`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md)
factorizes it where the argument `sparse` takes the sparse path, and
[`vcov()`](https://dmodverse.github.io/dMod2/reference/vcov.objlist.md)
solves it for the columns of `which`.

`hessian = "expected"` returns, in place of the Gauss-Newton Hessian of
the forward sweep, the expected (Fisher) information of the likelihood:
`2 / sigma^2` in the prediction and `4 / sigma^2` in sigma per data
point, positive semidefinite also where an error model depends on the
parameters. The Gauss-Newton Hessian keeps the exact curvature of the
loss in the prediction and sigma, which is indefinite there and shows
the optimizer a saddle.

With `deriv = FALSE` it returns the value alone.

## See also

[+.fn](https://dmodverse.github.io/dMod2/reference/plus-.fn.md),
[`P()`](https://dmodverse.github.io/dMod2/reference/P.md)

## Examples

``` r
## A prediction function interpolating parameter values piecewise linearly
times <- 0:5
grid <- data.frame(name = "A", time = times, row.names = paste0("p", times))
x <- Xd(grid)

## An observation function
g <- Y(c(Aobs = "s*A"), states = "A", parameters = "s",
       modelname = "prediction_obs", outdir = tempdir())
#> Overwriting: /tmp/RtmpNgFK2N/prediction_obs.cpp

## One parameter transformation per condition
innerpars <- c(getParameters(x), getParameters(g))
trafo <- setNames(innerpars, innerpars)
p <- P(list(C1 = replaceSymbols(innerpars, paste0(innerpars, "_C1"), trafo),
            C2 = replaceSymbols(innerpars, paste0(innerpars, "_C2"), trafo)),
       modelname = "prediction_p", outdir = tempdir())

# \donttest{
## Compile the observation function and the transformations together
compile(g, p, output = "prediction_ex")
#> Warning: A shared library named 'prediction_ex' is already loaded; overwriting it is not portable. Using 'prediction_ex_3' instead.
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC]
#> reusing 3 unchanged object(s)

outerpars <- getParameters(p)
pars <- setNames(seq(0.1, 1, length.out = length(outerpars)), outerpars)

## Unobserved states, then observables
plot((x * p)(times, pars))

plot((g * x * p)(times, pars))

# }
```
