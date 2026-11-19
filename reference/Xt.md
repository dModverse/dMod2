# Prediction Function Returning Time

Builds a prediction function whose only state is `time`, so that an
observation function can define a model as an explicit function of time,
without an ODE.

## Usage

``` r
Xt(condition = NULL)
```

## Arguments

- condition:

  `NULL` for a prediction valid in every condition, or the name of the
  condition it belongs to.

## Value

A [prdfn](https://dmodverse.github.io/dMod2/reference/prdfn.md), called
as `x(times, pars, ...)`, see section Calling a dMod function. The
derivatives of the prediction are zero.

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

With `deriv = FALSE` it returns the value alone.

## See also

[`Y()`](https://dmodverse.github.io/dMod2/reference/Y.md)

## Examples

``` r
# \donttest{
x <- Xt()
g <- Y(c(y = "a*time^2+b"), f = NULL, parameters = c("a", "b"),
       compile = TRUE, modelname = "Xt_example_obs", outdir = tempdir())
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC -fopenmp]

times <- seq(-1, 1, by = .05)
pars <- c(a = .1, b = 1)

plot((g*x)(times, pars))

# }
```
