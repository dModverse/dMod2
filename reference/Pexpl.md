# Explicit Parameter Transformation

Builds a [parfn](https://dmodverse.github.io/dMod2/reference/parfn.md)
that evaluates inner parameters as expressions in outer parameters, with
their derivatives. It is evaluable only after compilation.

## Usage

``` r
Pexpl(
  trafo,
  parameters = NULL,
  attachInput = FALSE,
  condition = NULL,
  compile = FALSE,
  modelname = NULL,
  verbose = FALSE,
  deriv = TRUE,
  deriv2 = FALSE,
  derivMode = "forward",
  outdir = .dmodOutdir()
)
```

## Arguments

- trafo:

  Named character or
  [eqnvec](https://dmodverse.github.io/dMod2/reference/eqnvec.md); names
  are inner parameters, values are expressions in the outer parameters,
  in the language of
  [cppDE::expressions](https://dmodverse.github.io/cppDE/reference/expressions.html).

- parameters:

  Outer parameters; default `NULL`, the symbols of `trafo`. Names given
  here that `trafo` does not map pass through as identity.

- attachInput:

  Logical. Append the inputs that `trafo` does not map; varying ones
  with identity derivatives, fixed ones as fixed. Default `FALSE`. Can
  be changed with
  [`controls()`](https://dmodverse.github.io/dMod2/reference/controls.md).

- condition:

  `NULL` (default) for a transformation valid in every condition, or the
  name of the condition it belongs to.

- compile:

  Logical. `TRUE` compiles now, `FALSE` (default) leaves it to
  [`compile()`](https://dmodverse.github.io/dMod2/reference/compile.md).

- modelname:

  Base name of the generated files, default `NULL`, which is
  `"expl_parfn"`, followed by `_<condition>`.

- verbose:

  Logical. Print the compiler output, default `FALSE`.

- deriv, deriv2:

  Build the first derivatives, `attr(, "deriv")`, `[inner, outer]`, and
  the second, `attr(, "deriv2")`, `[inner, outer, outer]`. Default
  `TRUE` and `FALSE`; `deriv2` needs `deriv = TRUE`.

- derivMode:

  Which derivative code to build, any of `"forward"` (default),
  `"reverse"` and `"forward-reverse"`. `"reverse"` is needed for
  `sweep = "reverse"`, `"forward-reverse"` for `sweep = "reverse"` with
  `deriv2 = TRUE`.

- outdir:

  Directory for the generated source and shared object, default
  `getOption("dMod2.outdir")`, else the working directory.

## Value

A [parfn](https://dmodverse.github.io/dMod2/reference/parfn.md), called
as `p(pars, ...)`, see section Calling a dMod function.

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
[`vcov()`](https://rdrr.io/r/stats/vcov.html) solves it for the columns
of `which`.

With `deriv = FALSE` it returns the value alone.

## See also

[`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md),
[`P()`](https://dmodverse.github.io/dMod2/reference/P.md)

## Examples

``` r
# \donttest{
trafo <- eqnvec(k1 = "exp(log_k1)", k2 = "exp(log_k2)", A = "A0")
p <- Pexpl(trafo, condition = "C1", compile = TRUE,
           modelname = "Pexpl_example", outdir = tempdir())
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC -fopenmp]
out <- p(c(log_k1 = 0, log_k2 = -1, A0 = 2))
out$C1
#> Parameter vector:
#>   k1 :       1
#>   k2 :  0.367879
#>    A :       2
#> 
#> Attributes:
#>   deriv  : 3 x 3 matrix
#>   deriv2 : <none>
#>   fixed  : <none>
attr(out$C1, "deriv")
#>    log_k1    log_k2 A0
#> k1      1 0.0000000  0
#> k2      0 0.3678794  0
#> A       0 0.0000000  1
# }
```
