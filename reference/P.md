# Parameter Transformation

Builds a [parfn](https://dmodverse.github.io/dMod2/reference/parfn.md)
with [`Pexpl()`](https://dmodverse.github.io/dMod2/reference/Pexpl.md)
or [`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md).

## Usage

``` r
P(
  trafo = NULL,
  parameters = NULL,
  condition = NULL,
  compile = FALSE,
  modelname = NULL,
  method = NULL,
  cores = NULL,
  verbose = FALSE,
  deriv = TRUE,
  deriv2 = FALSE,
  ...
)
```

## Arguments

- trafo:

  An [eqnvec](https://dmodverse.github.io/dMod2/reference/eqnvec.md),
  named character or
  [eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md), or
  a list of those named by condition, which builds one transformation
  per condition. `NULL` (default) returns `NULL`. Expressions follow
  [cppDE::expressions](https://rdrr.io/pkg/cppDE/man/expressions.html).

- parameters:

  Outer parameter names. Default `NULL`, the symbols of `trafo`.

- condition:

  Condition label of a single `trafo`, default `NULL`; a list takes the
  list names instead.

- compile:

  Logical. `TRUE` compiles all conditions into one shared object named
  `modelname`, `FALSE` (default) leaves it to
  [`compile()`](https://dmodverse.github.io/dMod2/reference/compile.md).

- modelname:

  Base name of the generated files, default `NULL`, which is
  `"expl_parfn"` or `"impl_parfn"`, followed by `_<condition>`.

- method:

  `"explicit"`, `"implicit"`, or `NULL` (default) for `"implicit"` if
  `trafo` is an
  [eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md) and
  `"explicit"` otherwise.

- cores:

  Number of parallel jobs generating and compiling the sources of a
  list; on Windows the sources are generated one after the other. `NULL`
  (default) uses
  [`detectFreeCores()`](https://dmodverse.github.io/dMod2/reference/detectFreeCores.md).

- verbose:

  Logical. Print the compiler output, default `FALSE`.

- deriv, deriv2:

  Build first and second derivatives, default `TRUE` and `FALSE`.
  `deriv2` requires `deriv = TRUE`.

- ...:

  Further arguments of
  [`Pexpl()`](https://dmodverse.github.io/dMod2/reference/Pexpl.md)
  (`attachInput`, `derivMode`, `outdir`) or
  [`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md)
  (`forcings`, `keepRoot`, `flow`, `controlsPTC`, `outdir`), as `method`
  selects.

## Value

A [parfn](https://dmodverse.github.io/dMod2/reference/parfn.md), called
as `p(pars, ...)`, see section Calling a dMod function. The random
effects of a transformation built by
[`insert()`](https://dmodverse.github.io/dMod2/reference/define.md) stay
attached, so
[`randomEffects()`](https://dmodverse.github.io/dMod2/reference/randomEffects.md)
reads them from the result.

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
  [`Y()`](https://dmodverse.github.io/dMod2/reference/Y.md) and `P()`.

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

[`Pexpl()`](https://dmodverse.github.io/dMod2/reference/Pexpl.md),
[`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md),
[`define()`](https://dmodverse.github.io/dMod2/reference/define.md)

## Examples

``` r
# \donttest{
## One transformation per condition, all compiled into one shared object
trafo <- list(
  ctrl  = eqnvec(k1 = "exp(log_k1)", k2 = "exp(log_k2)", A = "A0"),
  treat = eqnvec(k1 = "exp(log_k1 + log_fold)", k2 = "exp(log_k2)", A = "A0"))
p <- P(trafo, compile = TRUE, modelname = "P_example", outdir = tempdir())
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC]
getParameters(p)
#> [1] "log_k1"   "log_k2"   "A0"       "log_fold"
p(c(log_k1 = 0, log_k2 = -1, A0 = 2, log_fold = log(3)))
#> $ctrl
#> Parameter vector:
#>   k1 :       1
#>   k2 :  0.367879
#>    A :       2
#> 
#> Attributes:
#>   deriv  : 3 x 3 matrix
#>   deriv2 : <none>
#>   fixed  : <none>
#> 
#> $treat
#> Parameter vector:
#>   k1 :       3
#>   k2 :  0.367879
#>    A :       2
#> 
#> Attributes:
#>   deriv  : 3 x 4 matrix
#>   deriv2 : <none>
#>   fixed  : <none>
#> 
# }
```
