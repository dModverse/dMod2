# Observation Function

Builds an [obsfn](https://dmodverse.github.io/dMod2/reference/obsfn.md)
that evaluates observables, and their first and second derivatives, on
the output of a prediction function such as
[`Xs()`](https://dmodverse.github.io/dMod2/reference/Xs.md).

## Usage

``` r
Y(
  g,
  f = NULL,
  states = NULL,
  parameters = NULL,
  condition = NULL,
  attachInput = FALSE,
  compile = FALSE,
  modelname = NULL,
  verbose = FALSE,
  cores = NULL,
  deriv = TRUE,
  deriv2 = FALSE,
  derivMode = "forward",
  outdir = .dmodOutdir()
)
```

## Arguments

- g:

  Named character vector or
  [eqnvec](https://dmodverse.github.io/dMod2/reference/eqnvec.md)
  defining the observation function, or a list of those named by
  condition, which builds one observation function per condition in a
  single call. Expressions follow
  [cppDE::expressions](https://dmodverse.github.io/cppDE/reference/expressions.html).

- f:

  Equations, anything
  [`as.eqnvec()`](https://dmodverse.github.io/dMod2/reference/as.eqnvec.md)
  accepts, or a prediction function. The states and parameters are read
  off `f`; for a prediction function the forcings count as states.

- states:

  Character vector of states, added to those read off `f`.

- parameters:

  Character vector of parameters, added to the symbols of `g` and `f`
  that are neither states nor `time`.

- condition:

  `NULL` for an observation function valid in every condition, or the
  name of the condition it belongs to.

- attachInput:

  Logical, append the incoming states to the output. Default `FALSE`;
  set `TRUE` to plot states alongside observables. Can be changed with
  [`controls()`](https://dmodverse.github.io/dMod2/reference/controls.md).

- compile:

  Logical. `TRUE` compiles the generated code now, `FALSE` (default)
  leaves it to
  [`compile()`](https://dmodverse.github.io/dMod2/reference/compile.md).
  The function is evaluable only compiled.

- modelname:

  Character, the base name of the generated C++ file and shared object.
  Default `NULL`, which is `"obsfn"`, followed by `_<condition>`.

- verbose:

  Logical, print compiler output to the R console.

- cores:

  Number of parallel jobs used to generate and compile the sources when
  `g` is a list, on Windows only to compile them; `NULL` auto-detects.
  Ignored for a single observation function, which is one source either
  way.

- deriv:

  Logical. If `TRUE` (default), build the first derivatives, returned as
  `attr(, "deriv")`, `[time, observable, parameter]`, when the call has
  `deriv = TRUE` and `derivMode` contains `"forward"`. `FALSE` returns
  none, whatever the call asks.

- deriv2:

  Logical. If `TRUE`, build the second derivatives, returned as
  `attr(, "deriv2")`, `[time, observable, parameter, parameter]`, when
  the call has `deriv2 = TRUE`; on a function built without them that
  call is an error. Requires `deriv = TRUE`. Default `FALSE`.

- derivMode:

  Which derivative code to build, any of `"forward"` (default),
  `"reverse"` and `"forward-reverse"`.

  - `"forward"`: the Jacobian; without it the function returns no
    `deriv` attribute.

  - `"reverse"`: needed for `sweep = "reverse"`.

  - `"forward-reverse"`: needed for `sweep = "reverse"` with
    `deriv2 = TRUE`.

- outdir:

  Character. Directory for the generated C++ source and the shared
  object, default `getOption("dMod2.outdir")`, else the working
  directory.

## Value

An [obsfn](https://dmodverse.github.io/dMod2/reference/obsfn.md), called
as `g(out, pars, ...)` with a
[prdlist](https://dmodverse.github.io/dMod2/reference/prdlist.md) `out`,
see section Calling a dMod function, or composed as `g * x`. It returns
a [prdlist](https://dmodverse.github.io/dMod2/reference/prdlist.md) of
the observables and their derivatives.

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
  `Y()` and [`P()`](https://dmodverse.github.io/dMod2/reference/P.md).

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

[`compile()`](https://dmodverse.github.io/dMod2/reference/compile.md),
[`controls()`](https://dmodverse.github.io/dMod2/reference/controls.md)

## Examples

``` r
## A prediction function interpolating parameter values piecewise linearly
times <- 0:5
grid <- data.frame(name = "A", time = times, row.names = paste0("p", times))
x <- Xd(grid)

## An observation function
g <- Y(c(Aobs = "s*A"), states = "A", parameters = "s",
       modelname = "prediction_obs", outdir = tempdir())

## One parameter transformation per condition
innerpars <- c(getParameters(x), getParameters(g))
trafo <- setNames(innerpars, innerpars)
p <- P(list(C1 = replaceSymbols(innerpars, paste0(innerpars, "_C1"), trafo),
            C2 = replaceSymbols(innerpars, paste0(innerpars, "_C2"), trafo)),
       modelname = "prediction_p", outdir = tempdir())

# \donttest{
## Compile the observation function and the transformations together
compile(g, p, output = "prediction_ex")
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC]

outerpars <- getParameters(p)
pars <- setNames(seq(0.1, 1, length.out = length(outerpars)), outerpars)

## Unobserved states, then observables
plot((x * p)(times, pars))

plot((g * x * p)(times, pars))

# }
```
