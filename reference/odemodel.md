# Compiled ODE Model

Generates and compiles the code of an ODE system and, as requested, of
the first- and second-order sensitivities and the adjoint.
[`Xs()`](https://dmodverse.github.io/dMod2/reference/Xs.md) or
[`Xf()`](https://dmodverse.github.io/dMod2/reference/Xf.md) turn the
result into a prediction function.

## Usage

``` r
odemodel(
  f,
  deriv = TRUE,
  deriv2 = FALSE,
  derivMode = "forward",
  forcings = NULL,
  events = NULL,
  fixed = NULL,
  modelname = "odemodel",
  backend = c("cppDE", "Sundials", "deSolve"),
  verbose = FALSE,
  outdir = .dmodOutdir(),
  ...
)
```

## Arguments

- f:

  Named character vector,
  [eqnvec](https://dmodverse.github.io/dMod2/reference/eqnvec.md) or
  [eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md) with
  the right-hand sides; the names are the states. The functions and
  operators an equation may use are listed in
  [cppDE::expressions](https://dmodverse.github.io/cppDE/reference/expressions.html).

- deriv:

  logical, build first-order sensitivities. Default `TRUE`. `FALSE`
  builds the states alone, for
  [`Xf()`](https://dmodverse.github.io/dMod2/reference/Xf.md), and
  `derivMode` has no effect.

- deriv2:

  logical, build second-order sensitivities. Default `FALSE`. `TRUE` is
  the same as adding `"forward-forward"` to `derivMode`, needs
  `backend = "cppDE"` and implies `deriv = TRUE`.

- derivMode:

  character, the derivative directions to build, one or more of:

  - `"forward"` (default): sensitivity equations.

  - `"reverse"`: adjoint gradient, needed for
    `obj(..., sweep = "reverse")`. Not available on
    `backend = "deSolve"`.

  - `"forward-forward"`: second-order sensitivities.

  - `"forward-reverse"`: adjoint Hessian, needed for
    `obj(..., sweep = "reverse", deriv2 = TRUE)`. `backend = "cppDE"`
    only.

  Each direction is compiled separately.

- forcings:

  character, the names of the forcings in `f`. Default `NULL`. Their
  data are given to
  [`Xs()`](https://dmodverse.github.io/dMod2/reference/Xs.md) or
  [`Xf()`](https://dmodverse.github.io/dMod2/reference/Xf.md).

- events:

  an
  [eventlist](https://dmodverse.github.io/dMod2/reference/eventlist.md),
  or a `data.frame` accepted by
  [`as.eventlist()`](https://dmodverse.github.io/dMod2/reference/as.eventlist.md).
  Default `NULL`. Events defined here also act on the sensitivities and
  the adjoint.

- fixed:

  character, the initial values and parameters without sensitivities.
  Default `NULL`.

- modelname:

  character, the base name of the generated files and symbols. Default
  `"odemodel"`. Files of the same name in `outdir` are overwritten.

- backend:

  `"cppDE"` (default), `"Sundials"` or `"deSolve"`: the integrators of
  cppDE, CVODE(S) through cppDE, or deSolve through cOde.

- verbose:

  logical, print the compiler output. Default `FALSE`.

- outdir:

  directory for the generated sources and shared objects. Default
  `getOption("dMod2.outdir")`, else the working directory.
  `backend = "deSolve"` writes to the working directory only.

- ...:

  Passed to
  [`cppDE::cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.html),
  [`cppDE::cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.html)
  or [`cOde::funC()`](https://rdrr.io/pkg/cOde/man/funC.html), according
  to `backend`. Among them:

  - `method`: the integration method, `"bdf"` (default), `"adams"`,
    `"rb4"` or `"tsit5"` on `cppDE`; `"bdf"` or `"adams"` on `Sundials`.
    On `deSolve` it is `options$method` of
    [`Xs()`](https://dmodverse.github.io/dMod2/reference/Xs.md) instead.

  - `rootfunc`: `"equilibrate"` or expressions whose roots end the
    integration (`cppDE`, `Sundials`).

  - `includeTimeZero`: integrate from 0 rather than from the first
    requested time, default `TRUE`.

  - `sparse`: sparse or dense linear algebra, `NULL` chooses from the
    Jacobian (`cppDE`, `Sundials`).

  - `stepTrace`: record per-step diagnostics, see `traceFile` in
    [`Xs()`](https://dmodverse.github.io/dMod2/reference/Xs.md)
    (`cppDE`, `Sundials`).

  - `compile`: default `TRUE`. `FALSE` writes the sources only;
    [`compile()`](https://dmodverse.github.io/dMod2/reference/compile.md)
    builds them later.

  - `estimate`, `outputs`, `gridpoints`: the parameters with
    sensitivities, additional outputs and spline nodes of the forcings
    (`deSolve`).

## Value

An object of class `odemodel`: a list with the compiled objects `func`
(states), `extended` (first-order sensitivities), and on `cppDE` /
`Sundials` `extended2` (second order) and `reversed`, `reversed2`
(adjoints), each `NULL` where not built. The attribute `"compileInfo"`
holds the sources and compiler flags
[`compile()`](https://dmodverse.github.io/dMod2/reference/compile.md)
uses.

## See also

[`Xs()`](https://dmodverse.github.io/dMod2/reference/Xs.md),
[`Xf()`](https://dmodverse.github.io/dMod2/reference/Xf.md),
[`compile()`](https://dmodverse.github.io/dMod2/reference/compile.md),
[`eventlist()`](https://dmodverse.github.io/dMod2/reference/eventlist.md),
[cppDE::expressions](https://dmodverse.github.io/cppDE/reference/expressions.html),
[`cppDE::cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.html),
[`cppDE::cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.html),
[`cOde::funC()`](https://rdrr.io/pkg/cOde/man/funC.html)

## Examples

``` r
# \donttest{
## A model with a forcing F, no sensitivities for "switch"
f <- eqnvec(A = "-k*A + switch*F")
model <- odemodel(f, forcings = "F", fixed = "switch", modelname = "odemodel_F",
                  outdir = tempdir(), compile = FALSE)
print(model)
#> dMod odemodel
#>   Backend: cppDE [bdf]
#>   Model:   odemodel_F
#>   Sens1:   odemodel_F_s
#>   Sens2:   not compiled (deriv2 = FALSE)
#> 
#> Equations:
#> Idx Inner <- Outer
#>   1     A <- -k*A+switch*F
#> 
#> States:
#> [1] "A"
#> 
#> Parameters:
#> [1] "k"      "switch"
#> 
#> Forcings:
#> [1] "F"

## The same model from an equation list
f2 <- addReaction(NULL, from = "", to = "A", rate = "switch*F", description = "production")
f2 <- addReaction(f2, from = "A", to = "", rate = "k*A", description = "degradation")
print(f2)
#> 
#>   Check Educt -> Product     Rate Description
#> 1             ->       A switch*F  production
#> 2           A ->              k*A degradation
model2 <- odemodel(f2, forcings = "F", fixed = "switch", modelname = "odemodel_F2",
                   outdir = tempdir(), compile = FALSE)

## One prediction per forcing; the data of each forcing cover the time range
tF <- seq(0, 5, 0.1)
x <- Xs(model, data.frame(name = "F", time = tF, value = sin(tF)), condition = "sin") +
  Xs(model, data.frame(name = "F", time = tF, value = exp(-tF)), condition = "exp")
x2 <- Xs(model2, data.frame(name = "F", time = c(0, 5), value = 0.1), condition = "const")

## Compile both models into one shared object
compile(x, x2, output = "odemodel_ex")
#> using C++ compiler: g++
#>   every source : -O2 -DNDEBUG -w -fPIC
#>   4 of 4 also  : -fopenmp

pars <- c(A = 1, k = 0.5, switch = 1)
plot((x + x2)(seq(0, 5, 0.01), pars))

# }
```
