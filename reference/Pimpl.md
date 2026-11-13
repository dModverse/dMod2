# Implicit Parameter Transformation

Solves \\f(x, p) = 0\\ for the states \\x\\ by pseudo-transient
continuation ([`cppDE::ptc()`](https://rdrr.io/pkg/cppDE/man/ptc.html))
and returns a
[parfn](https://dmodverse.github.io/dMod2/reference/parfn.md) with
sensitivities from the implicit function theorem.

## Usage

``` r
Pimpl(
  trafo,
  parameters = NULL,
  forcings = NULL,
  condition = NULL,
  keepRoot = TRUE,
  flow = inherits(trafo, "eqnlist"),
  compile = FALSE,
  modelname = NULL,
  verbose = FALSE,
  deriv = TRUE,
  deriv2 = FALSE,
  controlsPTC = list(),
  outdir = .dmodOutdir()
)
```

## Arguments

- trafo:

  Named character,
  [eqnvec](https://dmodverse.github.io/dMod2/reference/eqnvec.md) or
  [eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md)
  defining \\f\\, in the language of
  [cppDE::expressions](https://rdrr.io/pkg/cppDE/man/expressions.html).

- parameters:

  Outer parameter names. For an
  [eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md) the
  totals \\T\\ of the conserved quantities are added. States named here
  are held at their value instead of solved for, and a state whose
  right-hand side is `"0"` becomes a parameter.

- forcings:

  Forcing names, set to 0. Default `NULL`.

- condition:

  `NULL` (default) for a transformation valid in every condition, or the
  name of the condition it belongs to.

- keepRoot:

  If `TRUE` (default), roots are kept and repeated calls with the same
  values are answered from memory.

- flow:

  If `TRUE`, \\\dot{x} = f(x, p)\\ and a stable steady state is
  returned. If `FALSE`, any regular root. Default `TRUE` for an
  [eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md),
  else `FALSE`.

- compile:

  Logical. `TRUE` compiles now, `FALSE` (default) leaves it to
  [`compile()`](https://dmodverse.github.io/dMod2/reference/compile.md).

- modelname:

  Base name of the generated files, default `NULL`, which is
  `"impl_parfn"`, followed by `_<condition>`.

- verbose:

  Logical. Print the compiler output, default `FALSE`.

- deriv, deriv2:

  Build first and second derivatives, default `TRUE` and `FALSE`.
  `deriv2` requires `deriv = TRUE`.

- controlsPTC:

  Named list of solver controls; an entry not given keeps the default,
  an unknown entry is an error:

  `reltol`, `abstol`

  :   Convergence tolerances, default `1e-10`, `1e-14`.

  `flowTol`

  :   Relative local error per pseudo-time step, default `1`.

  `maxit`

  :   Iterations per start, default `ceiling(70 * log(n + 1))` for `n`
      states.

  `dtInit`

  :   Initial pseudo-time step relative to the fastest rate, default
      `1e-2`.

  `positive`

  :   Keep states positive, default `TRUE`.

  `stability`

  :   Require a stable root if `flow = TRUE`, default `TRUE`; where no
      start reaches one, the evaluation fails. With `FALSE` an unstable
      root is returned, on which a
      [`constraintStability()`](https://dmodverse.github.io/dMod2/reference/constraintStability.md)
      term rises.

  `stabilityTol`

  :   With `stability`, a root counts as stable while the real part of
      the leading eigenvalue on the manifold of the conserved quantities
      is at most `stabilityTol`, in the units of the rates; default `0`.
      A small positive value admits a resting state that departs from
      rest more slowly than the experiment lasts, as near a Hopf
      bifurcation.

  `nStarts`

  :   Random starts after the deterministic ones, default `20`.

  `startRange`

  :   Range of a random start, drawn uniformly per state on
      `startScale`, default `c(-5, 5)`.

  `startScale`

  :   `"log10"` or `"linear"`, default `"log10"` if `positive` and
      `"linear"` otherwise.

  `seed`

  :   Seed of the random starts, default `1`. The starts depend on the
      seed and the values; the global RNG is not used.

- outdir:

  Directory for the generated files, default
  `getOption("dMod2.outdir")`, else the working directory.

## Value

A [parfn](https://dmodverse.github.io/dMod2/reference/parfn.md), called
as `p(pars, ...)`, see section Calling a dMod function.

## Details

Conserved quantities of an
[eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md) enter
as \\C x = T\\. States without influx are 0. A state whose reactions all
vanish with them rests: it takes the value from `pars`, 0 where `pars`
has none, and is not among the parameters. The root is solved from the
states in `pars`, missing ones at 1, then from all states at 1, then
from `nStarts` random starts; with `flow = TRUE` the stable root the
flow leads to. If all fail, it is an error. A root depends on the
parameters, the states in `pars` and the controls only, not on the calls
before it. Identical values are solved once, also across conditions, and
conditions are solved in parallel with the cores of the evaluation.
`keepRoot` and `controlsPTC` can be changed with
[`controls()`](https://dmodverse.github.io/dMod2/reference/controls.md).

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

[`Pexpl()`](https://dmodverse.github.io/dMod2/reference/Pexpl.md),
[`P()`](https://dmodverse.github.io/dMod2/reference/P.md); the guide
chapter [Setting up
models](https://dmodverse.github.io/dMod2/articles/guide-model.html).

## Examples

``` r
# \donttest{
## Steady state of A <-> B; the total of A and B enters as a parameter
r <- addReaction(NULL, "A", "B", "k1*A")
r <- addReaction(r, "B", "A", "k2*B")
p <- Pimpl(r, condition = "C1", compile = TRUE,
           modelname = "Pimpl_example", outdir = tempdir())
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC -fopenmp]
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC -fopenmp]
getParameters(p)
#> [1] "k1"      "k2"      "total_1"
pars <- setNames(c(1, 0.5, 3), getParameters(p))
out <- p(pars)
out$C1
#> Parameter vector:
#>         A :       1
#>         B :       2
#>        k1 :       1
#>        k2 :     0.5
#>   total_1 :       3
#> 
#> Attributes:
#>   deriv  : 5 x 3 matrix
#>   deriv2 : <none>
#>   fixed  : <none>
attr(out$C1, "deriv")
#>                 k1        k2   total_1
#> A       -0.6666667  1.333333 0.3333333
#> B        0.6666667 -1.333333 0.6666667
#> k1       1.0000000  0.000000 0.0000000
#> k2       0.0000000  1.000000 0.0000000
#> total_1  0.0000000  0.000000 1.0000000
# }
```
