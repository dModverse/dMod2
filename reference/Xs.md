# Prediction Function of an ODE Model with Sensitivities

Turns an
[`odemodel()`](https://dmodverse.github.io/dMod2/reference/odemodel.md)
into a prediction function that integrates the states and, with
`deriv = TRUE`, their sensitivities.

## Usage

``` r
Xs(odemodel, ...)

# S3 method for class 'deSolve'
Xs(
  odemodel,
  forcings = NULL,
  events = NULL,
  names = NULL,
  condition = NULL,
  options = list(),
  optionsSens = list(),
  fcontrol = NULL,
  ...
)

# S3 method for class 'cppDE'
Xs(
  odemodel,
  forcings = NULL,
  events = NULL,
  names = NULL,
  condition = NULL,
  options = list(),
  optionsSens = list(),
  optionsReverse = NULL,
  ...
)
```

## Arguments

- odemodel:

  An
  [`odemodel()`](https://dmodverse.github.io/dMod2/reference/odemodel.md).
  The backend of the model decides which method runs.

- ...:

  Not used; an argument not named here is an error.

- forcings:

  `data.frame` with columns `name` (character or factor), `time` and
  `value`, one row per data point of a forcing declared in
  `odemodel(forcings = )`. See section Forcings.

- events:

  `deSolve` backend only. An
  [eventlist](https://dmodverse.github.io/dMod2/reference/eventlist.md),
  or a `data.frame`
  [`as.eventlist()`](https://dmodverse.github.io/dMod2/reference/as.eventlist.md)
  accepts. Applied to the states only, the sensitivities are not
  corrected, so define events in
  [`odemodel()`](https://dmodverse.github.io/dMod2/reference/odemodel.md)
  unless the prediction is used for simulation alone. The other backends
  take events only from
  [`odemodel()`](https://dmodverse.github.io/dMod2/reference/odemodel.md).

- names:

  Character vector, the states and forcings to return. `NULL` returns
  all states followed by the forcings. On the `cppDE` and `Sundials`
  backends a name that is no state or forcing is dropped, and only an
  empty result is an error; on the `deSolve` backend such a name is an
  error.

- condition:

  `NULL` for a prediction valid in every condition, or the name of the
  condition it belongs to.

- options:

  Named list of solver options for every solve. See section Solver
  options.

- optionsSens:

  Named list of solver options for the solves with sensitivities,
  overriding the entries of `options` it names. See section Solver
  options.

- fcontrol:

  `deSolve` backend only. List with the interpolation settings of the
  forcings, `method`, `rule`, `f` and `ties` as in
  [`stats::approxfun()`](https://rdrr.io/r/stats/approxfun.html), passed
  to [`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html) as
  argument `fcontrol`.

- optionsReverse:

  Named list controlling the backward pass of `sweep = "reverse"`. Needs
  `backend = "cppDE"` and a model built with `derivMode = "reverse"`.
  `NULL`, the default, sweeps the step grid of the value pass. An
  unknown entry is an error. Entries:

  - `refine`: `TRUE` checks each step of the sweep against the
    tolerances of the solve and sweeps a failing step again in substeps,
    see
    [`cppDE::adjointControl()`](https://dmodverse.github.io/cppDE/reference/adjointControl.html).
    Default `FALSE`. The value and the `chi2` attribute of the reverse
    objective of
    [`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md)
    then include the first-order change from the computed trajectory to
    the one the sweep differentiates, so value and gradient belong to
    one function, which need not be closer to the exact value; a call
    with `deriv = FALSE` sweeps as well. The change is the attribute
    `correction` of the result. Forward values, the sweeps of
    [`em()`](https://dmodverse.github.io/dMod2/reference/em.md) and
    [`msem()`](https://dmodverse.github.io/dMod2/reference/msem.md)
    (forward over reverse) and multiple shooting run without refine;
    [`em()`](https://dmodverse.github.io/dMod2/reference/em.md) refuses
    a prediction with it.

  - `gradtol`: under `refine`, the absolute tolerance on each step's
    share of the gradient. Default `NULL`, no check on the gradient.

## Value

A [prdfn](https://dmodverse.github.io/dMod2/reference/prdfn.md), called
as `x(times, pars, ...)`, see section Calling a dMod function. It
returns a
[prdlist](https://dmodverse.github.io/dMod2/reference/prdlist.md) with
one [prdframe](https://dmodverse.github.io/dMod2/reference/prdframe.md)
per condition: the states over `times`, with the sensitivities in
`attr(, "deriv")`, `[time, state, parameter]`, when `deriv = TRUE`, and
the second-order sensitivities in `attr(, "deriv2")`,
`[time, state, parameter, parameter]`, when `deriv2 = TRUE` (`cppDE`
backend, model built with `deriv2 = TRUE`).

## Solver options

On the `cppDE` and `Sundials` backends the lists take the entries below
and pass them to
[`cppDE::solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.html).
An entry not given keeps the default; an unknown entry is an error, also
when set later with
[`controls()`](https://dmodverse.github.io/dMod2/reference/controls.md).

- `atol`, `rtol`: absolute and relative error tolerance, default `1e-6`.

- `maxsteps`: largest number of steps of one solve, default `1e6`.

- `maxattempts`: consecutive rejected steps before the solve fails,
  default `50`.

- `hini`: first step size, default `0`, which estimates it.

- `roottol`: tolerance of root-triggered events and of
  `rootfunc = "equilibrate"`, default `1e-6`.

- `maxroot`: number of times a root event may fire, default `1`.

- `onFailure`: `"stop"` (default), `"warn"` or `"silent"` on a failed
  solve. With `"warn"` and `"silent"` the result ends at the time
  reached, so an objective sees fewer data points rather than a failure.

- `traceFile`: CSV file the per-step trace is written to, default
  `NULL`; needs `odemodel(..., stepTrace = TRUE)`.

- `maxOrder`: upper bound on the order of the multistep methods, 1 to 5
  for `"bdf"` and 1 to 12 for `"adams"`; default `NULL`, the method's
  own.

- `stald`: stability limit detection for `"bdf"`, default `FALSE`; see
  [`cppDE::solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.html).

- `sensErrCon`: `optionsSens` only. `FALSE` takes the sensitivities out
  of the error test, which is cheaper and gives less accurate
  sensitivities; default `TRUE`.

`options` applies to every solve. The solves with sensitivities, those
of `deriv = TRUE` and of `sweep = "reverse"` with `deriv2 = TRUE`, take
`optionsSens` merged over `options`. A fit with the default forward
sweep therefore runs under `optionsSens`, and the value and backward
pass of `sweep = "reverse"` under `options` alone. The integration
method is fixed when the model is built, `odemodel(..., method = )`.

On the `deSolve` backend the entries of both lists are arguments of
[`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html), among them
`method`, `atol`, `rtol`, `maxsteps` and `hini`, merged the same way.
`method` defaults to `"lsoda"` for the solves without sensitivities and
to `"lsodes"` for those with, unless `options` or `optionsSens` names
it.

## Forcings

On the `cppDE` and `Sundials` backends each forcing is the monotone
cubic Hermite interpolant (PCHIP) of the forcing data points, see
[`cppDE::forcingValues()`](https://dmodverse.github.io/cppDE/reference/forcingValues.html).
Outside the points it holds the value of the first and the last one; a
single point gives a constant. The integrator does not stop at the
points, so a jump is best given by an event or by two close points.

On the `deSolve` backend the forcings are interpolated linearly and held
constant outside their data, unless `fcontrol` says otherwise.

The forcings are returned after the states, without derivatives, so
observables of [`Y()`](https://dmodverse.github.io/dMod2/reference/Y.md)
can contain them. Different forcings per condition come from one `Xs()`
per condition, joined by `+`.

## Controls

`forcings`, `names` and the solver options are kept, as given, in the
controls of the returned function and read at every solve, so
[`controls()`](https://dmodverse.github.io/dMod2/reference/controls.md)
can change them later; on the `deSolve` backend `events` and `fcontrol`
as well. The controls are named after the arguments, `options`,
`optionsSens` and `optionsReverse`. A replaced list is merged over the
defaults as described in section Solver options, not over the list it
replaces, so an entry it leaves out returns to the default.

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

[`Xf()`](https://dmodverse.github.io/dMod2/reference/Xf.md) for
predictions without sensitivities,
[`controls()`](https://dmodverse.github.io/dMod2/reference/controls.md),
[`cppDE::solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.html);
the guide chapters [Setting up
models](https://dmodverse.github.io/dMod2/articles/guide-model.html) and
[Sensitivities](https://dmodverse.github.io/dMod2/articles/guide-sensitivities.html),
the methods article [Sensitivities and
Hessians](https://dmodverse.github.io/dMod2/articles/methods-sensitivities.html)
and the cppDE article [Forward
sensitivities](https://dmodverse.github.io/cppDE/articles/methods-forward.html).

## Examples

``` r
# \donttest{
f <- eqnvec(A = "-k1*A", B = "k1*A - k2*B")
model <- odemodel(f, modelname = "Xs_AB", outdir = tempdir())
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC -fopenmp]
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC -fopenmp]

## Tighter tolerances for the solves with sensitivities, as used in a fit
x <- Xs(model, optionsSens = list(atol = 1e-8, rtol = 1e-8))
times <- seq(0, 10, 0.1)
pars <- c(A = 1, B = 0, k1 = 0.5, k2 = 0.2)
pred <- x(times, pars)
plot(pred)

dim(attr(pred[[1]], "deriv"))
#> [1] 101   2   4

## Change the options of the existing prediction function
controls(x, NULL, "options") <- list(atol = 1e-10, rtol = 1e-10)
head(x(times, pars, deriv = FALSE)[[1]])
#>      time         A          B
#> [1,]  0.0 1.0000000 0.00000000
#> [2,]  0.1 0.9512294 0.04828208
#> [3,]  0.2 0.9048374 0.09325337
#> [4,]  0.3 0.8607080 0.13509426
#> [5,]  0.4 0.8187308 0.17397599
#> [6,]  0.5 0.7788008 0.21006106
# }
```
