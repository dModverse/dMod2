# Non-Linear Optimization by a Trust-Region Method

`trust` minimizes (or maximizes) a smooth objective function for which
value, gradient and Hessian are available. It solves a Moré-Sorensen
trust-region subproblem exactly, via the full eigendecomposition of the
Hessian.

## Usage

``` r
trust(
  objfun,
  parinit,
  rinit = 0.1,
  rmax = 10,
  iterlim = 1000L,
  hessianMethod = c("gn", "bfgs", "sr1", "exact"),
  phases = NULL,
  switches = NULL,
  parupper = NULL,
  parlower = NULL,
  boundary = c("reflective", "clip"),
  tolControl = NULL,
  qnControl = NULL,
  stepControl = NULL,
  minimize = TRUE,
  blather = FALSE,
  printIter = FALSE,
  traceFile = NULL,
  sparse = "auto",
  ...
)

trustL1(
  objfun,
  parinit,
  mu = 0 * parinit,
  one.sided = FALSE,
  lambda = 1,
  rinit = 0.1,
  rmax = 10,
  iterlim = 1000L,
  parupper = NULL,
  parlower = NULL,
  boundary = c("reflective", "clip"),
  tolControl = NULL,
  stepControl = NULL,
  minimize = TRUE,
  blather = FALSE,
  printIter = FALSE,
  traceFile = NULL,
  gate = NULL,
  fuse = NULL,
  ...
)
```

## Arguments

- objfun:

  Function whose first argument is a numeric vector of parameters,
  returning a list with components `value`, `gradient` and `hessian`. It
  is called with the logical arguments `deriv` and `hessian`, and with
  `deriv2 = TRUE` when an exact Hessian is requested, so it must accept
  them or `...`. `trustL1` calls it with the parameters and `...` only.

- parinit:

  Named numeric starting vector. Must be finite. Values outside
  `[parlower, parupper]` are clipped with a warning; with
  `boundary = "reflective"` they are moved just inside the box.

- rinit:

  Initial trust-region radius. Default `0.1`.

- rmax:

  Maximum trust-region radius. Default `10`.

- iterlim:

  Maximum number of trust-region iterations over all phases. Default
  `1000`.

- hessianMethod:

  Source of the model Hessian of a run without `phases`: `"gn"`
  (default) the Gauss-Newton Hessian of the objective, `"bfgs"` or
  `"sr1"` a quasi-Newton update seeded as `qnControl$hessianInit` says,
  `"exact"` the true Hessian at every iterate, which needs an objective
  that takes `deriv2`. `"bfgs"`, `"sr1"` and `"exact"` require
  `boundary = "reflective"`. Also selects the control defaults below.
  With `phases`, each phase names a method, and this argument is an
  error.

- phases:

  List of phases, each a list of settings, run one after the other and
  from the first again after the last. The run hands over to the next
  phase at a value, model, step or stagnation stop and when a phase uses
  up the `iterlim` of that phase. A handover restarts the trust region
  at `rinit` and discards the quasi-Newton pairs. Default `NULL`: one
  phase from `hessianMethod` and `qnControl`, evaluated with the `sweep`
  in `...`. Names of the list label the phases in the result, the
  Hessian source does otherwise. Not available with multiple shooting
  and `boundary = "clip"`.

  `hessianMethod`

  :   Hessian source of the phase, as `hessianMethod`. Required.

  `iterlim`

  :   Maximum number of iterations of one run of the phase. Default
      `NULL`: only the total `iterlim`.

  `sweep`

  :   `"forward"` or `"reverse"`, the `sweep` the objective is evaluated
      with in the phase; a `sweep` in `...` is an error with `phases`.
      Default `NULL`: for `"bfgs"` and `"sr1"` reverse if the objective
      can be evaluated backwards at `parinit`, as with
      `derivMode = c("forward", "reverse")` throughout the chain,
      forward otherwise. A Gauss-Newton Hessian is always computed
      forward. A handover between phases of different sweep evaluates
      the current iterate once more. The result reports the sweeps as
      `phaseSweep`.

  `qnControl`

  :   Quasi-Newton settings of the phase, merged over `qnControl`; see
      there. `hessianInit` defaults to `"gn"` in every phase but the
      first: a quasi-Newton phase starts from the Gauss-Newton Hessian
      at the iterate of the handover.

- switches:

  Number of handovers allowed. Default `NULL`: one fewer than there are
  phases, so each runs once. A higher value cycles through the phases
  again, `0` runs the first phase only. The result reports the count as
  `nSwitch`.

- parupper, parlower:

  Upper and lower bounds. Unnamed: the first element applies to all
  parameters. Named: entries are matched by name, missing parameters are
  unbounded. Default `NULL`, unbounded.

- boundary:

  Box-bound handling, `"reflective"` (default) or `"clip"`. See section
  Box bounds.

- tolControl:

  List of convergence thresholds, merged into the defaults.

  `ftol`

  :   Change in objective value. Default `1e-6`.

  `mtol`

  :   Predicted model decrease. Default `1e-6`.

  `gtol`

  :   First-order optimality measure. Default `1e-6`.

  `xtol`

  :   Step norm. Default `0`, disabled.

  `rmin`

  :   Lower limit on the trust radius; falling below it stops the run
      with `converged = FALSE`. Default `0`, disabled.

  `ctol`

  :   Multiple shooting only: scaled gap size below which the trajectory
      counts as continuous; no other test stops the run before. Default
      `1e-6`.

  `trustL1` reads `ftol`, `mtol`, `gtol`, `xtol` and `rmin`.

- qnControl:

  List of quasi-Newton settings, merged into the defaults and read in a
  `"bfgs"` or `"sr1"` phase; `hessianInit` also in a first phase of
  source `"gn"`.

  `hessianInit`

  :   Seed of the approximation: `"gn"` (default) the Gauss-Newton
      Hessian at `parinit`; `"identity"` the identity matrix, after
      which the objective is asked for gradients only until a reseed;
      `"exact"` the true Hessian at `parinit`, then gradients only.
      `"exact"` needs `deriv2`, see `derivMode` in
      [`odemodel()`](https://dmodverse.github.io/dMod2/reference/odemodel.md).
      In a first phase of source `"gn"`, `"exact"` makes the first step
      use the exact Hessian at `parinit`.

  `hessianReseed`

  :   `"never"` (default): a stalled quasi-Newton phase, or one whose
      trust radius collapses, stops. `"stall"`: it continues from a
      fresh Hessian at the current iterate, the exact one with
      `hessianInit = "exact"` and the Gauss-Newton one otherwise,
      provided the iterate moved since the last reseed. The result
      reports the count as `nReseed`.

  `qnMemory`

  :   Number of `(s, y)` pairs kept. Default `0`: every update
      accumulates onto the seed. A positive value rebuilds the
      approximation in each iteration from the last `qnMemory` pairs,
      without the seed.

  `qnCautious`

  :   A BFGS pair is used only when \\s^\top y \> c\\\|s\|\\\|y\|\\ with
      \\c\\ = `qnCautious`. Default `1e-8`, `0` disables the test;
      `"sr1"` is exempt. The result reports the rejected pairs as
      `qnSkipped`.

  `qnRejected`

  :   Whether `"sr1"` also updates from a rejected trial point. Default
      `TRUE` for a phase whose source is `"sr1"`, `FALSE` otherwise.

- stepControl:

  List of step and acceptance settings, merged into the defaults.

  `thetaMax`

  :   Smallest fraction of the distance to a bound a step may use. The
      fraction is `max(thetaMax, 1 - m)`, below 1, with `m` the
      first-order optimality measure. Default `0.99995`.
      `boundary = "reflective"` only.

  `nonmonotone`

  :   Zhang-Hager relaxation in `[0, 1)`: a step is rated against a
      weighted average of past objective values and may increase the
      objective. Default `0`, monotone.

  `acceptance`

  :   Multiple shooting only: `"filter"` (default), the filter of
      Fletcher and Leyffer on objective and gap size; `"merit"`, the l2
      exact penalty; `"natural"`, Bock's natural level function, which
      needs `hessianMethod = "gn"`, a fixed sigma and no error model.

  `twoPhase`

  :   Multiple shooting only: if `TRUE`, the parameters are first fitted
      with the node values held fixed, then the nodes are laid out
      afresh and released. Default `FALSE`.

  `restore`

  :   Multiple shooting only: if `TRUE`, a run whose steps keep being
      rejected with open gaps enters a restoration phase, at most five
      times. Default `FALSE`.

  `regularize`

  :   Multiple shooting only: directions of the condensed problem with a
      singular value below `regularize` times the largest are held fixed
      in a step. Default `0`, off. With `acceptance = "natural"`
      directions below `1e-7` are always held fixed.

  `soc`

  :   Multiple shooting only: if `TRUE`, a rejected step whose gaps did
      not close as predicted is retried once with a second-order
      correction. Default `TRUE`.

  `anneal`

  :   Multiple shooting only: if `TRUE`, the run first minimizes the
      objective plus a rising penalty on the gaps, then enforces
      continuity. Needs forward sensitivities and is skipped without
      them. Default `TRUE`.

  `trustL1` takes `thetaMax` only.

- minimize:

  If `TRUE` (default) minimize, if `FALSE` maximize.

- blather:

  If `TRUE`, the result also contains the per-iteration trace listed
  under Value. Default `FALSE`.

- printIter:

  If `TRUE`, print the number of the evaluation and the objective value
  at each evaluation. Default `FALSE`.

- traceFile:

  Path of a CSV file written with one row per evaluation: the number of
  the evaluation (column `Iteration`), objective value and parameters.
  Default `NULL`, no file.

- sparse:

  `"auto"` (default), `TRUE` or `FALSE`. With `TRUE`, `objfun` is called
  with `hessian = "sparse"`, and the subproblem of every step is solved
  on sparse Cholesky factors of the Hessian it returns (a `Matrix`), a
  low-rank part `lowrank` (see
  [`objlist()`](https://dmodverse.github.io/dMod2/reference/objlist.md))
  added by Woodbury, where the dense path takes the eigendecomposition:
  the same steps to rounding, with every setting of `trust()`; the
  Hessian comes from the objective, `hessianMethod` `"gn"` or `"exact"`,
  and `boundary = "reflective"`. `"auto"` takes the sparse path where
  that holds, outside multiple shooting, when `parinit` has at least
  `parameters` entries and a step on the sparse Hessian at `parinit`
  costs fewer flops than the dense one, by the symbolic factorization;
  `parameters` (default 100) comes from `getOption("dMod2.sparseAuto")`,
  where a `fill` replaces the cost by the rule that the sparse Hessian
  fills at most `fill` of the matrix.

- ...:

  Further named arguments passed to `objfun`.

- mu:

  Named numeric vector of reference values for the L1-penalized
  parameters. Names must be a subset of `names(parinit)`; only the named
  parameters receive a penalty. With `boundary = "reflective"` each must
  lie strictly inside `[parlower, parupper]`. Defaults to a zero vector
  covering all of `parinit`.

- one.sided:

  Logical. If `TRUE`, the penalty is one-sided,
  `lambda * max(0, mu - p)`: only values below `mu` are penalized.
  Otherwise it is the symmetric `lambda * |p - mu|`.

- lambda:

  Strength of the L1 penalty. Either a scalar (broadcast to all entries
  of `mu`) or a named numeric aligned with `mu`.

- gate:

  Character, names in `mu` whose penalty is the gate `lambda * (p - mu)`
  with `p >= mu`: the parameter can sit exactly on `mu` and leaves it
  only upward. Starts below `mu` are moved onto it. Needs
  `boundary = "reflective"`.

- fuse:

  List of fusion blocks. A block penalizes
  `lambda * sum_{i<j} w_ij |p_i - p_j|` over the parameters of the
  block, plus `lambda * w_ia |p_i - a|` toward an optional anchor `a`.
  Parameters with equal values move as one and split again when the data
  pull them apart. A block is a character vector of parameter names or a
  list with `pars`, and optionally `lambda` (default: the scalar
  `lambda`), `anchor` (default none) and `weights`, an `(m+1) x (m+1)`
  matrix whose last row and column refer to the anchor (default all 1).
  Fused parameters must not appear in `mu`; if `fuse` is given and `mu`
  is not, `mu` is empty. Needs `boundary = "reflective"`.

## Value

A list of class `trustfit` with components

- `argument`:

  Named numeric, the best iterate.

- `value`, `gradient`, `hessian`:

  Objective value, gradient and model Hessian at `argument`.

- `iterations`:

  Number of iterations.

- `neval`:

  Number of objective evaluations.

- `qnEval`:

  Evaluations spent in a quasi-Newton phase.

- `qnSkipped`:

  Pairs rejected by the cautious test.

- `nSwitch`:

  Number of handovers between phases.

- `phaseIterations`:

  Named integer, the iterations of each phase run, in run order and
  named after the phase. The sum is `iterations` without the iteration
  that ends the run on failed evaluations and without the iterations
  that reseed instead of taking a step.

- `phaseSweep`:

  Named character, the sweep of each element of `phases`, present when
  `phases` is set.

- `nReseed`:

  Number of reseeds by `qnControl$hessianReseed = "stall"`, which allows
  each run of a phase the same number.

- `evalBySource`:

  Named integer, evaluations per Hessian source.

- `converged`:

  Logical.

- `atBound`:

  Named logical, the parameters held by a bound.

- `stopReason`:

  The test that ended the run, see section Convergence.

With `blather = TRUE` the list also contains the matrices `argpath`
(accepted iterates) and `argtry` (trial points), one row per iteration,
and the per-iteration vectors `steptype`, `stepback`, `accept`, `r`
(radius), `rho` (ratio of actual to predicted reduction), `valpath`,
`valtry`, `preddiff`, `stepnorm` and, with `boundary = "reflective"`,
`hessianSource`.

`boundary = "clip"` returns only `argument`, `value`, `gradient`,
`hessian`, `iterations`, `converged`, `atBound` and `stopReason`, plus
the trace.

On the sparse path, see `sparse`, `hessian` is a sparse `Matrix`, and
`lowrank` where the objective returns one: the Hessian is then `hessian`
plus that part.

`trustL1` returns a list without class with `argument`, `value`,
`gradient`, `hessian`, `iterations`, `converged`, `atBound` and
`stopReason`, plus the trace without `hessianSource`. `value` and
`gradient` include the penalty, whose gradient counts zero on a kink;
`hessian` is that of `objfun`.

A multiple-shooting run returns `argument`, `value`, `gradient`,
`hessian`, `iterations`, `neval`, `qnSkipped`, `converged`,
`stopReason`, `atBound` and a list `multipleShooting` with the final
`nodes`, the `gaps` per condition, their scaled size `violation`, the
counts `nRelaxed` (iterations with a node step below the full one),
`nSOC` (accepted second-order corrections), `nSplit` (segments cut),
`nRestore` (restoration phases) and `nStages` (annealing stages), the
`cuts`, the `acceptance` and `sweep` used, `mu` (merit) or `filterSize`
(filter), and `trace` when `blather = TRUE`.

## Box bounds

`boundary = "reflective"` (default) uses the interior
trust-region-reflective scheme of Coleman and Li. Iterates stay strictly
inside the box and never land exactly on a bound; `atBound` in the
result reports which parameters a bound holds.

`boundary = "clip"`: active-set reduction with componentwise clipping,
with a single `"gn"` phase. The run stops when a trial point, accepted
or not, changes the value by less than `ftol` or has a predicted
reduction below `mtol`, both reported as `stopReason = "fvalue"`.
`gtol`, `xtol`, `rmin` and stagnation are not tested; three failed
evaluations in a row end the run with `stopReason = "iterlim"`.

## Convergence

With `boundary = "reflective"` the run stops as soon as one of the
following holds. The gradient test is evaluated at the current iterate;
the value and step tests count only on an accepted step, the model test
also on a rejected Newton step from inside the trust region.

- `gtol`: the first-order optimality measure \\\max_i \|v_i g_i\|\\ is
  at most `gtol`, with \\v_i\\ the distance to the bound in the descent
  direction \\-g_i\\, 1 if that bound is infinite. Without bounds it is
  \\\max_i \|g_i\|\\.

- `ftol`: \\\|f - f\_{\mathrm{try}}\| \<\\ `ftol`.

- `mtol`: the predicted reduction \\\|g^\top p + \frac{1}{2} p^\top H
  p\| \<\\ `mtol`.

- `xtol`: the step norm is below `xtol`.

- stagnation: five consecutive rejected steps left the objective flat
  within `ftol`.

- `rmin`: the trust radius is below `rmin`, or a step no longer changes
  the parameters in floating point. Reports `converged = FALSE`.

`stopReason` in the result names the test that fired: `"gradient"`,
`"fvalue"`, `"preddiff"`, `"step"`, `"stagnation"`, `"radius"`,
`"objfun"` (three failed evaluations in a row, ten with multiple
shooting) or `"iterlim"`. Stagnation reports `converged = TRUE`.

Only `gtol` tests first-order optimality. While `switches` remain, see
`phases`, an `"fvalue"`, `"preddiff"`, `"step"` or `"stagnation"` stop
hands over to the next phase instead of ending the run. In a
quasi-Newton phase `"fvalue"`, `"preddiff"` and `"step"` do not end the
run. A phase that reaches the `iterlim` of that phase ends like a soft
stop, with `stopReason = "iterlim"` and `converged = FALSE` when no
handover is left. The total `iterlim` ends the run in any phase, with
`stopReason = "iterlim"` and `converged = FALSE`.

`argument`, `value`, `gradient`, `hessian` and `atBound` in the result
belong to the best iterate visited, which is the last accepted one when
`stepControl$nonmonotone = 0`.

## Multiple shooting

An objective built by
[`normL2`](https://dmodverse.github.io/dMod2/reference/normL2.md)`(..., multipleShootingControl = )`,
alone or in a sum, is optimized by multiple shooting with the same
controls. Not available then: `hessianMethod = "exact"`, more than one
phase, `switches`, `qnControl$hessianInit = "exact"`,
`qnControl$hessianReseed`, `qnControl$qnMemory`,
`stepControl$nonmonotone`, `boundary = "clip"`, `minimize = FALSE` and
`traceFile`. The run reports `converged = TRUE` only when the scaled
gaps between the segments are below `tolControl$ctol`. The method is
described in the methods article [Multiple
shooting](https://dmodverse.github.io/dMod2/articles/methods-multiple-shooting.html).

## References

Coleman, T. F. and Li, Y. (1996). An interior trust region approach for
nonlinear minimization subject to bounds. *SIAM Journal on Optimization*
6(2), 418-445.

Fröhlich, F. and Sorger, P. K. (2022). Fides: Reliable trust-region
optimization for parameter estimation of ODE models. *PLoS Computational
Biology* 18(7), e1010322.

Bock, H. G. (1983). Recent advances in parameter identification
techniques for ODE. In *Numerical Treatment of Inverse Problems in
Differential and Integral Equations*, 95-121. Birkhäuser.

Fletcher, R. and Leyffer, S. (2002). Nonlinear programming without a
penalty function. *Mathematical Programming* 91, 239-269.

Horbelt, W., Timmer, J. and Voss, H. U. (2002). Parameter estimation in
nonlinear delayed feedback systems from noisy data. *Physics Letters A*
299, 513-521.

Voss, H. U., Timmer, J. and Kurths, J. (2004). Nonlinear dynamical
system identification from uncertain and indirect measurements.
*International Journal of Bifurcation and Chaos* 14(6), 1905-1933.

Peifer, M. and Timmer, J. (2007). Parameter estimation in ordinary
differential equations for biochemical processes using the method of
multiple shooting. *IET Systems Biology* 1(2), 78-88.

## See also

[`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md)
for multi-start fits,
[profile()](https://dmodverse.github.io/dMod2/reference/profile.md) for
profile likelihoods,
[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md) for
the objective, the guide chapter
[Optimization](https://dmodverse.github.io/dMod2/articles/guide-optimization.html),
the methods article [Trust-region
optimization](https://dmodverse.github.io/dMod2/articles/methods-trust.html).

## Examples

``` r
rosenbrock <- function(x, ...) {
  a <- x[["a"]]
  b <- x[["b"]]
  list(value = 100 * (b - a^2)^2 + (1 - a)^2,
       gradient = c(a = -400 * a * (b - a^2) - 2 * (1 - a),
                    b = 200 * (b - a^2)),
       hessian = matrix(c(1200 * a^2 - 400 * b + 2, -400 * a, -400 * a, 200),
                        2, 2, dimnames = list(c("a", "b"), c("a", "b"))))
}

fit <- trust(rosenbrock, c(a = -1.2, b = 1), iterlim = 200)
fit$argument
#> a b 
#> 1 1 
fit$stopReason
#> [1] "fvalue"

# With an upper bound on a
fit <- trust(rosenbrock, c(a = -1.2, b = 1), parupper = c(a = 0.5),
             iterlim = 200)
fit$atBound
#>     a     b 
#>  TRUE FALSE 

# Gauss-Newton, then BFGS with a separate budget, then Gauss-Newton again
fit <- trust(rosenbrock, c(a = -1.2, b = 1), iterlim = 200,
             phases = list(list(hessianMethod = "gn"),
                           list(hessianMethod = "bfgs", iterlim = 50)),
             switches = 2)
fit$phaseIterations
#>   gn bfgs 
#>   26    0 
fit$evalBySource
#>    gn  bfgs   sr1 exact 
#>    27     0     0     0 
```
