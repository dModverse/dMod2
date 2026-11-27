# Multiple shooting

Single shooting integrates from the initial state over the whole time
axis. For oscillating or spiking systems the objective is rugged, with
many local minima and narrow valleys: a parameter change shifts every
later period, and a peak moved by more than the peak width leaves the
residuals without direction, so few starts reach the narrow valleys.
[`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md) offers
multiple shooting instead ([trust
region](https://dmodverse.github.io/dMod2/articles/methods-trust.md)).

## Usage

Multiple shooting is a control of the objective:

``` r

obj <- normL2(data, g * x * p, multipleShootingControl = TRUE)
fit <- trust(obj, start)
fit$multipleShooting$violation
controls(obj, "multipleShootingControl") <- NULL   # single shooting again
```

The chain must be `g * x * p` with one
[`Xs()`](https://dmodverse.github.io/dMod2/reference/Xs.md) built on
`odemodel(..., includeTimeZero = FALSE)`. Not available with it:
`hessianMethod = "exact"`, `qnControl$hessianInit = "exact"`, several
phases, `switches`, `hessianReseed`, `qnMemory`, `nonmonotone`,
`boundary = "clip"`, `minimize = FALSE`, `traceFile`.

| entry of `multipleShootingControl` | default | controls |
|----|----|----|
| `nodes` | `"auto"` | `"auto"`, `"transitions"`, or node times |
| `growth`, `misfit` | `10`, `3` | thresholds at which `"auto"` cuts a segment |
| `minPoints` | states + 1 | data points per half of a cut |
| `charts` | all linear | `"log10"`, `"linear"` or `"angle"` per state |
| `scale` | from the data | scale of the gaps and node steps |
| `init` | `"data"` | node values: `"data"`, `"spline"`, `"simulation"` or a list |
| `breaks` | none | node times without continuity |

| entry of `stepControl` | default | controls |
|----|----|----|
| `acceptance` | `"filter"` | `"filter"`, `"merit"` or `"natural"` |
| `soc` | `TRUE` | second-order correction of a rejected step |
| `anneal` | `TRUE` | penalty on the gaps raised in stages |
| `twoPhase`, `restore`, `regularize` | off | further globalizations |

`converged = TRUE` requires the scaled gaps below `tolControl$ctol`. The
result has a list `multipleShooting` with `nodes`, `gaps`, `violation`
and the counts of the globalization.

## Segments and gaps

Multiple shooting [\[1\]](#ref-bock1983recent),
[\[2\]](#ref-baake1992fitting) cuts the time axis at nodes \\\tau_1 \<
\dots \< \tau_M\\ into segments and integrates segment \\j\\ from node
values \\s_j\\, the first from the initial state of the parameter
transformation. The node values are fitted with the parameters. The gap
\\c_j\\ is the difference between the end state of segment \\j\\ and the
node value of the next segment, and the continuity constraint demands
every gap to be zero:

\\c_j = x(\tau\_{j+1};\\ \tau_j, s_j, \theta) - s\_{j+1} = 0 . \tag{1}\\

Observed states start at the data, so every segment stays close to the
data it covers for any parameters. At convergence the trajectory is
continuous and the result is a stationary point of the single-shooting
objective.

## Condensed step

Bock’s generalized Gauss-Newton method, Gauss-Newton for least squares
with equality constraints, linearizes the residuals and the constraints
(1) and takes the step from the linear problem. With the end-state
sensitivities \\G_j = \partial x_j / \partial s_j\\ and \\P_j = \partial
x_j / \partial \theta\\,

\\\Delta s\_{j+1} = \alpha\\ c_j + G_j\\ \Delta s_j + P_j\\ \Delta
\theta , \tag{2}\\

and the recursion gives \\\Delta s = \alpha z + S\\ \Delta\theta\\.
Substituted into the linearized residuals it leaves a Gauss-Newton model
in \\\Delta\theta\\ alone, of single-shooting size, so subproblem solver
and bounds of
[`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md) apply
unchanged. This elimination of the node steps is condensing. The nodes
follow the step by the linearization, not by simulation, so the gaps
close only at convergence.

The trust region bounds parameters and nodes together,

\\\lVert D^{-1} \Delta\theta \rVert^2 + \lVert V (\alpha z + S\\
\Delta\theta) \rVert^2 \le r^2 , \tag{3}\\

with \\D = \operatorname{diag}(\lvert v \rvert^{1/2})\\ for the
distances \\v\\ to the bounds of [Bound
constraints](https://dmodverse.github.io/dMod2/articles/methods-trust.html#bound-constraints),
the identity without bounds, since a region on \\\theta\\ alone leaves
\\S\\ \Delta\theta\\ unbounded, which is huge across a spike. \\V\\ is
one over the state scale (`multipleShootingControl$scale`): data range
for observed states, simulation range otherwise, a decade for log10 node
values. The scaled gap size \\h = \lVert V c \rVert\\ is the Euclidean
norm of the gaps in units of the state scale. \\\alpha = 1\\ closes the
linearized gaps (Bock’s step); under an active region \\\alpha \le 1\\
is the share of the node part of the full step that the radius allows. A
link of the condensing chain that would propagate a gap through
expanding segments beyond a bound is dropped for the step. Other
summands, a prior or a single-shooting
[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md),
enter with value, gradient and Gauss-Newton Hessian.

| control | effect |
|----|----|
| `stepControl$acceptance = "filter"` (default) | filter of Fletcher and Leyffer [\[3\]](#ref-fletcher2002nonlinear) on the objective function \\\ell(\theta)\\, \\-2 \log L\\ or \\-2 \log\\ of the posterior, and the scaled gap size \\h\\: a trial point is accepted if it lowers \\h\\ or \\\ell\\ by a margin against the current point and every pair \\(h, \ell)\\ stored from earlier iterates |
| `stepControl$acceptance = "merit"` | exact penalty \\\ell + \mu h\\ |
| `stepControl$soc` (on) | one second-order correction of a rejected step whose gaps did not close as predicted: the step repeated with the gaps at the trial point closed by the current linearization |
| `tolControl$ctol` | `converged = TRUE` only for \\h\\ below it; stagnation, `rmin`, ten failed solves in a row, `iterlim`, and under `acceptance = "natural"` a zero full step or one failed solve end the run with `converged = FALSE` above it |
| `hessianMethod = "gn"` | condensed Gauss-Newton blocks |
| `hessianMethod = "bfgs"`, `"sr1"` | one quasi-Newton block per segment Lagrangian, from gradients alone, so usable with `sweep = "reverse"` |

With open gaps the nodes absorb what the parameters do not explain, so
the data term does not measure the fit.

## Continuity annealing

From a distant start the condensed step can be forced through a regime
the data never show, such as a spike threshold. `stepControl$anneal`
(on) therefore first minimizes

\\\ell + w \lVert V c \rVert^2 \tag{4}\\

over parameters and nodes on the sparse Gauss-Newton matrix, \\w\\
rising tenfold per stage from one, then hands over to the exact method:
precision annealing [\[4\]](#ref-ye2015improved) with a fixed schedule.
A stage ends when a step gains less than a thousandth of the objective;
annealing ends early where a whole stage no longer closes the gaps. It
needs forward sensitivities.

## Node placement

`multipleShootingControl$nodes = "auto"` (default) starts from ten
segments per condition, or one per half oscillation of the data if more,
and halves a segment where

- the spectral radius of \\G_j\\ exceeds `growth` (default 10), the
  expansion of the linearized flow in the coordinates of the nodes;
- the RMS of the weighted residuals exceeds `misfit` (default 3) and
  twice the condition median;
- or the solve fails.

Cuts happen before the first step and at every iterate, never below half
the finest data spacing, and keep `minPoints` (default: states plus one)
data points per half, so no node fits a segment alone.

Over a period an oscillator has growth near one; what grows is the phase
sensitivity to the frequency, linearly in the periods, until the segment
fits the data a period off. One segment per half oscillation, counted at
turning points of the data, bounds that slip. A misfit on every segment,
from wrong parameters, cuts nothing, since no node removes it.

| setting | node values |
|----|----|
| `charts = c(phi = "angle")` | gap counted modulo \\2\pi\\ |
| `charts` entry `"log10"` | node values of a positive state on the log scale |
| default | observed states at the data; others from a run with observed states reset to the data at each data time, replacement synchronization [\[5\]](#ref-pecora1990synchronization) |
| `init = "spline"` | smoothing spline, roughness by generalized cross-validation [\[6\]](#ref-horbelt2002parameter) |
| `init = "simulation"` | simulation at the start parameters |
| list as `fit$multipleShooting$nodes` | given values, node times in `nodes` |

Under replacement synchronization, driven states such as a voltage-gated
variable forget their start and follow the measured trajectory.

`nodes = "transitions"` puts a fixed node just before every fast change
(a spike, the jump of a relaxation oscillator), a run of samples whose
rate lies more than eight times
[`mad()`](https://rdrr.io/r/stats/mad.html), the median absolute
deviation scaled to a normal standard deviation, off the median rate,
keeps the even nodes of ten segments where no change is near, and does
not cut segments during the fit. Without a fast change in the data it is
`"auto"`. Each change then starts in the segment where the data have it,
and a wrong timing shows as a gap before the next node. Numeric node
times, for all conditions or as a list by condition, fix the layout.

## Further step controls

`stepControl$acceptance = "natural"` selects Bock’s natural level
function [\[7\]](#ref-peifer2007parameter): the full Gauss-Newton step
damped by a factor from a curvature estimate, cut back by a corrector,
and accepted when the simplified step at the trial point (current
Jacobians, trial residuals and gaps) has shrunk. The filter takes over
where the damping finds no step. It needs the Gauss-Newton Hessian, a
fixed sigma and no error model.

Off by default, since each helps on some problems and costs on others:

| control | effect |
|----|----|
| `stepControl$twoPhase` | fit parameters with nodes held at the start [\[6\]](#ref-horbelt2002parameter), then lay out nodes at the result and release them; helps with noisy node values, as on spiking traces, at several times the cost elsewhere |
| `stepControl$regularize` | fix condensed directions with singular value below that fraction of the largest [\[7\]](#ref-peifer2007parameter) |
| `multipleShootingControl$breaks` | nodes without continuity constraint [\[8\]](#ref-voss2004nonlinear); Gauss-Newton Hessian only |
| `stepControl$restore` | filter restoration, steps that only reduce the gaps where the filter accepts no step, at most five times: close the linearized gaps, damped if needed, else put the nodes on the current trajectory; on noisy spike trains this makes the layout cut more segments and fewer starts succeed |

## FitzHugh-Nagumo oscillator

The FitzHugh-Nagumo oscillator [\[9\]](#ref-fitzhugh1961impulses),

\\\dot V = c\\(V - V^3/3 + R), \qquad \dot R = -(V - a + b R)/c ,
\tag{5}\\

is fitted to the voltage \\V\\ alone, 81 points on \\\[0, 20\]\\ with
noise of 0.1.

``` r

fhnDir <- file.path(tempdir(), "fitzhughNagumo")
dir.create(fhnDir, showWarnings = FALSE)
fhn <- eqnvec(V = "c*(V - V^3/3 + R)", R = "-(V - a + b*R)/c")
fhnModel <- odemodel(fhn, modelname = "fhn_ode", includeTimeZero = FALSE, compile = FALSE,
                     outdir = fhnDir)
fx <- Xs(fhnModel, options = list(atol = 1e-10, rtol = 1e-10))
fg <- Y(c(y = "V"), f = fx, modelname = "fhn_obs", compile = FALSE, outdir = fhnDir)
fp <- P(eqnvec(V = "V", R = "R", a = "a", b = "b", c = "exp(lc)"), condition = "C",
        modelname = "fhn_p", compile = FALSE, outdir = fhnDir)
compile(fx, fg, fp, output = "fhn", cores = 6)
fhnPrd <- fg * fx * fp

fhnTruth <- c(V = -1, R = 1, a = 0.2, b = 0.2, lc = log(3))
tt <- seq(0, 20, by = 0.25)
set.seed(1)
yy <- fhnPrd(tt, fhnTruth, deriv = FALSE)[[1]][, "y"]
fhnData <- as.datalist(data.frame(name = "y", time = tt, value = yy + rnorm(length(tt), 0, 0.1),
                                  sigma = 0.1, condition = "C"))
fhnSS <- normL2(fhnData, fhnPrd)
fhnMS <- normL2(fhnData, fhnPrd, multipleShootingControl = TRUE)
```

Forty starts perturb the initial states by a standard deviation of 1.5
and \\a\\, \\b\\, \\\log c\\ by 1.

``` r

set.seed(7)
fhnStarts <- parframe(do.call(rbind, lapply(1:40, function(i)
  fhnTruth + c(rnorm(2, 0, 1.5), rnorm(3, 0, 1)))))
fhnFit <- function(obj, ...)
  as.parframe(mstrust(obj, center = fhnStarts, cores = 4, iterlim = 300, ...))
fhnRuns <- list("single shooting" = fhnFit(fhnSS),
                "multiple shooting, filter" = fhnFit(fhnMS),
                "multiple shooting, merit" = fhnFit(fhnMS, stepControl = list(acceptance = "merit")))
```

A hit has a single-shooting objective within 0.1 of the best, evaluated
anew because open gaps can report a data term below the optimum.

``` r

fhnRuns <- readRDS("results/fitzhughNagumo.rds")
fhnValue <- function(f) apply(as.matrix(f[, names(fhnTruth)]), 1, function(th)
  fhnSS(th, deriv = FALSE)$value)
fhnBest <- min(unlist(lapply(fhnRuns, fhnValue)))
do.call(rbind, lapply(names(fhnRuns), function(method) {
  f <- fhnRuns[[method]]
  hits <- sum(fhnValue(f) <= fhnBest + 0.1)
  data.frame(method, converged = sum(f$converged), hits, medianEval = stats::median(f$neval),
             evalPerHit = if (hits > 0) round(sum(f$neval) / hits) else NA_real_)
}))
#>                      method converged hits medianEval evalPerHit
#> 1           single shooting        19   17        301        424
#> 2 multiple shooting, filter        37   37         38         74
#> 3  multiple shooting, merit        36   36         38         88
```

Multiple shooting reaches the optimum from more than twice as many
starts, at a sixth of the evaluations per hit. Most missed
single-shooting starts hit the iteration limit; the three missed
multiple-shooting starts take almost a third of the multiple-shooting
evaluations, annealing included. An evaluation covers the time axis
once, as in single shooting, plus the node-value sensitivities. Filter
and merit function reach the optimum from nearly as many starts at a
similar cost.

## Stored results

- `vignettes/articles/results/fitzhughNagumo.rds`: the FitzHugh-Nagumo
  fits, from the script `fitzhughNagumo.R` beside it.

Reference:
[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md),
[`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md),
[`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md),
[`controls()`](https://dmodverse.github.io/dMod2/reference/controls.md)

## References

\[1\]

H. G. Bock, “Recent advances in parameter identification techniques for
O.D.E.” in *Numerical treatment of inverse problems in differential and
integral equations*, vol. 2, P. Deuflhard and E. Hairer, Eds., in
Progress in scientific computing, vol. 2., Boston: Birkhäuser, 1983, pp.
95–121. doi:
[10.1007/978-1-4684-7324-7_7](https://doi.org/10.1007/978-1-4684-7324-7_7).

\[2\]

E. Baake, M. Baake, H. G. Bock, and K. M. Briggs, “Fitting ordinary
differential equations to chaotic data,” *Physical Review A*, vol. 45,
no. 8, pp. 5524–5529, 1992, doi:
[10.1103/PhysRevA.45.5524](https://doi.org/10.1103/PhysRevA.45.5524).

\[3\]

R. Fletcher and S. Leyffer, “Nonlinear programming without a penalty
function,” *Mathematical Programming*, vol. 91, no. 2, pp. 239–269,
2002, doi:
[10.1007/s101070100244](https://doi.org/10.1007/s101070100244).

\[4\]

J. Ye, N. Kadakia, P. J. Rozdeba, H. D. I. Abarbanel, and J. C. Quinn,
“Improved variational methods in statistical data assimilation,”
*Nonlinear Processes in Geophysics*, vol. 22, no. 2, pp. 205–213, 2015,
doi: [10.5194/npg-22-205-2015](https://doi.org/10.5194/npg-22-205-2015).

\[5\]

L. M. Pecora and T. L. Carroll, “Synchronization in chaotic systems,”
*Physical Review Letters*, vol. 64, no. 8, pp. 821–824, 1990, doi:
[10.1103/PhysRevLett.64.821](https://doi.org/10.1103/PhysRevLett.64.821).

\[6\]

W. Horbelt, J. Timmer, and H. U. Voss, “Parameter estimation in
nonlinear delayed feedback systems from noisy data,” *Physics Letters
A*, vol. 299, no. 5–6, pp. 513–521, 2002, doi:
[10.1016/S0375-9601(02)00748-X](https://doi.org/10.1016/S0375-9601(02)00748-X).

\[7\]

M. Peifer and J. Timmer, “Parameter estimation in ordinary differential
equations for biochemical processes using the method of multiple
shooting,” *IET Systems Biology*, vol. 1, no. 2, pp. 78–88, 2007, doi:
[10.1049/iet-syb:20060067](https://doi.org/10.1049/iet-syb:20060067).

\[8\]

H. U. Voss, J. Timmer, and J. Kurths, “Nonlinear dynamical system
identification from uncertain and indirect measurements,” *International
Journal of Bifurcation and Chaos*, vol. 14, no. 6, pp. 1905–1933, 2004,
doi:
[10.1142/S0218127404010345](https://doi.org/10.1142/S0218127404010345).

\[9\]

R. FitzHugh, “Impulses and physiological states in theoretical models of
nerve membrane,” *Biophysical Journal*, vol. 1, no. 6, pp. 445–466,
1961, doi:
[10.1016/S0006-3495(61)86902-6](https://doi.org/10.1016/S0006-3495(61)86902-6).
