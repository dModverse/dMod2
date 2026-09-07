# Trust-region optimiser: Stage 1b

Version 0.7.1 at the time of writing, so pre-1.0 and one breaking interface change is
affordable.

## Context

Stage 1 landed in `e1a88bc`: `trust()` takes its model Hessian from an
interchangeable source (`gn`, `bfgs`, `sr1`, `hybrid`), the quasi-Newton updates are Powell-damped
with a Li-Fukushima cautious test, and the Moré-Sorensen subproblem already handles the indefinite
matrices SR1 produces. That is fides' Hessian-source axis, reproduced.

Reading it against the current literature turned up three defects that are verified in the code
rather than suspected, plus three methods that are standard practice and absent. None needs the
adjoint; all raise the quality of the optimiser the adjoint would eventually feed. The ordering
below is cheap before expensive, and the measurement harness comes before both, because each stage
is kept or reverted on a number rather than on an argument.

**The three defects.**

1. **SR1 never learns from a rejected step.** The quasi-Newton update sits inside
   `if (accept && eval_ok)` (`src/trust_kernel.cpp:349`). Nocedal and Wright, sec. 6.2, perform the
   SR1 update after a rejected trial point as well: `(s, y)` there is valid curvature, and that is
   SR1's advantage inside a trust region. The trial gradient is already computed and held in
   `grad_try`; today it is discarded.
2. **The identity seed is never sized.** `qn_gamma` is computed in the dense branch
   (`src/trust_kernel.cpp:373-377`) and never used, because `qn_update()` takes no `gamma`; it is
   live only in the limited-memory branch through `qn_assemble()`. The missing step is
   Oren-Luenberger sizing `B <- (y'y / s'y) I` after the first accepted pair, which matters for
   `hessianInit = "identity"`, the seed the kernel's own comment calls "all a gradient-only
   derivative scheme can supply".
3. **The exact Hessian exists and is not selectable.** `src/residual_kernel.cpp` computes it under
   `use_deriv2_exact`, second-order sigma terms included. A caller reaches it only by smuggling
   `deriv2 = TRUE` through `trust()`'s `...`, which works because `R/trust.R:189` splats `dots` into
   every call, but fixes it for the whole run. It cannot be switched per phase, `.trustControl`
   rejects it in a `control =` list, and in a quasi-Newton phase it is silently a no-op because the
   kernel forces `hessian = FALSE` (`trust_kernel.cpp:307`) and `R/objClass.R:209` then drops
   `deriv2`.

**One unrelated bug found on the way.** `R/trust.R:259-264` builds `trustL1`'s callback as a
one-argument function, but `eval_objfun` (`src/trust_driver.h:181`) always calls
`objfun(x, build_hessian)`. Any `trustL1(..., extra_arg = ...)` therefore raises "unused argument",
which `catch (...)` swallows into `eval_ok = false` and reports as
`"parinit not feasible: objfun failed"`. The no-dots branch survives only because it passes `objfun`
itself, whose `...` absorbs the second positional. No test covers it;
`tests/testthat/test-trustL1.R` is entirely dots-free. Stage 0 fixes it, because stage 6 widens that
callback again.

## Scope

**In.** SR1 on rejected steps; sizing of the quasi-Newton seed; non-monotone acceptance; a switching
rule that can also switch back; structured secant; `exact` as a first-class Hessian source; and the
interface redesign that stops the string list from exploding.

**Out, with reasons.** Recorded so they are not re-proposed later.

- **Hierarchical optimisation** (Loos et al. 2018, analytic elimination of scalings, offsets and
  sigma). Excluded by decision: not always available in closed form.
- **Adaptive trust-region scaling**, a per-iteration `D_k = diag(max(D_{k-1}, sqrt(diag H)))`, and
  the associated initial-radius heuristic. Excluded by decision. `parscale` stays the constant
  change of variables it is today, baked in once at `src/trust_kernel.cpp:155-168`.
- **Scale corrections the log parametrisation already makes.** `parscale` is off by default,
  Oren-Luenberger sizing of a quasi-Newton seed was removed rather than switched off, and adaptive
  trust-region scaling stays out. Two reasons: the units are already uniform on a log scale, and
  each of these estimates a scale from far too little information. Sizing takes one `(s, y)` pair
  from the first step of a random start and stamps its curvature on every direction, which a
  spectrum spread over decades cannot support; with `qnMemory = 0` the updates accumulate onto the
  seed, so the estimate is never forgotten. Measured on Boehm over 100 shared starts, sizing costs
  an identity seed almost every optimum it would otherwise reach.
- **Multi-secant and block quasi-Newton** (Schnabel; Gao and Goldfarb). Enforcing `B S = Y` over a
  window instead of one pair per iteration is cheap at `n <= 150`, but a symmetric `B` cannot
  generally satisfy it, `S'Y` would have to be symmetric, so a correct version is a least-change
  projection and a stage of its own. `qn_assemble()` already rebuilds from a window and that arm
  (`bfgs_id_m10`) has never beaten the dense one. Re-propose only against it as the incumbent.
- **Sampled quasi-Newton** (Berahas, Jahani, Takac). Extra gradients along random directions per
  iteration, `(S, Y)` built from them: gradients bought, iterations saved. Under forward
  sensitivities an extra gradient costs what a Jacobian costs, so it is dominated by `gn`
  outright. It belongs to the reverse-mode work, conditional on the adjoint, not to this plan.
- **Hessian-vector products by forward-over-reverse** (Steihaug-Toint, GLTR). Not the same object
  as the `J^T(Jv)` ruled out above, since the full Hessian is indefinite and a trust region wants
  that. Ruled out anyway by the second argument there: a sloppy spectrum has no effective rank to
  exploit, so the Krylov count lands on the parameter count.
- **DFP, Broyden class, Barzilai-Borwein.** At `n <= 150` with a dense approximation the margin over
  BFGS and SR1 does not repay a second family to maintain and test.
- **Non-linear CG, spectral gradient.** Not competitive against a dense quasi-Newton at this size.
- **Matrix-free or Krylov Gauss-Newton.** Ruled out on structural grounds: at a few hundred
  parameters forming and factorising `J` densely costs microseconds next to one ODE solve, and a
  sloppy spectrum has no effective rank a Krylov method could exploit.
- **Sketched Gauss-Newton.** A small effective rank makes the sketch cheap
  and the step wrong in the directions carrying the decrease. If revisited it belongs to the
  reverse-AD track as a seed or preconditioner, not as a subproblem Hessian.

**`trustL1` is declared frozen.** `src/trustL1_kernel.cpp` carries its own driver
(`ref_init`, `ref_propose`, `ref_accept`, `ref_result`) with no `qn_*` machinery, and its
`trustL1_impl` takes no `hessianMethod`. Adding Hessian sources to `trust()` needs no change there.
Non-monotone acceptance would have to be mirrored into `ref_accept` by hand. The decision is not to,
and to say so in the roxygen. Record it at stage 3 rather than rediscover it later.

## Interface

Done once, up front, so every later stage slots in instead of appending another argument. Grouped
into thematic control lists on the `profile()` pattern, because the flat list had reached 26
arguments.

    trust(objfun, parinit, rinit = 0.1, rmax = 10, iterlim = 100L,
          hessianMethod = c("gn", "bfgs", "sr1", "hybrid"),
          parscale = NULL, parupper = NULL, parlower = NULL,
          tolControl  = NULL,   # ftol, mtol, gtol, xtol, rmin
          qnControl   = NULL,   # hessianInit, qnMemory, qnCautious, qnRejected, qnSizing
          stepControl = NULL,   # boundary, theta.max, nonmonotone
          minimize = TRUE, blather = FALSE, printIter = FALSE, traceFile = NULL,
          ...)

`.trustDefaults(hessianMethod)` holds the per-method defaults and is the single source of truth for
what each group accepts; `.mergeControl()` merges a user list into it and rejects unknown names.
Stages 4 to 6 add `hessianFallback` and `fallbackLimit` as formals beside `hessianMethod`, and
`ssm` and `exact` to its choices.

Only `qnRejected` genuinely varies by method: `TRUE` for `sr1` (Nocedal and Wright, sec. 6.2),
`FALSE` otherwise, where it is inert anyway. The other defaults are the literature values already in
use: `theta.max = 0.99995` (Coleman-Li, as in fides), `qnCautious = 1e-8` (Li-Fukushima),
`qnMemory = 0` (dense is right below a few hundred parameters), `qnSizing = TRUE` for an identity
seed (Oren-Luenberger, N and W sec. 6.1). `nonmonotone` defaults to `0`: Zhang and Hager report
`eta = 0.85` as generally best, but every production trust region (fides, TRON, MINPACK) is monotone
by default, so the value is documented rather than switched on.

**The flat arguments are gone, not deprecated.** `.trustMoved` maps each old name to its new
`group$member` and `.trustRejectMoved()` raises that mapping as an error from both `trust()` and
`mstrust()`. The `mstrust()` guard is the load-bearing one: it forwards on
`intersect(names(formals(trust)), ...)`, so without it a stale `qnMemory = 10L` would be routed
silently to the objective. `fterm` and `mterm` survive on `trustL1`, which keeps its flat
tolerances.

The trigger stays what it is: the first soft stop (`fvalue`, `preddiff`, `step`) or stagnation, at
`src/trust_kernel.cpp:401-421`. The existing single-shot switch there generalises into "switch
source, reset counters and radius, increment switch count".

Touchpoints: `R/trust.R` (formals, roxygen, the control helpers), the `mstrust` guard in
`R/statistics.R`, the shared catalogue `inst/benchmarks/hessianSourceSettings.R`, the seeding chunk
of `dev/optimisation/Optimisation.Rmd`, `profile()`'s `optControl` users such as
`inst/examples/BA_transport.R`, and `tests/testthat/test-trust.R`. The kernel is untouched: the
regrouping is entirely R-level, and `trust_impl` keeps its positional signature.

## Stages

Each stage says what is measured and what decides whether it stays. A stage that does not clear its
gate is reverted, not kept because it is theoretically right.

**Status, 2026-09-07.** Measured on Boehm, 100 shared starts, `set.seed(20260905)`, one objective
and one start set for every arm.

| stage | outcome |
|---|---|
| 0 callback and harness | landed |
| 1 SR1 from rejected steps | landed, gate cleared: `sr1` with an identity seed goes from 18 to 34 starts reaching the optimum, and from 973 to 580 evaluations per one |
| 2 sizing the seed | **removed, gate failed.** Sizing costs an identity seed almost every optimum, 13 to 1 for `bfgs` and 34 to 1 for `sr1`, with most runs stagnating at a common plateau. Removed rather than defaulted off, see "Out, with reasons" |
| 3 non-monotone acceptance | landed, inert at the default `0`. The `eta` gate has not been run |
| 4 selectable fallback and switch back | landed, gate cleared narrowly: `gn -> sr1` raises gradient stops from 13 to 28 without losing a start, but wins no start that `gn` alone does not reach |
| 5 structured secant | deferred to a follow-up branch |
| 6 `exact` | deferred to a follow-up branch |

One result the stage list did not anticipate: the separation on this model is not between Hessian
sources inside a run but between primary methods. `sr1` from an identity seed reaches the optimum
from several times as many starts as `gn`, at a fraction of the cost per optimum, and matches its
typical value.

**Bachmann answers that, and reverses it.** 113 parameters, 100 shared starts on Helix, three arms
of 50 array tasks. Measured against the plateau every converging Gauss-Newton start lands on:

| arm | best | median | reaching the plateau | median evaluations | evaluations per plateau start |
|---|---|---|---|---|---|
| `gn` | -389.79 | -389.13 | 50 of 100 | 128 | 276 |
| `gn -> sr1` | -389.79 | -389.13 | 50 of 100 | 142 | 308 |
| `sr1`, identity | **-393.96** | -387.71 | 30 of 100 | 463 | 1590 |

So the identity-seeded arm is worse on every routine measure at this size, which is the predicted
narrowing gone past narrowing, and at the same time it is the only arm that leaves the plateau at
all: two of its starts end below -390 and the best at -393.96, with none of its 113 parameters at a
bound. The default stays `gn`. `sr1` from an identity seed is a way to find a better optimum at
several times the price, not a replacement.

The stage-4 gate reads the same at both sizes: the handover matches `gn` on quality and costs about
a tenth more, and it moves the run off the model tests onto `stagnation`. It is kept as the
interface that generalises the hard-wired hybrid, not as a better default.

**Why the handover does not inherit what `sr1` finds.** Matched start by start, `gn` and
`gn -> sr1` differ by a median of 3e-6 and are bit-identical on 21 of 100; the quasi-Newton phase is
14 of 142 evaluations. On the two starts where `sr1` from an identity seed reaches -393.4 and
-393.96, both Gauss-Newton arms return exactly -389.789. So `gn` has already chosen the basin by the
time the fallback fires, and the fallback only polishes inside it. What finds the deeper optimum is
`sr1` running the descent, not `sr1` refining a Gauss-Newton solution. That is the experiment the
reversed handover `sr1 -> gn` tests.

**The prior these numbers were measured under was wrong, and the run was repeated.** Centred at a
flat -1 for a model whose published parameters span twelve decades, it cost 108.97 at the published
optimum, against a 4.17 gap between the local optima it was meant to help resolve, and it moved the
optimum roughly 86 units in the data term. Recentring it on the published order of magnitude for
every parameter more than a decade from -1 leaves the centre within one decade everywhere and brings
the same `sigma = 4` down to 1.04 at the published optimum. That version, with 200 starts per arm and
five arms, is the one the shipped numbers come from.

A defect found while measuring, and fixed before any of the gates above were read: the best-iterate
snapshot introduced with stage 3 ran only at the head of an iteration, so a run ending on `fvalue`,
`preddiff`, `step` or the iteration limit reported the second to last accepted point. Every gn
number in the previous vignette was produced by that build.

### Stage 0. Make the harness answer the question, and fix the callback

No optimiser change.

- Fix the `trustL1` one-argument callback bug described above, with a test that passes an extra
  argument through `...`.
- The metric already exists in one of the two harnesses: `evalPerHit = sum(neval) / hits`
  (`dev/optimisation/bachmann_hessianSource.R:151`), fides' `phi` inverted. Add it to the vignette's
  `summary()` (`Optimisation.Rmd:349-355`) so both report the same number, and freeze the success
  threshold, `best + 0.1` at `:343-347`, as the definition.
- Make the variant list data rather than three near-copies: one `settings` definition shared by
  `Optimisation.Rmd`, `bachmann_hessianSource.R` and `inst/benchmarks/bench_hessianSource.R`.
- Add the counters the later gates need: number of source switches, and gradient evaluations per
  source. `qnEval` and `qnSkipped` exist in the result list; per-phase counts and `nSwitch` do not.
- The vignette cache is memoised by name (`Optimisation.Rmd:53-61`, key `"multistart"`), so every
  stage that changes accepted trajectories needs `dev/render-optimisation.R --refit`, not just a
  re-render.

Boehm (9 parameters, 100 shared starts, `set.seed(20260905)`) is the gate for every stage. Bachmann
(113 parameters, SLURM through `distributedComputing`) confirms only what has already cleared Boehm.

### Stage 1. SR1 update after rejected steps

At `src/trust_kernel.cpp:349`, form `(s, y)` from the trial point whenever `eval_ok`, and apply the
SR1 update irrespective of `accept`. BFGS keeps updating on accepted steps only; its positive
definiteness rests on the curvature condition holding along an accepted step. Gated by `qnRejected`.

Gate: on `.quadratic_objfn` (`tests/testthat/test-trust.R:21-29`) SR1 must reproduce the exact
Hessian within `n+1` accepted steps, which it cannot do today. On Boehm, `sr1` must not regress.

### Stage 2. Size the quasi-Newton seed

Use the `qn_gamma` already computed: after the first accepted pair, scale the seed by `y'y / s'y`
before the first update, when `hessianInit = "identity"`. Leave a `gn` seed alone, since it is
already in the right units and rescaling discards what it was chosen for. Gated by `qnSizing`.

Gate: `bfgs` and `sr1` with `hessianInit = "identity"` improve on Boehm; `gn`-seeded runs stay
bit-identical, which is a regression test rather than a measurement.

### Stage 3. Non-monotone acceptance

Zhang-Hager: carry `Q_{k+1} = eta Q_k + 1`, `C_{k+1} = (eta Q_k C_k + f_{k+1}) / Q_{k+1}`, and score
`rho` (`trust_kernel.cpp:326`) against `C_k` instead of `f_used`. `nonmonotone = 0` reproduces
today's behaviour exactly.

Two things this must not break:

- The convergence tests (`ftol`, `mtol`, stagnation, `trust_kernel.cpp:401-410`) keep reading the
  true objective. Acceptance loosens; termination does not.
- The final iterate need not be the best one visited. Track the best
  `(value, argument, gradient, hessian)` and return that from `:437-463`.

Gate: `evalPerHit` on Boehm at `eta` in {0.5, 0.85} against `eta = 0`, for `gn` and for the hybrid.
Stagnation stops should become rarer, which is the mechanism under test. Needs a valley objective;
the suite has no Rosenbrock, so add `.rosenbrock_objfn` beside the three existing local helpers.

### Stage 4. Switching back

`fallbackLimit > 1` lets a stalled fallback phase hand control back to `hessianMethod`, re-seeding
the quasi-Newton approximation from a fresh `H_GN` at the current point, one Jacobian per switch.
This is the existing switch path at `:416-421` plus a counter, so it is small once the interface
work is done. While there, verify with a test what happens when the switch fires on a stagnation
stop: `accept` is false, so the next iteration skips the recompute block at `:234` and re-solves
with the cached `Bhat` at the reset radius. That is probably correct, since `H_full` is deliberately
not reset, but nothing currently asserts it.

Gate: does it raise the fraction of evaluations spent in the quasi-Newton phase without lowering the
success rate? That fraction is what decides whether an adjoint gradient would have a caller at all.

### Stage 5. Structured secant (Dennis-Gay-Welsch)

The Gauss-Newton track's answer to a stalled `gn`, and the one that keeps the exact part exact.
Write `H = H_GN + A` and learn only `A`:

    y_sharp = y - H_GN(x_{k+1}) s          # y = g_{k+1} - g_k
    A_{k+1} = update(A_k, s, y_sharp)      # SR1 or BFGS form, A_0 = 0

This needs nothing new from the objective. The kernel already reads the Gauss-Newton Hessian at
every trial point (`Htry_mat`, `trust_kernel.cpp:312-316`), so `H_GN(x_{k+1}) s` is on hand, and
stating the update through the objective's own `H_GN` keeps it free of any convention about `J` and
`r`, which the kernel never sees separately. `qn_update_B()` (`:54`) is reused unchanged on
`(s, y_sharp)`. The one structural change: `want_h` (`:307`) must stay `true` in this mode, unlike
the quasi-Newton phases.

Because `H_GN` is needed every iteration, this is a Gauss-Newton-track method. It cannot run on
adjoint gradients and is not a candidate consumer for the reverse sweep. Say so in the roxygen so
the two tracks are not conflated later.

Settle two choices against the references before implementing, not by taste: the exact form of
`y_sharp`, where Nocedal and Wright, sec. 10.3, and the original DGW paper differ in which Jacobian
is used, and whether `A` is sized by a scalar before updating, as DGW do.

Gate: `ssm` against `gn` and the hybrid on Boehm, then Bachmann. This is the stage most likely to
win on large-residual problems, and the one whose novelty claim most needs checking against what
fides already ships.

### Stage 6. `exact` as a Hessian source

Last, because it is the most expensive and the narrowest.

`eval_objfun` (`src/trust_driver.h:178-188`) passes exactly two positional arguments. Widen it to
three so the kernel can ask for `deriv2` per evaluation, and update every call site plus both R
wrappers in `R/trust.R`; the `trustL1` one is the bug from stage 0. The subproblem needs no change,
since the exact Hessian is indefinite away from the optimum and the hard-hard branch of
Moré-Sorensen handles that natively. That is why this is cheap here and expensive elsewhere.

**Know what it costs before building it.** `deriv2 = TRUE` requires a model compiled with
`odemodel(..., deriv2 = TRUE)`; `Xs.cppDE` dispatches to a third compiled object and hard-errors
otherwise (`R/prediction.R:301-312`), and `Xs.deSolve` and `Xf` refuse outright. The solve then
carries `p(p+1)/2` extra sensitivity blocks, about 45 for Boehm and about 6400 for Bachmann. This is
a small-model method, and the plan says so rather than discovering it on the cluster.

**It must not degrade silently.** If the `deriv2` attribute is missing,
`normL2_kernel.cpp:164-190` leaves `has_d2pred` false and the kernel falls back to Gauss-Newton
without a word. `hessianMethod = "exact"` has to error at the first evaluation instead.

With the interface above, `hessianMethod = "gn"` with `hessianFallback = "exact"` is the alternative
to the fides hybrid: when Gauss-Newton stalls, `S = sum_i r_i grad^2 r_i` is what it stalled on, so
compute `S` exactly for the remaining iterations rather than approximate it. That is the comparison
fides could not run, because both of their phases cost the same.

Gate: `gn -> exact` against `gn -> bfgs` and against `ssm` on Boehm, reported in `evalPerHit` and in
wall clock. Bachmann only if Boehm is favourable; on that model the per-evaluation price is expected
to be decisive on its own.

## Verification

Toolchain per `CLAUDE.md`. Call the binaries explicitly and run tests serially.

    $env:TESTTHAT_PARALLEL="false"; $env:NOT_CRAN="true"
    & "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" -e 'devtools::test(filter = "trust")'

- **Unit tests** extend `tests/testthat/test-trust.R`, the only file mentioning `hessianMethod`
  anywhere in `tests/`. Follow its conventions: `## ---- Title ----` banners, lowercase declarative
  test names describing behaviour, `info = hm` inside method loops. Stages 1 to 4 are optimiser-only
  and belong on the closed-form helpers (`.quadratic_objfn:21`, `.flat_objfn:492`,
  `.misscaled_objfn:499`, plus a new `.rosenbrock_objfn`), so a failure implicates the optimiser and
  not the integrator. Stages 5 and 6 need a compiled model; reuse `fx_decay_compiled()`
  (`helper-fixtures.R:54`) behind `skip_if_no_compile()`.
- **Regression per stage.** Each stage ships a test that its off setting (`qnRejected = FALSE`,
  `qnSizing = FALSE`, `nonmonotone = 0`, `fallbackLimit = 1`) reproduces the pre-stage trajectory
  bit for bit. That is what makes the gates interpretable: benchmark movement is then caused by the
  flag and not by drift.
- **Interface.** A test that `hessianMethod = "hybrid"` still maps onto `gn` with
  `hessianFallback = "bfgs"`, and that the invariant asserted at `test-trust.R:610-629`, a
  quasi-Newton run ending only on `gradient`, `stagnation`, `radius`, `objfun` or `iterlim`,
  survives the rewrite.
- **Benchmark per stage.** `inst/benchmarks/bench_hessianSource.R` on Boehm, same start set,
  reporting `evalPerHit` plus wall clock. Bachmann through
  `dev/optimisation/bachmann_hessianSource.R` (`BACHMANN_SUBMIT`, `BACHMANN_COLLECT`) only for
  stages that clear Boehm.
- **Documentation, once at the end.** `devtools::document()`, then re-render the optimisation
  vignette with `dev/render-optimisation.R --refit` and commit the PDF alongside the source change.
  Never `devtools::build_vignettes()`, which deletes the only copy.

## Risks

- **The vignette's prose, not just its numbers.** `results.rds` is a cache that reproduces the
  shipped PDF. Stages 1 to 4 change accepted trajectories, so the text around the tables has to be
  re-read after a refit, not only re-rendered.
- **`deriv2` cost is super-linear in `n_theta`.** Stage 6 can look excellent on Boehm and be
  unusable on Bachmann. Its wall-clock gate is not optional, and the roxygen has to state the
  `odemodel(deriv2 = TRUE)` precondition.
- **Stage 5 may already exist upstream.** fides' Hessian-approximation list contains structured
  variants. Check what the paper actually benchmarks before describing anything here as new.
- **Two drivers drifting.** Freezing `trustL1` is a decision, not an accident, but it means
  `trust()` and `trustL1()` stop sharing acceptance semantics from stage 3 onward. It belongs in the
  roxygen of both.
