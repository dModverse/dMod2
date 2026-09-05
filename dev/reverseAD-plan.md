# Reverse-mode AD for dMod2 and cppDE

State of the question as of 2026-09-04. Written to be re-argued from scratch, so it separates what
is structure, what is external evidence, and what is still unknown. Nothing here rests on an
in-house measurement; the numbers that would decide it are listed at the end and none of them exist
yet.

## The idea

Forward sensitivities propagate tangent directions through the ODE, one per parameter direction
seeded at `t0`. A reverse sweep seeded from the scalar objective costs one backward pass regardless
of that count. For many-parameter models that should be the cheaper way to a gradient.

## Why the first-order adjoint has no caller today

This is the structural fact the whole question turns on, and it needs no measurement.

In dMod2 the gradient and the Gauss-Newton Hessian come out of the same object:

    g = J^T r        H = J^T J        J = dr/dtheta

So the Hessian is free **given** the gradient, and the gradient is expensive **because** it comes
from J. A pure Gauss-Newton trust region therefore gains nothing from an adjoint: it needs J in
every iteration anyway. Any adjoint scheme first needs an optimiser step that consumes the gradient
without consuming J.

BFGS and SR1 are exactly that. They are the missing caller.

## What decides it

Let

    factor = t(obj with derivatives) / t(obj value only)
    F      = evaluations the gradient-consuming scheme needs / evaluations trust() needs
    c      = cost of one adjoint gradient, in value-solves

then

    gain = factor / (c * F)

`c` is not 2. One adjoint gradient is a forward solve plus a backward solve, the adjoint system of a
stiff problem is itself stiff and needs implicit steps, and checkpointing adds restore or
re-integration. Budget 2 to 4 and treat every tenth as coming straight off the gain.

**There is no ceiling of the form `factor / c`.** Such a bound would need `F >= 1`, i.e. that no
scheme reaches the optimum in fewer evaluations than today's Gauss-Newton trust region. Froehlich
and Sorger report a scheme that does, against pure Gauss-Newton as the baseline, on the PEtab
collection, so `F < 1` is reachable and the expression is open at the top.

## External evidence

**Froehlich and Sorger 2022 (Fides), PLoS Comput Biol.** A trust region with interchangeable Hessian
approximations, benchmarked on the PEtab collection, already cited in `R/trust.R`.

- Metric is `phi = gamma / n_grad`, successes per gradient evaluation. Not wall clock.
- Pure Gauss-Newton: strong global convergence, slow local convergence.
- Pure BFGS and SR1: fast local, poor global. BFGS fails on 1 problem, SR1 on 2.
- Winner is a hybrid: Gauss-Newton until the trust radius fails to update for `n_hybrid` consecutive
  iterations, then BFGS. `n_hybrid = 50` best of 25 / 50 / 75 / 100. Average 1.51-fold improvement
  over pure Gauss-Newton, range 0.56 to 6.34, best method on 7 of 13 problems.
- **Both phases used the same derivative computation.** The reported gain is purely a count effect,
  no phase was cheaper per evaluation.
- Their own caveat: differences between implementations of the same algorithm can exceed differences
  between algorithms, pseudo-inverse against Cholesky being their example.

The last two points are why this matters here. Fides left the price effect on the table because in
their setup both phases cost the same. With a cheap gradient in the BFGS phase, dMod2 could take a
second, multiplicative effect that does not exist in their numbers.

**Stapor, Froehlich, Hasenauer 2018, Bioinformatics.** Second-order adjoint sensitivity analysis.
Established methods scale linearly in the number of states and **quadratically** in the number of
parameters; the second-order adjoint scales linearly in both. That settles what `deriv2` costs, and
it is the strongest single case for reverse mode.

**Sloppiness and Krylov methods.** Sloppy models have eigenvalues spread roughly uniformly in log
over many decades (Gutenkunst et al. 2007). That is the absence of clustering, and clustering is
what Krylov convergence needs. Krylov resolves dominant eigendirections first, while a trust-region
step draws its decrease from the soft tail. This is the worst available pairing.

## Matrix-free Gauss-Newton: ruled out by structure

`H v = J^T (J v)` needs one tangent and one adjoint sweep, so a Steihaug-Toint trust region would
need no explicit Jacobian. Break-even is `2k < n_theta` for `k` Krylov iterations.

The case fails for two independent reasons, neither of which is a measurement.

1. **Matrix-free solves a memory problem dMod2 does not have.** It is standard where the parameter
   dimension makes forming `J` impossible, `10^6` and upward in PDE-constrained inverse problems.
   The condition under which it wins there is stated as: the number of Hessian-vector products must
   be of the order of the *effective rank*, not the parameter dimension. At `n_theta` in the
   hundreds, forming `J` and factorising densely costs microseconds next to one ODE solve.
2. **A sloppy spectrum has no effective rank worth exploiting.** Uniform-in-log decay over the
   observed range puts the effective rank at a substantial fraction of `n_theta`, which lands
   break-even and the iteration count on top of each other. There is no regime where it wins by a
   margin.

If it is ever revisited, two things have to be right. The operator has to be `J^T(J v)`, positive
semi-definite by construction, not the full Hessian, which carries indefinite error-model terms. And
it has to be preconditioned, because unpreconditioned Krylov at this conditioning measures the
conditioning and not the method. LSQR or LSMR is the right vehicle rather than CG on the normal
equations: the rate is the same but the stability is not.

## The plan: consumer before producer

The dependency runs from the consumer to the producer, not the other way. Every stage below is a
prerequisite for the next, and each one can end the programme cheaply.

### Stage 1: interchangeable Hessian source in `trust()`

No derivative work at all. This is the stage that decides whether anything after it has a caller.

`trust_subproblem.h` takes the Hessian as an eigendecomposition and minimises `g^T p + p^T H p / 2`
subject to `||p|| <= r`. It does not care where the matrix came from, and it handles indefinite `H`
natively through the hard-hard branch. That matters: SR1 produces indefinite matrices, and
Moré-Sorensen solves them exactly where a truncated-CG scheme could only stop. The exact subproblem
solver is an advantage here, not something to be preserved reluctantly.

**The trigger already exists.** `trust_driver.h` documents the stagnation state, and
`trust_kernel.cpp:230` and `:261` maintain and act on it:

    if (accept || !eval_ok || dval >= ftol) n_stall = 0;
    if (n_stall >= kStallLimit) { converged = true; stop_reason = "stagnation"; break; }

That is in substance the Fides switching condition. Today it ends the run; Fides switches the
Hessian source and continues. `accept` is decided a few lines above, and `Blather` already records
`accept`, `r` and `rho` per iteration. The driver comment already reads "Mirrors fides'
make_non_degenerate", so this is continuation rather than graft.

The seam is two lines, where the Hessian is read out of the objective's return:

    trust_kernel.cpp:94    Hmat0    = as<NumericMatrix>(out_init["hessian"]);
    trust_kernel.cpp:200   Htry_mat = as<NumericMatrix>(out_try["hessian"]);

A small Hessian source slots in there with `gn` as pass-through, `bfgs` and `sr1` maintaining their
own update, and `hybrid` switching on the existing stall counter. `Hmat0` doubles as the seed, so
the quasi-Newton approximation starts from an exact `J^T J` rather than from the identity, which on
this conditioning is the difference between a usable and a useless start. After a switch the counter
resets and the second stagnation ends the run.

`trust_subproblem.h` is not touched. The numerical core stays Rcpp-free and testable, which is the
whole point of that split.

This is C++ in `src/`, not R.

**Why the update cannot live in the objective.** It was worth checking, because it would have kept
the change out of the kernel entirely. It does not work, for four reasons, the first of which is
decisive:

1. The objective is called at every trial point, accepted or not. A quasi-Newton pair `(s, y)` must
   come from accepted iterates, and `accept` is decided after the call returns.
2. Objectives compose. The quasi-Newton approximation of a sum is not the sum of the approximations,
   so the state cannot sit inside `normL2`, which is below the composition.
3. Multistart and `cores` would carry mutable state across starts and across threads, breaking the
   property that an objective is a function of its parameters.
4. The hybrid trigger needs the trust radius, which exists only in the kernel.

### Stage 2: `buildHessian` on the objective call

A call-time argument next to `deriv` in `objClass.R`, not a constructor argument, because the phase
changes inside a run. `trust()` passes it through the objfun callback, and the kernel tolerates a
`NULL` hessian in the quasi-Newton phase.

Today it saves the `J^T J` contraction only. J is still needed for the gradient, so the ODE work is
unchanged and the saving is dense BLAS, real but small. Its value is that it is the same seam the
adjoint uses later, when `buildHessian = FALSE` comes to mean "do not form J at all, take the
reverse sweep". Pulling the seam now means stage 3 changes an implementation, not an interface.

Open: whether this is a third flag beside `deriv` and `deriv2`, or whether `deriv` becomes ordinal.
The constructors already reject invalid combinations of the existing two, and a third multiplies
that matrix. Pre-1.0 is the moment to decide; after the CRAN release it is a breaking change.

### Stage 3: the reverse sweep, if stage 1 gives it a caller

Only worth starting if stage 1 shows that a useful fraction of evaluations lands in the
quasi-Newton phase. That fraction is the entire target, since it is the only phase the adjoint makes
cheaper.

**It must run at objective level**, through `normL2(dataL, g * x * p)`, or the seed never becomes a
scalar:

    normL2   emits  w = (w_pred, w_sigma)   per condition, per observation row
      g      w_y = (dg/dy)^T w_pred      w_phi += (dg/dphi)^T w_pred
      x      w_phi += adjoint solve seeded with w_y at the observation times    <- cppDE
      p      w_theta = (dphi/dtheta)^T w_phi

The seed has two parts because `normL2`'s gradient runs through the prediction (`dwr`) and through
the error model (`dsigma`); see `src/residual_kernel.cpp`.

**Split across the repos.** Branch `devel-reverseAD` in both, cppDE first.

cppDE: seeded adjoint solve for `x`; `derivMode = c("forward", "reverse")` plus `seed`, with
`sens1ini` and `sens2ini` folded into `seed` in the same change, because the generated `.Call` seam
is positional over a fixed 15-argument signature and should be touched once; codegen emits the
analytic vjp routines, `J^T lambda` for the state part and `(df/dp)^T lambda` for the parameter
part, straight from SymPy, no tape; checkpointing for the replay.

dMod2: analytic vjps for `g` and `p`; the `normL2` seed; vjp composition in `*`; and
`Remotes: dModverse/cppDE@devel-reverseAD` on the dMod2 branch until cppDE's lands, else CI pulls
cppDE master and fails.

The two directions are symmetric, so one argument pair fits `solveODE` and `cppFUN`:

| | seed | enters at | yields |
|---|---|---|---|
| forward | k columns of `dphi/dtheta` | t0, the input end | k **columns** of the Jacobian |
| reverse | k rows of `dL/dy` | the observation times, the output end | k **rows** of the Jacobian |

`seed` is optional for forward, where its absence means the identity over all directions, and
mandatory for reverse, where the identity over the outputs is a large loss on any model with more
observation rows than parameters.

**Main design risk: the control flow inverts.** Today `prd(times, pars, deriv = TRUE)` computes
everything and `normL2` contracts it. In reverse, `normL2` needs the value first, forms the seed,
then asks for `vjp(w)`. That is a second evaluation path through `g * x * p`, not a flag on the
existing one, and `*` has to compose vjps in the opposite associativity.

**Main implementation risk, and a fork the old plan did not mark.** The stepper is a variable-order
Nordsieck multistepper, so a step's dependence spans a history window plus order and step-size
changes. There are two different things one can build here and they are not interchangeable:

- a **discrete** adjoint, which differentiates the actual steps and therefore agrees with the forward
  gradient to machine precision, and whose derivation for this stepper is the hardest single piece
  of work in either repo;
- a **continuous** adjoint with checkpointing and dense output, CVODES style, which is a well-trodden
  path but agrees only to solver tolerance.

The choice is not cosmetic in the scheme this plan proposes. BFGS and SR1 build curvature from
gradient differences `y = g_{k+1} - g_k`, and differences are more noise-sensitive than the
quantities themselves. The switch happens by construction in the stagnation phase, where steps are
small and those differences are worst conditioned. A gradient that is only tolerance-accurate is
therefore least trustworthy exactly where the scheme is supposed to earn its keep, and the curvature
condition `s^T y > 0` will start skipping updates, degrading BFGS toward gradient descent.

Derive and verify for fixed order and fixed step size first, then order changes, then step-size
changes. cppDE also has one-step methods with dense output; a discrete adjoint for Rosenbrock4 is
tractable and is the right place to establish the verification harness before touching the
multistepper.

**Verification before any timing.** The forward path is the oracle. The reverse gradient must equal
the `deriv = TRUE` gradient on a small model, and later the second-order path must equal `deriv2`.
Belongs in `tests/testthat`, not in `bench/`. Note that "to machine precision" is only an achievable
bar for the discrete adjoint; for the continuous one the bar is the solver tolerance, and the test
has to say which is being asserted.

## What would decide it

Fix these before running anything, so the result is not interpreted afterwards.

1. **Does the hybrid reduce evaluations at equal success rate, on dMod2's implementation?** Run
   multistart over the same start set for `gn`, `bfgs`, `hybrid`, and optionally `sr1`, with one
   trust-region code path, identical stopping rule, identical solver tolerances, and a success
   defined in advance against the known best value. Count gradient evaluations, not seconds.
   If the hybrid does not beat `gn`, stages 2 and 3 have no caller and the programme ends here.
2. **What fraction of evaluations lands in the quasi-Newton phase?** One extra column in the trace.
   This is the sizing number for the adjoint: a small fraction makes the adjoint pointless even if
   the hybrid wins.
3. **What is `factor`, honestly?** Across several model sizes, `min` over repetitions, in a script
   that exists. Nothing in the repo currently measures it.
4. **Does `factor` bend upward with model size?** If it grows linearly, `F` plausibly grows with it
   and the ratio cancels at every size. Only superlinear growth from the cache cliff flips that.

## Open, and cheap, and independent of all of the above

**Seed column pruning.** `R/prediction.R:257` builds `sens1ini` as `[n_phi, n_theta]` and zero-fills
it, so every condition integrates every tangent direction, including those whose column of
`dphi/dtheta` is entirely zero. Such a column has zero forcing and zero initial value, so its tangent
is identically zero for all time. In a multi-condition PEtab problem many parameters are scalings,
offsets and noise terms that enter through `g` rather than through `x`, and those columns are
structurally zero in every condition.

This is a hypothesis, not a measurement. The test needs no fit: import a problem and count zero
columns per condition. If it holds, it lowers `factor` directly, with no adjoint and no correctness
risk, through the same `seed` API that stage 3 wants anyway.

It also cuts the other way and that ordering matters: lowering `factor` lowers the adjoint's payoff,
because the adjoint would then be winning against an already cheaper forward path. **Measure pruning
before deciding on the adjoint, not after**, or the most expensive work in either repo gets
justified with a `factor` that the cheap change then halves.

**`parscale` is `NULL` by default.** The trust region ball assumes comparable coordinates. Two
different `D` exist in `trust()`: the Coleman-Li `diag(|v|^{1/2})`, which is bound distance and runs
always, and `parscale`, which is conditioning and is off. The mechanism is built and unused.
