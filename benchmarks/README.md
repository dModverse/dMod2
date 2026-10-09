# dMod2 benchmark suite

The objective of PEtab benchmark problems, imported with `importPEtab()` and
evaluated at their published best fit: the value, the gradient by forward
sensitivities, the gradient by the reverse sweep with and without
`optionsReverse = list(refine = TRUE)`, and optionally by CVODES adjoint
sensitivities on the Sundials backend. Where cppDE's own suite
(`cppDE/benchmarks`) times the solver, this one times the whole chain
`normL2(data, g * x * p, e)` that a fit evaluates, R included.

## Quick start

```sh
source ~/dModverse/env.sh                         # OMP_NUM_THREADS=1, libSBML
Rscript benchmarks/fetch-models.R                 # clone the collection once
Rscript benchmarks/run-benchmarks.R --tier tiny   # a few minutes
```

`--petab-root` points at an existing clone of the
[Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab)
collection instead; both its layouts are read. `--help` lists the options.

## Tiers

| tier | problems |
|---|---|
| tiny | Boehm, Raia, Lucarelli |
| medium | + Bachmann, Isensee |
| full | + Chen, Lang |

A tier is an explicit list, so runs compare across commits. Problems are
chosen by what they exercise, not by size, and a run prints how many of its
problems have each trait, with a warning for a trait none covers:

| trait | what it exercises |
|---|---|
| conditions | many experimental conditions, the batched solves |
| events | events and switching times |
| preeq | preequilibration through a `Pimpl` steady state |
| parameters | many estimated parameters, where the adjoint pays off |
| klu | the sparse linear solver |
| errormodel | estimated error model parameters |

## What is measured

Per problem, on one core (`options(dMod.cores = 1)`):

- `sec`, the minimum over `--nrep` bursts of at least `--target` seconds;
- `cost`, the time of an arm in value evaluations;
- `err_rel`, `max |g - g_ref| / max |g_ref|` against the forward gradient at
  atol 1e-12, rtol 1e-10, and `err_tau`, the same difference against
  `rtol |g_ref| + gradtol` with `gradtol = rtol 1e-3 max |g_ref|`, which is at
  most 1 within tolerance; `refine` runs under that `gradtol`. At the best
  fit the gradient is close to zero, so both measure against a small
  reference;
- `solve_frac`, the share of a call spent inside cppDE, the rest being dMod2's
  R code.

Results land in `benchmarks/results/<tag>/`: `results.csv` with one row per
problem and arm, `README.md` with the traits and the summary table,
`run-info.txt`, and the figures `01-cost.png` and `02-error.png`. `cache/` and
`results/` are not tracked.

`compare-results.R base.csv new.csv` puts two runs side by side: the
geometric mean of new / base per arm, then each problem.

The minimum is reported because the machines scatter. Run on an idle machine
with `OMP_NUM_THREADS=1`; a fit beside it makes every number meaningless.

## Scripts

`scripts/` holds the linear scripts stepped through in RStudio from the
package root: `bench_gradientCost.R` (the cost curve over the number of
estimated parameters), `bench_adjointPerCondition.R` (the discrete adjoint
against CVODES ASA per condition) and `bench_hessianSource.R` with its shared
`hessianSourceSettings.R` (the Hessian sources of `trust()`).
