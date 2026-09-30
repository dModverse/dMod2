# symmetryDetection at scale

Target: phosphoproteomic ODE models with about 10^3 states, 10^3 to 10^4 parameters, 10^2
conditions and 10^2 to 10^3 observables. Acceptance: a generated model of that size is detected,
with its rank certified, in under 10 min on 10 cores and under 16 GB.

## Benchmark

`bench/symmetry_tgfb_phospho.R`: the canonical Smad core of M017 plus ERK, JNK/p38 (TRAF6-TAK1),
NF-kB, PI3K/AKT/mTOR and Rho/Par6 with phosphosites, Smad3 linker crosstalk and 30 genes with the
pathway feedbacks (DUSP6, IkBa, A20, Smad7, PMEPA1, FB3/FB4, autocrine TGFB1). 131 states, 284
free parameters, 74 observables (targeted proteomics and phosphoproteomics in absolute units,
mRNA with one scale per gene), 16 conditions: control with and without TGF-b, ten inhibitors,
two knockdowns, two overexpressions. The joint system is about 2400 columns wide.

The goal is this model under `equilibrate = TRUE` with its rank certified and its directions
reconstructed. The family (S6) gives the growth curves on the way there.

## Status

Benchmark at 4.3 min on 10 cores, about 11 GB: rank 266 / 305, 39 directions, all exact
scalings, rank proven. Lie order 10, gap order 1.

| Stage | Done | Open |
|---|---|---|
| S1 | continuation first, direct solve only where it fails | continuation without event gaps; direct fallback leaves eliminated pivots and denominators out of the saturation |
| S2 | saturation on condition-local elimination, series rows with per-row precision; gap order and Lie-order certificate stop at the bound of the exact scalings | the reference reduction for the reconstruction is still wide; pivot rows not yet kept to lift directions |
| S3 | exhausted rows leave the series pivot search; sparse pivot rows, threaded updates | blocked elimination with delayed reduction; incremental echelon |
| S4 | jets in column blocks under a memory budget, nested threads | sparse lanes; Krylov rows at resting points; reverse mode |
| S5 | | reconstruction in parameter space; general directions along a continuation (reported as not reconstructed, rank then not proven) |
| S6 | benchmark model | generator with a size knob; bench scripts per stage with min over repetitions |
| S7 | | supplementary notes for continuation genericity and condition-local elimination |

Acceptance runs still due: M017 with autocrine ligand (a general direction is expected), the
full test suite on the final state.

## Where the cost went before S1 to S4

Measured on M017 with autocrine ligand (5 distinct resting models, 49 states, 140 parameters,
38 conditions), on a small resting-state family and on the benchmark.

| Stage | Growth | Observed |
|---|---|---|
| Direct resting-state solve (elimination, reduction, msolve) | exponential in the coupled core | 15-state core of degree 6: msolve > 600 s |
| Root of the solve in GF(p) | about 0.63 per condition, 0.63^K per joint point | 2 of 12 family solves without a root |
| Saturation before msolve | eliminated pivots and denominators are not saturated | spurious positive-dimensional components |
| Joint stacking | width nShared + K nStates, dense elimination over it | 1772 columns; the series rank over them did not finish in 10 min |
| Dual lanes of the jets | width nz + 1 per condition, dense | dominant after the solves |
| Modular elimination | scalar `mulmod`, series path single-threaded | 6x from sparsity and threads alone |

## Principles

1. No nonlinear solve on the main path.
2. Work per condition; conditions meet only in parameter space.
3. Linear algebra is the cost model: blocked, delayed reduction, parallel.
4. Cost follows structure: sparse lanes, Krylov rows at resting points.

## Stages

### S1 Continuation as the main path

A resting model that differs from the reference by free parameters is continued from the
forward-seeded reference state: Newton steps over GF(p)[[eps]], one sparse linear solve per order.
Resting models carry no switches (forcings are held at 0), so conditions differ by parameters or
by held states only. Since a continued parameter takes a random value anyway, eps = 1 is a generic
point and the rank over GF(p)((eps)) is the generic rank.

- Direct solves remain for small cores and as fallback when the resting Jacobian is singular.
- No root lottery, no boundary components, polynomial cost.
- Done when: every family member with K conditions solves every point without redraw.

### S2 Condition-local elimination

The joint kernel in parameter space is the intersection over conditions of
{dθ : O_k^θ dθ ∈ range O_k^x}. Each condition eliminates its own state columns (pivots on state
columns first); the rows without a state pivot are its constraints on dθ. Only these are
stacked, over nShared columns. The pivot rows are kept to lift a parameter direction back to its
state parts.

- Width drops from nShared + K nStates to nShared; conditions run independently and in parallel.
- The series axis stays inside one condition.
- Done when: M017 joint rank matches the wide stack, and time grows linearly in K.

### S3 Fast modular linear algebra

- Blocked elimination with delayed reduction (64-bit accumulation, one reduction per block).
- Incremental echelon: new jet orders are reduced against the current basis, not re-eliminated.
- Series elimination on the same kernel, per condition.
- Done when: rank of a dense 5000 x 5000 matrix mod p in seconds.

### S4 Jets at scale

- Sparse dual lanes: a lane is carried only where its column has reached the tape.
- At a resting point the rows are C J^k (lemma `lem:stop-rest`): sparse products with the resting
  Jacobian, no Taylor tape.
- Reverse mode where rows are far fewer than columns.
- Done when: jets of a 1000-state model to order 20 in seconds.

### S5 Reconstruction at scale

- Scalings from sparse integer elimination on the monomial exponents.
- General directions only on the residual support; inner detections on the subnetwork of that
  support.
- Done when: a 1000-state model with 100 scalings and a few general directions reconstructs in
  under a minute.

### S6 Model family and benches

`bench/` generator for phosphoproteomic networks: multisite phosphorylation, kinase and
phosphatase cascades, feedback, inhibitor conditions, fold-change observables, with a size knob
from 50 to 5000 states. One bench script per stage, reporting min over repetitions. Regression
coverage of every stage in `tests/testthat`.

### S7 Theory

Supplementary notes, Lean where cheap:

- continuation rank equals the generic rank for continued free parameters;
- condition-local elimination gives the joint kernel.

## Order

| Milestone | Stages | Acceptance |
|---|---|---|
| M1 | S1, S2, S6 | benchmark: every point solved without redraw, rank certified (reached) |
| M2 | S3, S4 | benchmark rank in under 10 min on 10 cores (reached, S3 and S4 in part) |
| M3 | S5, S6, S7 | benchmark and M017 autocrine reconstructed; family at 1000 states |

## Risks

- Held states that differ between conditions (overexpression) start the continuation at 0.
- Sparsely observed modules need high jet orders; the Lie-order certificate bounds them but
  costs its own rank per segment.
- Row memory at 10^3 observables times order 20 per condition; the incremental echelon keeps
  only the basis.
