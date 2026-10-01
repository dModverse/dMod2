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

Phospho benchmark (linear transcription) at 6.4 min on 23 cores (knecht6): rank 266 / 305, 39
directions, all exact scalings, rank proven. Lie order 10, gap order 1.

Ion-channel benchmark (S8), knecht6, 23 cores, `gaugePreference = NULL`:

| `.nDend` | states | parameters | rank / dim | directions | detection | reduction |
|---|---|---|---|---|---|---|
| 1 | 32 | 100 | 126 / 132 | 2 scalings, 4 general | 9.5 min -> 42 s | 1.2 min -> 5 s |
| 4 | 68 | 157 | 216 / 225 | 2 scalings, 7 general | 12.2 min -> 118 s | 28.6 min -> 6 s |
| 12 | 164 | 309 | 456 / 473 | 2 scalings, 15 general | (over 36 min in the block scans alone) -> 7.1 min | 11 s |

At 12 dendrites the 7.1 min split into the first step (stacked rank and scalings, 87 s) and
the second: scalings 20 s, narrow reconstruction of 15 directions 120 s, exact scalings 22 s,
field proof 37 s. The phospho benchmark is unchanged at 6.2 min.

Every rank proven by invariant fields, every direction reduced (the cAMP translation through
the section cAMP*k_cAMP = 1). The stacked rank is taken at Lie order 9 instead of the block
orders 42-56 (4 observables per condition).

Family (`bench/symmetry_family.R`, equilibrate, knecht5, 11 cores): 2 / 4 / 8 modules
(34 / 68 / 136 states) in 26 s / 104 s / 920 s, every direction reduced in under 4 s. The
global-symmetry fields (one analysis per regime) dominate from 8 modules on.

M017 (TGFbModelling, autocrine): small gene pool 3.2 min, rank 120 / 140; full gene pool
8.0 min, rank 214 / 263 (knecht5, 10 cores). All scalings, rank proven, reduced in under a
second. No general direction in the current M017.

| Stage | Done | Open |
|---|---|---|
| S1 | continuation first, direct solve only where it fails; forward seed matched by augmenting paths and state-solved balances, so saturable kinetics keep it | continuation without event gaps; conditions that differ by numbers (knockdowns, overexpression) still solve directly: with Michaelis-Menten cycles (phospho, 339 parameters) msolve runs for hours |
| S2 | saturation on condition-local elimination, series rows with per-row precision; gap order and Lie-order certificate stop at the bound of the exact scalings | the reference reduction for the reconstruction is still wide; pivot rows not yet kept to lift directions |
| S3 | exhausted rows leave the series pivot search; sparse pivot rows, threaded updates; rank profile per order from one jet | blocked elimination with delayed reduction |
| S4 | jets in column blocks on threads under a memory budget; directional jets for verification | sparse lanes; Krylov rows at resting points; reverse mode |
| S5 | without a resting state: stacked rank, narrow reconstruction at the order its support needs, proof by the closed forms as fields | general directions along a continuation (reported as not reconstructed, rank then not proven); field analyses per regime under equilibrate |
| S6 | benchmark models; family generator with a size knob | family at 1000 states |
| S7 | | supplementary notes for continuation genericity, condition-local elimination, and the rank proof by own fields |
| S8 | ion-channel benchmark detected, proven and reduced up to 309 parameters | the narrow reconstruction per direction (relevance probe now by group testing) |
| S9 | gauge by the sparsity of the remaining kernel, cross-block verification by symbol overlap | independent blocks in parallel |
| S10 | split orbits and the product section, radical sign rules, memoised sign certificates, root-free charts first | a time budget per block |

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

### S8 Transient models with exponentials: the ion-channel benchmark

`bench/symmetry_ionchannel.R`: a multicompartment neuron (soma, AIS, a dendritic cable of
`.nDend` compartments) with Boltzmann gates, GHK calcium flux, calcium pools with buffers, SK
channels, and two neuromodulators (cAMP shifts the HCN activation, PIP2 gates KM). States
enter through exp(), so `equilibrate` does not apply: initial values are free, all conditions
start from the same resting state. Somatic patch, voltage-sensitive dye and a calcium dye with
one loading scale per region; blockers, forskolin, carbachol and two current steps make the
conditions. About 300 parameters at `.nDend = 16`. General directions by construction: the
cAMP translation (cAMP + c against Vh_HCN and k_cAMP) and the PIP2 turnover; the scalings are
gauged by `gaugePreference = NULL`.

- Exponential recast in polynomial time: the independence test of the exponents over Q and the
  lift of the start by cached derivatives, not by sympy ranks and fresh parses per point.
- Sparse observation of a large transient model needs Lie orders near the rank (four
  observables and rank 104 of one condition: order 55). Jet cost per order and the number of
  orders are the cost model, not the solves.
- Done when: the benchmark detects, reconstructs its general directions and reduces, and the
  reduced model is re-detected identifiable, within the budget of the phospho benchmark.

### S9 Reduction at scale

`symmetryReduction()` on a large model with several general directions:

- independent curved blocks solved in parallel;
- cross-block verification only against the blocks whose generators move a symbol of the
  invariant, not against every block;
- sympy work in the smallest symbol set of a block (compressions and module reduction first,
  as today), with a time bound per block that leaves a block `invariantOnly` instead of
  stalling the whole reduction;
- the scaling stage on the sparse integer weight matrix of all scalings at once.
- Done when: the ion-channel benchmark with its general directions reduces in under a minute,
  and a model with ten independent general directions (bench family) in time linear in their
  number.

### S10 Chart search

The chart of a general block decides how fast and whether `symmetryReduction()` succeeds.
Measured on the ion-channel benchmark: the four-dendrite PIP2 block took 28 min under a
gauge that couples every pool through the scaling and 1 min under one that does not; the
cAMP translation spends 60 s of sympy (cancel, fraction, solve) before it is reported
invariant-only.

- Gauge by the sparsity of the remaining kernel: the general directions keep their minimal
  supports (step 1 of the gauged detection, `summary()` suggests it, a preferred gauge that
  couples warns).
- Split orbits: where every face a block reaches is positive only under a sign condition on
  one invariant, and the conditions of two faces are complementary, the orbits leave the
  orthant through different faces. No face and no rational pin covers them; the block is
  reported with the split and the face of each sign, without the balance search.
- Memoised sign certificates and positivity forms per expression, candidates ordered by
  cost, a time budget per block that leaves it invariant-only instead of stalling.
- Independent blocks in parallel.
- Done when: every block of the ion-channel benchmarks is reduced or reported in seconds,
  and every motif of `bench/motifs` stays reduced.

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
| M4 | S8, S9, S10 | ion-channel benchmark detected, reconstructed and reduced; reduction linear in the blocks |

## Risks

- Held states that differ between conditions (overexpression) start the continuation at 0.
- Sparsely observed modules need high jet orders; the Lie-order certificate bounds them but
  costs its own rank per segment.
- Row memory at 10^3 observables times order 20 per condition; the incremental echelon keeps
  only the basis.
