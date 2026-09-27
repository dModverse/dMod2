# symmetryDetection / symmetryReduction: Lie structure, speed and reach

Two benchmarks drive this plan:

- **M011** (TGFbModelling `M011_H1975TGFb_small`): 38 conditions, events at -2880 / -30 / 0,
  138 coordinates, rank 120, 17 scalings, 1 general direction. Detection-bound.
- **`bench/motifs/`**: small motifs, one difficulty each. Reduction-bound. `run_all.sh` runs
  them in parallel; each run re-detects the reduced model and must find it identifiable.

## Theory the plan builds on

| Fact | Consequence |
|---|---|
| Every jet L_f^k h is invariant under every generator: X(L_f^k h) = 0 (the rows of the observability matrix) | Generators and invariants can be read off the jets instead of searched for |
| A residual direction in the free-column gauge is the only kernel vector supported on its support S with entry 1 at the free column; rows restricted to S have rank \|S\| - 1 | Samples need duals on S only (width \|S\|+1, not nz+1): narrow kernel |
| Same fact, symbolically: X_S ⟂ ∇_S Φ_j for m-1 jets with independent S-gradients | X_i = (-1)^i det(∂Φ/∂z_S without column i); for \|S\| = 2, X ∝ (∂_b Φ, -∂_a Φ), the Hamiltonian field of a jet |
| A generator is defined up to a factor g(z) (a line field); the free-column gauge divides by N_f | Wide entries are often the normaliser's variables, not the direction's |
| Invariants of a block = identifiable combinations; a section Σ meeting every orbit once is a global chart | Charts need a section, not only invariants |
| A face {Z = 0} of the positive orthant is such a section when every orbit reaches it with the rest positive | Face section: solve I(z') = I(z), z'_Z = 0, certify z'(z) > 0 on the orthant; the reduced model switches rates off |
| Y with [X, Y] ∈ span(X) on the support of X gives X(det(X, Y_1..Y_{m-1})) = div X · det(...) | 1/det is a Jacobi multiplier: integrating factors from symmetries of X, no ansatz |
| An entry invariant under the translation ∂_a - ∂_b (every entry, jointly) depends on a, b only through a + b | Translation groups: fit in one variable, fold the sum back; in the reduction, search invariants in the sum |
| Noether needs a symplectic structure; the parameter space has none | Noether itself adds no invariant search; its adjoint form (λ·X conserved) is a cheap validator only |

## Motifs (bench/motifs)

| Motif | Difficulty | Baseline (0.8.1) | Now |
|---|---|---|---|
| cat3 / cat4 / cat5 | catenary compartments, leaks everywhere: N-1 curved directions, invariants of degree up to N | reduction fails (300 s) / timeout / timeout | face section, 2.7 s / 7 s / 26 s, reduced model identifiable |
| mam3 / mam4 | mammillary compartments | fails (245 s) / timeout | face section, 2.8 s / 6.6 s |
| tworoute | two unobserved routes, invariants k3·ka, k4·kb, ka + kb | invariants found, not reduced (greedy carrier bug) | section ka = kb, 2.2 s |
| autocrine_free | autocrine loop, free initial values: scalings entangled with a curved direction | reduction timeout | scalings pinned outside the curved support, face kdg = 0, 1 s |
| route3 / route4 | n routes: invariants k_i·r_i and Σ r_i, no face | not reduced / not reduced | equal-share section r_i = q/n, 3.2 s / 4.1 s |
| autocrine_ss | autocrine loop from steadyStates(), dose event | fine | face kdg = 0 instead of a balance |
| wide4 / wide12 / wide30 | M011's pattern: degradation kdg·P/(Km + R1 + … + Rn), entry ksec·(Km + ΣR)² | reduction 185 s / 1 dir open / 2 dirs open | 34 s / detect 8 s + reduce 6 s / detect 9 s + reduce 47 s, all identifiable after reduction |
| receptor_kd | knockdown switch long before the stimulus; species named Ci | steadyStates() crashed (Ci) | reduced, 1.5 s |
| wideprod8 / wideprod30 | degradation kdg·P·Π/(Km + Π), Π = R1⋯Rn | open directions | detect 19 s / 19 s, reduce 10 s / 18 s, identifiable after reduction |
| TGFbModelling M001–M010 (battery, `work/tgfb/battery`) | the real models, 38 conditions, gauges fixed by `fixed` | crash (named `g`) | all analysed, every reduced model identifiable |
| TGFbModelling M011 full gene pool | 283 coordinates, 49 scalings + 1 general direction | (not run) | general direction on 2 of 234 columns with the translation group Km_Smad7 + four receptor pools, 109 s; 50 of 50 reduced in 6 s |
| switch_late, hill, enzyme, cat2 | events after t0, Hill production, identifiable enzyme | fine | fine (faster) |

## Done

| Step | Effect |
|---|---|
| Lie order per condition (`NtChain`), capped at the block's own saturation | M011 kernel call in a worker 43.5 s → 21.0 s |
| Joint saturation: lower the deepest blocks while the stacked rank stays | M011 all 38 conditions at order 13 instead of 36; worker call 7.0 s |
| `symObsNullChainPointBatch`: (point, chain) jets in one OpenMP batch on the plain gap path | no fork |
| Sparse Laurent and general rational entries sample a whole stage per batch | were serial `kcall` loops |
| `gaugePreference` (default `FALSE`) | implemented |
| **Carrier matching** (`.symRedMatchCarriers`): admissible carriers per invariant, matched by backtracking | tworoute reduced; the greedy pick gave up with an empty reason |
| **Face sections** (`.symRedFaceSection`), tried before the balance/pin search | compartment motifs reduced in seconds; `info$modelExprs` keeps denominators out of the zero set |
| **Narrow kernel** (`narrowOne` in the observability engine): relevance probe, sampling and interpolation on dual columns S only; verified on the full kernel | plain and plain-gap paths; `DMOD_SYM_NONARROW` switches it off |
| **Outside gauges**: a scaling overlapping a curved direction is pinned outside every curved support | autocrine_free: reduction timeout → 1 s |
| **Equal-share sections** for r > 2 gauge coordinates | route4 reduced in 4 s |
| **Translation groups** in detection (`.symTranslationGroups`, `.symShiftGroups`) and reduction (`.symRedTranslationCompress`) | wide30: 2 open directions → closed in 2 s and reduced |
| **Multiplicative groups** (members held at 1 during the fit) and **unmoved monomials** in the reduction (`.symRedUnmovedCompress`), both before the module reduction | wideprod30 (entry ksec·(Km + R1⋯R30)²): open direction → detect 19 s, reduce 18 s |
| **Jet closed form** (`jetGenerator`, Python): X_S = generalised cross product of the S-gradients of \|S\|-1 jets; tried first for \|S\| ≤ 6, time-bounded | wideprod8/30: the \|S\| = 2 direction from the order-1 jets in about 1 s |
| **Named per-condition `g`** unnamed on input | TGFbModelling M001–M009 crashed in the log chart; all now analysed and reduced |
| **Sympy hygiene**: no factor(cancel()) of wide expressions in `.symTidy`; wide directions classified through factor lists (`_classify_factored`) | wide12 finalisation 520 s → 0.1 s; wide30 classification 12+ min → 0.6 s |
| **steadyStates()**: names sympy resolves (Ci, S, E, Q, gamma) aliased | receptor motif with species Ci solved |

Switches: `DMOD_SYM_NOCAP`, `DMOD_SYM_NOJOINTCAP`, `DMOD_SYM_NOPOINTBATCH`, `DMOD_SYM_NONARROW`, `DMOD_SYM_NOGROUP`.

## Work packages

### WP1 reduction: charts and speed
1. ✅ Carrier matching.
2. ✅ Face sections; zero sets pre-filtered (an invariant must stay finite and nonzero on the face).
   ✅ Outside gauges for overlapping scalings; ✅ equal-share sections; ✅ translation groups.
3. Time the chart search: autocrine and every block that reaches the balance search; memoise
   the sympy sign tests, cap candidates by cost, never spend minutes on a block that a face or
   a single balance closes.
4. Zero compatibility through faces: report the faces the orbit reaches (a face section is
   such a certificate) instead of "unknown".

### WP2 detection: jets instead of full kernels
1. ✅ Narrow kernel for the minimal-support / free-column residual directions; ✅ translation
   groups on it.
2. ✅ Closed form from jets (`jetGenerator`): m-1 lowest-order jets with independent
   S-gradients, their minors, verified exactly. Tried first, time-bounded.
3. Narrow kernel on the equilibrate/joint and recast paths (the steady-state seed has full
   width; restrict its columns).
4. ✅ M011: `gaugePreference` end to end, small and full gene pool: the general direction
   closes and reduces (k_dg_TGFB1 = 0), the reduced model re-detected identifiable.

### WP3 invariants: Lie structure in the cascade
1. Jets as invariant candidates: first-segment output Taylor coefficients (identifiable by
   construction) seed the invariant search before the degree ladder.
2. Lie multiplier stage: scaling symmetries of X itself (integer kernel over the exponents of
   its components) and generators with [X, Y] ∈ span(X) give 1/det(X, Y) as the integrating
   factor of Stage 5, and invariants Y(I) from a found I.

### WP4 evidence and hygiene
1. Motif table above kept current; every reduced motif re-detected identifiable.
2. Tests for carrier matching, face sections, narrow kernel, jet generator.
3. NEWS, reference docs, vignette section.

### Carried over
- steadyStates 1.4: solves M011 with and without the autocrine ligand in ~20 s; leaves some
  states free that 1.3 solved (C234, Smad7). Check against every TGFbModelling model before it
  becomes the default.
