# symmetryDetection: speed and gauge

Measured on TGFbModelling `M011_H1975TGFb_small`: 38 conditions, events at -2880 / -30 / 0,
138 coordinates, rank 120, 17 scalings, 1 general direction.

## Done

| Step | Effect on M011 |
|---|---|
| Lie order per condition (`NtChain`), capped at the block's own saturation | kernel call in a worker 43.5 s to 21.0 s |
| Joint saturation: lower the deepest blocks while the stacked rank stays | all 38 conditions at order 13 instead of 36; worker call 7.0 s |
| `symObsNullChainPointBatch`: (point, chain) jets in one OpenMP batch on the plain gap path | no fork; sparse entries sample through it |
| Sparse Laurent and general rational entries sample a whole stage per batch | were serial `kcall` loops |
| `gaugePreference` (default `FALSE`) | implemented, not yet run |

Switches: `DMOD_SYM_NOCAP`, `DMOD_SYM_NOJOINTCAP`, `DMOD_SYM_NOPOINTBATCH`.

## gaugePreference

`FALSE` one-step analysis as before. `NULL` any gauge. A character vector ranks the
coordinates to fix, `*` as wildcard (`"scale_*"`).

1. `reconstruct = FALSE`: rank and exact scalings.
2. Gauge: a basis of the scaling weights, greedy in order: preference rank, then coordinates
   outside the general directions' support, then the rest. A coordinate enters only if it
   raises the weight rank, so the gauge is always valid.
3. `fixed = gauge, reconstruct = TRUE`: only the general directions, fewer variables, in the
   gauge a fit uses.
4. Result: scalings as found, general directions in the gauge, `$gauge`; `print`, `summary`
   and `symmetryReduction` default `fixed` to it.

Which coordinates are fixed does not change rank or the number of general directions. It does
change their form: a fixed coordinate inside a direction's support adds scaling corrections to
it. Hence the preference for coordinates outside the support.

## Open

- Run `gaugePreference` on M011 and a model without events; add tests.
- The reconstruction of one general direction still dominates (M011: support {TGFB1,
  k_dg_TGFB1}, not finished after 9 min serial). Candidates:
  - ansatz on the relevant parameters (low degree, one batch of samples, exact check) in place
    of the sparse term-count search;
  - a higher `relevanceCap` so more entries take the dense, batched fit;
  - `reconstControl(timeout = )` returns the support; enough when the direction is read by hand.
- The relevance probe (one kernel per coordinate) is ~90 s on M011. The support from step 1 does
  not bound it: entries depend on more parameters than they are supported on.

## steadyStates 1.4 (on master)

- Solves M011 with and without the autocrine ligand in ~20 s, positive, mod p test passed;
  1.3 did not finish.
- Order: isolated states, rate constants, coupled states; targeted restart on dead ends.
  Leaves some states free that 1.3 solved (C234, Smad7 in M011): check before it becomes the
  default, against every model in TGFbModelling.
