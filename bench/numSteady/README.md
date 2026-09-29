# numSteady

Numerical steady states with `Pimpl()` on a real model and on the
Benchmark-Models-PEtab problems.

## Model

`model.R`: TGF-beta receptor model reduced from TGFbModelling
M016_H1975TGFb, with H1975 data (`data.csv`, 933 points, 29 conditions).
Receptor transcripts with siRNA and plasmid copies, receptor trafficking,
ligand binding, a Smad core with three conserved totals, negative feedback by
Smad7 and positive feedback by FB3. The steady state differs between the five
arms Ctrl, R1Knd, R2Knd, OE_R1 and OE_R2; FB3 makes some parameter sets
bistable.

## Scripts

Run from the package root.

| Script | Content |
|---|---|
| `battery.R [ndraw] [nwalk]` | Cold solves at random parameters in all arms, warm solves along random walks; residual, conservation, stability, IFT sensitivities against finite differences |
| `exact.R [n]` | Steady states from `steadyStates()` at random values of its free symbols; `Pimpl()` has to return them |
| `fit.R [nfits] [cores]` | One `trust()` fit and a multistart with the `Pimpl()` statistics |
| `petab.R <dir> <out.rds> [problems] [nevals]` | Benchmark problems: objective, predictions against `simulatedData`, objective and gradient at perturbed parameters |

## Results

dMod2 0.9.0, cppDE 0.11.0.

- `battery.R 60 10`: 300/300 cold and 250/250 warm solves, median 0.008 s,
  4.9 iterations per solve. Conservation error 6.5e-16, all roots stable,
  IFT sensitivities within 2.4e-3 of central differences.
- `exact.R 40`: of 190 stable exact states, 185 are returned from a cold
  start (median relative error at most 2e-12) and 5 lie in a second stable
  basin; none wrong, none failed. Started at the exact state, `Pimpl()`
  returns it unchanged.
- `fit.R 16 8`: in the `trust()` fit 95% of the calls are answered from
  memory, every solve starts from the kept root of its condition, 4.0
  iterations per solve, no failure. Multistart: 12 of 16 fits converge.
- `petab.R` on 10 problems of the collection, against dMod2 0.8.2 with
  `Pequil()`:

| Problem | 0.8.2 | 0.9.0 |
|---|---|---|
| Brannmark_JBC2010 | 283.7741, predictions to 3.4e-3 | 283.7784, predictions to 1.4e-5 |
| Isensee_JCB2018 | 4 of 20 evaluations fail | none fail, predictions to 6.1e-10 |
| Raimundez_PCB2020 | 2002.0190, 9.9 s per evaluation | 2002.0103, 8.9 s per evaluation |
| Zheng_PNAS2012 | 0.027 s per evaluation | 0.018 s per evaluation |
| other six | unchanged | unchanged |

Values are -2 log L. Brannmark matches the AMICI reference log L = -141.8891.
