# Symmetry motifs

Small models on which `symmetryDetection()` / `symmetryReduction()` of dMod2 0.8.1 were
slow or failed, one difficulty each (see `motifs.R`). Results and history: `dev/planSymmetry.md`.

| Script | Does |
|---|---|
| `motifs.R` | `.motifs$<name>()` returns the `symmetryDetection()` arguments |
| `run_motif.R <name> [outdir] [extra args]` | detect with reconstruction, reduce, re-detect the reduced model; writes `<outdir>/<name>.rds` and one summary line |
| `run_all.sh [names]` | all motifs (or the named ones) in parallel, each capped at `TMO` seconds (default 1200) |

Environment: `DMOD_LOADER` (a script attaching dMod2, default `library(dMod2)`), `OUT`
(result directory), `TMO`, `EXTRA` (extra `symmetryDetection()` arguments as R code),
`DMOD_MOTIF_CACHE` (where steady states are cached; runs sharing it must not start
together, or the first two may race on the cache file).

```bash
OUT=results TMO=600 ./run_all.sh cat3 route4 autocrine_ss
```
