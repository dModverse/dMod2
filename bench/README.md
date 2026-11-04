# bench

Scripts outside the built package, run from the package root against the installed
package. The scripts are linear and meant to be stepped through in RStudio.

| folder | what |
|---|---|
| `harness/` | the PEtab benchmark suite: value and gradient of the whole objective chain, see its README |
| `gradient/` | gradient cost over the number of parameters, forward against reverse, the discrete adjoint against CVODES ASA |
| `hessian/` | the Hessian sources of `trust()` on multi-start fits |
| `symmetry/` | `symmetryDetection()` and `symmetryReduction()` on models from several fields |
| `multipleShooting/` | multiple shooting fits on oscillating and switching models |
| `numSteady/` | numerical steady states with `Pimpl()` on a real model and on PEtab problems |
| `nlme/` | the methods of `em()` against nlmixr2 as multistarts from the same starts, one core per fit: theophylline, warfarin PK/PD, mavoglurant with occasions, and TMDD over simulated data sets (bias, RMSE, coverage); for warfarin, mavoglurant and TMDD also single fits on 10 cores from starts halfway out |
| `sparse/` | sparse against dense Hessians on simulated alignments: `trust()` over the number of experiments, `em()` over the number of groups, one `align()` of 1.65e6 readings |
| `structureSelection/` | cell-line specific parameters: `dirichletProcess()` against the L_q scan of Hauber et al. (2023), with known truth |
