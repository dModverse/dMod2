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
