dMod2 - Dynamic Modeling and Parameter Estimation in R
================

<!-- badges: start -->

[![R-CMD-check](https://github.com/dModverse/dMod2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dModverse/dMod2/actions/workflows/R-CMD-check.yaml)
[![docs](https://img.shields.io/badge/docs-dmodverse.github.io%2FdMod2-blue)](https://dmodverse.github.io/dMod2/)

<!-- badges: end -->

dMod2 builds ODE models of reaction networks, generates and compiles C++
code for them, and estimates their parameters from data. Prediction,
observation, parameter transformation and objective functions return
gradient and Hessian with their value, propagated forward or in one
reverse sweep.

The [website](https://dmodverse.github.io/dMod2/) has a [Get
started](https://dmodverse.github.io/dMod2/articles/dMod2.html) page, a
[guide](https://dmodverse.github.io/dMod2/articles/guide-model.html)
from the setup of a model to fits, profile likelihoods and mixed
effects, the derivations under Methods, the [example
projects](https://dmodverse.github.io/dMod2/articles/index.html) and the
reference.

## System requirements

Model code is generated and compiled at run time, so a **C++17
compiler** and **Python with SymPy** are required; reticulate provides
the Python packages on first use. On Linux the compilers are usually
present, Windows users need
[Rtools](https://cran.r-project.org/bin/windows/Rtools/). [Installation
and
troubleshooting](https://dmodverse.github.io/dMod2/articles/installation.html)
lists the optional parts and the errors of a first compile. The default
backend is [cppDE](https://github.com/dModverse/cppDE), which generates
C++ with forward and reverse derivatives to second order;
[cOde](https://github.com/dkaschek/cOde) with deSolve is available as an
alternative.

## Installation

``` r
remotes::install_github("dModverse/dMod2")
```

The backend [cppDE](https://github.com/dModverse/cppDE) comes along as a
`Remotes:` dependency and does not have to be installed separately.
Working on the package from a checkout is described in
[CONTRIBUTING.md](.github/CONTRIBUTING.md).

## Building a model

Models are assembled from four kinds of function objects that compose
with `*` (concatenation, chain rule) and `+` (direct sum over
experimental conditions):

| Object | Constructor | Role |
|----|----|----|
| `odemodel`, `prdfn` | `odemodel()`, `Xs()` | the compiled ODE with sensitivities |
| `obsfn` | `Y()` | observation and error functions |
| `parfn` | `P()` | parameter transformations |
| `objfn` | `normL2()`, `constraintL2()` | objective functions |

The [example
projects](https://dmodverse.github.io/dMod2/articles/index.html) in
[`examples/`](https://github.com/dModverse/dMod2/tree/master/examples)
work through published models of signaling pathways and pharmacometrics,
each with data, scripts and saved results.

## Mixed-effects models

Random effects are written as effect formulas on the parameter
transformation, in the notation of lme4, and `em()` fits them by
marginal maximum likelihood: FOCEI, adaptive sparse-grid quadrature,
Laplace EM or SAEM. Besides normal random effects, Laplace,
spike-and-slab and Dirichlet-process ones select which groups differ,
with their hyperparameters estimated in the same fit. The guide covers
[mixed-effects
models](https://dmodverse.github.io/dMod2/articles/guide-nlme.html) and
[group-specific
effects](https://dmodverse.github.io/dMod2/articles/guide-group-effects.html).

## Alignment

`align()` brings measurements in relative units from separate
experiments, such as Western blots and qPCR runs, to common values, with
the scaling factors of each experiment written as free expressions,
optionally with random effects per group
([guide](https://dmodverse.github.io/dMod2/articles/guide-alignment.html)).

## Structural identifiability

`symmetryDetection()` finds the directions in parameters and initial
values along which the observations of a model do not change, and
`symmetryReduction()` turns them into a parameter transformation that
removes them. Both run in the Python package
[symident](https://github.com/dModverse/symident), which dMod2 installs
from PyPI on first use. If `RETICULATE_PYTHON` points to your own
Python, install it there with `pip install symident`. A steady-state
start on coupled steady states needs msolve, which `symmetryDetection()`
offers to build when it is missing, or which `SYMIDENT_MSOLVE` points
to.

## Distributed computing

`distributedComputing()` runs fits, profiles or any other R code as
tasks in the background: on the local machine, on machines reached by
ssh with the cores each may use, or as an array job of a SLURM cluster
through a job script template. The
[guide](https://dmodverse.github.io/dMod2/articles/guide-remote-computing.html)
shows how.

## Citation

There is no publication for dMod2 yet. The paper below describes the
predecessor [dMod](https://github.com/JetiLab/dMod) and is the reference
for the modeling framework:

Kaschek D, Mader W, Fehling-Kaschek M, Rosenblatt M, Timmer J (2019).
Dynamic Modeling, Parameter Estimation, and Uncertainty Analysis in R.
*Journal of Statistical Software*, 88(10), 1-32.
[doi:10.18637/jss.v088.i10](https://doi.org/10.18637/jss.v088.i10)
