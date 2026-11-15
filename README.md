dMod2: Parameter Estimation in Dynamical Models
================

<!-- badges: start -->

[![R-CMD-check](https://github.com/dModverse/dMod2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dModverse/dMod2/actions/workflows/R-CMD-check.yaml)
[![docs](https://img.shields.io/badge/docs-dmodverse.github.io%2FdMod2-blue)](https://dmodverse.github.io/dMod2/)

<!-- badges: end -->

dMod2 builds ODE models of reaction networks, compiles them to C/C++
with exact derivatives, and estimates their parameters from data. It
succeeds [dMod](https://github.com/JetiLab/dMod); [Migrating from
dMod](https://dmodverse.github.io/dMod2/articles/migrating.html) lists
what a dMod script has to change.

## Installation

``` r
remotes::install_github("dModverse/dMod2")
```

The ODE backend [cppDE](https://github.com/dModverse/cppDE) is installed
with it. Models compile at run time and need a C++17 compiler and Python
with SymPy, see [Installation and
troubleshooting](https://dmodverse.github.io/dMod2/articles/installation.html).
Working on the package from a checkout is described in
[CONTRIBUTING.md](https://github.com/dModverse/dMod2/blob/master/.github/CONTRIBUTING.md).

## Documentation

The [website](https://dmodverse.github.io/dMod2/) starts with [Get
started](https://dmodverse.github.io/dMod2/articles/dMod2.html). The
guide covers

- [model
  setup](https://dmodverse.github.io/dMod2/articles/guide-model.html)
  and [structural
  identifiability](https://dmodverse.github.io/dMod2/articles/guide-symmetries.html),
- [sensitivities](https://dmodverse.github.io/dMod2/articles/guide-sensitivities.html),
  [optimization](https://dmodverse.github.io/dMod2/articles/guide-optimization.html),
  [profile
  likelihood](https://dmodverse.github.io/dMod2/articles/guide-profiles.html)
  and [prediction
  uncertainty](https://dmodverse.github.io/dMod2/articles/guide-prediction-uncertainty.html),
- [mixed-effects
  models](https://dmodverse.github.io/dMod2/articles/guide-nlme.html)
  and [sparse
  deviations](https://dmodverse.github.io/dMod2/articles/guide-sparse-deviations.html),
- [alignment of relative
  data](https://dmodverse.github.io/dMod2/articles/guide-alignment.html),
  [PEtab and
  SBML](https://dmodverse.github.io/dMod2/articles/guide-petab.html) and
  [remote
  computing](https://dmodverse.github.io/dMod2/articles/guide-remote-computing.html).

Methods, example projects and the reference are listed in the
[articles](https://dmodverse.github.io/dMod2/articles/index.html) and
the [reference
index](https://dmodverse.github.io/dMod2/reference/index.html).

## Citation

There is no publication for dMod2 yet. The paper below describes the
predecessor dMod and is the reference for the modeling framework:

Kaschek D, Mader W, Fehling-Kaschek M, Rosenblatt M, Timmer J (2019).
Dynamic Modeling, Parameter Estimation, and Uncertainty Analysis in R.
*Journal of Statistical Software*, 88(10), 1-32.
[doi:10.18637/jss.v088.i10](https://doi.org/10.18637/jss.v088.i10)
