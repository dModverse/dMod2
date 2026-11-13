# dMod2

dMod2 builds ODE models of reaction networks, compiles them to C/C++
with exact derivatives, and estimates their parameters from data.

``` r

remotes::install_github("dModverse/dMod2")
```

## Where to go

- [Get started](https://dmodverse.github.io/dMod2/articles/dMod2.md): a
  two-state model from reactions to profile likelihoods.
- [Installation and
  troubleshooting](https://dmodverse.github.io/dMod2/articles/installation.md):
  compilers, Python, SBML, symident.
- [Migrating from
  dMod](https://dmodverse.github.io/dMod2/articles/migrating.md): what
  changed for users of dMod.
- [Guide](https://dmodverse.github.io/dMod2/articles/guide-model.md):
  the NFκB project of Oppelt et al. (2018), chapter by chapter, from the
  setup of the model to fits, profiles and mixed effects.
- [Methods](https://dmodverse.github.io/dMod2/articles/index.md): the
  derivations behind the guide.
- [Examples](https://dmodverse.github.io/dMod2/articles/index.md):
  projects with data, scripts and saved results.
- [Reference](https://dmodverse.github.io/dMod2/reference/index.md):
  every exported function.

The ODE backend is [cppDE](https://dmodverse.github.io/cppDE/),
structural identifiability runs in the Python package
[symident](https://github.com/dModverse/symident).

## Citation

There is no publication for dMod2 yet. The paper below describes the
predecessor [dMod](https://github.com/JetiLab/dMod) and is the reference
for the modeling framework:

Kaschek D, Mader W, Fehling-Kaschek M, Rosenblatt M, Timmer J (2019).
Dynamic Modeling, Parameter Estimation, and Uncertainty Analysis in R.
*Journal of Statistical Software*, 88(10), 1-32.
[doi:10.18637/jss.v088.i10](https://doi.org/10.18637/jss.v088.i10)
