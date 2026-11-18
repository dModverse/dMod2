# dMod2

dMod2 builds ODE models of reaction networks, compiles them to C/C++ with exact derivatives,
and estimates their parameters from data.

```r
remotes::install_github("dModverse/dMod2")
```

## Where to go

- [Get started](articles/dMod2.html): a model of bile acid transport from reactions to profile
  likelihoods.
- [Installation and troubleshooting](articles/installation.html): compilers, Python, SBML, symident.
- [Migrating from dMod](articles/migrating.html): what changed for users of dMod.
- [Guide](articles/guide-model.html): the NFκB project of Oppelt et al. (2018), chapter by chapter,
  from the setup of the model to fits, profiles and prediction bands, then mixed effects, sparse
  deviations, alignment, PEtab and remote computing.
- [Methods](articles/index.html): the derivations behind the guide.
- [Examples](articles/index.html): projects with data, scripts and saved results.
- [Reference](reference/index.html): every exported function.

The ODE backend is [cppDE](https://dmodverse.github.io/cppDE/), structural identifiability runs in
the Python package [symident](https://github.com/dModverse/symident).

## Citation

There is no publication for dMod2 yet. The paper below describes the predecessor [dMod](https://github.com/JetiLab/dMod) and is the reference for the modeling framework:

Kaschek D, Mader W, Fehling-Kaschek M, Rosenblatt M, Timmer J (2019). Dynamic Modeling, Parameter Estimation, and Uncertainty Analysis in R. *Journal of Statistical Software*, 88(10), 1-32. [doi:10.18637/jss.v088.i10](https://doi.org/10.18637/jss.v088.i10)
