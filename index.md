# dMod2

dMod2 builds ODE models of reaction networks, compiles them to C/C++
with exact derivatives, and estimates their parameters from data.

``` r

remotes::install_github("dModverse/dMod2")
```

## Documentation

- [Get started](https://dmodverse.github.io/dMod2/articles/dMod2.md): a
  model of bile acid transport from reactions to profile likelihoods.
- [Installation and
  troubleshooting](https://dmodverse.github.io/dMod2/articles/installation.md):
  compilers, Python, SBML, symident.
- [Migrating from
  dMod](https://dmodverse.github.io/dMod2/articles/migrating.md): what
  changed for users of dMod.
- [Guide](https://dmodverse.github.io/dMod2/articles/guide-model.md):
  the NFκB project of Oppelt et al. (2018), chapter by chapter, from the
  setup of the model to fits, profiles and prediction bands, then mixed
  effects, sparse deviations, alignment, PEtab and remote computing.
- [Methods](https://dmodverse.github.io/dMod2/articles/index.html#methods):
  the derivations behind the guide.
- [Examples](https://dmodverse.github.io/dMod2/articles/index.html#examples-signaling-models):
  projects with data, scripts and saved results.
- [Reference](https://dmodverse.github.io/dMod2/reference/index.md):
  every exported function.

The ODE backend is [cppDE](https://dmodverse.github.io/cppDE/),
structural identifiability runs in the Python package
[symident](https://dmodverse.github.io/symident/).
