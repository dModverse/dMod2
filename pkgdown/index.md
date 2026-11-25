# dMod2

dMod2 builds ODE models of reaction networks, compiles them to C/C++ with exact derivatives,
and estimates their parameters from data.

```r
remotes::install_github("dModverse/dMod2")
```

## Documentation

- [Get started](articles/dMod2.html): a model of bile acid transport from reactions to profile
  likelihoods.
- [Installation and troubleshooting](articles/installation.html): compilers, Python, SBML, symident.
- [Migrating from dMod](articles/migrating.html): what changed for users of dMod.
- [Guide](articles/guide-model.html): the NFκB project of Oppelt et al. (2018), chapter by chapter,
  from the setup of the model to fits, profiles and prediction bands, then mixed effects, sparse
  deviations, alignment, PEtab and remote computing.
- [Methods](articles/index.html#methods): the derivations behind the guide.
- [Examples](articles/index.html#examples-signaling-models): projects with data, scripts and saved results.
- [Reference](reference/index.html): every exported function.

The ODE backend is [cppDE](https://dmodverse.github.io/cppDE/), structural identifiability runs in
the Python package [symident](https://dmodverse.github.io/symident/).
