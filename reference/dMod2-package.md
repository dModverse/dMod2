# dMod2: Parameter Estimation in Dynamical Models

Models are built from reactions or equations, compiled to C++ with
sensitivities, and composed from observation, prediction and parameter
transformation functions as `g * x * p`. Objective functions over
several experimental conditions are minimized by trust-region
optimization, and parameter uncertainty is assessed by profile
likelihoods.

## Options

- `dMod2.cores`:

  Cores over the experimental conditions inside one evaluation of a
  prediction, observation or objective, the default of their `cores`
  argument. Default 1.

- `dMod2.outdir`:

  Directory the model constructors such as
  [`odemodel()`](https://dmodverse.github.io/dMod2/reference/odemodel.md),
  [`Y()`](https://dmodverse.github.io/dMod2/reference/Y.md) and
  [`P()`](https://dmodverse.github.io/dMod2/reference/P.md) write their
  sources and shared objects to, and where
  [`distributedComputing()`](https://dmodverse.github.io/dMod2/reference/distributedComputing.md)
  looks for model files when the input objects name none. Default the
  working directory.

- `dMod2.evaluation`:

  `"C++"` evaluates the objectives of
  [`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md),
  with their predictions and the random-effects fits of
  [`em()`](https://dmodverse.github.io/dMod2/reference/em.md) built on
  them, in C++ per condition where the model allows it (see Details of
  [`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md));
  `"R"` evaluates every model through the function objects in R. Default
  `"C++"`.

- `dMod2.compile.cmdlimit`:

  Longest command line in characters that
  [`compile()`](https://dmodverse.github.io/dMod2/reference/compile.md)
  passes to the compiler and archiver. Default 96000, 8000 on Windows.

- `dMod2.sparseAuto`:

  Thresholds of `sparse = "auto"` in
  [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md) and
  [`em()`](https://dmodverse.github.io/dMod2/reference/em.md), a list of
  `parameters` and `fill`, see
  [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md).
  Default `list(parameters = 500, fill = 0.05)`.

## Environment variables

- `DMOD_LIBSBML_PYTHON`:

  Python interpreter with libsbml for the SBML interface.

- `SYMIDENT_<NAME>`:

  Switch `<name>` of symident for
  [`symmetryDetection()`](https://dmodverse.github.io/dMod2/reference/symmetryDetection.md)
  and
  [`symmetryReduction()`](https://dmodverse.github.io/dMod2/reference/symmetryReduction.md),
  see the Details there.

## See also

Models:
[`eqnlist()`](https://dmodverse.github.io/dMod2/reference/eqnlist.md),
[`addReaction()`](https://dmodverse.github.io/dMod2/reference/addReaction.md),
[`odemodel()`](https://dmodverse.github.io/dMod2/reference/odemodel.md),
[`Xs()`](https://dmodverse.github.io/dMod2/reference/Xs.md),
[`Y()`](https://dmodverse.github.io/dMod2/reference/Y.md),
[`P()`](https://dmodverse.github.io/dMod2/reference/P.md). Data and
objectives:
[`as.datalist()`](https://dmodverse.github.io/dMod2/reference/datalist.md),
[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md),
[`constraintL2()`](https://dmodverse.github.io/dMod2/reference/constraintL2.md).
Fitting and uncertainty:
[`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md),
[`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md),
[profile()](https://dmodverse.github.io/dMod2/reference/profile.md),
[`plotProfile()`](https://dmodverse.github.io/dMod2/reference/plotProfile.md).
Exchange formats:
[`importPEtab()`](https://dmodverse.github.io/dMod2/reference/importPEtab.md),
[`exportPEtab()`](https://dmodverse.github.io/dMod2/reference/exportPEtab.md).
Identifiability:
[`symmetryDetection()`](https://dmodverse.github.io/dMod2/reference/symmetryDetection.md).
Articles: [Get
started](https://dmodverse.github.io/dMod2/articles/dMod2.html),
[Installation and
troubleshooting](https://dmodverse.github.io/dMod2/articles/installation.html).

## Author

**Maintainer**: Simon Beyer <simon.beyer@fdm.uni-freiburg.de>

Authors:

- Simon Beyer <simon.beyer@fdm.uni-freiburg.de>

- Daniel Kaschek <daniel.kaschek@gmail.com> (Created dMod 1.0; original
  design, philosophy)
