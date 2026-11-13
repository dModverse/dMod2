# Import a PEtab Problem

Reads a PEtab v1 or v2 problem, a YAML file with an SBML model and TSV
tables, and builds the dMod problem: prediction, observation and error
functions, parameter transformation, data and objective. The SBML model
is read by
[`importSbml()`](https://dmodverse.github.io/dMod2/reference/importSbml.md).

## Usage

``` r
importPEtab(
  yamlPath,
  backend,
  compile = TRUE,
  cores = 1L,
  modelname = NULL,
  deriv = TRUE,
  derivMode = "forward",
  options = NULL,
  optionsSens = NULL,
  sparse = NULL,
  outdir = .dmodOutdir()
)
```

## Arguments

- yamlPath:

  Path to the PEtab YAML manifest.

- backend:

  Required, one of `"deSolve"`, `"cppDE"` or `"Sundials"`, passed to
  [`odemodel()`](https://dmodverse.github.io/dMod2/reference/odemodel.md).

- compile:

  Logical. If `TRUE` (default), the generated sources are compiled.
  `FALSE` leaves them uncompiled for inspection.

- cores:

  Number of parallel compilation jobs, passed to
  [`compile()`](https://dmodverse.github.io/dMod2/reference/compile.md).
  Default `1`.

- modelname:

  Base name of the generated files. Default `NULL`, the base name of the
  YAML file.

- deriv:

  Logical. If `FALSE`, the model is built without sensitivities; the
  objective is then evaluated with `deriv = FALSE`. Default `TRUE`.

- derivMode:

  Derivative directions to build, passed to
  [`odemodel()`](https://dmodverse.github.io/dMod2/reference/odemodel.md):
  `"forward"` (default), `"reverse"` or `c("forward", "reverse")`.
  `"reverse"` needs `backend = "cppDE"` or `"Sundials"`. First order
  only.

- options, optionsSens:

  Lists passed to
  [`Xs()`](https://dmodverse.github.io/dMod2/reference/Xs.md): the
  solver options of every solve and their overrides for the solves with
  sensitivities. Default `NULL`, the backend's defaults.

- sparse:

  `NULL` (default) lets the cppDE or Sundials backend choose a sparse
  (KLU) or dense linear solver; `TRUE` or `FALSE` sets it.

- outdir:

  Directory of the generated sources and shared objects, one source per
  condition. Default `getOption("dMod2.outdir")`, else the working
  directory.

## Value

A list of class `"petabproblem"` with components `dataList`,
`reactions`, `odemodel`, `g`, `x`, `p`, `e`, `prd` (the composite
`g * x * p`), `obj`, `bestfit`, `parlower` and `parupper`. `obj` has the
PEtab fixed parameters set; `bestfit` is the PEtab `nominalValue` of
each estimated parameter, with the PEtab scales in the attribute
`"petab_scales"`. The attribute `"petab_meta"` of the list holds what
[`exportPEtabObject()`](https://dmodverse.github.io/dMod2/reference/exportPEtabObject.md)
needs, among it the fixed parameters as `fixed`.

## Details

Supported: one problem per YAML file with any number of v2 model files;
preequilibration (at most two periods per v2 experiment); observable
transformations `lin`, `log` and `log10`; the noise distribution
`normal`, and `log-normal` in v2; a v2 mapping table. Priors of the
distributions `uniform`, `normal`, `log-normal`, `cauchy`, `chisquare`,
`exponential`, `gamma`, `laplace`, `log-laplace`, `log-uniform` and
`rayleigh`, also with the `parameterScale` prefix, are added to `obj`,
truncated to the parameter bounds. Not supported: model languages other
than SBML, other noise distributions, and a plain `normal` or `laplace`
prior on a parameter whose scale is not `lin`.

## See also

[`exportPEtabObject()`](https://dmodverse.github.io/dMod2/reference/exportPEtabObject.md),
[`readPEtabTables()`](https://dmodverse.github.io/dMod2/reference/readPEtabTables.md),
[`importSbml()`](https://dmodverse.github.io/dMod2/reference/importSbml.md);
the guide chapter
[PEtab](https://dmodverse.github.io/dMod2/articles/guide-petab.html).

## Examples

``` r
# \donttest{
yaml <- system.file("extdata/petab_boehm/Boehm.yaml", package = "dMod2")
petab <- importPEtab(yaml, backend = "cppDE", cores = 1, outdir = tempdir())
#> Sparse Jacobian detected (8x8, 21 nnz, 67.188% sparse)
#> Sparse Jacobian detected (8x8, 21 nnz, 67.188% sparse)
#> using C++ compiler: g++
#>   every source : -O2 -DNDEBUG -w -fPIC -DKLU -I/usr/include/suitesparse
#>   2 of 5 also  : -DKLUBTF=0 -DKLUAMD=0 -fopenmp
petab
#> <PEtab problem Boehm_2>
#>   source:        /home/runner/work/_temp/Library/dMod2/extdata/petab_boehm/Boehm.yaml
#>   conditions:    1 (1 sim, 0 preeq)
#>   observables:   pSTAT5A_rel, pSTAT5B_rel, rSTAT5A_rel
#>   measurements:  48
#>   bestfit (n=9): Epo_degradation_BaF3, k_exp_hetero, k_exp_homo, k_imp_hetero, k_imp_homo, k_phos, sd_pSTAT5A_rel, sd_pSTAT5B_rel, sd_rSTAT5A_rel
#>   fixed   (n=4): ratio, specC17, cyt, nuc

# The objective at the nominal values
petab$obj(petab$bestfit)$value
#> [1] 276.4441

# The prediction needs the fixed parameters
fixed <- attr(petab, "petab_meta")$fixed
times <- seq(0, 240, length.out = 61)
plot(petab$prd(times, c(petab$bestfit, fixed)), petab$dataList)


# Write the problem back to PEtab
yamlOut <- exportPEtabObject(petab, file.path(tempdir(), "boehm"),
                             formatVersion = "1")
readPEtabTables(yamlOut)$parameters
#>             parameterId parameterScale lowerBound upperBound nominalValue
#> 1  Epo_degradation_BaF3          log10      1e-05      1e+05 2.698251e-02
#> 2          k_exp_hetero          log10      1e-05      1e+05 1.000680e-05
#> 3            k_exp_homo          log10      1e-05      1e+05 6.170228e-03
#> 4          k_imp_hetero          log10      1e-05      1e+05 1.636792e-02
#> 5            k_imp_homo          log10      1e-05      1e+05 9.774938e+04
#> 6                k_phos          log10      1e-05      1e+05 1.576651e+04
#> 7        sd_pSTAT5A_rel          log10      1e-05      1e+05 3.852612e+00
#> 8        sd_pSTAT5B_rel          log10      1e-05      1e+05 6.591478e+00
#> 9        sd_rSTAT5A_rel          log10      1e-05      1e+05 3.152713e+00
#> 10                ratio            lin       -Inf        Inf 6.930000e-01
#> 11              specC17            lin       -Inf        Inf 1.070000e-01
#>    estimate
#> 1         1
#> 2         1
#> 3         1
#> 4         1
#> 5         1
#> 6         1
#> 7         1
#> 8         1
#> 9         1
#> 10        0
#> 11        0
# }
```
