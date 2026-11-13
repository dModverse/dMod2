# Export an Imported PEtab Problem

Writes a problem returned by
[`importPEtab()`](https://dmodverse.github.io/dMod2/reference/importPEtab.md)
back to PEtab: the parameter, observable, condition and measurement
tables, for v2 also the experiment table, the SBML model and a YAML
file. Conditions the importer split are merged again. For a problem
built in dMod, use
[`exportPEtab()`](https://dmodverse.github.io/dMod2/reference/exportPEtab.md).

## Usage

``` r
exportPEtabObject(
  petab,
  dir,
  modelID = NULL,
  formatVersion = "2.0.0",
  overwrite = FALSE
)
```

## Arguments

- petab:

  A `petabproblem` returned by
  [`importPEtab()`](https://dmodverse.github.io/dMod2/reference/importPEtab.md),
  or a list with the same components.

- dir:

  Output directory, created if missing.

- modelID:

  SBML model identifier. Default `NULL`: the model id of the imported
  problem, else `"dMod_export"`.

- formatVersion:

  `"2.0.0"` (default) or `"1"`.

- overwrite:

  Logical. If `FALSE` (default), existing files in `dir` are an error.

## Value

Path to the written YAML file, invisibly.

## Details

Limits:

- Parameter scales are read from `attr(petab$bestfit, "petab_scales")`;
  without it they are `"lin"`.

- v2 has no parameter scales: values and bounds are written on the
  linear scale, with a warning.

- v2 has no `log10` observable transformation: it is written as
  `"log-normal"` noise, with a warning.

- Initial values set per condition only in the parameter transformation
  are not written to a v2 condition table; use `formatVersion = "1"` to
  keep them.

- SBML `<algebraicRule>` elements are not imported and therefore not
  written.

## See also

[`importPEtab()`](https://dmodverse.github.io/dMod2/reference/importPEtab.md),
[`exportPEtab()`](https://dmodverse.github.io/dMod2/reference/exportPEtab.md)

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
#> <PEtab problem Boehm>
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
