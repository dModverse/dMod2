# Read the Tables of a PEtab Problem

Reads the TSV tables named in a PEtab YAML file into data frames. The
SBML model is not read.

## Usage

``` r
readPEtabTables(yamlPath)
```

## Arguments

- yamlPath:

  Path to the PEtab YAML file.

## Value

A list with the data frames `parameters`, `conditions`, `measurements`,
`observables`, `experiments` and `mapping` (each `NULL` when the problem
has no such table), `sbmlPath`, the path of the first model,
`sbmlPaths`, the SBML path of every model id, and `formatVersion`
(integer).

## See also

[`readPEtabYaml()`](https://dmodverse.github.io/dMod2/reference/readPEtabYaml.md),
[`importPEtab()`](https://dmodverse.github.io/dMod2/reference/importPEtab.md)

## Examples

``` r
yaml <- system.file("extdata/petab_boehm/Boehm.yaml", package = "dMod2")
tables <- readPEtabTables(yaml)
head(tables$measurements)
#>   observableId preequilibrationConditionId simulationConditionId measurement
#> 1  pSTAT5A_rel                          NA          model1_data1    7.901073
#> 2  pSTAT5A_rel                          NA          model1_data1   66.363494
#> 3  pSTAT5A_rel                          NA          model1_data1   81.171324
#> 4  pSTAT5A_rel                          NA          model1_data1   94.730308
#> 5  pSTAT5A_rel                          NA          model1_data1   95.116483
#> 6  pSTAT5A_rel                          NA          model1_data1   91.441717
#>   time observableParameters noiseParameters                datasetId
#> 1  0.0                   NA  sd_pSTAT5A_rel model1_data1_pSTAT5A_rel
#> 2  2.5                   NA  sd_pSTAT5A_rel model1_data1_pSTAT5A_rel
#> 3  5.0                   NA  sd_pSTAT5A_rel model1_data1_pSTAT5A_rel
#> 4 10.0                   NA  sd_pSTAT5A_rel model1_data1_pSTAT5A_rel
#> 5 15.0                   NA  sd_pSTAT5A_rel model1_data1_pSTAT5A_rel
#> 6 20.0                   NA  sd_pSTAT5A_rel model1_data1_pSTAT5A_rel
tables$parameters
#>             parameterId          parameterName parameterScale lowerBound
#> 1  Epo_degradation_BaF3 EPO_{degradation,BaF3}          log10      1e-05
#> 2          k_exp_hetero         k_{exp,hetero}          log10      1e-05
#> 3            k_exp_homo           k_{exp,homo}          log10      1e-05
#> 4          k_imp_hetero         k_{imp,hetero}          log10      1e-05
#> 5            k_imp_homo           k_{imp,homo}          log10      1e-05
#> 6                k_phos               k_{phos}          log10      1e-05
#> 7                 ratio                  ratio            lin     -5e+00
#> 8        sd_pSTAT5A_rel  \\sigma_{pSTAT5A,rel}          log10      1e-05
#> 9        sd_pSTAT5B_rel  \\sigma_{pSTAT5B,rel}          log10      1e-05
#> 10       sd_rSTAT5A_rel  \\sigma_{rSTAT5A,rel}          log10      1e-05
#> 11              specC17                specC17            lin     -5e+00
#>    upperBound nominalValue estimate
#> 1      100000 2.698251e-02        1
#> 2      100000 1.000680e-05        1
#> 3      100000 6.170228e-03        1
#> 4      100000 1.636792e-02        1
#> 5      100000 9.774938e+04        1
#> 6      100000 1.576651e+04        1
#> 7           5 6.930000e-01        0
#> 8      100000 3.852612e+00        1
#> 9      100000 6.591478e+00        1
#> 10     100000 3.152713e+00        1
#> 11          5 1.070000e-01        0
```
