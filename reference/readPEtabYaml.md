# Read a PEtab YAML File

Reads a PEtab v1 or v2 YAML file, as selected by the `format_version`
entry, and resolves the paths of the files it names. Only one problem
per file is supported.

## Usage

``` r
readPEtabYaml(yamlPath)
```

## Arguments

- yamlPath:

  Path to the PEtab YAML file.

## Value

A list with `baseDir`, `formatVersion` (integer major version),
`parameterFile`, and a one-element list `problems` whose entry holds the
absolute paths `sbmlFile`, `conditionFile`, `measurementFile` and
`observableFile`. For v2 the entry also holds `experimentFile` and
`mappingFile` (possibly `NULL`), `modelID`, the first model id, and
`models`, the SBML path of every model id.

## See also

[`readPEtabTables()`](https://dmodverse.github.io/dMod2/reference/readPEtabTables.md),
[`importPEtab()`](https://dmodverse.github.io/dMod2/reference/importPEtab.md)

## Examples

``` r
yaml <- system.file("extdata/petab_boehm/Boehm.yaml", package = "dMod2")
str(readPEtabYaml(yaml))
#> List of 4
#>  $ baseDir      : chr "/home/runner/work/_temp/Library/dMod2/extdata/petab_boehm"
#>  $ formatVersion: int 1
#>  $ parameterFile: chr "/home/runner/work/_temp/Library/dMod2/extdata/petab_boehm/parameters_Boehm_JProteomeRes2014.tsv"
#>  $ problems     :List of 1
#>   ..$ :List of 8
#>   .. ..$ sbmlFile       : chr "/home/runner/work/_temp/Library/dMod2/extdata/petab_boehm/model_Boehm_JProteomeRes2014.xml"
#>   .. ..$ conditionFile  : chr "/home/runner/work/_temp/Library/dMod2/extdata/petab_boehm/experimentalCondition_Boehm_JProteomeRes2014.tsv"
#>   .. ..$ measurementFile: chr "/home/runner/work/_temp/Library/dMod2/extdata/petab_boehm/measurementData_Boehm_JProteomeRes2014.tsv"
#>   .. ..$ observableFile : chr "/home/runner/work/_temp/Library/dMod2/extdata/petab_boehm/observables_Boehm_JProteomeRes2014.tsv"
#>   .. ..$ experimentFile : NULL
#>   .. ..$ mappingFile    : NULL
#>   .. ..$ modelID        : chr NA
#>   .. ..$ models         : Named chr "/home/runner/work/_temp/Library/dMod2/extdata/petab_boehm/model_Boehm_JProteomeRes2014.xml"
#>   .. .. ..- attr(*, "names")= chr "model"
```
