# Loaded Shared Objects in the Working Directory

Only the working directory is searched, not `getOption("dMod2.outdir")`
or the directories the objects were generated in.

## Usage

``` r
getLocalDLLs()
```

## Value

Character vector with the names of the loaded shared objects whose path
lies in the working directory.

## See also

[`loadDLL()`](https://dmodverse.github.io/dMod2/reference/loadDLL.md)

## Examples

``` r
## Names of the loaded shared objects below the working directory
getLocalDLLs()
#> NULL
```
