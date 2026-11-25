# Load Fits Written by mstrust

Reads the per-fit files that
[`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md)
writes with `output = TRUE`, for example after an aborted run.

## Usage

``` r
load.parlist(folder)
```

## Arguments

- folder:

  The folder `<resultPath>/<name>/trial-<k>-<date>/interRes` of the run,
  which holds one `fit-<i>.Rda` per completed fit.

## Value

A parlist with one entry per file, in the order of the fit numbers.

## See also

[`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md)

## Examples

``` r
folder <- file.path(tempdir(), "interRes")
dir.create(folder, showWarnings = FALSE)
for (i in 1:3)
  saveRDS(list(value = 10 / i, converged = TRUE, iterations = 5L, argument = c(a = i)),
          file.path(folder, paste0("fit-", i, ".Rda")))
fits <- load.parlist(folder)
as.parframe(fits)
#>   index     value converged iterations a
#> 1     3  3.333333      TRUE          5 3
#> 2     2  5.000000      TRUE          5 2
#> 3     1 10.000000      TRUE          5 1
unlink(folder, recursive = TRUE)
```
