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

A parlist.

## See also

[`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md)
