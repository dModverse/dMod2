# Jobs of distributedComputing()

Lists the jobs submitted with
[`distributedComputing()`](https://dmodverse.github.io/dMod2/reference/distributedComputing.md)
that were not deleted with `purge(local = TRUE)`, from every working
directory. Their functions come with
`distributedComputing(NULL, jobname = , recover = TRUE)`.

## Usage

``` r
dMod2Jobs()
```

## Value

A data frame with one row per job: `jobname`, `on` (where it runs),
`tasks`, `jobid` (of SLURM), `submitted`, `wd` (the working directory it
was submitted from) and `local` (whether the local job folder is still
there).

## See also

[`distributedComputing()`](https://dmodverse.github.io/dMod2/reference/distributedComputing.md)

## Examples

``` r
dMod2Jobs()
#> [1] jobname   on        tasks     jobid     submitted wd        local    
#> <0 rows> (or 0-length row.names)
```
