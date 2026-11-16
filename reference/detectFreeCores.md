# Detect Number of Free Cores

Estimates the free cores as the number of cores minus the 1-minute load
average, at least 1. With `machine = NULL` inside a SLURM job, they are
the cores the job was given (`SLURM_CPUS_ON_NODE`) and the job processes
are bound to. Works on Linux and macOS, locally or on remote machines
via `ssh`. On the local Windows machine the load is the processor load
since the previous call, over at least 0.25 s and at most a minute.
Under `R CMD check` with `_R_CHECK_LIMIT_CORES_` set, the result is at
most 2.

## Usage

``` r
detectFreeCores(machine = NULL)
```

## Arguments

- machine:

  Character vector of machines: `"localhost"` is the local machine, any
  other entry an `ssh` destination such as `"host"` or `"user@host"`.
  `NULL` (default) for the local machine.

## Value

Numeric vector of free cores, one per machine, named by `machine` when
it is given, with attributes `"ncores"` (number of cores) and `"used"`
(load average).

## Examples

``` r
detectFreeCores()
#> [1] 3
#> attr(,"ncores")
#> [1] 4
#> attr(,"used")
#> [1] 1.47
```
