# Non-Linear Optimization, Multi-Start

Runs [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md)
from several starting points: `center` plus random offsets drawn by
`samplefun`, or the rows of a parframe.

## Usage

``` r
mstrust(
  objfun,
  center,
  rinit = 0.1,
  rmax = 10,
  fits = 20,
  cores = 1,
  optmethod = "trust",
  samplefun = "rnorm",
  resultPath = ".",
  name = "mstrust",
  stats = FALSE,
  output = FALSE,
  cautiousMode = FALSE,
  startFromCenter = FALSE,
  retry = TRUE,
  nTries = 10L,
  ...
)
```

## Arguments

- objfun:

  Objective function, see
  [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md).

- center:

  Named numeric parameter center, or a parframe whose rows are used as
  starting points; `fits` is then the number of rows. The random starts,
  those of retries included, are drawn before the fits are spread over
  `cores`, so [`set.seed()`](https://rdrr.io/r/base/Random.html)
  reproduces a run.

- rinit, rmax:

  Initial and maximum trust-region radius, see
  [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md).
  Default `0.1` and `10`.

- fits:

  Number of fits. Default `20`.

- cores:

  Number of parallel fits, or `c(fits = , conditions = )` to also give
  each fit a condition axis, which selects a PSOCK backend. Keep the
  product below the number of cores. Default `1`.

- optmethod:

  Optimizer, a function or a function name, called via
  `do.call(optmethod, ...)`. Default `"trust"`. Arguments in `...` are
  routed by the formals of this optimizer.

- samplefun:

  Sampler for the random offsets, a function or a function name, with an
  argument `n`, the number of parameters. Default `"rnorm"`.

- resultPath:

  Folder under which files are written. Default `"."`, the working
  directory.

- name:

  Folder name under `resultPath`. Default `"mstrust"`.

- stats:

  Logical, print the summary to the console. Default `FALSE`.

- output:

  Logical, write every fit and the final parlist to disk. Default
  `FALSE`.

- cautiousMode:

  Logical. With `output = TRUE`, also write every fit deparsed to an
  `.R` file in `interRes/`, which is kept after the run while the `.Rda`
  files in it are removed. Default `FALSE`.

- startFromCenter:

  Logical. If `TRUE`, the first fit starts at `center` without offset.
  Default `FALSE`.

- retry:

  Logical. If `TRUE` (default), a fit that fails with an error is
  repeated from a new random start, up to `nTries` attempts. Not used
  when `center` is a parframe.

- nTries:

  Maximum number of attempts per fit, including the first. Default `10`.

- ...:

  Arguments passed by name to `optmethod`, `samplefun` or `objfun`.
  Names matching neither the optimizer nor the sampler go to `objfun`.

## Value

A parlist with one entry per fit, failed fits included. Each entry is
the result of the optimizer without `hessian`, plus `parinit`, the
starting point, and `elapsed`, the wall time in seconds including
retries; a fit that failed is a list with element `error`. For the
covariance at a fit, use `vcov(objfun(fit$argument))`. The attribute
`"run"` holds the start time, total wall time, number of fits, cores,
optimizer and sampler with their arguments, and the versions of R, dMod2
and cppDE.
[`as.parframe()`](https://dmodverse.github.io/dMod2/reference/as.parframe.md)
turns the parlist into a table.

## Details

Nothing is written to disk unless `output`, `cautiousMode` or a
`traceFile` argument for
[`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md) asks
for it. Files then go to `<resultPath>/<name>/trial-<k>-<date>/`: the
per-fit files in the subfolder `interRes/`, which is removed on
completion unless `cautiousMode = TRUE`, the final `parameterList.Rda`
(with `output`) and the log `mstrust.log`.
[`load.parlist()`](https://dmodverse.github.io/dMod2/reference/load.parlist.md)
reads the fits from `interRes/` of a run that is still going or was
aborted.

## See also

[`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md),
[`msParframe()`](https://dmodverse.github.io/dMod2/reference/msParframe.md),
[`as.parframe()`](https://dmodverse.github.io/dMod2/reference/as.parframe.md),
[`load.parlist()`](https://dmodverse.github.io/dMod2/reference/load.parlist.md);
the guide chapters
[Optimization](https://dmodverse.github.io/dMod2/articles/guide-optimization.html)
and [Remote
computing](https://dmodverse.github.io/dMod2/articles/guide-remote-computing.html),
which explains `cores`.

## Examples

``` r
rosenbrock <- function(x, ...) {
  a <- x[["a"]]
  b <- x[["b"]]
  list(value = 100 * (b - a^2)^2 + (1 - a)^2,
       gradient = c(a = -400 * a * (b - a^2) - 2 * (1 - a),
                    b = 200 * (b - a^2)),
       hessian = matrix(c(1200 * a^2 - 400 * b + 2, -400 * a, -400 * a, 200),
                        2, 2, dimnames = list(c("a", "b"), c("a", "b"))))
}

# Random starts around center, sd goes to rnorm()
fits <- mstrust(rosenbrock, center = c(a = 0, b = 0), fits = 5, sd = 2,
                iterlim = 200)
as.parframe(fits)
#>   index        value converged iterations stopReason neval qnEval qnSkipped
#> 1     3 5.703357e-18      TRUE         19     fvalue    20      0         0
#> 2     1 9.263308e-17      TRUE         26     fvalue    27      0         0
#> 3     2 4.196884e-16      TRUE         19     fvalue    20      0         0
#> 4     4 1.311916e-15      TRUE         31     fvalue    32      0         0
#> 5     5 5.083307e-12      TRUE         34     fvalue    35      0         0
#>   nSwitch elapsed         a         b
#> 1       0   0.001 1.0000000 1.0000000
#> 2       0   0.024 1.0000000 1.0000000
#> 3       0   0.001 1.0000000 1.0000000
#> 4       0   0.001 1.0000000 0.9999999
#> 5       0   0.002 0.9999978 0.9999955

# Reproducible starts from a parframe
starts <- msParframe(c(a = 0, b = 0), n = 5, sd = 2)
fits <- mstrust(rosenbrock, center = starts, iterlim = 200)
as.parframe(fits)
#>   index        value converged iterations stopReason neval qnEval qnSkipped
#> 1     2 1.729742e-18      TRUE          7     fvalue     8      0         0
#> 2     4 3.350785e-15      TRUE         18     fvalue    19      0         0
#> 3     5 1.829942e-14      TRUE         12     fvalue    13      0         0
#> 4     1 2.178106e-14      TRUE         15     fvalue    16      0         0
#> 5     3 2.523921e-14      TRUE         20     fvalue    21      0         0
#>   nSwitch elapsed         a         b
#> 1       0   0.001 1.0000000 1.0000000
#> 2       0   0.001 1.0000000 1.0000000
#> 3       0   0.001 0.9999999 0.9999997
#> 4       0   0.001 0.9999999 0.9999998
#> 5       0   0.000 0.9999999 0.9999997
```
