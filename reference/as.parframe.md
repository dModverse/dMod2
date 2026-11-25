# Coerce to a Parameter Frame

Coerce to a Parameter Frame

## Usage

``` r
# S3 method for class 'parframe'
as.parframe(x, ...)

# S3 method for class 'parlist'
as.parframe(x, sortBy = "value", ...)

as.parframe(x, ...)
```

## Arguments

- x:

  object to be coerced, e.g. a `parlist`.

- ...:

  not used.

- sortBy:

  character, the column by which the returned parameter frame is sorted.
  Default `"value"`.

## Value

Object of class
[parframe](https://dmodverse.github.io/dMod2/reference/parframe.md). For
a `parlist`, the meta-information columns are `index`, `value`,
`converged` and `iterations`, plus the termination and evaluation counts
and the wall time `elapsed` the fits report. Fits of
[`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md) with
several `phases` add one column `iterations.<phase>` per phase, the
iterations summed over the runs of that phase. The attribute `"run"` of
the parlist is kept.

## See also

[`parframe()`](https://dmodverse.github.io/dMod2/reference/parframe.md),
[`as.parvec.parframe()`](https://dmodverse.github.io/dMod2/reference/as.parvec.parframe.md),
and
[`parlist()`](https://dmodverse.github.io/dMod2/reference/parlist.md)
for an example.

## Examples

``` r
# \donttest{
## Observation of a regression curve
regfn <- c(y = "sin(a*time)")
g <- Y(regfn, parameters = "a", compile = TRUE, modelname = "parlist_obs",
       outdir = tempdir())
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC -fopenmp]
x <- Xt(condition = "C1")

## Data
set.seed(1)
data <- datalist(
  C1 = data.frame(
    name = "y",
    time = 1:5,
    value = sin(1:5) + rnorm(5, 0, .1),
    sigma = .1
  )
)

pars <- c(a = 1)
times <- seq(0, 5, .1)
plot((g*x)(times, pars), data)


## Fits from random starts, stored in a parlist
obj <- normL2(data, g*x)
out <- as.parlist(lapply(1:5, function(i) {
  trust(obj, pars + rnorm(length(pars), 0, 1), rinit = 1, rmax = 10)
}))
summary(out)
#> Best and worst fit
#>       index     value converged iterations stopReason neval qnEval qnSkipped
#> best      2 -10.24463      TRUE          8     fvalue     9      0         0
#> worst     3 440.78695      TRUE         16     fvalue    17      0         4
#>       nSwitch iterations.gn iterations.bfgs        a
#> best        0             8               0 0.989118
#> worst       0            16               0 2.363339
#> 
#> Statistics of fit outcome
#> Fits aborted:       0
#> Fits not converged: 0
#> Fits converged:     5
#> Fits total:         5
#> 
#> 
#> Termination reason
#> fvalue              5

## Parameter frame and best fit
parframe <- as.parframe(out)
plotValues(parframe)

unique(parframe)
#>   index     value converged iterations stopReason neval qnEval qnSkipped
#> 1     2 -10.24463      TRUE          8     fvalue     9      0         0
#> 2     1 223.37336      TRUE          6     fvalue     7      0         0
#> 3     3 440.78695      TRUE         16     fvalue    17      0         4
#>   nSwitch iterations.gn iterations.bfgs           a
#> 1       0             8               0  0.98911802
#> 2       0             6               0 -0.07418182
#> 3       0            16               0  2.36333881
bestfit <- as.parvec(parframe)
plot((g*x)(times, bestfit), data)

# }
```
