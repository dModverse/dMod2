# Parameter List

A parameter list holds the results of several optimization runs, as
returned by
[`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md).
[`c()`](https://rdrr.io/r/base/c.html) concatenates parameter lists.

## Usage

``` r
as.parlist(x = NULL)

# S3 method for class 'parlist'
print(x, ...)

# S3 method for class 'parlist'
summary(object, ...)

# S3 method for class 'parlist'
c(...)

parlist(...)
```

## Arguments

- x:

  list of fit results, each as returned by
  [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md).
  Default `NULL` returns `NULL`.

- ...:

  for `parlist()`, fit results, one entry each; for
  [`c()`](https://rdrr.io/r/base/c.html), fit results or parameter
  lists, whose entries are concatenated; for
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html), not used.

- object:

  object of class `parlist`.

## Value

Object of class `parlist`.
[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) return their argument
invisibly.

## See also

[`load.parlist()`](https://dmodverse.github.io/dMod2/reference/load.parlist.md),
[`plot.parlist()`](https://dmodverse.github.io/dMod2/reference/plot.parlist.md),
[`as.parframe()`](https://dmodverse.github.io/dMod2/reference/as.parframe.md)

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
#> worst     3 440.78695      TRUE         16     fvalue    17      0         0
#>       nSwitch        a
#> best        0 0.989118
#> worst       0 2.363339
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

bestfit <- as.parvec(parframe)
plot((g*x)(times, bestfit), data)

# }
```
