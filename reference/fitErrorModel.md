# Fit an Error Model to Replicate Data

Fits the variance of replicate measurements as a function of their mean
by maximum likelihood, assuming that the sample variance follows a
scaled \\\chi^2\\ distribution with \\n - 1\\ degrees of freedom.
Requires the package optimx.

## Usage

``` r
fitErrorModel(
  data,
  factors,
  errorModel = "exp(s0)+exp(srel)*x^2",
  par = c(s0 = 1, srel = 0.1),
  lower = NULL,
  upper = NULL,
  plotting = FALSE,
  blather = FALSE,
  ...
)
```

## Arguments

- data:

  A data frame as returned by
  [`reduceReplicates()`](https://dmodverse.github.io/dMod2/reference/reduceReplicates.md),
  with columns `value` (mean of the replicates), `sigma` (standard error
  of the mean) and `n` (number of replicates).

- factors:

  Character vector, the columns of `data` whose value combinations
  define groups. The model is fitted separately per group.

- errorModel:

  Character, the variance of a single measurement as an expression in
  the mean `x` and the parameters. Default `"exp(s0)+exp(srel)*x^2"`.

- par:

  Named numeric vector of starting values for the parameters of
  `errorModel`. Default `c(s0 = 1, srel = 0.1)`.

- lower, upper:

  Bounds on the parameters, named or in the order of `par`. Default
  `NULL`, unbounded.

- plotting:

  Logical. If `TRUE`, plot the variances, the fitted model with 68 and
  95 percent bands per group. Default `FALSE`.

- blather:

  Logical. If `TRUE`, return the extended data frame described under
  Value. Default `FALSE`.

- ...:

  Further arguments passed to
  [`optimx::optimr()`](https://rdrr.io/pkg/optimx/man/optimr.html),
  which runs `"L-BFGS-B"`.

## Value

`data` with `sigma` replaced by the standard error of the mean that the
fitted model gives.

With `blather = TRUE`, a data frame with the columns of `data`, `sigma`
replaced as above, one column per parameter of `errorModel` with the
fitted value, the variance bands `cbLower68`, `cbUpper68`, `cbLower95`
and `cbUpper95`, the group label `condidnt` and the input `sigma` as
`sigmaLS`. The attribute `"errorModel"` holds `errorModel`.

## See also

[`reduceReplicates()`](https://dmodverse.github.io/dMod2/reference/reduceReplicates.md)

## Examples

``` r
set.seed(1)
mu <- rep(c(0.1, 1, 10, 100), each = 4)
data <- data.frame(name = "y", time = rep(1:4, each = 4), condition = "ctrl",
                   value = mu + rnorm(16, sd = sqrt(0.01 + 0.04 * mu^2)))
reduced <- reduceReplicates(data)
#> Dropped columns: 
fitErrorModel(reduced, factors = "name", par = c(s0 = -4, srel = -3))
#>    time      value      sigma n name condition
#> 1     1  0.1080779 0.05321124 4    y      ctrl
#> 5     2  1.0410762 0.12172064 4    y      ctrl
#> 9     3 11.0873653 1.17336938 4    y      ctrl
#> 13    4 91.2201745 9.64444405 4    y      ctrl
```
