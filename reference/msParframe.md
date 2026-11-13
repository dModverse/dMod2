# Reproducible Random Starting Points

Draws a parframe of starting points for the `center` argument of
[`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md)
from the seed `seed`. The state of the global random number generator is
left unchanged.

## Usage

``` r
msParframe(
  pars,
  n = 20,
  seed = 12345,
  samplefun = stats::rnorm,
  keepfirst = TRUE,
  ...
)
```

## Arguments

- pars:

  Named numeric vector. If `samplefun` has an argument `mean`, the draws
  are centered on `pars`.

- n:

  Integer, number of rows. Default `20`.

- seed:

  Seed for the random number generator. Default `12345`.

- samplefun:

  Random number generator such as
  [`rnorm()`](https://rdrr.io/r/stats/Normal.html) or
  [`runif()`](https://rdrr.io/r/stats/Uniform.html). Default
  [`stats::rnorm()`](https://rdrr.io/r/stats/Normal.html).

- keepfirst:

  Logical. If `TRUE` (default), the first row is `pars`.

- ...:

  Arguments passed to `samplefun`.

## Value

A parframe without metanames, one column per parameter.

## See also

[`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md),
[`parframe()`](https://dmodverse.github.io/dMod2/reference/parframe.md)

## Examples

``` r
# Five starts around pars, the first one at pars
msParframe(c(k1 = 1, k2 = 10), n = 5)
#>          k1        k2
#> 1 1.0000000 10.000000
#> 2 1.5855288 10.709466
#> 3 0.8906967  9.546503
#> 4 1.6058875  8.182044
#> 5 1.6300986  9.723816

# A standard deviation per parameter
msParframe(c(k1 = 1, k2 = 10), n = 5, sd = c(0.1, 2))
#>          k1        k2
#> 1 1.0000000 10.000000
#> 2 1.0585529 11.418932
#> 3 0.9890697  9.093006
#> 4 1.0605887  6.364088
#> 5 1.0630099  9.447632
```
