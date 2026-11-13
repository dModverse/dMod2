# Soft Rayleigh Constraint on Positive Parameters

Soft Rayleigh Constraint on Positive Parameters

## Usage

``` r
constraintRayleigh(sigma, attrName = "prior", condition = NULL)
```

## Arguments

- sigma:

  Named numeric vector of scale parameters. The names select the
  constrained parameters.

- attrName:

  Character. Name of the attribute holding the constraint value, default
  `"prior"`.

- condition:

  Character vector, the conditions of a sum of objectives in which the
  term is evaluated. `NULL` (default) evaluates it in every one.

## Value

An objective function of class `objfn`, called as `obj(pars, ...)`, see
[dModfn](https://dmodverse.github.io/dMod2/reference/dModfn.md). The
objective value is the full `-2 log` density.

## Details

Computes the Rayleigh prior on the `-2 log` scale, \$\$-2\log p +
4\log\sigma + p^2/\sigma^2.\$\$ Unlike the exponential it vanishes at 0,
so it keeps a parameter away from both 0 and large values. The value is
`Inf` for a non-positive parameter.

## See also

[constraintGamma](https://dmodverse.github.io/dMod2/reference/constraintGamma.md),
[constraintExponential](https://dmodverse.github.io/dMod2/reference/constraintExponential.md)

## Examples

``` r
prior <- constraintRayleigh(sigma = c(k1 = 3))
prior(pars = c(k1 = 5))$value
#> [1] 3.953351
```
