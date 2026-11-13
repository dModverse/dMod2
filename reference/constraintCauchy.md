# Soft Cauchy Constraint on Parameters

Soft Cauchy Constraint on Parameters

## Usage

``` r
constraintCauchy(mu, sigma = 1, attrName = "prior", condition = NULL)
```

## Arguments

- mu:

  Named numeric vector of prior locations. The names select the
  constrained parameters.

- sigma:

  Numeric, scalar or named and aligned with `mu`, default `1`. The
  Cauchy scale.

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

Computes the Cauchy prior on the `-2 log` scale, \$\$2\log(\pi\sigma) +
2\log(1 + ((p-\mu)/\sigma)^2).\$\$ The heavy tails make it the tolerant
alternative to
[constraintL2](https://dmodverse.github.io/dMod2/reference/constraintL2.md):
a parameter far from `mu` is pulled far more weakly than a Gaussian
would pull it.

## See also

[constraintL2](https://dmodverse.github.io/dMod2/reference/constraintL2.md),
[constraintL1](https://dmodverse.github.io/dMod2/reference/constraintL1.md)

## Examples

``` r
prior <- constraintCauchy(mu = c(k1 = 0), sigma = 2)
prior(pars = c(k1 = 1))$value
#> [1] 4.122041
```
