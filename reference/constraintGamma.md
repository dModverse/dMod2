# Soft Gamma Constraint on Positive Parameters

Soft Gamma Constraint on Positive Parameters

## Usage

``` r
constraintGamma(shape, scale = 1, attrName = "prior", condition = NULL)
```

## Arguments

- shape:

  Named numeric vector of shape parameters. The names select the
  constrained parameters.

- scale:

  Numeric, scalar or named and aligned with `shape`, default `1`. The
  gamma scale, not the rate.

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

Computes the gamma prior on the `-2 log` scale, \$\$-2(a-1)\log p +
2p/s + 2\log\Gamma(a) + 2a\log s\$\$ with shape \\a\\ and scale \\s\\.
The value is `Inf` for a non-positive parameter.

## See also

[constraintExponential](https://dmodverse.github.io/dMod2/reference/constraintExponential.md),
[constraintChisq](https://dmodverse.github.io/dMod2/reference/constraintChisq.md)

## Examples

``` r
prior <- constraintGamma(shape = c(k1 = 3), scale = 5)
prior(pars = c(k1 = 5))$value
#> [1] 6.60517
```
