# Soft L1 Constraint on Parameters

Builds a Laplace prior on parameters as an objective function.

## Usage

``` r
constraintL1(mu, ...)

# Default S3 method
constraintL1(mu, sigma = 1, attrName = "prior", condition = NULL, ...)
```

## Arguments

- mu:

  Named numeric vector of prior locations. The names select the
  constrained parameters.

- ...:

  Passed to the method.

- sigma:

  Numeric, scalar or named and aligned with `mu`, default `1`. The
  Laplace scale, the reciprocal penalty strength.

- attrName:

  Character. Name of the attribute holding the constraint value, default
  `"prior"`.

- condition:

  Character vector, the conditions of a sum of objectives in which the
  term is evaluated. `NULL` (default) evaluates it in every one.

## Value

An objective function of class `objfn`, called as `obj(pars, ...)`, see
[dModfn](https://dmodverse.github.io/dMod2/reference/dModfn.md). The
objective value is the full `-2 log` density, while
[`constraintL2()`](https://dmodverse.github.io/dMod2/reference/constraintL2.md)
returns the penalty form without the constant.

## Details

Computes the Laplace prior on the `-2 log` scale, \$\$2\log(2\sigma) +
2\|p-\mu\|/\sigma,\$\$ the L1 counterpart of
[`constraintL2()`](https://dmodverse.github.io/dMod2/reference/constraintL2.md).
At \\p = \mu\\ the term is not differentiable and the gradient is
returned as 0.

## See also

[`constraintL2()`](https://dmodverse.github.io/dMod2/reference/constraintL2.md)

## Examples

``` r
prior <- constraintL1(mu = c(k1 = 0), sigma = 2)
prior(pars = c(k1 = 1))$value
#> [1] 3.772589
```
