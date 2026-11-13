# Soft Chi-Squared Constraint on Positive Parameters

Soft Chi-Squared Constraint on Positive Parameters

## Usage

``` r
constraintChisq(df, attrName = "prior", condition = NULL)
```

## Arguments

- df:

  Named numeric vector of degrees of freedom. The names select the
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

Computes the chi-squared prior on the `-2 log` scale, \$\$-(k-2)\log p +
p + k\log 2 + 2\log\Gamma(k/2)\$\$ with `k` degrees of freedom. The
value is `Inf` for a non-positive parameter.

## See also

[constraintGamma](https://dmodverse.github.io/dMod2/reference/constraintGamma.md)

## Examples

``` r
prior <- constraintChisq(df = c(k1 = 4))
prior(pars = c(k1 = 5))$value
#> [1] 4.553713
```
