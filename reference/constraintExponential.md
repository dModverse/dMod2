# Soft Exponential Constraint on Positive Parameters

Soft Exponential Constraint on Positive Parameters

## Usage

``` r
constraintExponential(scale, attrName = "prior", condition = NULL)
```

## Arguments

- scale:

  Named numeric vector of scale parameters, not rates. The names select
  the constrained parameters.

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

Computes the exponential prior on the `-2 log` scale, \$\$2p/s + 2\log
s.\$\$ It is the gamma prior at shape 1 and pulls a parameter toward 0
with a constant force, which makes it the smooth one-sided counterpart
of an L1 penalty at `mu = 0`. The value is `Inf` for a negative
parameter.

## See also

[constraintGamma](https://dmodverse.github.io/dMod2/reference/constraintGamma.md),
[constraintL1](https://dmodverse.github.io/dMod2/reference/constraintL1.md)

## Examples

``` r
prior <- constraintExponential(scale = c(k1 = 3))
prior(pars = c(k1 = 5))$value
#> [1] 5.530558
```
