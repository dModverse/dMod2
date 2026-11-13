# Threshold for a Profile Likelihood

How far a profile has to rise for a parameter value to leave the
confidence region. `"chisq"` is the asymptotic \\\chi^2_1\\ quantile.
`"F"` is the finite-sample threshold \\n\log(1 +
F\_{1,\nu,1-\alpha}/\nu)\\, \\\nu = n - p\\, exact for Gaussian errors
with a single sigma profiled out and the generalization of the classical
t interval.

## Usage

``` r
profileThreshold(level = 0.95, method = c("chisq", "F"), n = NULL, p = NULL)
```

## Arguments

- level:

  the confidence level required.

- method:

  `"chisq"` (default) or `"F"`.

- n:

  number of data points. Only for `method = "F"`.

- p:

  number of estimated mean parameters, that is without the error-model
  parameters. Only for `method = "F"`.
  [remlLeverage](https://dmodverse.github.io/dMod2/reference/remlLeverage.md)
  reports the effective count, which is what to pass when directions are
  non-identifiable.

## Value

A single number.

## See also

[profile()](https://dmodverse.github.io/dMod2/reference/profile.md),
[`confint.parframe()`](https://dmodverse.github.io/dMod2/reference/confint.parframe.md)
