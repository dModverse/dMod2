# Confidence Intervals from Profiles

Reads confidence intervals off profile likelihoods. Where a profile does
not reach the threshold, the last part of it is extrapolated linearly.

## Usage

``` r
# S3 method for class 'parframe'
confint(
  object,
  parm = NULL,
  level = 0.95,
  ...,
  val.column = "data",
  method = c("chisq", "F"),
  n = NULL,
  p = NULL
)
```

## Arguments

- object:

  A `parframe` returned by
  [profile()](https://dmodverse.github.io/dMod2/reference/profile.md).

- parm:

  Names of the parameters to give intervals for. Default `NULL`, all
  profiled parameters.

- level:

  The confidence level. Default `0.95`.

- ...:

  Not used.

- val.column:

  Column of the profile the threshold applies to, for example `"value"`
  or `"data"`. Default `"data"`, or `"value"` for a profile without a
  `data` column.

- method:

  Calibration of the threshold, `"chisq"` (default) or `"F"`, see
  [`profileThreshold()`](https://dmodverse.github.io/dMod2/reference/profileThreshold.md).

- n:

  number of data points. Only for `method = "F"`.

- p:

  number of estimated mean parameters, that is without the error-model
  parameters. Only for `method = "F"`.
  [remlLeverage](https://dmodverse.github.io/dMod2/reference/remlLeverage.md)
  reports the effective count, which is what to pass when directions are
  non-identifiable.

## Value

A data frame with one row per parameter and columns `name`, `value` (the
value at the optimum), `lower` and `upper`. A bound is `NA` or infinite
when the profile on that side is too short or does not rise toward the
threshold.

## See also

[profile()](https://dmodverse.github.io/dMod2/reference/profile.md),
[`profileThreshold()`](https://dmodverse.github.io/dMod2/reference/profileThreshold.md)
