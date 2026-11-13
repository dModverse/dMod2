# Diagnostic plots of an alignment

`plotAligned()` shows the aligned values with their standard errors,
`plotScaled()` the measured readings mapped onto the values next to them
(a censored reading has no value and is left out),
[`plotResiduals()`](https://dmodverse.github.io/dMod2/reference/plotResiduals.md)
the standardized residuals against the prediction.

## Usage

``` r
plotAligned(x, ...)

plotScaled(x, ...)

# S3 method for class 'alignment'
plotResiduals(parframe, ...)
```

## Arguments

- x:

  An `alignment` from
  [`align()`](https://dmodverse.github.io/dMod2/reference/align.md).

- ...:

  Ignored.

- parframe:

  An `alignment`, for
  [`plotResiduals()`](https://dmodverse.github.io/dMod2/reference/plotResiduals.md).

## Value

A ggplot.
