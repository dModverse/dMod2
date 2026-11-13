# Residuals Between Data and Model Prediction

Matches data to the prediction by time and name and computes the
residuals and their derivatives. Values below `lloq` are set to `lloq`.

## Usage

``` r
res(data, out, err = NULL)
```

## Arguments

- data:

  `data.frame` with columns `time`, `name`, `value`, `sigma` and `lloq`,
  one row per data point. Rows with `sigma = NA` take it from `err`.

- out:

  Prediction matrix, a
  [prdframe](https://dmodverse.github.io/dMod2/reference/prdframe.md):
  column `time`, then one column per observable. Optional attributes
  `"deriv"`, `[time, name, parameter]`, and `"deriv2"`,
  `[time, name, parameter, parameter]`.

- err:

  Optional error-model matrix, in the layout of `out`. Default `NULL`.

## Value

An
[`objframe()`](https://dmodverse.github.io/dMod2/reference/objframe.md)
with one row per data point and columns `time`, `name`, `value`,
`prediction`, `sigma`, `residual`, `weighted.residual`, `bloq` and
`weighted.0`. Attributes `"deriv"` and `"deriv.err"` are the derivatives
of the prediction and of `sigma`, `[data point, parameter]`; `"deriv2"`
and `"deriv2.err"`, `[data point, parameter, parameter]`, the second
derivatives when `out` or `err` have them. Each is `NULL` when missing
on the input.

## See also

[`objframe()`](https://dmodverse.github.io/dMod2/reference/objframe.md),
[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md)
