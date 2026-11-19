# Objective Frame

An objective frame stores residuals and their derivatives with respect
to parameters, as
[`res()`](https://dmodverse.github.io/dMod2/reference/res.md) returns
them.

## Usage

``` r
objframe(
  mydata,
  deriv = NULL,
  deriv.err = NULL,
  deriv2 = NULL,
  deriv2.err = NULL
)
```

## Arguments

- mydata:

  Data frame with columns `time`, `name`, `value`, `prediction`,
  `sigma`, `residual`, `weighted.residual`, `bloq` and `weighted.0`;
  further columns are dropped.

- deriv:

  numeric matrix `[n_residuals, p]` of first-order derivatives of the
  residuals

- deriv.err:

  numeric matrix `[n_residuals, p]` of first-order derivatives of the
  error model

- deriv2:

  numeric 3D array `[n_residuals, p, p]` of second-order derivatives of
  residuals with respect to parameters. Optional.

- deriv2.err:

  numeric 3D array `[n_residuals, p, p]` of second-order derivatives of
  the error model. Optional.

## Value

An object of class `"objframe"` (data.table) with attributes `"deriv"`
and `"deriv.err"`. These arrays have the same parameter axes as those
returned by
[prdframe](https://dmodverse.github.io/dMod2/reference/prdframe.md) and
[res](https://dmodverse.github.io/dMod2/reference/res.md). When
`deriv2`/`deriv2.err` are supplied, the corresponding 3D arrays are
attached as `"deriv2"` / `"deriv2.err"`.
