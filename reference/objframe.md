# Objective Frame

An objective frame stores residuals and their derivatives with respect
to parameters. It is typically created by
[res](https://dmodverse.github.io/dMod2/reference/res.md) and used
internally in objective functions.

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

  data.table produced by
  [res](https://dmodverse.github.io/dMod2/reference/res.md)

- deriv:

  numeric matrix of first-order derivatives of residuals (Jacobian)

- deriv.err:

  numeric matrix of first-order derivatives of the error model

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
