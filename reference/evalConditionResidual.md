# L2 Objective of One Condition

Computes the contribution of one condition to
[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md).

## Usage

``` r
evalConditionResidual(
  dataI,
  predictionI,
  errfn = NULL,
  fixed = NULL,
  cn = NULL,
  eCondNames = NULL,
  deriv = TRUE,
  deriv2 = FALSE,
  optBLOQ = c("M3", "M1", "M4NM", "M4BEAL"),
  df = Inf
)
```

## Arguments

- dataI:

  `data.frame` of one condition of a
  [datalist](https://dmodverse.github.io/dMod2/reference/datalist.md),
  with columns `name`, `time`, `value` and `sigma`.

- predictionI:

  [prdframe](https://dmodverse.github.io/dMod2/reference/prdframe.md) of
  that condition.

- errfn:

  Optional
  [obsfn](https://dmodverse.github.io/dMod2/reference/obsfn.md), the
  error model. Default `NULL`.

- fixed:

  Optional named numeric vector of fixed parameters, passed to `errfn`.
  Default `NULL`.

- cn:

  Character, the condition name. Required when `errfn` is set.

- eCondNames:

  Character vector, the conditions `errfn` is defined for. Default
  `NULL`, all.

- deriv:

  Logical. `TRUE` (default) returns the gradient and the Gauss-Newton
  Hessian.

- deriv2:

  Logical. `TRUE` returns the exact Hessian, from the second derivatives
  of the prediction. Default `FALSE`.

- optBLOQ:

  Character, the treatment of data below the limit of quantification, as
  in
  [`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md).
  Default `"M3"`.

- df:

  Degrees of freedom of a Student-t likelihood as in
  [`normT()`](https://dmodverse.github.io/dMod2/reference/normT.md);
  `Inf` (default) is the normal one.

## Value

An [objlist](https://dmodverse.github.io/dMod2/reference/objlist.md).

## See also

[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md)
