# Prediction Frame

A prediction frame stores model predictions as a matrix with the column
`"time"` and one column per variable. The attributes are `"deriv"`, the
first-order sensitivities with respect to the outer parameters (see
[`P()`](https://dmodverse.github.io/dMod2/reference/P.md)), `"deriv2"`,
the second-order sensitivities, and `"parameters"`, the inner parameters
used for the prediction.

Prediction lists
([prdlist](https://dmodverse.github.io/dMod2/reference/prdlist.md))
returned by [`Xs()`](https://dmodverse.github.io/dMod2/reference/Xs.md),
[`Xd()`](https://dmodverse.github.io/dMod2/reference/Xd.md) or
[`Xf()`](https://dmodverse.github.io/dMod2/reference/Xf.md) consist of
prediction frames. A custom prediction function (see `P2X` in
[`prdfn()`](https://dmodverse.github.io/dMod2/reference/prdfn.md))
returns one.

## Usage

``` r
prdframe(prediction = NULL, deriv = NULL, deriv2 = NULL, parameters = NULL)
```

## Arguments

- prediction:

  numeric matrix of model predictions. Default `NULL` gives an empty
  matrix.

- deriv:

  3D numeric array `[time, variable, parameter]` of first-order
  sensitivities. Default `NULL`.

- deriv2:

  4D numeric array `[time, variable, parameter, parameter]` of
  second-order sensitivities, symmetric in the last two axes. Default
  `NULL`.

- parameters:

  named numeric vector of the inner parameters. Default `NULL`.

## Value

Object of class `prdframe`, a matrix with the attributes above.

## See also

[`getDerivs()`](https://dmodverse.github.io/dMod2/reference/getDerivs.md),
[`getDerivs2()`](https://dmodverse.github.io/dMod2/reference/getDerivs2.md)
