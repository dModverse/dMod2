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

## Examples

``` r
## Prediction frames of a decay at two rates, without derivatives
times <- 0:4
frame1 <- prdframe(cbind(time = times, A = exp(-times)), parameters = c(k = 1))
frame2 <- prdframe(cbind(time = times, A = exp(-2 * times)), parameters = c(k = 2))
frame1
#>      time          A
#> [1,]    0 1.00000000
#> [2,]    1 0.36787944
#> [3,]    2 0.13533528
#> [4,]    3 0.04978707
#> [5,]    4 0.01831564
#> 
#> The prediction contains 1st-order derivatives: no
#> The prediction contains 2nd-order derivatives: no

## One frame per condition, in long format
pred <- prdlist(slow = frame1, fast = frame2)
as.data.frame(pred)
#>    time name        value condition sigma
#> 1     0    A 1.0000000000      slow   NaN
#> 2     1    A 0.3678794412      slow   NaN
#> 3     2    A 0.1353352832      slow   NaN
#> 4     3    A 0.0497870684      slow   NaN
#> 5     4    A 0.0183156389      slow   NaN
#> 6     0    A 1.0000000000      fast   NaN
#> 7     1    A 0.1353352832      fast   NaN
#> 8     2    A 0.0183156389      fast   NaN
#> 9     3    A 0.0024787522      fast   NaN
#> 10    4    A 0.0003354626      fast   NaN
```
