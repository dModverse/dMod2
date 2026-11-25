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

## Examples

``` r
## Prediction with the derivative of the observable y in the parameter k
times <- 0:2
J <- array(-times * exp(-times), c(3, 1, 1), list(NULL, "y", "k"))
out <- prdframe(cbind(time = times, y = exp(-times)), deriv = J, parameters = c(k = 1))

data <- data.frame(time = times, name = "y", value = c(1.05, 0.35, 0.15),
                   sigma = 0.1, lloq = -Inf)
r <- res(data, out)
r
#>     time   name value prediction sigma    residual weighted.residual   bloq
#>    <int> <char> <num>      <num> <num>       <num>             <num> <lgcl>
#> 1:     0      y  1.05  1.0000000   0.1 -0.05000000        -0.5000000  FALSE
#> 2:     1      y  0.35  0.3678794   0.1  0.01787944         0.1787944  FALSE
#> 3:     2      y  0.15  0.1353353   0.1 -0.01466472        -0.1466472  FALSE
#>    weighted.0
#>         <num>
#> 1:  10.000000
#> 2:   3.678794
#> 3:   1.353353
attr(r, "deriv")
#>               k
#> [1,]  0.0000000
#> [2,] -0.3678794
#> [3,] -0.2706706
```
