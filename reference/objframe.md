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
