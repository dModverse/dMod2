# L2 Objective of a Validation Data Point

Builds an objective function that compares the prediction of one
quantity at one time with a data value held by a parameter.

## Usage

``` r
datapointL2(
  name,
  time,
  parameter,
  sigma = 1,
  attrName = "validation",
  condition
)
```

## Arguments

- name:

  Character, the name of the predicted quantity, a state or an
  observable.

- time:

  Numeric of length 1, the time of the data point.

- parameter:

  Character of length 1, the name of the parameter holding the data
  value.

- sigma:

  Numeric of length 1, the standard deviation of the data point, default
  `1`.

- attrName:

  Character. Name of the attribute holding the value, default
  `"validation"`. The value is also returned as `chi2` contribution, see
  `attrName` in
  [`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md).

- condition:

  Character, the condition of the prediction. No default.

## Value

An objective function of class `objfn`, called as `obj(pars, ...)`, see
section Calling a dMod function. It returns an
[objlist](https://dmodverse.github.io/dMod2/reference/objlist.md) with
the prediction as attribute `prediction`.

## Details

The value is \$\$\left(\frac{x(t) - v}{\sigma}\right)^2\$\$ with the
prediction \\x(t)\\ and the value of `parameter`, \\v\\. The prediction
is read from `env`, where
[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md)
stores it, so the objective is evaluated as summand of an objective that
predicts, or with the `env` of such an objective. `time`, `sigma` and
`attrName` can be changed with
[`controls()`](https://dmodverse.github.io/dMod2/reference/controls.md).

## Calling a dMod function

The leading arguments depend on the class and are matched by position or
by name:

- prediction function: `x(times, pars, ...)`, returning a
  [prdlist](https://dmodverse.github.io/dMod2/reference/prdlist.md);

- observation function: `g(out, pars, ...)` with a
  [prdlist](https://dmodverse.github.io/dMod2/reference/prdlist.md)
  `out`, returning a
  [prdlist](https://dmodverse.github.io/dMod2/reference/prdlist.md);

- parameter transformation: `p(pars, ...)`, returning a list with one
  [parvec](https://dmodverse.github.io/dMod2/reference/parvec.md) per
  condition;

- objective function: `obj(pars, ...)`, returning an
  [objlist](https://dmodverse.github.io/dMod2/reference/objlist.md).

`pars` is a named numeric vector. If it has a `deriv` attribute, as the
output of a parameter transformation does, the derivatives are taken
with respect to the parameters behind it. A composition is called
according to the resulting class: `g * x * p` is a prediction function,
`(g * x * p)(times, pars)`.

The other arguments are named:

- `fixed`: named numeric vector of parameters that enter the evaluation
  without derivatives. Default `NULL`.

- `deriv`: `TRUE` (default) returns first derivatives with respect to
  `pars`. Every function of the chain has to be built with
  `derivMode = "forward"`, see
  [`odemodel()`](https://dmodverse.github.io/dMod2/reference/odemodel.md),
  [`Y()`](https://dmodverse.github.io/dMod2/reference/Y.md) and
  [`P()`](https://dmodverse.github.io/dMod2/reference/P.md).

- `deriv2`: `TRUE` returns second derivatives as well. Default `FALSE`.
  Every function of the chain has to be built with
  `derivMode = "forward-forward"`.

- `hessian`: objective functions only, see the table below. Default
  `NULL`.

- `sweep`: objective functions only. `"forward"` (default) propagates
  the sensitivities through the chain; `"reverse"` computes the gradient
  in one backward pass, whose cost does not grow with the number of
  parameters. `"reverse"` needs every function of the chain built with
  `derivMode = "reverse"`, and with `"forward-reverse"` for
  `deriv2 = TRUE`.

- `conditions`: character vector, the conditions to evaluate. Default
  all conditions of the function.

- `env`: environment for intermediate results. Default `NULL`, a new
  one. An objective function stores the prediction there as `prediction`
  and returns the environment as attribute `env`.

- `cores`: number of threads for the ODE solves and the residual kernel.
  Default `getOption("dMod2.cores", 1L)`.

Prediction, observation and parameter functions accept `hessian` and
`sweep` and ignore them. With `deriv = TRUE` an objective function
returns this Hessian:

|             |                |                      |                      |
|-------------|----------------|----------------------|----------------------|
| `sweep`     | `hessian`      | `deriv2 = FALSE`     | `deriv2 = TRUE`      |
| `"forward"` | `NULL`, `TRUE` | Gauss-Newton         | exact                |
| `"forward"` | `FALSE`        | none                 | none, with a warning |
| `"reverse"` | `NULL`         | none                 | exact                |
| `"reverse"` | `TRUE`         | none, with a warning | exact                |
| `"reverse"` | `FALSE`        | none                 | none, with a warning |

`hessian = "sparse"` returns the Gauss-Newton Hessian of the forward
sweep as a sparse `Matrix`, assembled from the blocks of the conditions;
[`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md)
factorizes it where the argument `sparse` takes the sparse path, and
[`vcov()`](https://dmodverse.github.io/dMod2/reference/vcov.objlist.md)
solves it for the columns of `which`.

`hessian = "expected"` returns, in place of the Gauss-Newton Hessian of
the forward sweep, the expected (Fisher) information of the likelihood:
`2 / sigma^2` in the prediction and `4 / sigma^2` in sigma per data
point, positive semidefinite also where an error model depends on the
parameters. The Gauss-Newton Hessian keeps the exact curvature of the
loss in the prediction and sigma, which is indefinite there and shows
the optimizer a saddle.

With `deriv = FALSE` it returns the value alone.

## See also

[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md),
[`constraintL2()`](https://dmodverse.github.io/dMod2/reference/constraintL2.md)

## Examples

``` r
## A prediction without an ODE: the values of A at the times 0 to 5 are parameters
grid <- data.frame(name = "A", time = 0:5, row.names = paste0("p", 0:5))
x <- Xd(grid, condition = "C1")
pars <- setNames(c(1, 0.8, 0.6, 0.5, 0.4, 0.3), row.names(grid))
data <- as.datalist(data.frame(name = "A", time = c(0, 2, 4), value = c(1, 0.6, 0.4),
                               sigma = 0.1, condition = "C1"))

## The data term also predicts at t = 2.5, where the data point v reads the prediction
obj <- normL2(data, x, times = 2.5) +
  datapointL2(name = "A", time = 2.5, parameter = "v", sigma = 0.05, condition = "C1")
obj(c(pars, v = 0.6))
#> value
#> ==================
#> -7.301879
#> gradient[1:7] (full length = 7)
#> ==================
#>  p0  p1  p2  p3  p4  p5   v 
#>   0   0 -20 -20   0   0  40 
#> 
#> hessian[1:6,1:6]
#> ==================
#>     p0 p1  p2  p3  p4 p5
#> p0 200  0   0   0   0  0
#> p1   0  0   0   0   0  0
#> p2   0  0 400 200   0  0
#> p3   0  0 200 200   0  0
#> p4   0  0   0   0 200  0
#> p5   0  0   0   0   0  0
#> 
#> 
#> attributes
#> ==================
#>  data             -8.301879
#>  validation       1
#>  prediction       0.55
#>  chi2_data        0
#>  chi2_validation  1
#>  env              <environment>
```
