# Soft L2 Constraint on Parameters

Builds a Gaussian prior on parameters as an objective function, to be
added to a data term such as
[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md).

## Usage

``` r
constraintL2(mu, ...)

# Default S3 method
constraintL2(mu, sigma = 1, attrName = "prior", condition = NULL, ...)
```

## Arguments

- mu:

  Named numeric vector of prior means. The names select the constrained
  parameters.

- ...:

  Passed to the method.

- sigma:

  Numeric standard deviations, scalar or named and aligned with `mu`,
  default `1`; or a character vector of parameter names that estimate
  them on log scale, \\\sigma = \exp(s)\\. One kind for all entries.

- attrName:

  Character. Name of the attribute holding the constraint value, default
  `"prior"`.

- condition:

  Character vector, the conditions of a sum of objectives in which the
  term is evaluated. `NULL` (default) evaluates it in every one.

## Value

An objective function of class `objfn`, called as `obj(pars, ...)`, see
section Calling a dMod function. The objective value is the penalty
form, without the constant of the Gaussian density; the other
`constraint*()` functions return the full `-2 log` density.

## Details

The value is the sum over the constrained parameters of \$\$(p-\mu)^2 /
\sigma^2,\$\$ plus \\2\log\sigma\\ per parameter if `sigma` is
estimated.

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
  `pars`.

- `deriv2`: `TRUE` returns second derivatives as well. Default `FALSE`.
  Every function of the chain has to be built for them, see the
  arguments `deriv2` and `derivMode` of
  [`odemodel()`](https://dmodverse.github.io/dMod2/reference/odemodel.md),
  [`Y()`](https://dmodverse.github.io/dMod2/reference/Y.md) and
  [`P()`](https://dmodverse.github.io/dMod2/reference/P.md).

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

With `deriv = FALSE` it returns the value alone.

## See also

[`constraintL1()`](https://dmodverse.github.io/dMod2/reference/constraintL1.md),
[`constraintCauchy()`](https://dmodverse.github.io/dMod2/reference/constraintCauchy.md),
[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md)

## Examples

``` r
prior <- constraintL2(mu = c(k1 = 0, k2 = 0), sigma = 2)
prior(pars = c(k1 = 1, k2 = -0.5))
#> value
#> ==================
#> 0.3125
#> gradient[1:2] (full length = 2)
#> ==================
#>    k1    k2 
#>  0.50 -0.25 
#> 
#> hessian[1:2,1:2]
#> ==================
#>     k1  k2
#> k1 0.5 0.0
#> k2 0.0 0.5
#> 
#> 
#> attributes
#> ==================
#>  prior  0.3125

## A common sigma, estimated on log scale by the parameter s
prior_s <- constraintL2(mu = c(k1 = 0, k2 = 0), sigma = "s")
prior_s(pars = c(k1 = 1, k2 = -0.5, s = 0))$gradient
#>   k1   k2    s 
#>  2.0 -1.0  1.5 
```
