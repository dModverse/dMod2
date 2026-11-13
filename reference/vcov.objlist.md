# Variance-Covariance Matrix of a Fit

Methods of [`stats::vcov()`](https://rdrr.io/r/stats/vcov.html) for an
evaluation of an objective function and for the result of
[`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md).

## Usage

``` r
# S3 method for class 'objlist'
vcov(object, fixed = NULL, which = NULL, ...)

# S3 method for class 'trustfit'
vcov(object, parupper = NULL, parlower = NULL, which = NULL, ...)
```

## Arguments

- object:

  An evaluation of an objective function, class `objlist`, such as
  `obj(pars)`, or a result of
  [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md),
  class `trustfit`.

- fixed:

  Names of parameters held fixed, for example at a bound. Default
  `NULL`.

- which:

  Names of the parameters whose block of the covariance is returned.
  Default `NULL`, all of them. With a sparse Hessian (a `Matrix`, from
  `hessian = "sparse"`) only these columns are solved for.

- ...:

  Not used.

- parupper, parlower:

  Named upper and lower bounds. Parameters at or beyond a bound count as
  fixed, as do those that `object$atBound` marks. Default `NULL`.

## Value

Square matrix with the parameter names as dimnames: the inverse of half
the Hessian of `object` (pseudo-inverse when singular), with rows and
columns of fixed parameters set to zero. A sparse Hessian is factorized
by sparse Cholesky; when singular, a pseudo-inverse is taken up to 5000
free parameters and an error raised beyond. For a `trustfit`, `NULL`
when it has no `argument`, and an error when it has no Hessian; the fits
of [`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md)
have none, use `vcov(obj(fit$argument))` for them.

## See also

[`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md),
[`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md)

## Examples

``` r
obj <- function(x, ...) {
  list(value = sum(x^2 / c(1, 4)),
       gradient = 2 * x / c(1, 4),
       hessian = diag(2 / c(1, 4)))
}
fit <- trust(obj, c(a = 1, b = 1))
vcov(fit)
#>   a b
#> a 1 0
#> b 0 4

# The same from an evaluation of the objective
ev <- obj(fit$argument)
vcov(objlist(ev$value, ev$gradient, ev$hessian))
#>   a b
#> a 1 0
#> b 0 4
```
