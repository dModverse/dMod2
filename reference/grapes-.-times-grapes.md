# Multiplication of Objective Functions with Scalars

`x1 %.*% x2` multiplies an object of class `objfn` or `objlist` by a
scalar.

## Usage

``` r
x1 %.*% x2
```

## Arguments

- x1:

  numeric of length one.

- x2:

  object of class `objfn` or `objlist`.

## Value

For an `objfn`, an `objfn` whose value, gradient and Hessian are scaled
by `x1`. For an `objlist`, the `objlist` with every numeric entry and
attribute scaled. Otherwise `x1 * x2`.

## See also

[+.objfn](https://dmodverse.github.io/dMod2/reference/plus-.objfn.md)

## Examples

``` r
obj <- constraintL2(mu = c(a = 0, b = 0), sigma = 1)
obj2 <- 2 %.*% obj
obj2(c(a = 1, b = 2))$value
#> [1] 10
```
