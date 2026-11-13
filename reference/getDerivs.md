# Extract the First Derivatives of an Object

Extract the First Derivatives of an Object

## Usage

``` r
getDerivs(x, ...)

# S3 method for class 'parvec'
getDerivs(x, ...)

# S3 method for class 'prdframe'
getDerivs(x, ...)

# S3 method for class 'prdlist'
getDerivs(x, ...)

# S3 method for class 'list'
getDerivs(x, ...)

# S3 method for class 'objlist'
getDerivs(x, ...)
```

## Arguments

- x:

  object of class `parvec`, `prdframe`, `prdlist`, `objlist`, or a list
  of such objects.

- ...:

  not used.

## Value

Depends on the class of `x`:

- `parvec`: the Jacobian matrix.

- `prdframe`: a `prdframe` with the column `time` and one column per
  pair of variable and parameter.

- `prdlist`: a `prdlist` of such `prdframe`s.

- `objlist`: the gradient, a named numeric vector.

- `list`: a list with the result for each element.

An error is raised if `x` has no first derivatives.

## See also

[`getDerivs2()`](https://dmodverse.github.io/dMod2/reference/getDerivs2.md),
[`parvec()`](https://dmodverse.github.io/dMod2/reference/parvec.md)

## Examples

``` r
J <- matrix(c(1, 0, 0, 2), 2, dimnames = list(c("a", "b"), c("x", "y")))
p <- as.parvec(c(a = 1, b = 2), deriv = J)
getDerivs(p)
#>   x y
#> a 1 0
#> b 0 2
```
