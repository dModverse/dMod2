# Extract the Second Derivatives of an Object

Extract the Second Derivatives of an Object

## Usage

``` r
getDerivs2(x, ...)

# S3 method for class 'parvec'
getDerivs2(x, ...)

# S3 method for class 'prdframe'
getDerivs2(x, ...)

# S3 method for class 'prdlist'
getDerivs2(x, ...)

# S3 method for class 'list'
getDerivs2(x, ...)

# S3 method for class 'objlist'
getDerivs2(x, ...)
```

## Arguments

- x:

  object of class `parvec`, `prdframe`, `prdlist`, `objlist`, or a list
  of such objects.

- ...:

  not used.

## Value

Depends on the class of `x`:

- `parvec`: a 3D array `[parameter, theta, theta]`.

- `prdframe`: a `prdframe` with the column `time` and one column per
  variable and pair of parameters, for the upper triangle including the
  diagonal.

- `prdlist`: a `prdlist` of such `prdframe`s.

- `objlist`: the Hessian matrix.

- `list`: a list with the result for each element.

An error is raised if `x` has no second derivatives.

## See also

[`getDerivs()`](https://dmodverse.github.io/dMod2/reference/getDerivs.md)

## Examples

``` r
H <- array(0, dim = c(1, 2, 2), dimnames = list("a", c("x", "y"), c("x", "y")))
H["a", "x", "y"] <- H["a", "y", "x"] <- 1
p <- as.parvec(c(a = 1), deriv2 = H)
getDerivs2(p)
#> , , x
#> 
#>   x y
#> a 0 1
#> 
#> , , y
#> 
#>   x y
#> a 1 0
#> 
```
