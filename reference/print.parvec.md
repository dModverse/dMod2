# Print Parameter Vectors

Prints a parameter vector with information about derivatives and fixed
parameters.

## Usage

``` r
# S3 method for class 'parvec'
print(x, ...)
```

## Arguments

- x:

  object of class `parvec`.

- ...:

  not used.

## Value

`x`, invisibly.

## Examples

``` r
## Parameter vector without derivatives
p <- as.parvec(c(a = 1, b = 2, c = 3))
p
#> Parameter vector:
#>   a :       1
#>   b :       2
#>   c :       3
#> 
#> Attributes:
#>   deriv  : <none>
#>   deriv2 : <none>
#>   fixed  : <none>

## Jacobian and second derivatives
J <- matrix(c(1, 0, 0,
              0, 2, 0,
              0, 0, 3),
            nrow = 3, byrow = TRUE,
            dimnames = list(c("a", "b", "c"), c("x", "y", "z")))

H <- array(0, dim = c(3, 3, 3),
           dimnames = list(c("a", "b", "c"),
                           c("x", "y", "z"),
                           c("x", "y", "z")))
H["a", "x", "x"] <- 1
H["b", "y", "y"] <- 2
H["c", "z", "z"] <- 3

p2 <- as.parvec(c(a = 10, b = 20, c = 30), deriv = J, deriv2 = H)
p2
#> Parameter vector:
#>   a :      10
#>   b :      20
#>   c :      30
#> 
#> Attributes:
#>   deriv  : 3 x 3 matrix
#>   deriv2 : 3 x 3 x 3 array
#>   fixed  : <none>
getDerivs(p2)
#>   x y z
#> a 1 0 0
#> b 0 2 0
#> c 0 0 3
getDerivs2(p2)
#> , , x
#> 
#>   x y z
#> a 1 0 0
#> b 0 0 0
#> c 0 0 0
#> 
#> , , y
#> 
#>   x y z
#> a 0 0 0
#> b 0 2 0
#> c 0 0 0
#> 
#> , , z
#> 
#>   x y z
#> a 0 0 0
#> b 0 0 0
#> c 0 0 3
#> 

## Subsetting keeps the matching derivatives
p_sub <- p2[c("a", "b")]
getDerivs(p_sub)
#>   x y z
#> a 1 0 0
#> b 0 2 0

## Concatenation
p3 <- as.parvec(c(d = 4, e = 5))
c(p2, p3)
#> Parameter vector:
#>   a :      10
#>   b :      20
#>   c :      30
#>   d :       4
#>   e :       5
#> 
#> Attributes:
#>   deriv  : 3 x 3 matrix
#>   deriv2 : 5 x 3 x 3 array
#>   fixed  : d, e

## From numbers
parvec(a = 1, b = 2)
#> Parameter vector:
#>   a :       1
#>   b :       2
#> 
#> Attributes:
#>   deriv  : <none>
#>   deriv2 : <none>
#>   fixed  : <none>
```
