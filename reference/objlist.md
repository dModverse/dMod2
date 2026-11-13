# Objective List

An objective list holds an objective value with gradient and Hessian, as
objective functions such as
[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md),
[`constraintL2()`](https://dmodverse.github.io/dMod2/reference/constraintL2.md)
and
[`datapointL2()`](https://dmodverse.github.io/dMod2/reference/datapointL2.md)
return them. Further numeric attributes are added along when two
objective lists are added with `+`, see
[sumobjlist](https://dmodverse.github.io/dMod2/reference/plus-.objlist.md).

## Usage

``` r
objlist(value, gradient, hessian)
```

## Arguments

- value:

  Numeric of length 1.

- gradient:

  Named numeric vector, or `NULL`.

- hessian:

  Matrix with row and column names those of `gradient`, or `NULL`.

## Value

An object of class `objlist`.

## See also

[`as.objlist()`](https://dmodverse.github.io/dMod2/reference/as.objlist.md)

## Examples

``` r
a <- objlist(1, c(a = 1, b = 2),
             matrix(2, nrow = 2, ncol = 2,
                    dimnames = list(c("a", "b"), c("a", "b"))))
a + a
#> value
#> ==================
#> 2
#> gradient[1:2] (full length = 2)
#> ==================
#> a b 
#> 2 4 
#> 
#> hessian[1:2,1:2]
#> ==================
#>   a b
#> a 4 4
#> b 4 4
#> 
#> 
#> attributes
#> ==================
```
