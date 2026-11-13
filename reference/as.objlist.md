# Zero Objective List

Builds an
[objlist](https://dmodverse.github.io/dMod2/reference/objlist.md) with
value, gradient and Hessian zero for the parameters of `p`.

## Usage

``` r
as.objlist(p)
```

## Arguments

- p:

  Named numeric vector.

## Value

An [objlist](https://dmodverse.github.io/dMod2/reference/objlist.md)
with value `0`, gradient `rep(0, length(p))` and Hessian
`matrix(0, length(p), length(p))`, named by `p`.

## Examples

``` r
p <- c(A = 1, B = 2)
as.objlist(p)
#> value
#> ==================
#> 0
#> gradient[1:2] (full length = 2)
#> ==================
#> A B 
#> 0 0 
#> 
#> hessian[1:2,1:2]
#> ==================
#>   A B
#> A 0 0
#> B 0 0
#> 
#> 
#> attributes
#> ==================
```
