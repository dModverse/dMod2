# Identity Function

`Id() * f` and `f * Id()` return `f`.

## Usage

``` r
Id()
```

## Value

A function of class `idfn`.

## Examples

``` r
x <- Xt()
identical(Id() * x, x)
#> [1] TRUE
identical(x * Id(), x)
#> [1] TRUE
(Id() * x)(0:3, c(a = 1))
#> NULL:
#>      time
#> [1,]    0
#> [2,]    1
#> [3,]    2
#> [4,]    3
#> 
#> The prediction contains 1st-order derivatives: yes [4 x 1 x 1]
#> The prediction contains 2nd-order derivatives: no
```
