# Coordinate Transformation for Data Frames

Applies a symbolically defined transformation to the `value` column of a
data frame. A `sigma` column is transformed by first-order error
propagation.

## Usage

``` r
coordTransform(data, transformations)
```

## Arguments

- data:

  Data frame with columns `name` (character) and `value` (numeric),
  optionally `sigma` (numeric).

- transformations:

  Character, one transformation for all names, or a named list of
  characters, one per entry of the `name` column. Each transformation is
  an expression in a single symbol, e.g. `"log(x)"`.

## Value

The data frame with transformed `value` and `sigma`.

## Examples

``` r
d <- data.frame(name = rep(c("u", "v"), each = 3), time = 0:2, value = 1:6,
                sigma = 0.2)
coordTransform(d, "log(y)")
#>   name time     value      sigma
#> 1    u    0 0.0000000 0.20000000
#> 2    u    1 0.6931472 0.10000000
#> 3    u    2 1.0986123 0.06666667
#> 4    v    0 1.3862944 0.05000000
#> 5    v    1 1.6094379 0.04000000
#> 6    v    2 1.7917595 0.03333333
coordTransform(d, list(u = "y^2", v = "exp(-y)"))
#>   name time       value        sigma
#> 1    u    0 1.000000000 0.4000000000
#> 2    u    1 4.000000000 0.8000000000
#> 3    u    2 9.000000000 1.2000000000
#> 4    v    0 0.018315639 0.0036631278
#> 5    v    1 0.006737947 0.0013475894
#> 6    v    2 0.002478752 0.0004957504
```
