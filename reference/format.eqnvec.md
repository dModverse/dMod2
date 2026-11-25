# Format Equation Vector with Spaces

Format Equation Vector with Spaces

## Usage

``` r
# S3 method for class 'eqnvec'
format(x, ...)
```

## Arguments

- x:

  object of class
  [eqnvec](https://dmodverse.github.io/dMod2/reference/eqnvec.md), or a
  named parsable character vector.

- ...:

  not used.

## Value

Named character vector with spaces around `+`, `-`, `*` and `/`.

## Examples

``` r
v <- eqnvec(y = "2*x + offset")
print(v)
#> Idx Inner <- Outer
#>   1     y <- 2*x+offset
is.eqnvec(v)
#> [1] TRUE

format(eqnvec(y = "2*x+offset"))
#>                y 
#> "2 * x + offset" 
summary(v)
#> y = f( x, offset)
```
