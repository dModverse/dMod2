# Summary of an Equation Vector

Prints, for each equation, the symbols it depends on.

## Usage

``` r
# S3 method for class 'eqnvec'
summary(object, ...)
```

## Arguments

- object:

  object of class
  [eqnvec](https://dmodverse.github.io/dMod2/reference/eqnvec.md).

- ...:

  not used.

## Value

`object`, invisibly.

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
