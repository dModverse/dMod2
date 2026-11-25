# Print Equation Vector

Print Equation Vector

## Usage

``` r
# S3 method for class 'eqnvec'
print(x, width = 140, pander = FALSE, ...)
```

## Arguments

- x:

  object of class
  [eqnvec](https://dmodverse.github.io/dMod2/reference/eqnvec.md).

- width:

  numeric, width of the print-out. Default 140.

- pander:

  logical, use pander for output (R Markdown). Default `FALSE`.

- ...:

  not used.

## Value

`x`, invisibly.

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
