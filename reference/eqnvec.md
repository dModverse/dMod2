# Generate Equation Vector Object

The eqnvec object stores explicit algebraic equations, like the
right-hand sides of an ODE, observation functions or parameter
transformations as named character vectors.

## Usage

``` r
is.eqnvec(x)

eqnvec(...)
```

## Arguments

- x:

  object of any class

- ...:

  mathematical expressions as characters to be coerced, the right-hand
  sides of the equations

## Value

Object of class `eqnvec`, a named character vector.

## See also

[eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md)

## Examples

``` r
v <- eqnvec(y = "2*x + offset")
print(v)
#> Idx Inner <- Outer
#>   1     y <- 2*x+offset
is.eqnvec(v)
#> [1] TRUE
```
