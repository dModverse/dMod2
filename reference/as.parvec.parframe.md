# Select a Parameter Vector from a Parameter Frame

Select a Parameter Vector from a Parameter Frame

## Usage

``` r
# S3 method for class 'parframe'
as.parvec(x, index = 1, ...)

# S3 method for class 'parlist'
as.parvec(x, index = 1, ...)
```

## Arguments

- x:

  object of class `parframe`, e.g. from
  [`as.parframe()`](https://dmodverse.github.io/dMod2/reference/as.parframe.md),
  or the `parlist` of
  [`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md).

- index:

  integer, select the row with the `index`-th lowest objective value.
  Default 1, the best fit. Without a `value` column, the `index`-th row.

- ...:

  not used.

## Value

Object of class `parvec` with the parameter columns of the selected row.
A warning is issued if that fit did not converge.

## See also

[`as.parframe()`](https://dmodverse.github.io/dMod2/reference/as.parframe.md),
[`parvec()`](https://dmodverse.github.io/dMod2/reference/parvec.md)
