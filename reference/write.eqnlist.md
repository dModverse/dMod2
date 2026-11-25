# Write Equation List into a CSV File

Compartments and reaction frames are not written; a warning is issued if
the equation list has any.

## Usage

``` r
write.eqnlist(eqnlist, ...)
```

## Arguments

- eqnlist:

  object of class
  [eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md).

- ...:

  arguments going to
  [`utils::write.csv()`](https://rdrr.io/r/utils/write.table.html).
  `row.names` defaults to `FALSE` and `na` to `""`.

## Value

`NULL`, invisibly. Called for the file it writes.

## See also

[`as.data.frame.eqnlist()`](https://dmodverse.github.io/dMod2/reference/as.data.frame.eqnlist.md)

## Examples

``` r
el <- addReaction(eqnlist(), "A", "B", "k1*A", description = "conversion")
file <- tempfile(fileext = ".csv")
write.eqnlist(el, file = file)
read.csv(file)
#>   Description Rate  A B
#> 1  conversion k1*A -1 1
```
