# Combine Data Frames or Matrices by Rows

Binds data frames or matrices by rows after adding the columns each one
lacks: `NA` for data frames, 0 for matrices.

## Usage

``` r
combine(...)
```

## Arguments

- ...:

  Data frames or matrices, with column names that need not overlap.
  `NULL` entries are dropped.

## Value

A data frame or matrix with the union of the column names.

## Examples

``` r
data1 <- data.frame(Description = "reaction 1", Rate = "k1*A", A = -1, B = 1)
data2 <- data.frame(Description = "reaction 2", Rate = "k2*B", B = -1, C = 1)
combine(data1, data2)
#>   Description Rate  A  B  C
#> 1  reaction 1 k1*A -1  1 NA
#> 2  reaction 2 k2*B NA -1  1
```
