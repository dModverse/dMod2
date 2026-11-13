# Bind a Named List of Data Frames into One Data Frame

Adds to each data frame a column `condition` with the name of the data
frame in the list and binds them by rows.

## Usage

``` r
lbind(mylist)
```

## Arguments

- mylist:

  A named list of data frames; a column missing in one of them is `NA`
  there.

## Value

A data frame with the original columns and `condition`.

## Examples

``` r
lbind(list(C1 = data.frame(x = 1:2), C2 = data.frame(x = 3)))
#>   x condition
#> 1 1        C1
#> 2 2        C1
#> 3 3        C2
```
