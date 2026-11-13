# Direct Sum of Datasets

Combines two datalists.

## Usage

``` r
# S3 method for class 'datalist'
e1 + e2
```

## Arguments

- e1:

  object of class `datalist`.

- e2:

  object of class `datalist`.

## Value

Object of class `datalist` for the union of conditions.

## Details

A condition present in both datalists is taken from `e2`, with a
warning. The condition grids are combined.

## See also

[`as.datalist()`](https://dmodverse.github.io/dMod2/reference/datalist.md)

## Examples

``` r
mydata1 <- data.frame(
  name = "A",
  time = 0:1,
  value = 1:2,
  sigma = .1,
  compound = c("DEM", "APAP"),
  dose = "0.1"
)

mydata2 <- data.frame(
  name = "A",
  time = 0:1,
  value = 3:4,
  sigma = .1,
  compound = c("APAP", "DCF"),
  dose = "0.1"
)

data1 <- as.datalist(mydata1, splitBy = c("compound", "dose"))
data2 <- as.datalist(mydata2, splitBy = c("compound", "dose"))

## Condition APAP_0.1 is taken from data2, with a warning
data <- data1 + data2
#> Warning: Condition APAP_0.1 existed and has been overwritten.
data
#> DEM_0.1:
#>   name time value sigma lloq
#> 1    A    0     1   0.1 -Inf
#> APAP_0.1:
#>   name time value sigma lloq
#> 1    A    0     3   0.1 -Inf
#> DCF_0.1:
#>   name time value sigma lloq
#> 2    A    1     4   0.1 -Inf
attr(data, "condition.grid")
#>          compound dose
#> DEM_0.1       DEM  0.1
#> APAP_0.1     APAP  0.1
#> DCF_0.1       DCF  0.1
```
