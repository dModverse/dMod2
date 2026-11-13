# Access the Covariates in the Data

Access the Covariates in the Data

## Usage

``` r
covariates(x, ...)

# S3 method for class 'datalist'
covariates(x, ...)

# S3 method for class 'data.frame'
covariates(x, ...)
```

## Arguments

- x:

  a [datalist](https://dmodverse.github.io/dMod2/reference/datalist.md)
  or a `data.frame` with the columns
  `c("name", "time", "value", "sigma", "lloq")`.

- ...:

  not used.

## Value

A data frame with one row per condition: the `condition.grid` of a
datalist, or the unique combinations of the non-standard columns of a
data frame.

## See also

[`as.datalist()`](https://dmodverse.github.io/dMod2/reference/datalist.md)
