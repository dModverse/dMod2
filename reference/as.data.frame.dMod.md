# Coerce to a Data Frame

Coerce to a Data Frame

## Usage

``` r
# S3 method for class 'datalist'
as.data.frame(x, ...)

# S3 method for class 'prdlist'
as.data.frame(x, ..., data = NULL, errfn = NULL)
```

## Arguments

- x:

  any R object

- ...:

  not used.

- data:

  object of class `datalist`. Default `NULL`. The condition grid of the
  data adds covariate columns to the output.

- errfn:

  object of class `obsfn`, the error model that predicts the column
  `sigma`. Default `NULL` sets `sigma` to `NaN`.

## Value

a data frame
