# Subset a Parameter Vector

Subsets a `parvec` together with the derivatives.

## Usage

``` r
# S3 method for class 'parvec'
x[..., drop = FALSE]
```

## Arguments

- x:

  object of class `parvec`.

- ...:

  subsetting indices.

- drop:

  logical, drop derivative columns that are zero after subsetting.
  Default `FALSE`.

## Value

Object of class `parvec`.
