# Rows of a Parameter Frame with Unique Values

Keeps one row per distinct objective value.

## Usage

``` r
# S3 method for class 'parframe'
unique(x, incomparables = FALSE, tol = 1, ...)
```

## Arguments

- x:

  parameter frame

- incomparables:

  not used. Argument exists for compatibility with S3 generic.

- tol:

  tolerance to decide when values are assumed to be equal, see
  [`plotValues()`](https://dmodverse.github.io/dMod2/reference/plotValues.md).

- ...:

  Not used. Select rows first with
  [`subset()`](https://rdrr.io/r/base/subset.html), as in
  `unique(subset(x, converged))`.

## Value

The rows of `x` with distinct values, sorted by value.
