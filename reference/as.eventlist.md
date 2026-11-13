# Coerce to Eventlist

Coerce to Eventlist

## Usage

``` r
as.eventlist(x, ...)

# S3 method for class 'list'
as.eventlist(x, ...)

# S3 method for class 'data.frame'
as.eventlist(x, ...)
```

## Arguments

- x:

  list or data frame with the entries `var`, `time`, `value` and
  `method`, and optionally `root` (default `NA`).

- ...:

  not used.

## Value

Object of class `eventlist`.

## See also

[`eventlist()`](https://dmodverse.github.io/dMod2/reference/eventlist.md)
