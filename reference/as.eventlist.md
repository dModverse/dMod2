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

## Examples

``` r
events <- data.frame(var = c("A", "B"), time = c(5, 10), value = c(1, 0),
                     method = c("add", "replace"))
as.eventlist(events)
#> Eventlist with 2 events:
#>   var time value root  method
#> 1   A    5     1 <NA>     add
#> 2   B   10     0 <NA> replace
```
