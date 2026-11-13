# Aligned values as dMod2 data

Aligned values as dMod2 data

## Usage

``` r
# S3 method for class 'alignment'
as.datalist(x, time = "time", name = "symbol", condition = NULL, ...)
```

## Arguments

- x:

  An `alignment` from
  [`align()`](https://dmodverse.github.io/dMod2/reference/align.md).

- time, name, condition:

  Columns of the aligned values that become `time`, `name` and
  `condition`; by default `"time"`, `"symbol"` (the biological symbol)
  and the remaining columns pasted.

- ...:

  Ignored.

## Value

A [datalist](https://dmodverse.github.io/dMod2/reference/datalist.md)
with `value` and `sigma` on the scale of the biological symbol.
