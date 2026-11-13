# Coerce to a Parameter Frame

Coerce to a Parameter Frame

## Usage

``` r
# S3 method for class 'parframe'
as.parframe(x, ...)

# S3 method for class 'parlist'
as.parframe(x, sortBy = "value", ...)

as.parframe(x, ...)
```

## Arguments

- x:

  object to be coerced, e.g. a `parlist`.

- ...:

  not used.

- sortBy:

  character, the column by which the returned parameter frame is sorted.
  Default `"value"`.

## Value

Object of class
[parframe](https://dmodverse.github.io/dMod2/reference/parframe.md). For
a `parlist`, the meta-information columns are `index`, `value`,
`converged` and `iterations`, plus the termination and evaluation counts
and the wall time `elapsed` the fits report. Fits of
[`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md) with
several `phases` add one column `iterations.<phase>` per phase, the
iterations summed over the runs of that phase. The attribute `"run"` of
the parlist is kept.

## See also

[`parframe()`](https://dmodverse.github.io/dMod2/reference/parframe.md),
[`as.parvec.parframe()`](https://dmodverse.github.io/dMod2/reference/as.parvec.parframe.md),
and
[`parlist()`](https://dmodverse.github.io/dMod2/reference/parlist.md)
for an example.
