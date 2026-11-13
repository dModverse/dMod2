# Coerce to an Equation Vector

Coerce to an Equation Vector

## Usage

``` r
as.eqnvec(x, ...)

# S3 method for class 'character'
as.eqnvec(x = NULL, names = NULL, ...)

# S3 method for class 'eqnlist'
as.eqnvec(x, ...)
```

## Arguments

- x:

  object of class `character` or `eqnlist`.

- ...:

  for the `eqnlist` method, arguments going to
  [`getFluxes()`](https://dmodverse.github.io/dMod2/reference/getFluxes.md),
  e.g. `type`; otherwise not used.

- names:

  character, the left-hand sides of the equations. Default `NULL` uses
  `names(x)`.

## Value

Object of class
[eqnvec](https://dmodverse.github.io/dMod2/reference/eqnvec.md).

## Details

For an `eqnlist`, the right-hand sides of the ODE are assembled from
[`getFluxes()`](https://dmodverse.github.io/dMod2/reference/getFluxes.md).
