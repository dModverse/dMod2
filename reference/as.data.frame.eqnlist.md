# Coerce Equation List into a Data Frame

Coerce Equation List into a Data Frame

## Usage

``` r
# S3 method for class 'eqnlist'
as.data.frame(x, ...)
```

## Arguments

- x:

  object of class
  [eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md).

- ...:

  not used.

## Value

A `data.frame` with the columns "Description" (character), "Rate"
(character) and one column per state, the columns of the stoichiometric
matrix. The compartment layout is kept as attributes and read back by
[`as.eqnlist()`](https://dmodverse.github.io/dMod2/reference/eqnlist.md).
`NULL` for an empty equation list.
