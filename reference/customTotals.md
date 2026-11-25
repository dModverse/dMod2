# Set or Reset User-Defined Conserved Quantities

Attaches a named list of conservation expressions to the eqnlist. They
replace the detected basis, also as parameters of
[`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md). Each
expression must be linear in the states and conserved by the reactions,
and together they must span as many quantities as
[`conservedQuantities()`](https://dmodverse.github.io/dMod2/reference/conservedQuantities.md)
finds.

## Usage

``` r
customTotals(eqnlist, totals)
```

## Arguments

- eqnlist:

  object of class
  [eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md).

- totals:

  named list of expressions, or `NULL` or
  [`list()`](https://rdrr.io/r/base/list.html) to return to the detected
  basis.

## Value

The eqnlist with `$totals` updated.

## See also

[`conservedQuantities()`](https://dmodverse.github.io/dMod2/reference/conservedQuantities.md)

## Examples

``` r
## A reversible conversion conserves A + B
el <- addReaction(eqnlist(), "A", "B", "k1*A")
el <- addReaction(el, "B", "A", "k2*B")
conservedQuantities(el)
#> $total_1
#> [1] "1*A+1*B"
#> 

## Name the conserved quantity, then return to the detected basis
el <- customTotals(el, list(total = "A + B"))
el$totals
#> $total
#> [1] "A + B"
#> 
#> attr(,"custom")
#> [1] TRUE
el <- customTotals(el, NULL)
```
