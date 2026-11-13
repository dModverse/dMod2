# Subset of an Equation List

Subset of an Equation List

## Usage

``` r
# S3 method for class 'eqnlist'
subset(x, ...)
```

## Arguments

- x:

  object of class
  [eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md).

- ...:

  logical expression for subsetting.

## Value

Object of class
[eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md).

## Details

The expression in `...` can refer to `Educt`, `Product`, `Rate` and
`Description`. Within it, `%in%` searches the educts and products of
each reaction (see examples).

## See also

[`addReaction()`](https://dmodverse.github.io/dMod2/reference/addReaction.md),
[`as.eqnlist()`](https://dmodverse.github.io/dMod2/reference/eqnlist.md)

## Examples

``` r
reactions <- data.frame(Description = c("Activation", "Deactivation"), 
                        Rate = c("act*A", "deact*pA"), A=c(-1,1), pA=c(1, -1) )
f <- as.eqnlist(reactions)
subset(f, "A" %in% Educt)
#> Conserved quantities:
#>   total_1 = 1*A+1*pA
#> 
#>   Check Educt -> Product  Rate Description
#> 1           A ->      pA act*A  Activation
subset(f, "pA" %in% Product)
#> Conserved quantities:
#>   total_1 = 1*A+1*pA
#> 
#>   Check Educt -> Product  Rate Description
#> 1           A ->      pA act*A  Activation
subset(f, grepl("act", Rate))
#> Conserved quantities:
#>   total_1 = 1*A+1*pA
#> 
#>   Check Educt -> Product     Rate  Description
#> 1           A ->      pA    act*A   Activation
#> 2          pA ->       A deact*pA Deactivation
```
