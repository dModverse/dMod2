# Time-Course Data for the Bile Acid Demonstration Model

Bile acid amounts in the buffer (`buffer`) and in the cells (`cellular`)
between 0.1 and 41 time units, in two conditions, `"closed"` and
`"open"`.

## Format

A data frame with 32 rows and the columns `name`, `time`, `value`,
`sigma`, `lloq` (lower limit of quantification, `-Inf` throughout) and
`condition`.

## See also

The example project [Bile acid
transport](https://dmodverse.github.io/dMod2/articles/bileAcidTransport.html)
builds and fits the model.

## Examples

``` r
data(badata)
plotData(as.datalist(badata))
```
