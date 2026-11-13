# Time-Course Data for the JAK-STAT Signaling Pathway

Phosphorylated Epo receptor (`pEpoR`), total phosphorylated STAT in the
cytoplasm (`tpSTAT`) and total STAT in the cytoplasm (`tSTAT`), measured
between minute 0 and 60 after stimulation with Epo.

## Format

A data frame with 47 rows and the columns `time` (minutes), `name`,
`value`, `sigma` and `condition` (`"Epo"`).

## Source

Swameye I, Mueller TG, Timmer J, Sandra O, Klingmueller U (2003),
Identification of nucleocytoplasmic cycling as a remote sensor in
cellular signaling by databased modeling. PNAS 100(3):1028-1033.

## See also

[`as.datalist()`](https://dmodverse.github.io/dMod2/reference/datalist.md)

## Examples

``` r
data(jakstat)
plotData(as.datalist(jakstat))
```
