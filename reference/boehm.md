# Time-Course Data for STAT5 Dimerization in BaF3-EpoR Cells

Relative amounts of phosphorylated STAT5A (`pSTAT5A_rel`),
phosphorylated STAT5B (`pSTAT5B_rel`) and of STAT5A within the total
STAT5 pool (`rSTAT5A_rel`), measured at 16 time points after stimulation
with Epo.

## Format

A data frame with 48 rows and the columns `time` (minutes), `name`,
`value`, `sigma` (`NA` throughout; the errors are estimated in the fit)
and `condition` (`"Boehm2014"`).

## Source

Boehm ME, Adlung L, Schilling M, Roth S, Klingmueller U, Lehmann WD
(2014), Identification of isoform-specific dynamics in
phosphorylation-dependent STAT5 dimerization by quantitative mass
spectrometry and mathematical modeling. J Proteome Res 13(12):5685-5694.
Values as distributed in the PEtab benchmark collection, model
`Boehm_JProteomeRes2014`.

## See also

The example project [STAT5 dimerization (Boehm
2014)](https://dmodverse.github.io/dMod2/articles/boehmJProteomeRes2014.html)
builds and fits the model.

## Examples

``` r
data(boehm)
plotData(as.datalist(boehm))
```
