# Time-Course Data for IFN-alpha Signaling in Huh7.5 Cells

1918 measurements on human hepatoma cells stimulated with IFN-alpha at 0
h, at 24 h or at both: phosphorylated and total STAT1 and STAT2 and IRF9
in the cytoplasm and the nucleus, USP18, seven target gene transcripts
and the molecules per cell of STAT1, STAT2, IRF9 and USP18. Two
experiments perturb the USP18 feedback, by knock-down and by
overexpression.

## Format

A data frame with 1918 rows and 8 columns:

- `time`:

  hours after the start of the experiment, when the first dose and the
  doxycycline are given.

- `name`:

  observable, 22 levels. Immunoblots (`*_Cyt`, `*_Nuc`) and qPCR
  (`*_mRNA`, fold change against untreated cells) are on the natural
  log, molecules per cell (`mc_*`) on the linear scale.

- `value`:

  the measurement, blots scaled across replicates.

- `sigma`:

  standard deviation from the combined scaling and error model.

- `condition`:

  one of 25 conditions, `<experiment>_<ifn1>_<ifn2>`, for example
  `"WT_10_5000"` for 10 U/ml at 0 h and 5000 U/ml at 24 h.

- `experiment`:

  `"WT"`, `"USP18kd"` (USP18 knock-down by siRNA) or `"USP18oe"` (USP18
  overexpression, induced by doxycycline from 0 h on).

- `ifn1`, `ifn2`:

  IFN-alpha dose at 0 h and at 24 h in U/ml.

## Source

Kok F, Rosenblatt M, Teusel M, et al. (2020), Disentangling molecular
mechanisms regulating sensitization of interferon alpha signal
transduction. Mol Syst Biol 16:e8955,
[doi:10.15252/msb.20198955](https://doi.org/10.15252/msb.20198955) .
Values as fitted for the Huh7.5 core model of the publication.

## See also

The example project [IFN-α sensitization (Kok
2020)](https://dmodverse.github.io/dMod2/articles/kokMSB2020.html)
builds the model.

## Examples

``` r
data(kok)
head(kok)
#>   time     name      value     sigma      condition experiment ifn1 ifn2
#> 1 24.0 IRF9_Cyt -3.1334241 0.2728416 USP18oe_0_5000    USP18oe    0 5000
#> 2 24.5 IRF9_Cyt -3.3229917 0.2728416 USP18oe_0_5000    USP18oe    0 5000
#> 3 25.0 IRF9_Cyt -2.9660025 0.2728416 USP18oe_0_5000    USP18oe    0 5000
#> 4 26.0 IRF9_Cyt -2.1692097 0.2728416 USP18oe_0_5000    USP18oe    0 5000
#> 5 27.0 IRF9_Cyt -1.7567120 0.2728416 USP18oe_0_5000    USP18oe    0 5000
#> 6 28.0 IRF9_Cyt -0.9608051 0.2728416 USP18oe_0_5000    USP18oe    0 5000
```
