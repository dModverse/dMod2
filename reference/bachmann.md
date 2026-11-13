# Time-Course Data for JAK2-STAT5 Signaling in CFU-E Cells

541 measurements from thirteen experiments on Epo stimulated erythroid
progenitor cells: phosphorylated JAK2 and Epo receptor, total and
phosphorylated STAT5, the feedback proteins CIS, SOCS3 and SHP1, their
transcripts, and dose responses at four times.

## Format

A data frame with 541 rows and 11 columns:

- `time`:

  minutes after stimulation.

- `name`:

  observable, 20 levels.

- `value`:

  signal divided by the maximum signal of the observable.

- `sigma`:

  `NA` throughout; the errors are estimated in the fit.

- `condition`:

  one of 36 experimental conditions.

- `experiment`:

  one of 13 experiments; selects the scale and offset parameters.

- `epo_level`:

  Epo dose.

- `ActD`:

  1 where actinomycin D was added, else 0. The PEtab copy of the problem
  codes it the other way round.

- `CISoe`, `SOCS3oe`, `SHP1oe`:

  1 where CIS, SOCS3 or SHP1 is overexpressed, else 0.

## Source

Bachmann J, Raue A, Schilling M, Boehm ME, Kreutz C, Kaschek D, Busch H,
Gretz N, Lehmann WD, Timmer J, Klingmueller U (2011), Division of labor
by dual feedback regulators controls JAK2/STAT5 signaling over broad
ligand range. Mol Syst Biol 7:516. Values as distributed in the PEtab
benchmark collection, model `Bachmann_MSB2011`, which is installed with
the package, see
`system.file("extdata", "petab_bachmann", package = "dMod2")`.

## See also

The example project [JAK2-STAT5 signaling (Bachmann
2011)](https://dmodverse.github.io/dMod2/articles/bachmannMSB2011.html)
builds and fits the model.

## Examples

``` r
data(bachmann)
head(bachmann)
#>   time   name      value sigma condition experiment epo_level ActD CISoe
#> 1    0 CIS_au 0.01700748    NA      long       long  1.25e-07    0     0
#> 2    5 CIS_au 0.03881933    NA      long       long  1.25e-07    0     0
#> 3   10 CIS_au 0.02134380    NA      long       long  1.25e-07    0     0
#> 4   20 CIS_au 0.03708581    NA      long       long  1.25e-07    0     0
#> 5   40 CIS_au 0.12873028    NA      long       long  1.25e-07    0     0
#> 6   60 CIS_au 0.28249349    NA      long       long  1.25e-07    0     0
#>   SOCS3oe SHP1oe
#> 1       0      0
#> 2       0      0
#> 3       0      0
#> 4       0      0
#> 5       0      0
#> 6       0      0
```
