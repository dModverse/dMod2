# Mavoglurant Pharmacokinetics

Plasma concentration of mavoglurant (`cp`) in 120 subjects after an
intravenous infusion of 10 min (20 min in one occasion); 78 subjects
were dosed on two occasions. Time counts from the start of the infusion
of each occasion.

## Format

A data frame with 2427 rows, one per reading, and 13 columns:

- `name`:

  `"cp"`.

- `time`:

  hours after the start of the infusion of the occasion.

- `value`:

  the concentration.

- `sigma`:

  `NA` throughout.

- `condition`:

  subject and occasion, `"<ID>_<OCC>"`.

- `ID`, `OCC`:

  subject and occasion.

- `dose`:

  the infused dose in mg, 25, 37.5 or 50.

- `duration`:

  the duration of the infusion in hours, `AMT / RATE` of the source.

- `AGE`, `SEX`, `WT`, `HT`:

  age (years), sex (1 or 2), weight (kg) and height (m).

## Source

Package nlmixr2data, data set `mavoglurant` (GPL \>= 3), the readings
without the rows of missing values.

Wendling T, Tsamandouras N, Dumitras S, Pigeolet E, Ogungbenro K, Aarons
L (2016), Reduction of a whole-body physiologically based
pharmacokinetic model to stabilise the Bayesian analysis of clinical
data. The AAPS Journal 18(1):196-209.

## See also

The example project [Theophylline, warfarin and
mavoglurant](https://dmodverse.github.io/dMod2/articles/pharmacometrics.html)
fits a two-compartment model with variability between occasions.

## Examples

``` r
data(mavoglurant)
head(mavoglurant)
#>   name  time value sigma condition  ID OCC dose  duration AGE SEX   WT       HT
#> 1   cp 0.200   491    NA     793_1 793   1   25 0.3333333  42   1 94.3 1.769997
#> 2   cp 0.250   605    NA     793_1 793   1   25 0.3333333  42   1 94.3 1.769997
#> 3   cp 0.367   556    NA     793_1 793   1   25 0.3333333  42   1 94.3 1.769997
#> 4   cp 0.533   310    NA     793_1 793   1   25 0.3333333  42   1 94.3 1.769997
#> 5   cp 0.700   237    NA     793_1 793   1   25 0.3333333  42   1 94.3 1.769997
#> 6   cp 1.200   147    NA     793_1 793   1   25 0.3333333  42   1 94.3 1.769997
```
