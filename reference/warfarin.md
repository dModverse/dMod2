# Warfarin Pharmacokinetics and Pharmacodynamics

Plasma concentration of warfarin (`cp`, mg/L) and prothrombin complex
activity (`pca`, percent of normal) in 32 subjects after one oral dose
at time 0. The four concentrations recorded as 0, at 0.5 h, are left
out; their limit of quantification is not known.

## Format

A data frame with 479 rows, one per reading, and 10 columns:

- `name`:

  `"cp"` or `"pca"`.

- `time`:

  hours after the dose.

- `value`:

  the reading.

- `sigma`:

  `NA` throughout.

- `condition`:

  the subject, `"ID1"` to `"ID32"`.

- `id`:

  the subject, 1 to 32.

- `dose`:

  the oral dose in mg.

- `wt`, `age`, `sex`:

  weight (kg), age (years) and sex.

## Source

Package nlmixr2data, data set `warfarin` (GPL \>= 3), after Funaki T,
Holford N, Fujita S (2018), Population PKPD analysis using nlmixr2 and
NONMEM, PAGJA 2018.

O'Reilly RA, Aggeler PM, Leong LS (1963), Studies on the coumarin
anticoagulant drugs: the pharmacodynamics of warfarin in man. Journal of
Clinical Investigation 42(10):1542-1551.

O'Reilly RA, Aggeler PM (1968), Studies on coumarin anticoagulant drugs:
initiation of warfarin therapy without a loading dose. Circulation
38:169-177.

## See also

The example project [Theophylline, warfarin and
mavoglurant](https://dmodverse.github.io/dMod2/articles/pharmacometrics.html)
fits a PK/PD model.

## Examples

``` r
data(warfarin)
head(warfarin)
#>   name time value sigma condition id dose   wt age  sex
#> 1   cp    1   1.9    NA       ID1  1  100 66.7  50 male
#> 2   cp    2   3.3    NA       ID1  1  100 66.7  50 male
#> 3   cp    3   6.6    NA       ID1  1  100 66.7  50 male
#> 4   cp    6   9.1    NA       ID1  1  100 66.7  50 male
#> 5   cp    9  10.8    NA       ID1  1  100 66.7  50 male
#> 6   cp   12   8.6    NA       ID1  1  100 66.7  50 male
```
