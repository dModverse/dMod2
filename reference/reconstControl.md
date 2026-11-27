# Settings of the Observability Engine

Saturation and closed-form reconstruction settings for
`symmetryDetection(control = reconstControl())`. Raise the caps to
reconstruct wide or high-degree directions, at the cost of more samples.

## Usage

``` r
reconstControl(
  relevanceCap = 6L,
  relevanceCapDir = 24L,
  relevanceCapSparse = 30L,
  degreeCap = 4L,
  sampleSlack = 5L,
  probeRetries = 8L,
  laurentDegNum = 4L,
  laurentDegDen = 2L,
  laurentCandCap = 200000L,
  termCap = 60L,
  generalDegNum = 4L,
  generalDegDen = 3L,
  gapOrderCap = 8L,
  minSupportCandCap = 20000L,
  perPrimeCap = 120L,
  perPrimeMinPrimes = 3L,
  homotopy = TRUE,
  homotopyDegree = NULL,
  degreeCapMax = NULL,
  exactProof = NULL,
  exactWorkCap = NULL,
  seed = NULL
)
```

## Arguments

- relevanceCap:

  Maximum number of coordinates in one entry of a direction for the
  dense fit; wider entries use sparse interpolation. Defaults to 6.

- relevanceCapDir:

  Maximum number of coordinates in one direction; a wider one is
  reported by the support. Defaults to 24.

- relevanceCapSparse:

  Maximum number of coordinates in one entry for the sparse
  (Ben-Or/Tiwari) fit; a wider entry is reported by the support.
  Defaults to 30.

- degreeCap:

  Total degree bound of the dense rational fit. Defaults to 4.

- sampleSlack:

  Number of samples beyond the minimum of the fit. Defaults to 5.

- probeRetries:

  Number of retries of the relevance probe when a perturbation changes
  the pivots. Defaults to 8.

- laurentDegNum, laurentDegDen:

  Numerator degree and monomial denominator degree bounds of the sparse
  Laurent fit. Default to 4 and 2.

- laurentCandCap:

  Maximum number of candidate monomials of the Laurent and general
  sparse fits. Defaults to 200000.

- termCap:

  Maximum number of terms of a sparse entry. Defaults to 60.

- generalDegNum, generalDegDen:

  Numerator and denominator degree bounds of the general sparse rational
  fit. Default to 4 and 3.

- gapOrderCap:

  Maximum order of the power series in the time between events. Curve
  samples of resting conditions start at this order and double it at
  most twice where a kernel at eps = 1 needs more. Defaults to 8.

- minSupportCandCap:

  Maximum number of column subsets searched for directions with small
  support. Defaults to 20000.

- perPrimeCap:

  Maximum number of samples per prime for the reconstruction under
  `equilibrate = TRUE`. Defaults to 120.

- perPrimeMinPrimes:

  Minimum number of primes with samples for a reconstruction under
  `equilibrate = TRUE`. Defaults to 3.

- homotopy:

  Whether directions still open under `equilibrate = TRUE` are fitted
  along lines in parameter space as a last route. Defaults to `TRUE`.

- homotopyDegree:

  Total degree bound of the entries the homotopy route fits.

- degreeCapMax:

  Degree up to which an entry past `degreeCap` is probed again, the
  bound doubled each time.

- exactProof:

  Whether a closed form that differs from the lift of the form is proven
  through the fields of each regime, one detection per regime.

- exactWorkCap:

  Work bound of that proof: regimes times coordinates times tape
  instructions; above it the proof is skipped and the limit noted.

- seed:

  Seed of the random choices.

## Value

A list of class `reconstcontrol` with the settings.

## Details

`homotopyDegree`, `degreeCapMax`, `exactProof`, `exactWorkCap` and
`seed` default to `NULL`, which keeps the value of symident (4, 8,
`TRUE`, 500000 and 0).

## See also

[`symmetryDetection()`](https://dmodverse.github.io/dMod2/reference/symmetryDetection.md)

## Examples

``` r
reconstControl(degreeCap = 6L)
#> $relevanceCap
#> [1] 6
#> 
#> $relevanceCapDir
#> [1] 24
#> 
#> $relevanceCapSparse
#> [1] 30
#> 
#> $degreeCap
#> [1] 6
#> 
#> $sampleSlack
#> [1] 5
#> 
#> $probeRetries
#> [1] 8
#> 
#> $laurentDegNum
#> [1] 4
#> 
#> $laurentDegDen
#> [1] 2
#> 
#> $laurentCandCap
#> [1] 200000
#> 
#> $termCap
#> [1] 60
#> 
#> $generalDegNum
#> [1] 4
#> 
#> $generalDegDen
#> [1] 3
#> 
#> $gapOrderCap
#> [1] 8
#> 
#> $minSupportCandCap
#> [1] 20000
#> 
#> $perPrimeCap
#> [1] 120
#> 
#> $perPrimeMinPrimes
#> [1] 3
#> 
#> $homotopy
#> [1] TRUE
#> 
#> $homotopyDegree
#> NULL
#> 
#> $degreeCapMax
#> NULL
#> 
#> $exactProof
#> NULL
#> 
#> $exactWorkCap
#> NULL
#> 
#> $seed
#> NULL
#> 
#> attr(,"class")
#> [1] "reconstcontrol" "list"          
```
