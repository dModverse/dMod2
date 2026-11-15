# Remove Detected Symmetries by Reparametrization

Turns the non-identifiable directions of a
[`symmetryDetection()`](https://dmodverse.github.io/dMod2/reference/symmetryDetection.md)
result into a parameter transformation that removes them. A scaling
fixes one coordinate to 1. A general direction is replaced by the
invariants of that direction, held by new parameters `q_<k>`, searched
up to the degree bounds `dPoly`, `dDarboux` and `dExp`. The computation
runs in the Python package `symident`.

## Usage

``` r
symmetryReduction(
  object,
  fixed = NULL,
  positive = NULL,
  dPoly = 3L,
  dDarboux = 2L,
  dExp = 2L,
  separable = TRUE,
  reportZeroCompatibility = FALSE,
  alternatives = FALSE,
  timeout = Inf,
  allowOpen = FALSE,
  allowChecked = FALSE,
  check = TRUE,
  workCap = 4000L,
  ...
)
```

## Arguments

- object:

  A `symmetrydetection` result from
  [`symmetryDetection()`](https://dmodverse.github.io/dMod2/reference/symmetryDetection.md).

- fixed:

  Coordinates with known values; scalings they remove are dropped.
  `NULL` (default) uses the gauge of `object`.

- positive:

  Coordinates known to be positive: `TRUE` (all), `FALSE` (none) or a
  character vector. `NULL` (default) takes the declaration of the
  detection in `object`.

- dPoly:

  Degree bound of polynomial and rational invariants. Defaults to 3.

- dDarboux:

  Degree bound of the Darboux polynomials. Defaults to 2.

- dExp:

  Numerator degree bound of exponential factors; `0` skips them.
  Defaults to 2.

- separable:

  Logical. Solve separable blocks by quadratures. Defaults to `TRUE`.

- reportZeroCompatibility:

  Logical. Report which coordinates each block can drive to 0. Defaults
  to `FALSE`.

- alternatives:

  Logical. Also return the charts of the other valid pins in
  `$alternatives`. Defaults to `FALSE`.

- timeout:

  Time limit in seconds for the section search of one block; a block not
  solved in time keeps the invariants (`"invariant_only"`). Unlike
  `workCap`, it depends on the load of the machine. Defaults to `Inf`.

- allowOpen:

  Logical. Reduce the directions with a closed form and leave those of
  level `"open"` in place; otherwise such a direction is an error. A
  direction of level `"verified"` or `"checked"`, and a section claim
  left unproven, give a warning. Defaults to `FALSE`.

- allowChecked:

  Logical. Reduce also on a direction of level `"checked"`, which holds
  at random points but has no proof, and on a result whose saturation
  check failed; otherwise either is an error. Defaults to `FALSE`.

- check:

  Logical. Run the detection again on the model with the trafo applied
  and record in `checks$loop` whether exactly the remaining directions
  are left. Defaults to `TRUE`; large models may turn it off.

- workCap:

  Number of candidate sections the search of one block may try; reaching
  it is an entry in `limits`. Defaults to 4000.

- ...:

  Not used.

## Value

An object of class `symmetryreduction`:

- `blocks`:

  one entry per set of coupled directions, with `invariants`, `section`,
  `newParameters` (the invariant each `q_<k>` holds), `solvedDomain`,
  `coverage` and `status` (`"fixed"`, `"reduced"`, `"invariant_only"` or
  `"unresolved"`), `level`, `generatorLevel` (the lowest level of the
  directions in the block), `unproven` (the claims of the section that
  are not proven) and, for a curved block, `count` (`"exact"` or
  `"upper bound"`), `injective` and `reach`. A section on a face
  (`face`) holds `faceWhy`, why the model is not shown continuous up to
  the face (`NULL` when it is), and `faceLevel`, `"checked"` when the
  resting state under `equilibrate` is checked at random points only.

- `level`:

  how far the reduction is established, the lowest over the reduced
  blocks: `"exact"`, `"verified"` or `"checked"` (a direction rests on a
  closed form verified with an error bound, or checked without one, or a
  face on a checked resting state) or `"open"` (a claim of a section is
  not proven, named under `unproven`).

- `assumptions`:

  what the reduction assumes: `domain` (the declared positive
  coordinates) and `starts` (the sides of switches declared at the
  start), with `errorBound`, the probability that a verified closed form
  is wrong.

- `limits`:

  the work budgets and time limits that stopped a section or invariant
  search.

- `checks`:

  `pins` and `invariants`, whether every pin holds under the trafo and
  every invariant of a block maps to the matching new parameter over Q,
  `failed`, and with `check = TRUE` `loop`: the `directions` the reduced
  model has, the number `expected`, `rank`, `dim` and `rankStatus`, and
  `ok`.

- `zeroCompatibility`:

  with `reportZeroCompatibility = TRUE`, the sets of coordinates that
  can vanish together and the `condition` for it.

- `trafo`:

  an [eqnvec](https://dmodverse.github.io/dMod2/reference/eqnvec.md)
  over all coordinates for
  [`P()`](https://dmodverse.github.io/dMod2/reference/P.md) or
  `symmetryDetection(trafo = )`; `NULL` if nothing was reduced.

- `family`:

  the valid gauges of each reduced block.

- `removed`, `remaining`:

  labels of the directions.

- `dependent`:

  labels of directions that depend on the others and need no reduction.

- `partial`:

  labels of removed directions whose chart covers only part of the
  positive orthant.

- `alternatives`:

  the other charts, each with `trafo`, `removed`, `remaining`,
  `rational` and `pins`.

- `coordinates`, `fixed`, `settings`, `call`:

  provenance.

[`print()`](https://rdrr.io/r/base/print.html) shows the verdict, the
trafo, the new parameters and the sections,
[`summary()`](https://rdrr.io/r/base/summary.html) also one line per
block. Both take `width`,
[`summary()`](https://rdrr.io/r/base/summary.html) also `verbose` and
`detailed`, which explains each section.

## Details

A chart is returned only if it is valid on the domain declared by
`positive`. A `q_<k>` that takes both signs is marked `[real-valued]`
and is estimated on a linear scale; a chart valid only for positive
values of the `q_<k>` has `coverage = "partial"`.
[`print()`](https://rdrr.io/r/base/print.html) lists the invariant
behind each `q_<k>`.

## See also

[`symmetryDetection()`](https://dmodverse.github.io/dMod2/reference/symmetryDetection.md),
[`P()`](https://dmodverse.github.io/dMod2/reference/P.md)

## Examples

``` r
f <- eqnvec(m = "ktx - dm*m", p = "ktl*m - dp*p")
g <- eqnvec(y = "p")
res <- symmetryDetection(f, g, reconstruct = TRUE, verbose = FALSE)
red <- symmetryReduction(res)
red
#> Reduced 1 of 1 direction.
#> 
#> Trafo (non-identity entries)
#>   ktl = 1
#> 
#> Sections:
#>   {X₁}  ktl = 1
#>   summary(x, detailed = TRUE) explains each section.
red$trafo
#> Idx Inner <- Outer
#>   1    dm <- dm
#>   2    dp <- dp
#>   3   ktl <- 1
#>   4   ktx <- ktx
#>   5     m <- m
#>   6     p <- p
symmetryDetection(f, g, trafo = red$trafo, verbose = FALSE)$identifiable
#> [1] TRUE
```
