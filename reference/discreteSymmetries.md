# Discrete Symmetries

Searches the discrete symmetries of a
[`symmetryDetection()`](https://dmodverse.github.io/dMod2/reference/symmetryDetection.md)
result: other points of the domain with the output of a given one. The
Taylor coefficients of the observables up to `order` are equated with
their values at a point; msolve counts the solutions over the algebraic
closure modulo a prime, and over the reals and in the declared domain at
`points` random rational points, where the number in the domain may
differ between regions. Where the other points are rational functions
\\\sigma\\ of the coordinates, \\\sigma\\ is reconstructed, checked on
the output at random points (`"checked"`) and, where every state has a
free initial value, proven by a map of the states (`"exact"`). A model
with directions is reduced by
[`symmetryReduction()`](https://dmodverse.github.io/dMod2/reference/symmetryReduction.md)
first and the reduced model is analyzed.

## Usage

``` r
discreteSymmetries(
  object,
  order = NULL,
  points = 5L,
  timeout = 600,
  budget = 300,
  cores = 1
)
```

## Arguments

- object:

  A `symmetrydetection` result from
  [`symmetryDetection()`](https://dmodverse.github.io/dMod2/reference/symmetryDetection.md)
  of a model with one condition and no events, switches or unknown
  inputs, detected without `equilibrate` and `reduceCQ`.

- order:

  Highest Taylor coefficient; `NULL` (default) twice the number of
  coordinates.

- points:

  Number of random points of the count in the domain, or a list of
  points, each a named list of the coordinates. Defaults to 5.

- timeout:

  Seconds one call of msolve may take. Defaults to 600.

- budget:

  Seconds the reconstruction of the maps may take; past it, the counts
  are reported without closed forms. Defaults to 300.

- cores:

  Threads of msolve. Defaults to 1.

## Value

`object` with `discreteSymmetries`, per \\\sigma\\ the `map` of the
coordinates it changes (an
[eqnvec](https://dmodverse.github.io/dMod2/reference/eqnvec.md)),
`changes`, `level`, `preservesDomain` (`TRUE`, `FALSE`, or `NULL` with
the `domainConditions` under which \\\sigma(\theta)\\ lies in the
domain) and `fixedSet`, the factors on which \\\sigma\\ fixes
\\\theta\\; and `info$discrete` with `order`, `count` (the solutions
over the algebraic closure at a random point modulo a prime), `counts`
per point (`complex`, `real`, `inDomain`, `undecided`), `closedForms`,
`countFinal` (every solution is the image of a proven \\\sigma\\),
`fundamentalDomain` (inequalities that keep one point of each orbit),
`limits` and, for a model with directions, `reduced`.

## Details

The search needs msolve;
[`symmetryDetection()`](https://dmodverse.github.io/dMod2/reference/symmetryDetection.md)
says where it is searched.

## References

Godfrey, K. (1983). Compartmental Models and Their Application. Academic
Press, London.

## See also

[`symmetryDetection()`](https://dmodverse.github.io/dMod2/reference/symmetryDetection.md),
[`symmetryReduction()`](https://dmodverse.github.io/dMod2/reference/symmetryReduction.md)

## Examples

``` r
if (FALSE) { # requireNamespace("reticulate", quietly = TRUE) && reticulate::py_module_available("symident.install") && !is.null(tryCatch(reticulate::import("symident.install")$find_msolve()[[1]], error = function(e) NULL))
# two compartments, an input u into the first, the second observed
f <- eqnvec(x1 = "-(k21 + k01)*x1 + k12*x2 + u", x2 = "k21*x1 - k12*x2")
res <- symmetryDetection(f, eqnvec(y = "x2"), knownInputs = "u", discrete = TRUE)
res
res$discreteSymmetries[[1]]$map
res$info$discrete$fundamentalDomain
}
```
