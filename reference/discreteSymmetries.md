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
first and the reduced model is analyzed. The search stops on work
limits, never on a clock, so the result is the same on every machine; a
count or proof that reaches a limit ends without a result and names the
limit in `limits` or in the certificate.

## Usage

``` r
discreteSymmetries(
  object,
  order = NULL,
  points = 5L,
  cores = 1,
  memory = NULL,
  experiments = NULL,
  msolveEntries = NULL,
  msolveTerms = NULL,
  searchEntries = NULL,
  eliminationTerms = NULL
)
```

## Arguments

- object:

  A `symmetrydetection` result from
  [`symmetryDetection()`](https://dmodverse.github.io/dMod2/reference/symmetryDetection.md).

- order:

  Highest Taylor coefficient; `NULL` (default) twice the number of
  coordinates.

- points:

  Number of random points of the count in the domain, or a list of
  points, each a named list of the coordinates. Defaults to 5.

- cores:

  Threads of msolve. Defaults to 1.

- memory:

  GiB of address space one call of msolve may take; past it the call
  fails. `NULL` (default) keeps the value of symident (8).

- experiments:

  Count on several experiments with shared parameters and free starts,
  which bounds the fiber of a model with free starts and no rational
  first integral: `NULL` (default) where a short count of one experiment
  finds none, `TRUE` first with the number of experiments from the rank
  criterion, a number `K` with `K` experiments, `FALSE` never.

- msolveEntries:

  Nonzero entries of the F4 matrices one call of msolve may reduce.

- msolveTerms:

  Terms of the largest system given to msolve.

- searchEntries:

  Entries of the matrices over finite fields the exact count may reduce.

- eliminationTerms:

  Terms of the largest polynomial of the elimination of the states in
  the exact count.

## Value

`object` with `discreteSymmetries`, per \\\sigma\\ the `map` of the
coordinates it changes (an
[eqnvec](https://dmodverse.github.io/dMod2/reference/eqnvec.md)),
`changes`, `level`, `preservesDomain` (`TRUE`, `FALSE`, or `NULL` with
the `domainConditions` under which \\\sigma(\theta)\\ lies in the
domain) and `fixedSet`, the factors on which \\\sigma\\ fixes
\\\theta\\; and `info$discrete` with `order`, `count` (the solutions
over the algebraic closure at a random point modulo a prime),
`countBound` (the probability that a generic point has more, `NULL`
where the count is exact or has no bound), `countDomain` (the points of
the fiber of a generic point inside the domain, where that number is
final), `counts` per point (`complex`, `real`, `inDomain`, `undecided`),
`closedForms`, `countFinal` (every solution is the image of a proven
\\\sigma\\), `fundamentalDomain` (inequalities that keep one point of
each orbit), `certificate` (what the result rests on, with `complete`:
`"exact"`, `"verified"`, `"axiom"`, `"no"` or `"open"`), `experiments`
(the count on several experiments, or why it does not apply), `limits`
and, where they apply, `reduced` (for a model with directions),
`eliminated` (under `reduceCQ`, the species solved through the totals),
`relaxed` (the model with non-rational terms of the states as new
states) and `unsupported`.

## Details

The search needs msolve;
[`symmetryDetection()`](https://dmodverse.github.io/dMod2/reference/symmetryDetection.md)
says where it is searched. The defaults of `msolveEntries`,
`msolveTerms`, `searchEntries` and `eliminationTerms` (`NULL`) keep the
values of symident.

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
res <- symmetryDetection(f, eqnvec(y = "x2"), knownInputs = "u", mode = "global")
res
res$discreteSymmetries[[1]]$map
res$info$discrete$fundamentalDomain
}
```
