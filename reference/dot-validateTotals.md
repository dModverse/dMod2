# Validate Conservation Expressions Against an Eqnlist

Each expression must be linear in `eqnlist$states` and conserve an
amount: the coefficients divided by the volumes of their states lie in
the null space of the stoichiometric matrix. Symbols other than states,
volumes among them, are checked at two sets of values. Together the
expressions must have the rank of the detected basis.

## Usage

``` r
.validateTotals(totals, eqnlist, tol = 1e-08)
```

## Arguments

- totals:

  Named list of conservation expressions.

- eqnlist:

  The [eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md)
  to validate against.

- tol:

  Numerical tolerance for the null-space check.

## Value

`TRUE` on success; otherwise
[`stop()`](https://rdrr.io/r/base/stop.html)s with a descriptive error.
