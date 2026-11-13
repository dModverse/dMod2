# Reparameterization

Replaces the left-hand side of `"lhs <- rhs"` by the right-hand side in
`trafo`, after the symbols of `expr` are replaced by the arguments in
`...`.
[`define()`](https://dmodverse.github.io/dMod2/reference/define.md) and
[`insert()`](https://dmodverse.github.io/dMod2/reference/define.md)
build on it.

## Usage

``` r
repar(trafo = NULL, expr, ..., reset = FALSE)
```

## Arguments

- trafo:

  Named character,
  [eqnvec](https://dmodverse.github.io/dMod2/reference/eqnvec.md) or a
  list of those. Default `NULL` starts from the identity of the
  left-hand sides, so the result maps each left-hand side to the
  matching right-hand side.

- expr:

  `"lhs <- rhs"` string. Inside a name, `_` separates symbols, so
  `"Delta_x"` with `x = "a"` becomes `"Delta_a"`. A `"lhs ~ rhs"` string
  or a formula is an error: `~` is reserved for model formulas.

- ...:

  Named character or numeric vectors of equal length, replacing the
  symbols of `expr` of the same name; element `i` of each gives the
  `i`-th substitution.

- reset:

  `FALSE` (default) replaces the left-hand side wherever it occurs in
  the expressions, `TRUE` sets the entry it names.

## Value

Same shape as `trafo`.

## See also

[`define()`](https://dmodverse.github.io/dMod2/reference/define.md),
[`insert()`](https://dmodverse.github.io/dMod2/reference/define.md)

## Examples

``` r
innerpars   <- letters[1:3]
constraints <- c(a = "b + c")
mycondition <- "cond1"
trafo <- repar(NULL,  "x <- x",      x = innerpars)
trafo <- repar(trafo, "x <- y",      x = names(constraints), y = constraints)
trafo <- repar(trafo, "x <- exp(x)", x = innerpars)
trafo <- repar(trafo, "x <- x + Delta_x_condition",
               x = innerpars, condition = mycondition)
```
