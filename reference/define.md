# Construct and Modify Parameter Transformations

Build, change and branch the symbolic transformations that
[`P()`](https://dmodverse.github.io/dMod2/reference/P.md) turns into a
[parfn](https://dmodverse.github.io/dMod2/reference/parfn.md).

- `define()` sets the entry named by the left-hand side of `expr` to the
  right-hand side, adding it if missing.

- `insert()` replaces the left-hand side by the right-hand side wherever
  it occurs in the expressions.

- `branch()` copies `trafo` once per condition and keeps `table` as
  attribute `tree`.

## Usage

``` r
define(trafo, expr, ..., conditionMatch = NULL)

insert(trafo, expr, ..., conditionMatch = NULL)

branch(
  trafo,
  table = NULL,
  conditions = rownames(table),
  apply = c("nothing", "insert", "define")
)
```

## Arguments

- trafo:

  Named character or
  [eqnvec](https://dmodverse.github.io/dMod2/reference/eqnvec.md), or a
  list of those named by condition as `branch()` returns it. `NULL`
  starts a new one.

- expr:

  `"lhs <- rhs"` string. The symbols in it are replaced by the arguments
  in `...` of the same name; see
  [`repar()`](https://dmodverse.github.io/dMod2/reference/repar.md). For
  `insert()` on a branched transformation also one or several effect
  formulas, e.g. `cl ~ wt + (1 | ID)`, which replace the symbol by a
  linear predictor over the columns of the branch table with fixed and
  random effects; see
  [`randomEffects()`](https://dmodverse.github.io/dMod2/reference/randomEffects.md).
  Named arguments in `...` are then placeholders of the formulas,
  `insert(trafo, x ~ 1 + (1 | ID), x = c("ka", "cl"))`.

- ...:

  Named vectors replacing the symbols of `expr` of the same name, one
  substitution per element. They are evaluated per condition, with the
  columns of `attr(trafo, "tree")`, `.currentTrafo` (the transformation
  of the condition) and `.currentSymbols` (the symbols of that
  transformation) in scope. An unnamed logical argument of `insert()`
  selects the conditions where it is `TRUE`. Next to effect formulas
  they are more formulas and placeholders, evaluated once in the calling
  environment.

- conditionMatch:

  Regular expression; only the conditions whose names match are changed.
  Default `NULL`, all.

- table:

  `data.frame` with one row per condition, row names the conditions.
  Kept as attribute `tree`; the columns can be read in `...` by name.

- conditions:

  Condition names; default `rownames(table)`.

- apply:

  One of `"nothing"` (default), `"insert"` and `"define"`: what is done
  with every non-missing entry `table[condition, column]`, with the
  column name as left-hand side and the entry as right-hand side.

## Value

Same shape as `trafo`, or a list named by condition for `branch()`.

## See also

[`repar()`](https://dmodverse.github.io/dMod2/reference/repar.md),
[`P()`](https://dmodverse.github.io/dMod2/reference/P.md)

## Examples

``` r
## Identity for four parameters, then log parameters
trafo <- define(NULL, "x <- x", x = c("A", "B", "k1", "k2"))
trafo <- insert(trafo, "x <- exp(x)", x = .currentSymbols)
trafo
#> Idx Inner <- Outer
#>   1     A <- exp(A)
#>   2     B <- exp(B)
#>   3    k1 <- exp(k1)
#>   4    k2 <- exp(k2)

## Replace k1 and k2 by expressions, fix B
trafo <- insert(trafo, "x <- y", x = c("k1", "k2"), y = c("q1 + q2", "q1 - q2"))
trafo <- define(trafo, "B <- 0")
trafo
#> Idx Inner <- Outer
#>   1     A <- exp(A)
#>   2     B <- 0
#>   3    k1 <- exp(q1+q2)
#>   4    k2 <- exp(q1-q2)

## One branch per condition; the non-missing table entries are inserted
conditions <- data.frame(dose = c(1, 10), q2 = c(NA, "q2_high"),
                         row.names = c("low", "high"))
trafoL <- branch(trafo, table = conditions, apply = "insert")

## Condition-specific changes, selected by name or read off the table
trafoL <- insert(trafoL, "x <- 2*x", x = "q1", conditionMatch = "^high")
trafoL <- define(trafoL, "A <- dose", dose = dose)
trafoL
#> $low
#> Idx Inner <- Outer
#>   1     A <- 1
#>   2     B <- 0
#>   3    k1 <- exp(q1+q2)
#>   4    k2 <- exp(q1-q2)
#> 
#> $high
#> Idx Inner <- Outer
#>   1     A <- 10
#>   2     B <- 0
#>   3    k1 <- exp((2*q1)+q2_high)
#>   4    k2 <- exp((2*q1)-q2_high)
#> 
#> attr(,"tree")
#>      dose      q2
#> low     1    <NA>
#> high   10 q2_high
```
