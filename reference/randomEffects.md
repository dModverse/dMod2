# Random effects of a parameter transformation

Returns the random-effect structure that
[`insert()`](https://dmodverse.github.io/dMod2/reference/define.md) with
effect formulas attached to a branched transformation: the random terms,
their groupings, distributions and correlation blocks, and the
parameters the formulas created.
[`em()`](https://dmodverse.github.io/dMod2/reference/em.md) takes it as
`random` and builds the prior of the random effects from it.
[`P()`](https://dmodverse.github.io/dMod2/reference/P.md) keeps the
structure, so `randomEffects()` takes the transformation or the
parameter function built from it.

An effect formula `cl ~ wt + (1 + dose | p | ID)` replaces the symbol
`cl` wherever it occurs in a branch by a linear predictor. The scale of
the symbol is the one the transformation gives it: after
`insert(trafo, "x <- exp(x)", x = "cl")` the random effects are normal
on the log scale. The parts of an effect formula:

- left-hand side:

  one or more symbols of the transformation, `ka + v ~ ...` applies the
  right-hand side to each.

- fixed effects:

  the terms outside the bars, coded as by
  [`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
  (treatment contrasts) over the columns of the branch table. The
  intercept keeps the name of the symbol (`cl`), the other columns are
  `cl_<column>` (`cl_wt`, `cl_log_wt`).

- random terms:

  `(1 | ID)`, `(1 + dose | ID)` with a random slope, `(1 | p | ID)` with
  an ID `p`: terms with the same ID and grouping are correlated across
  formulas, the columns of one term are correlated, other terms are
  independent. `gr(ID, dist = spikeSlab())` instead of `ID` sets the
  distribution, see
  [`normal()`](https://dmodverse.github.io/dMod2/reference/normal.md).
  The random effect of a group is `eta_<symbol>_<group>` for the
  intercept and `eta_<symbol>_<column>_<group>` for a slope; a symbol
  with random terms in a second grouping gets the grouping name in the
  stem, `eta_cl_center_<group>`.

## Usage

``` r
randomEffects(trafo)

# S3 method for class 'randomeffects'
print(x, ...)
```

## Arguments

- trafo:

  A branched transformation from
  [`insert()`](https://dmodverse.github.io/dMod2/reference/define.md)
  with effect formulas, or a parameter function
  [`P()`](https://dmodverse.github.io/dMod2/reference/P.md) built from
  one.

- x:

  A `randomeffects` object.

- ...:

  Ignored.

## Value

An object of class `randomeffects`, a list with `terms` (one row per
random effect column: `symbol`, `origin` (the formula), `term`,
`column`, `group`, `id`, `dist`, the fixed hyperparameters `prob`,
`concentration` and `tau` (`NA` where estimated) and `stem`, the name
stem of the etas), `blocks` (groups of rows: the normal rows of a block
are correlated, the rows of a spike-and-slab or Dirichlet-process block
share the hyperparameters), `subjects` (the groups of each grouping and
the conditions they own), `fixed` (the parameters the formulas created),
`etas` (eta names, one row per condition, one column per term row),
`groups` (per stem the group of each condition), `fixedDesign` (per
symbol the model matrix of the fixed effects, one row per condition) and
`etaNames` (`NULL` for random effects from
[`insert()`](https://dmodverse.github.io/dMod2/reference/define.md)).

## See also

[`insert()`](https://dmodverse.github.io/dMod2/reference/define.md),
[`em()`](https://dmodverse.github.io/dMod2/reference/em.md),
[`normal()`](https://dmodverse.github.io/dMod2/reference/normal.md)

## Examples

``` r
## Log parameters; clearance with a weight covariate and a random intercept and dose
## slope per subject, correlated with the random intercept of the absorption rate
cov <- data.frame(ID = c("s1", "s2", "s3"), wt = c(62, 75, 81), dose = c(100, 100, 200),
                  row.names = c("C1", "C2", "C3"))
trafo <- branch(define(NULL, "x <- x", x = c("ka", "cl", "v")), table = cov)
trafo <- insert(trafo, "x <- exp(x)", x = c("ka", "cl", "v"))
trafo <- insert(trafo,
                ka ~ 1 + (1 | p | ID),
                cl ~ wt + (1 + dose | p | ID),
                v  ~ 1 + (1 | gr(ID, dist = spikeSlab())))
trafo$C2
#> Idx Inner <- Outer
#>   1    cl <- exp((cl+cl_wt*75+eta_cl_s2+eta_cl_dose_s2*100))
#>   2    ka <- exp((ka+eta_ka_s2))
#>   3     v <- exp((v+eta_v_s2))
re <- randomEffects(trafo)
re
#> Random effects: 4 terms over ID
#>   block 1 (normal, ID): eta_ka, eta_cl, eta_cl_dose
#>   block 2 (spikeSlab, prob estimated, tau estimated, ID): eta_v
#>   fixed effects: ka, cl, cl_wt, v 
re$fixed
#> [1] "ka"    "cl"    "cl_wt" "v"    

## One random intercept for many parameters at once
randomEffects(insert(branch(define(NULL, "x <- x", x = c("ka", "cl", "v")), table = cov),
                     x ~ 1 + (1 | ID), x = c("ka", "cl", "v")))
#> Random effects: 3 terms over ID
#>   block 1 (normal, ID): eta_ka
#>   block 2 (normal, ID): eta_cl
#>   block 3 (normal, ID): eta_v
#>   fixed effects: ka, cl, v 
```
