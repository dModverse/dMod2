# Analytical Steady States

Computes symbolic steady-state expressions for use in a parameter
transformation, following the AlyssaPetit method \[1\], \[2\]. States
whose mass only leaks are set to zero first, with `zeroLoops = TRUE`
also states that are only produced through one another. States whose own
balance gives a ratio of positive sums are solved from it first. The
cycles of the graph of states left are cut then: at a state of a
conserved quantity, or at a state whose balance is solved for a rate
constant of one flux side instead. Each cut state stays a free
parameter. A state of a conserved quantity whose balance is left without
a positive unknown is solved in a further round, with another state of
the quantity free. The remaining states and the rate constants of the
cuts are then solved, each from a balance or a combination of balances
in amounts, such that the result is a ratio of sums of positive terms.
Only if that fails, a flux side with separate rate constants is cut as
well: the fluxes of the side keep free ratios, the new parameters
`r_<state>_<j>`. Compartment volumes are fixed. The solver is a Python
module called through reticulate; it needs numpy, sympy and scipy.

## Usage

``` r
steadyStates(
  model,
  file = NULL,
  forcings = NULL,
  givenCQs = NULL,
  neglect = NULL,
  testSteady = c("fast", "skip"),
  walltime = 0L,
  simplify = TRUE,
  solveQuadratic = FALSE,
  priority = NULL,
  zeroLoops = FALSE,
  resolve = TRUE,
  verbose = FALSE
)
```

## Arguments

- model:

  An [eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md),
  or the path of a csv file describing the model.

- file:

  Character, path of the RDS file the result is written to with
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html). For an `eqnlist`
  the model is also written to `<file>_model.csv` for the solver.
  Default `NULL` writes nothing; the solver then reads the model from a
  temporary file.

- forcings:

  Character vector, names of forcings. These states are held at zero and
  treated as exogenous. Default `NULL`.

- givenCQs:

  Unnamed character vector of conserved quantities, either as
  `c("A + pA = totA", "B + pB = totB")` or as `c("A + pA", "B + pB")`.
  Default `NULL`: a basis is derived, or taken from
  [`customTotals()`](https://dmodverse.github.io/dMod2/reference/customTotals.md)
  when the `eqnlist` has them. Each given quantity keeps one of the
  states in it free, the first one on a cycle of the graph of states.

- neglect:

  Character vector, states and rate parameters the solver must not solve
  for. Default `NULL`.

- testSteady:

  How the solution is verified: `"fast"` (default), a probabilistic
  check modulo a prime, or `"skip"`.

- walltime:

  Integer, time budget of the solver in seconds, `0` (default) for none.

- simplify:

  `TRUE` (default) simplifies every expression once with sympy, `FALSE`
  skips it, `"full"` simplifies further at a higher cost.

- solveQuadratic:

  Logical, whether the balance of a state, when quadratic in that state,
  may be solved for the positive root instead of for a rate constant,
  which introduces [`sqrt()`](https://rdrr.io/r/base/MathFun.html)
  terms. Default `FALSE`.

- priority:

  Character vector of state and rate-parameter names, most preferred
  first: named states are solved rather than kept free and solved
  earlier, named rate parameters are preferred as pivots. Default
  `NULL`.

- zeroLoops:

  Logical. `FALSE` (default) looks for a steady state with every state
  away from zero. `TRUE` sets to zero the largest set of states that are
  only produced through one another, such as a ligand that only the
  receptor complexes it forms activate: all fluxes of these states
  vanish where they are zero, and no conserved quantity lies on them.
  This is the steady state of a feedback loop that is off, as under
  `symmetryDetection(equilibrate = TRUE)`.

- resolve:

  Logical, whether the result is passed through
  [`resolveRecurrence()`](https://dmodverse.github.io/dMod2/reference/resolveRecurrence.md)
  so that no equation refers to another. Default `TRUE`.

- verbose:

  `TRUE` reports progress in a few lines and the result, `FALSE`
  (default) only the result, `"full"` every step of the solver.

## Value

Named character vector of steady-state equations in dMod format, or `0`
if no solution was found. An entry whose value equals the entry name is
a free parameter.

## References

\[1\] Rosenblatt, M., Timmer, J. and Kaschek, D. (2016). Customized
steady-state constraints for parameter estimation in non-linear ordinary
differential equation models. *Frontiers in Cell and Developmental
Biology* 4, 41.
[doi:10.3389/fcell.2016.00041](https://doi.org/10.3389/fcell.2016.00041)

\[2\] <https://github.com/marcusrosenblatt/AlyssaPetit>

## See also

[`define()`](https://dmodverse.github.io/dMod2/reference/define.md),
[`P()`](https://dmodverse.github.io/dMod2/reference/P.md)

## Examples

``` r
reactions <- eqnlist()
reactions <- addReaction(reactions, "Tca_buffer", "Tca_cyto",
                         "import_Tca*Tca_buffer", "Basolateral uptake")
reactions <- addReaction(reactions, "Tca_cyto", "Tca_buffer",
                         "export_Tca_baso*Tca_cyto", "Basolateral efflux")
reactions <- addReaction(reactions, "Tca_cyto", "Tca_canalicular",
                         "export_Tca_cana*Tca_cyto", "Canalicular efflux")
reactions <- addReaction(reactions, "Tca_canalicular", "Tca_buffer",
                         "transport_Tca*Tca_canalicular", "Transport bile")

steadies <- steadyStates(reactions)
steadies
#>                                                                                    Tca_cyto 
#>                                 "Tca_buffer*import_Tca/(export_Tca_baso + export_Tca_cana)" 
#>                                                                             Tca_canalicular 
#> "Tca_buffer*export_Tca_cana*import_Tca/(transport_Tca*(export_Tca_baso + export_Tca_cana))" 
#>                                                                                  Tca_buffer 
#>                                                                                "Tca_buffer" 

# Parameter transformation that puts the model in steady state
parameters <- getParameters(reactions)
trafo <- define(NULL, "x <- x", x = parameters)
trafo <- insert(trafo, "x <- y", x = names(steadies), y = steadies)
trafo
#> Idx           Inner <- Outer
#>   1      Tca_buffer <- Tca_buffer
#>   2 Tca_canalicular <- Tca_buffer*export_Tca_cana*import_Tca/(transport_Tca*(export_Tca_baso+export_Tca_cana))
#>   3        Tca_cyto <- Tca_buffer*import_Tca/(export_Tca_baso+export_Tca_cana)
#>   4 export_Tca_baso <- export_Tca_baso
#>   5 export_Tca_cana <- export_Tca_cana
#>   6      import_Tca <- import_Tca
#>   7   transport_Tca <- transport_Tca
```
