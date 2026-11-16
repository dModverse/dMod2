# Structural Non-Identifiabilities of an ODE Model

Finds the directions in parameters and initial values along which the
observations of a model do not change: the nullspace of the
observability matrix over finite fields, with exact closed forms and a
level that says whether each one is proven. The model `f`, the
observables `g` and an optional `trafo` are equations in any form
[`as.eqnvec()`](https://dmodverse.github.io/dMod2/reference/as.eqnvec.md)
accepts. The computation runs in the Python package `symident` through
`reticulate`, installed from PyPI on first use.

Right-hand sides, observables and initial values may contain
[`exp()`](https://rdrr.io/r/base/Log.html), `exp10()`, `b^x`, hyperbolic
and trigonometric functions, free exponents `x^n`, and
[`log()`](https://rdrr.io/r/base/Log.html), fractional powers and
[`abs()`](https://rdrr.io/r/base/MathFun.html) of `positive`
coordinates. `time` is the model time.

## Usage

``` r
symmetryDetection(
  f,
  g,
  trafo = NULL,
  parameters = NULL,
  fixed = NULL,
  gaugePreference = FALSE,
  forcings = NULL,
  unknownInputs = NULL,
  knownInputs = NULL,
  events = NULL,
  conditions = NULL,
  regimes = NULL,
  equilibrate = FALSE,
  reduceCQ = FALSE,
  freeInitial = NULL,
  reconstruct = FALSE,
  positive = TRUE,
  verify = TRUE,
  cores = 1,
  control = reconstControl(),
  scalingsOnly = FALSE,
  symEngine = c("modular", "symbolic"),
  discrete = FALSE,
  verbose = FALSE
)
```

## Arguments

- f:

  Right-hand sides: an
  [eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md), an
  [eqnvec](https://dmodverse.github.io/dMod2/reference/eqnvec.md) or a
  named character vector. `reduceCQ` and `equilibrate` use the reactions
  of an
  [eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md).

- g:

  Observables: one
  [eqnvec](https://dmodverse.github.io/dMod2/reference/eqnvec.md) or
  named character vector for all conditions, or a list with one per row
  of `conditions`. An observable `a*log(h) + c` with a number `a` is
  analyzed through `h`.

- trafo:

  Parameter transformation: one
  [eqnvec](https://dmodverse.github.io/dMod2/reference/eqnvec.md) for
  all conditions or a list with one per condition. An entry named like a
  parameter is substituted into `f`, `g`, event values and event times;
  an entry named like a state sets the initial value of the state. On
  the right-hand side a state name is the initial value of that state. A
  parameter that enters only through `k = "exp(lk)"`, `"exp10(lk)"` or
  `"10^lk"` is analyzed in `b^lk` and reported in `lk`.

- parameters:

  Character vector of additional symbols treated as parameters.

- fixed:

  Character vector of known symbols: a fixed parameter is a known
  constant, a fixed state has no unknown initial value.

- gaugePreference:

  `FALSE` (default) reconstructs every direction. Otherwise one
  coordinate per scaling is fixed first and only the general directions
  are reconstructed in that gauge. `NULL` chooses the gauge that leaves
  the general directions on the fewest coordinates; a character vector
  ranks the coordinates to fix (`*` as wildcard). The gauge is returned
  as `$gauge`.

- forcings:

  Character vector of input states, 0 at the start and at rest.

- unknownInputs:

  Character vector of unknown functions of time of no assumed form:
  symbols in `f` and `g` without an equation. A direction may change
  them along with the parameters; an exact one lists that change as
  `inputChange`. A character cell of `conditions` renames an unknown
  input in that condition.

- knownInputs:

  Character vector of known functions of time of no assumed form:
  symbols in `f` and `g` without an equation whose values and time
  derivatives are known.

- events:

  An
  [eventlist](https://dmodverse.github.io/dMod2/reference/eventlist.md)
  with event times as numbers or expressions in the parameters, in
  chronological order. Root events are not supported. An event value
  named like a column of `conditions` is read from it.

- conditions:

  Data frame with one row per condition and columns named by model
  symbols or event values. A numeric cell fixes the symbol, a character
  cell renames it. A column named like a dynamic state sets the initial
  value of the state. A direction is reported only if it is
  non-identifiable in every condition.

- regimes:

  For a model with switches, the regimes an experiment visits: a
  character vector of comparisons such as `"x > c"`, each selecting the
  regimes on that side of a switch, or a list of such vectors whose
  comparisons hold together; the first entry is the regime at the start.
  A list named by the rows of `conditions` gives them per condition.
  Without it every regime counts; a switch on parameters alone needs the
  side of the switch named here.

- equilibrate:

  Logical. Start at a steady state of `f` with the inputs at 0; the
  earliest events apply on top. Coupled steady states may need msolve;
  when it is missing, an interactive session asks whether to build it
  into the symident cache (on Windows with Rtools). Defaults to `FALSE`.

- reduceCQ:

  Logical,
  [eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md)
  only. Report the freedom of a conserved moiety on the total of the
  moiety (`TRUE`) or on the initial value of one species (`FALSE`,
  default).

- freeInitial:

  States, at most one per conserved quantity, holding the free resting
  value of their moiety under `equilibrate = TRUE, reduceCQ = FALSE`.

- reconstruct:

  Logical. Return general directions as exact rational functions instead
  of their support. Defaults to `FALSE`.

- positive:

  Coordinates known to be positive: `TRUE` (default, all), `FALSE`
  (none) or a character vector.

- verify:

  Logical. Check that the rank has saturated where the Lie order is not
  proven sufficient; the result is in `$info$verification`. Defaults to
  `TRUE`.

- cores:

  Number of threads and worker processes. Defaults to 1.

- control:

  A
  [`reconstControl()`](https://dmodverse.github.io/dMod2/reference/reconstControl.md)
  list.

- scalingsOnly:

  Logical. Only the scaling symmetries, from an exact integer kernel;
  `identifiable`, `rank` and `dim` are `NA`. Defaults to `FALSE`.

- symEngine:

  `"modular"` (default, finite fields) or `"symbolic"` (sympy, small
  models without `equilibrate` and later events).

- discrete:

  `FALSE` (default), `TRUE` or a list of arguments of
  [`discreteSymmetries()`](https://dmodverse.github.io/dMod2/reference/discreteSymmetries.md):
  also search the discrete symmetries.

- verbose:

  Logical. Print the result. Defaults to `FALSE`.

## Value

An object of class `symmetrydetection`:

- `identifiable`:

  `TRUE` or `FALSE`, `NA` with `scalingsOnly`; `FALSE` holds as far as
  `info$rankStatus` says.

- `rank`, `dim`:

  rank of the observability matrix and number of coordinates.

- `symmetries`:

  the directions, each a generator \\X = \sum_i \xi_i
  \partial\_{\theta_i}\\: `generator` (\\\xi\\ by coordinate), `weights`
  (a scaling's integer weights), `type` (`"scaling"` or `"general"`),
  `degree`, `support`, `explicit`, `level`, `stateComponents`, `reason`,
  `route`, `display`, `completeGenerator` and `factor`. `level` is
  `"exact"` for a direction proven over the rationals, `"verified"` for
  one checked at random points modulo a large prime with an error bound,
  `"checked"` for one checked without an error bound (an unproven Lie or
  gap order, or a model outside the scope of the bound) and `"open"`
  without closed form; `stateComponents` holds, for an exact direction,
  the components on the states of the field that proves it.

- `info`:

  `engine`, the Lie order and whether it is proven sufficient
  (`lieOrderUsed`, `lieOrderDriver`, `liePlateau`, `lieCertified`), the
  gap order of the events after the start and whether it is proven
  sufficient (`gapOrderUsed`, `gapCertified`, `NULL` without such
  events, and `gapCertificate`, why not), `rankProven`, `rankProof`,
  `rankFields` (a rank the invariant fields establish, checked at random
  points modulo a prime), `rankFieldsBound` (the probability that a
  wrong field passes), `rankStatus` (`status` `"proven"`, `"verified"`
  or `"lower bound"`, `why` and `probability`), `proof` (why each
  `"verified"` or `"checked"` direction is not proven), `assumptions`
  (`domain`, the declared positive coordinates, `starts`, the sides of
  switches declared at the start, and `errorBound`, the probability that
  a verified closed form or rank is wrong), `bound` (the probability
  that a wrong closed form passes the check, only where both orders are
  proven sufficient or the rank is settled), `limits` (what stopped part
  of the analysis), `errors` (exceptions turned into failed routes),
  `branches` (resting states with several roots modulo a prime, of which
  the first is taken), `conditions`, `segments`, `coordinates`,
  `settings`, `elapsed` and `verification`.

- `gauge`:

  the coordinates fixed under `gaugePreference`.

- `discreteSymmetries`:

  with `discrete`, see
  [`discreteSymmetries()`](https://dmodverse.github.io/dMod2/reference/discreteSymmetries.md).

- `call`:

  the matched call.

[`print()`](https://rdrr.io/r/base/print.html) shows the verdict and the
generators, [`summary()`](https://rdrr.io/r/base/summary.html) also the
computation. Both take `verbose` and `width`.

## Details

The analysis starts at the earliest event; later events split the time
line into segments. Directions are reported in the original coordinates.

Every environment variable `SYMIDENT_<NAME>` set in R, also with
[`Sys.setenv()`](https://rdrr.io/r/base/Sys.setenv.html) during the
session, is passed to symident for the call. The environment variables a
user may need are `SYMIDENT_MSOLVE`, the path of the msolve executable
(unset by default, msolve is then searched in the symident cache and on
the `PATH`), and `SYMIDENT_CACHE`, the directory of that cache (by
default `~/.cache/symident` on Linux). The other switches are
diagnostic, see the symident documentation at
<https://dmodverse.github.io/symident/>.

## See also

[`symmetryReduction()`](https://dmodverse.github.io/dMod2/reference/symmetryReduction.md),
[`reconstControl()`](https://dmodverse.github.io/dMod2/reference/reconstControl.md);
the guide chapter
[Symmetries](https://dmodverse.github.io/dMod2/articles/guide-symmetries.html).

## Examples

``` r
f <- eqnvec(m = "ktx - dm*m", p = "ktl*m - dp*p")
g <- eqnvec(y = "p")
res <- symmetryDetection(f, g, reconstruct = TRUE)
res
#> Result:  rank 5 / 6  |  1 non-identifiable direction
#> 
#> Generators  X = Σᵢ ξ(i) ∂ᵢ
#> 
#> Scalings:
#>   X₁  (exact)
#>       ξ(ktl) :  + ktl
#>       ξ(ktx) :  - ktx
#>       ξ(m)   :  - m
#> 
#> exact = proven over Q.
#> 
#> 
#> Discrete symmetries: 0
summary(res)
#> ------------------------------------------------------------
#> detect  |  observability   engine: modular (GF(p) + CRT)
#> ------------------------------------------------------------
#> Computation:
#>   Lie order 4 (gap order 0), set by condition 1
#>   rank: proven (1 directions, proven symmetries)
#>   settings: reduce_cq=False, equilibrate=False, reconstruct=True
#> 
#> Result:  rank 5 / 6  |  1 non-identifiable direction
#> 
#> Generators  X = Σᵢ ξ(i) ∂ᵢ
#> 
#> Scalings:
#>   X₁  (exact)
#>       ξ(ktl) :  + ktl
#>       ξ(ktx) :  - ktx
#>       ξ(m)   :  - m
#> 
#> exact = proven over Q.
#> 
#> Reduction:
#>   X₁ is a scaling: gauge it by holding any one of the coordinates fixed, at
#>   any value.
res$identifiable
#> [1] FALSE
```
