# Detect states that are structurally zero in steady state

Iterated three-layer test on the stoichiometric matrix of an `eqnlist`,
the one the steady-state solver of
[`steadyStates()`](https://dmodverse.github.io/dMod2/reference/steadyStates.md)
applies:

1.  *Neg-only column*: only outflux, must be zero at SS.

2.  *Pos-only column with single-state feeder*: each feeding flux
    vanishes at SS; if a feeder's rate involves exactly one state, that
    state must be zero.

3.  *Sink cluster (LP)*: subsets whose combined mass leaks monotonically
    (mass-balance LP via
    [`lpSolve::lp`](https://rdrr.io/pkg/lpSolve/man/lp.html)).

For each zero-state the column is dropped, the state symbol is
substituted by `"0"` in remaining rates, structurally-zero reactions are
removed, and detection re-runs (removals can expose new zero-states). A
state that some reaction changes but none is left to change, outside
every conserved moiety, rests: the steady state leaves the value of such
a state free.

## Usage

``` r
.zeroStatesFromSmatrix(eqnlist_obj)
```

## Arguments

- eqnlist_obj:

  An [eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md).

## Value

`list(zero_states, resting, eqnlist)` with the reduced system.
