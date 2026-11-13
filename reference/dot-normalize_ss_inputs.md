# Normalize steady-state inputs

Coerces `trafo` to an `eqnvec`, zeroes and removes forcings, removes
states that are structurally zero at steady state and promotes states
with a zero right-hand side to parameters; the resting ones among them
default to 0.

## Usage

``` r
.normalize_ss_inputs(trafo, parameters, forcings)
```

## Arguments

- trafo:

  An [eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md),
  [eqnvec](https://dmodverse.github.io/dMod2/reference/eqnvec.md) or
  named character vector.

- parameters:

  Outer parameters (may be `NULL`).

- forcings:

  Forcing names.

## Value

A named list with `trafo`, `states`, `zero_states`, `resting`,
`dependent`, `parameters`, `parms_all`, `totals` and, for an
[eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md),
`influx`, the sum of the producing fluxes per state.
