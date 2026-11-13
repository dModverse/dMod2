# Generate or Coerce to an Equation List

An eqnlist stores an ODE as a stoichiometric matrix, rate expressions,
state names and compartment information. `eqnlist()` builds one from
these parts, `as.eqnlist()` from a data frame of reactions.

## Usage

``` r
as.eqnlist(data, volumes, ...)

# S3 method for class 'data.frame'
as.eqnlist(
  data,
  volumes = NULL,
  compartments = NULL,
  compartmentOf = NULL,
  reactionCompartment = NULL,
  amountStates = NULL,
  ...
)

is.eqnlist(x)

eqnlist(
  smatrix = NULL,
  states = colnames(smatrix),
  rates = NULL,
  volumes = NULL,
  description = NULL,
  compartments = NULL,
  compartmentOf = NULL,
  reactionCompartment = NULL,
  amountStates = NULL,
  totals = NULL
)
```

## Arguments

- data:

  data frame with the columns "Description", "Rate" and one column per
  state.

- volumes:

  named character of state volumes. Without `compartments` and
  `compartmentOf`, each distinct expression becomes a compartment with
  ID `c1`, `c2`, ...

- ...:

  not used.

- compartments:

  named list keyed by compartment ID; each entry is a volume expression
  (character) or a list with fields `volume` and `rule` (`dV/dt`, see
  [`setCompartmentVolume()`](https://dmodverse.github.io/dMod2/reference/setCompartmentVolume.md)).

- compartmentOf:

  named character vector mapping states to compartment IDs. States not
  listed default to compartment `"defaultComp"` with volume `"1"`. May
  also name states that do not exist yet; see
  [`assignCompartment()`](https://dmodverse.github.io/dMod2/reference/assignCompartment.md).

- reactionCompartment:

  character of length `nrow(smatrix)`, the compartment in whose frame
  each rate is given; `NA` infers it from the educts. Required for
  reactions whose educts span several compartments. Default `NULL`.

- amountStates:

  character, the states given in amounts rather than concentrations
  (SBML's `hasOnlySubstanceUnits`). Their fluxes are not divided by a
  compartment volume. Default `NULL`.

- x:

  object to test.

- smatrix:

  numeric stoichiometric matrix, one row per reaction and one column per
  state. Default `NULL` gives an empty equation list.

- states:

  character, the state names. Default `colnames(smatrix)`.

- rates:

  character, the rate expressions, one per reaction.

- description:

  character, one description per reaction. Default `NULL` numbers the
  reactions.

- totals:

  named list of user-defined conservation expressions, see
  [`customTotals()`](https://dmodverse.github.io/dMod2/reference/customTotals.md).
  Default `NULL` uses the detected basis.

## Value

Object of class `eqnlist`, a list with the entries above and the
per-state `volumes`. `is.eqnlist()` returns a logical.

## Details

The data frame for `as.eqnlist()` contains the columns "Description"
(character), "Rate" (character) and one column per state, the columns of
the stoichiometric matrix. Layout attributes of the data frame, as set
by
[`as.data.frame.eqnlist()`](https://dmodverse.github.io/dMod2/reference/as.data.frame.eqnlist.md),
are used unless `volumes` is given.

## See also

[`addReaction()`](https://dmodverse.github.io/dMod2/reference/addReaction.md),
`as.eqnlist()`,
[`assignCompartment()`](https://dmodverse.github.io/dMod2/reference/assignCompartment.md),
[`setCompartmentVolume()`](https://dmodverse.github.io/dMod2/reference/setCompartmentVolume.md)

## Examples

``` r
# Generate eqnlist from the constructor
S <- matrix(c(-1, 1, 1, -1),
            nrow = 2, ncol = 2,
            dimnames = list(NULL, c("A", "B")))

rates <- c("k1*A", "k2*B")
description <- c("forward", "backward")

f <- eqnlist(smatrix = S, rates = rates, description = description)
print(f)
#> Conserved quantities:
#>   total_1 = 1*A+1*B
#> 
#>   Check Educt -> Product Rate Description
#> 1           A ->       B k1*A     forward
#> 2           B ->       A k2*B    backward

# Convert to data.frame
fdata <- as.data.frame(f)
print(fdata)
#>   Description Rate  A  B
#> 1     forward k1*A -1  1
#> 2    backward k2*B  1 -1

# Flat `volumes`: each distinct expression becomes a compartment c1, c2, ...
f <- as.eqnlist(fdata, volumes = c(A = "Vcyt", B = "Vnuc"))
print(f)
#> Conserved quantities:
#>   total_1 = 1*Vcyt*A+1*Vnuc*B
#> 
#> Compartments:
#>   c1 (V=Vcyt)
#>   c2 (V=Vnuc)
#> States by compartment:
#>   c1: A
#>   c2: B
#> 
#>   Check Educt -> Product Rate Description
#> 1           A ->       B k1*A     forward
#> 2           B ->       A k2*B    backward
print(as.eqnvec(f))
#> Idx Inner <- Outer
#>   1     A <- -1*(k1*A)+1*(k2*B)*(Vnuc/Vcyt)
#>   2     B <- 1*(k1*A)*(Vcyt/Vnuc)-1*(k2*B)
print(as.eqnvec(f, type = "amount"))
#> Idx Inner <- Outer
#>   1     A <- -1*(k1*A)*(Vcyt)+1*(k2*B)*(Vnuc)
#>   2     B <- 1*(k1*A)*(Vcyt)-1*(k2*B)*(Vnuc)

# Named compartments with explicit state assignment
f2 <- eqnlist(
  smatrix = S, rates = rates, description = description,
  compartments = list(cyt = "Vcyt", nuc = "Vnuc"),
  compartmentOf = c(A = "cyt", B = "nuc")
)
print(f2)
#> Conserved quantities:
#>   total_1 = 1*Vcyt*A+1*Vnuc*B
#> 
#> Compartments:
#>   cyt (V=Vcyt)
#>   nuc (V=Vnuc)
#> States by compartment:
#>   cyt: A
#>   nuc: B
#> 
#>   Check Educt -> Product Rate Description
#> 1           A ->       B k1*A     forward
#> 2           B ->       A k2*B    backward
print(as.eqnvec(f2))
#> Idx Inner <- Outer
#>   1     A <- -1*(k1*A)+1*(k2*B)*(Vnuc/Vcyt)
#>   2     B <- 1*(k1*A)*(Vcyt/Vnuc)-1*(k2*B)
```
