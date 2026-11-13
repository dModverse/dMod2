# Generate List of Fluxes from Equation List

Generate List of Fluxes from Equation List

## Usage

``` r
getFluxes(eqnlist, type = c("conc", "amount"))
```

## Arguments

- eqnlist:

  object of class
  [eqnlist](https://dmodverse.github.io/dMod2/reference/eqnlist.md).

- type:

  `"conc"` (default) or `"amount"` for fluxes in units of concentration
  or of amount.

## Value

Named list with one named character vector per state, the in- and
out-fluxes of the state.

## Examples

``` r
## Equation list
eq <- eqnlist()
eq <- addReaction(eq, "A", "pA", "act_A * A * stimulus", "Phosphorylation of A")
eq <- addReaction(eq, "pA", "A", "deact_A * pA", "Dephosphorylation of pA")
eq <- addReaction(eq, "2*pA", "pA_pA", "form_complex_pA * pA^2", "Complex formation of pA")
eq <- addReaction(eq, "B", "pB", "act_B * B * pA_pA", "Phosphorylation of B")
eq <- addReaction(eq, "pB", "B", "deact_B * pB", "Dephosphorylation of pB")

## Data frame of reactions
reactions <- getReactions(eq)
print(reactions)
#>   Check Educt -> Product                   Rate             Description
#> 1           A ->      pA   act_A * A * stimulus    Phosphorylation of A
#> 2          pA ->       A           deact_A * pA Dephosphorylation of pA
#> 3        2*pA ->   pA_pA form_complex_pA * pA^2 Complex formation of pA
#> 4           B ->      pB      act_B * B * pA_pA    Phosphorylation of B
#> 5          pB ->       B           deact_B * pB Dephosphorylation of pB

## Conserved quantities
print(conservedQuantities(eq))
#> $total_1
#> [1] "0.5*A+0.5*pA+1*pA_pA"
#> 
#> $total_2
#> [1] "1*B+1*pB"
#> 

## Fluxes
fluxes <- getFluxes(eq)
print(fluxes)
#> $A
#>        Phosphorylation of A     Dephosphorylation of pA 
#> "-1*(act_A * A * stimulus)"         "+1*(deact_A * pA)" 
#> 
#> $pA
#>          Phosphorylation of A       Dephosphorylation of pA 
#>    "1*(act_A * A * stimulus)"           "-1*(deact_A * pA)" 
#>       Complex formation of pA 
#> "-2*(form_complex_pA * pA^2)" 
#> 
#> $pA_pA
#>      Complex formation of pA 
#> "1*(form_complex_pA * pA^2)" 
#> 
#> $B
#>     Phosphorylation of B  Dephosphorylation of pB 
#> "-1*(act_B * B * pA_pA)"      "+1*(deact_B * pB)" 
#> 
#> $pB
#>    Phosphorylation of B Dephosphorylation of pB 
#> "1*(act_B * B * pA_pA)"     "-1*(deact_B * pB)" 
#> 

## Subsetting
subeq1 <- subset(eq, "pB" %in% Product)
print(subeq1)
#> Conserved quantities:
#>   pA_pA = 1*pA_pA
#>   total_2 = 1*B+1*pB
#> 
#>   Check Educt -> Product              Rate          Description
#> 1           B ->      pB act_B * B * pA_pA Phosphorylation of B
subeq2 <- subset(eq, grepl("not_available", Description))
print(subeq2)
#> 
#> NULL

## Equation vector of the ODE
f <- as.eqnvec(eq)
print(f)
#> Idx Inner <- Outer
#>   1     A <- -1*(act_A*A*stimulus)+1*(deact_A*pA)
#>   2     B <- -1*(act_B*B*pA_pA)+1*(deact_B*pB)
#>   3    pA <- 1*(act_A*A*stimulus)-1*(deact_A*pA)-2*(form_complex_pA*pA^2)
#>   4 pA_pA <- 1*(form_complex_pA*pA^2)
#>   5    pB <- 1*(act_B*B*pA_pA)-1*(deact_B*pB)
```
