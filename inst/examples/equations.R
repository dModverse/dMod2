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

## Conserved quantities
cq <- conservedQuantities(eq$smatrix)
print(cq)

## Fluxes
fluxes <- getFluxes(eq)
print(fluxes)

## Subsetting
subeq1 <- subset(eq, "pB" %in% Product)
print(subeq1)
subeq2 <- subset(eq, grepl("not_available", Description))
print(subeq2)

## Equation vector of the ODE
f <- as.eqnvec(eq)
print(f)
