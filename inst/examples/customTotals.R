## A reversible conversion conserves A + B
el <- addReaction(eqnlist(), "A", "B", "k1*A")
el <- addReaction(el, "B", "A", "k2*B")
conservedQuantities(el)

## Name the conserved quantity, then return to the detected basis
el <- customTotals(el, list(total = "A + B"))
el$totals
el <- customTotals(el, NULL)
