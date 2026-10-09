## Identity for four parameters, then log parameters
trafo <- define(NULL, "x ~ x", x = c("A", "B", "k1", "k2"))
trafo <- insert(trafo, "x ~ exp(x)", x = .currentSymbols)
trafo

## Replace k1 and k2 by expressions, fix B
trafo <- insert(trafo, "x ~ y", x = c("k1", "k2"), y = c("q1 + q2", "q1 - q2"))
trafo <- define(trafo, "B ~ 0")
trafo

## One branch per condition; the non-missing table entries are inserted
conditions <- data.frame(dose = c(1, 10), q2 = c(NA, "q2_high"),
                         row.names = c("low", "high"))
trafoL <- branch(trafo, table = conditions, apply = "insert")

## Condition-specific changes, selected by name or read off the table
trafoL <- insert(trafoL, "x ~ 2*x", x = "q1", conditionMatch = "^high")
trafoL <- define(trafoL, "A ~ dose", dose = dose)
trafoL
