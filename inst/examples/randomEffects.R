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
re <- randomEffects(trafo)
re
re$fixed

## One random intercept for many parameters at once
randomEffects(insert(branch(define(NULL, "x <- x", x = c("ka", "cl", "v")), table = cov),
                     x ~ 1 + (1 | ID), x = c("ka", "cl", "v")))
