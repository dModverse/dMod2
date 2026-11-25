## The distributions of a random term, written into gr() of an effect formula
normal()
spikeSlab()
spikeSlab(prob = 0.2, tau = 0.5)
spikeSlab(mean = log(0.5), spike = -Inf)
dirichletProcess()
dirichletProcess(concentration = 1)

## cell lines whose rates share the population value or a common deviation in clusters
cov <- data.frame(line = paste0("L", 1:4), row.names = paste0("C", 1:4))
trafo <- branch(define(NULL, "x <- x", x = c("kp", "kd")), table = cov)
trafo <- insert(trafo, "x <- exp(x)", x = c("kp", "kd"))
trafo <- insert(trafo, x ~ 1 + (1 | gr(line, dist = dirichletProcess())), x = c("kp", "kd"))
randomEffects(trafo)

## reaction rates switched off by a spike at -Inf, one group of all conditions
sel <- branch(define(NULL, "x <- x", x = c("kp", "kd")), table = cov)
sel <- insert(sel, "x <- exp(x)", x = c("kp", "kd"))
sel <- insert(sel, x ~ 0 + (1 | gr(dist = spikeSlab(mean = log(0.5), spike = -Inf))),
              x = c("kp", "kd"))
randomEffects(sel)
