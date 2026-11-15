## The distributions of a random term, written into gr() of an effect formula
normal()
spikeSlab()
spikeSlab(prob = 0.2, tau = 0.5)
dirichletProcess()
dirichletProcess(concentration = 1)

## cell lines whose rates share the population value or a common deviation in clusters
cov <- data.frame(line = paste0("L", 1:4), row.names = paste0("C", 1:4))
trafo <- branch(define(NULL, "x <- x", x = c("kp", "kd")), table = cov)
trafo <- insert(trafo, "x <- exp(x)", x = c("kp", "kd"))
trafo <- insert(trafo, x ~ 1 + (1 | gr(line, dist = dirichletProcess())), x = c("kp", "kd"))
randomEffects(trafo)
