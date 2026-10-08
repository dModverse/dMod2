\donttest{
## A model with a forcing F, no sensitivities for "switch"
f <- eqnvec(A = "-k*A + switch*F")
model <- odemodel(f, forcings = "F", fixed = "switch", modelname = "odemodel_F",
                  outdir = tempdir())
print(model)

## The same model from an equation list
f <- addReaction(NULL, from = "", to = "A", rate = "switch*F", description = "production")
f <- addReaction(f, from = "A", to = "", rate = "k*A", description = "degradation")
print(f)
model <- odemodel(f, forcings = "F", fixed = "switch", modelname = "odemodel_F2",
                  outdir = tempdir())

## One prediction per forcing; the data of each forcing cover the time range
tF <- seq(0, 5, 0.1)
x <- Xs(model, data.frame(name = "F", time = tF, value = sin(tF)), condition = "sin") +
  Xs(model, data.frame(name = "F", time = tF, value = exp(-tF)), condition = "exp") +
  Xs(model, data.frame(name = "F", time = c(0, 5), value = 0.1), condition = "const")

pars <- c(A = 1, k = 0.5, switch = 1)
plot(x(seq(0, 5, 0.01), pars))
}
