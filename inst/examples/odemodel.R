\donttest{
## A model with a forcing F, no sensitivities for "switch"
f <- eqnvec(A = "-k*A + switch*F")
model <- odemodel(f, forcings = "F", fixed = "switch", modelname = "odemodel_F",
                  outdir = tempdir(), compile = FALSE)
print(model)

## The same model from an equation list
f2 <- addReaction(NULL, from = "", to = "A", rate = "switch*F", description = "production")
f2 <- addReaction(f2, from = "A", to = "", rate = "k*A", description = "degradation")
print(f2)
model2 <- odemodel(f2, forcings = "F", fixed = "switch", modelname = "odemodel_F2",
                   outdir = tempdir(), compile = FALSE)

## One prediction per forcing; the data of each forcing cover the time range
tF <- seq(0, 5, 0.1)
x <- Xs(model, data.frame(name = "F", time = tF, value = sin(tF)), condition = "sin") +
  Xs(model, data.frame(name = "F", time = tF, value = exp(-tF)), condition = "exp")
x2 <- Xs(model2, data.frame(name = "F", time = c(0, 5), value = 0.1), condition = "const")

## Compile both models into one shared object
compile(x, x2, output = "odemodel_ex")

pars <- c(A = 1, k = 0.5, switch = 1)
plot((x + x2)(seq(0, 5, 0.01), pars))
}
