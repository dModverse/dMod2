## Accessors read the generated objects, no compiled code is needed
g <- Y(c(Aobs = "s*A"), states = "A", parameters = "s", compile = FALSE,
       modelname = "accessors_obs", outdir = tempdir())
p <- P(list(C1 = eqnvec(A = "A0", s = "s"), C2 = eqnvec(A = "2*A0", s = "s")),
       compile = FALSE, modelname = "accessors_p", outdir = tempdir())

getParameters(g)
getParameters(p)
getConditions(p)
getEquations(p, conditions = "C2")

modelname(g)
modelname(g, p)
modelname(p) <- "accessors_trafo"
modelname(p)
