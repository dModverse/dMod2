\donttest{
  # A prediction function interpolating parameter values piecewise linearly
  times <- 0:5
  grid <- data.frame(name = "A", time = times, row.names = paste0("p", times))
  x <- Xd(grid)

  # An observable and its observation function
  observables <- eqnvec(Aobs = "s*A")
  g <- Y(g = observables, f = NULL, states = "A", parameters = "s",
         compile = TRUE, modelname = "prediction_obs", outdir = tempdir())

  # One parameter transformation per experimental condition
  innerpars <- c(getParameters(x), getParameters(g))
  trafo <- structure(innerpars, names = innerpars)
  trafo_C1 <- cOde::replaceSymbols(innerpars, paste(innerpars, "C1", sep = "_"), trafo)
  trafo_C2 <- cOde::replaceSymbols(innerpars, paste(innerpars, "C2", sep = "_"), trafo)

  p <- P(trafo = trafo_C1, condition = "C1", compile = TRUE,
         modelname = "prediction_p1", outdir = tempdir()) +
    P(trafo = trafo_C2, condition = "C2", compile = TRUE,
      modelname = "prediction_p2", outdir = tempdir())

  # Outer parameters at random values
  outerpars <- getParameters(p)
  pars <- structure(runif(length(outerpars), 0, 1), names = outerpars)

  # Unobserved states, then observables
  plot((x*p)(times, pars))
  plot((g*x*p)(times, pars))
}
