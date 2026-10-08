\donttest{
old <- options(dMod.outdir = tempdir())

## Decay model, observation and parameter transformation, compiled together
f <- addReaction(eqnlist(), from = "A", to = "", rate = "k*A")
x <- Xs(odemodel(f, modelname = "profile_x", compile = FALSE))
g <- Y(c(y = "A"), f = x, attach.input = FALSE, modelname = "profile_g",
       compile = FALSE)
p <- P(eqnvec(A = "exp(logA)", k = "exp(logk)"), condition = "C1",
       modelname = "profile_p", compile = FALSE)
compile(x, g, p, output = "profile_example", cores = 1)

## Data, objective and fit
set.seed(1)
times <- seq(0, 5, by = 0.5)
data <- datalist(C1 = data.frame(
  name = "y", time = times, sigma = 0.05,
  value = 2 * exp(-0.5 * times) + rnorm(length(times), sd = 0.05)))
obj <- normL2(data, g * x * p)
myfit <- trust(obj, c(logA = 0, logk = 0), rinit = 1, rmax = 10)

## Profiles by integration and by repeated optimisation
profiles.approx <- profile(obj, myfit$argument, whichPar = c("logA", "logk"))
profiles.exact <- profile(obj, myfit$argument, whichPar = c("logA", "logk"),
                          method = "optimize")

## Plots and confidence intervals
plotProfile(list(approx = profiles.approx, exact = profiles.exact))
plotPaths(profiles.approx, whichPar = "logk")
confint(profiles.approx, val.column = "value")

options(old)
}
