\donttest{
## Prediction function
regfn <- c(y = "sin(a*time)")
g <- Y(regfn, parameters = "a", compile = TRUE, modelname = "parlist_obs",
       outdir = tempdir())
x <- Xt(condition = "C1")

## Data
set.seed(1)
data <- datalist(
  C1 = data.frame(
    name = "y",
    time = 1:5,
    value = sin(1:5) + rnorm(5, 0, .1),
    sigma = .1
  )
)

pars <- c(a = 1)
times <- seq(0, 5, .1)
plot((g*x)(times, pars), data)

## Fits from random starts, stored in a parlist
obj <- normL2(data, g*x)
out <- as.parlist(lapply(1:5, function(i) {
  trust(obj, pars + rnorm(length(pars), 0, 1), rinit = 1, rmax = 10)
}))
summary(out)

## Parameter frame and best fit
parframe <- as.parframe(out)
plotValues(parframe)
bestfit <- as.parvec(parframe)
plot((g*x)(times, bestfit), data)
}
