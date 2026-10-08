\donttest{
## Observation function and one parameter transformation per condition
fn <- eqnvec(
  sine = "1 + sin(6.28*omega*time)",
  cosine = "cos(6.28*omega*time)"
)
g <- Y(fn, parameters = "omega", outdir = tempdir())
x <- Xt()
trafo <- lapply(1:3, function(i) eqnvec(omega = paste0("omega_", i)))
names(trafo) <- paste0("frequency_", 1:3)
p <- P(trafo, outdir = tempdir())
compile(g, p, output = "plotting_example")

## Evaluate prediction
times <- seq(0, 1, .01)
pars <- structure(seq(1, 2, length.out = 3), names = getParameters(p))
prediction <- (g*x*p)(times, pars, deriv = FALSE)

## Plot prediction
plotPrediction(prediction)
plotPrediction(prediction, scales = "fixed")
plotPrediction(prediction, facet = "grid")
plotPrediction(prediction, scales = "fixed",
               transform = list(sine = "x^2", cosine = "x - 1"))

## Simulate data
set.seed(1)
dataset <- wide2long(prediction)
dataset <- dataset[seq(1, nrow(dataset), 5), ]
dataset$value <- dataset$value + rnorm(nrow(dataset), 0, .1)
dataset$sigma <- 0.1
data <- as.datalist(dataset, split.by = "condition")

## Plot data, and data with prediction
plotData(data)
plotCombined(prediction, data)
plotCombined(prediction, data, time <= 0.5 & condition == "frequency_1")
plotCombined(prediction, data, time <= 0.5 & condition != "frequency_1",
             facet = "grid")
plotCombined(prediction, data, aesthetics = list(linetype = "condition"))
}
