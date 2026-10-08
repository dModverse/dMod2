\donttest{
f <- eqnvec(A = "-k1*A", B = "k1*A - k2*B")
model <- odemodel(f, modelname = "Xs_AB", outdir = tempdir())

## Tighter tolerances for the solves with sensitivities, as used in a fit
x <- Xs(model, optionsSens = list(atol = 1e-8, rtol = 1e-8))
times <- seq(0, 10, 0.1)
pars <- c(A = 1, B = 0, k1 = 0.5, k2 = 0.2)
pred <- x(times, pars)
plot(pred)
dim(attr(pred[[1]], "deriv"))

## Change the options of the existing prediction function
controls(x, NULL, "options") <- list(atol = 1e-10, rtol = 1e-10)
head(x(times, pars, deriv = FALSE)[[1]])
}
