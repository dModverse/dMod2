\donttest{
f <- eqnvec(A = "-k1*A", B = "k1*A - k2*B")
model <- odemodel(f, modelname = "Xf_AB", outdir = tempdir())

## States alone, without sensitivities
x <- Xf(model)
pred <- x(seq(0, 10, 0.1), c(A = 1, B = 0, k1 = 0.5, k2 = 0.2))
plot(pred)
}
