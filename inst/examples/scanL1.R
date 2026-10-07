\dontrun{
## A decays into B, B decays. The model also carries a direct loss of A that
## the data do not need; its gate is pulled to exactly zero.
outdir <- tempdir()
f <- eqnvec(A = "-k1 * A - k2 * A", B = "k1 * A - k3 * B")
x <- Xs(odemodel(f, modelname = "scanL1_x", compile = FALSE, outdir = outdir),
        compile = FALSE)
g <- Y(eqnvec(obsA = "A", obsB = "B"), f, modelname = "scanL1_g",
       compile = FALSE, attach.input = FALSE, outdir = outdir)
trafo <- gateL1(eqnvec(A = "1", B = "0", k1 = "10^log10_k1",
                       k2 = "10^log10_k2", k3 = "10^log10_k3"),
                c("k1", "k2", "k3"))
p <- P(trafo, condition = "C1", modelname = "scanL1_p", compile = FALSE,
       outdir = outdir)
compile(g, x, p, output = "scanL1_example")
prd <- g * x * p

set.seed(1)
truth <- c(log10_k1 = log10(0.5), log10_k2 = -1, log10_k3 = log10(0.2),
           s_k1 = 1, s_k2 = 0, s_k3 = 1)
sim <- wide2long(prd(seq(0, 10, 0.5), truth, deriv = FALSE))
sim$value <- sim$value + rnorm(nrow(sim), 0, 0.02)
sim$sigma <- 0.02
obj <- normL2(as.datalist(sim), prd)

fit <- scanL1(obj, c(log10_k1 = 0, log10_k2 = -0.5, log10_k3 = 0),
              zero = attr(trafo, "gates"), lambda = 10^seq(-1, 4, length.out = 11),
              fits = 4)
fit
plot(fit)
plot(fit, type = "test")

## Spike-and-slab lasso: one multistart of the full model, then one chain over
## the spike strength; the structure is the one the chain settles on.
fitS <- scanL1(obj, c(log10_k1 = 0, log10_k2 = -0.5, log10_k3 = 0),
               zero = attr(trafo, "gates"), lambda = 10^seq(0, 4, length.out = 13),
               fits = 4, pathFits = 1, ssl = list(lambda1 = 1), select = "plateau")
fitS
plot(fitS, type = "inclusion")
plot(fitS, type = "waterfall")
}
