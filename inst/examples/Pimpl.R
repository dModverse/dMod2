\donttest{
## Steady state of A <-> B; the total of A and B enters as a parameter
r <- addReaction(NULL, "A", "B", "k1*A")
r <- addReaction(r, "B", "A", "k2*B")
p <- Pimpl(r, condition = "C1", compile = TRUE,
           modelname = "Pimpl_example", outdir = tempdir())
getParameters(p)
pars <- setNames(c(1, 0.5, 3), getParameters(p))
out <- p(pars)
out$C1
attr(out$C1, "deriv")
}
