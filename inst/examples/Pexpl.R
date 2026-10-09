\donttest{
trafo <- eqnvec(k1 = "exp(log_k1)", k2 = "exp(log_k2)", A = "A0")
p <- Pexpl(trafo, condition = "C1", compile = TRUE,
           modelname = "Pexpl_example", outdir = tempdir())
out <- p(c(log_k1 = 0, log_k2 = -1, A0 = 2))
out$C1
attr(out$C1, "deriv")
}
