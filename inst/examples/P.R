\donttest{
## One transformation per condition, all compiled into one shared object
trafo <- list(
  ctrl  = eqnvec(k1 = "exp(log_k1)", k2 = "exp(log_k2)", A = "A0"),
  treat = eqnvec(k1 = "exp(log_k1 + log_fold)", k2 = "exp(log_k2)", A = "A0"))
p <- P(trafo, compile = TRUE, modelname = "P_example", outdir = tempdir())
getParameters(p)
p(c(log_k1 = 0, log_k2 = -1, A0 = 2, log_fold = log(3)))
}
