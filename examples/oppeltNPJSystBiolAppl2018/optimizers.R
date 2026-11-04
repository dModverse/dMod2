# The Hessian sources of trust() on the full model: fits by Gauss-Newton, BFGS and SR1 with
# reverse gradients, the exact Hessian by forward over reverse, and Gauss-Newton handing over to
# BFGS, all from the same 100 starts. Saves results/optimizers.rds, which the guide reads.

# the chain with the reverse sweep, from sensitivities.R
if (!exists("objR")) source("sensitivities.R")

## ---- optimizers
# 100 starts around the means of the prior, as in fit.R. BFGS and SR1 take the reverse gradient
# and start from the identity; the BFGS phase after Gauss-Newton starts from its Hessian
free <- setdiff(names(mu), names(gauge))
starts <- msParframe(mu[free], n = 100, seed = 2, sd = 2, keepfirst = FALSE)
methods <- list(
  `Gauss-Newton` = list(hessianMethod = "gn"),
  BFGS = list(hessianMethod = "bfgs", sweep = "reverse", qnControl = list(hessianInit = "identity")),
  SR1 = list(hessianMethod = "sr1", sweep = "reverse", qnControl = list(hessianInit = "identity")),
  exact = list(hessianMethod = "exact", sweep = "reverse"),
  `GN, then BFGS` = list(phases = list(list(hessianMethod = "gn", iterlim = 20),
                                       list(hessianMethod = "bfgs"))))
# one fit of start `i` by method `m`, with its wall time
fitOne <- function(m, i) {
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(do.call(trust, c(list(objR, setNames(as.numeric(starts[i, free]), free),
                                        fixed = gauge, iterlim = 1000, tolControl = tolerances),
                                   methods[[m]])), error = function(e) NULL)
  data.frame(method = m, start = i, value = if (is.null(fit)) NA_real_ else fit$value,
             converged = !is.null(fit) && isTRUE(fit$converged),
             iterations = if (is.null(fit)) NA_integer_ else as.integer(fit$iterations),
             elapsed = proc.time()[["elapsed"]] - t0)
}
tasks <- expand.grid(method = names(methods), start = seq_len(nrow(starts)), stringsAsFactors = FALSE)
optimizers <- do.call(rbind, parallel::mclapply(seq_len(nrow(tasks)), function(k)
  fitOne(tasks$method[k], tasks$start[k]), mc.cores = getOption("dMod2.cores")))

dir.create("results", showWarnings = FALSE)
saveRDS(optimizers, "results/optimizers.rds")
