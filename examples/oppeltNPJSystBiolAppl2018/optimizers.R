# The Hessian sources of trust() on the full model: fits by Gauss-Newton, BFGS and SR1 with
# reverse gradients, the exact Hessian by forward over reverse, Gauss-Newton handing over to BFGS,
# and the structured secant updates with their Gauss-Newton hybrids, all from the same 100 starts.
# Saves results/optimizers.rds, which the guide reads.

# the chain with the reverse sweep, from sensitivities.R
if (!exists("objR")) source("sensitivities.R")

## ---- optimizers
# 100 starts drawn from the prior, as in fit.R; a start without a stable resting state is replaced
# by the next one of a second draw, as mstrust() draws anew. BFGS and SR1 take the reverse gradient
# and start from the identity; the BFGS phase after Gauss-Newton, after 20 iterations or once the
# radius stays unchanged, starts from the Gauss-Newton Hessian; the structured methods add a learned correction to the Gauss-Newton Hessian at every
# iterate, scaled with the residual norm or, on the Gauss-Newton branch of a hybrid, shifted by it
free <- setdiff(names(mu), names(gauge))
draw <- function(seed) msParframe(mu[free], n = 100, seed = seed, sd = sigma[free], keepfirst = FALSE)
feasible <- function(p) !inherits(try(objR(p, fixed = gauge, deriv = FALSE), silent = TRUE),
                                  "try-error")
starts <- draw(2)
spare <- draw(3)
spare <- spare[vapply(seq_len(nrow(spare)), function(i)
  feasible(setNames(as.numeric(spare[i, free]), free)), TRUE), ]
bad <- which(!vapply(seq_len(nrow(starts)), function(i)
  feasible(setNames(as.numeric(starts[i, free]), free)), TRUE))
starts[bad, ] <- spare[seq_along(bad), ]
methods <- list(
  `default (Fides hybrid)` = list(),
  `Gauss-Newton` = list(hessianMethod = "gn"),
  BFGS = list(hessianMethod = "bfgs", sweep = "reverse", qnControl = list(hessianInit = "identity")),
  SR1 = list(hessianMethod = "sr1", sweep = "reverse", qnControl = list(hessianInit = "identity")),
  exact = list(hessianMethod = "exact", sweep = "reverse"),
  `GN, then BFGS` = list(phases = list(list(hessianMethod = "gn", iterlim = 20),
                                       list(hessianMethod = "bfgs"))),
  `GN, then BFGS at stall 25` = list(phases = list(list(hessianMethod = "gn", radiusStall = 25),
                                                   list(hessianMethod = "bfgs"))),
  `GN, then BFGS at stall 50` = list(phases = list(list(hessianMethod = "gn", radiusStall = 50),
                                                   list(hessianMethod = "bfgs"))),
  SBFGS = list(hessianMethod = "sbfgs"),
  `SBFGS, scaled` = list(hessianMethod = "sbfgs", qnControl = list(qnScaling = "residual")),
  SSR1 = list(hessianMethod = "ssr1"),
  `SSR1, scaled` = list(hessianMethod = "ssr1", qnControl = list(qnScaling = "residual")),
  GNSBFGS = list(hessianMethod = "gnsbfgs"),
  `GNSBFGS, shifted` = list(hessianMethod = "gnsbfgs", qnControl = list(qnHybridShift = TRUE)),
  GNSSR1 = list(hessianMethod = "gnssr1"),
  `GNSSR1, shifted` = list(hessianMethod = "gnssr1", qnControl = list(qnHybridShift = TRUE)))
# one fit of start `i` by method `m`, with the wall time
fitOne <- function(m, i) {
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(do.call(trust, c(list(objR, setNames(as.numeric(starts[i, free]), free),
                                        fixed = gauge, iterlim = 1000, tolControl = tolerances),
                                   methods[[m]])), error = function(e) NULL)
  data.frame(method = m, start = i, value = if (is.null(fit)) NA_real_ else fit$value,
             converged = !is.null(fit) && isTRUE(fit$converged),
             iterations = if (is.null(fit)) NA_integer_ else as.integer(fit$iterations),
             stopReason = if (is.null(fit)) NA_character_ else fit$stopReason,
             elapsed = proc.time()[["elapsed"]] - t0)
}
tasks <- expand.grid(method = names(methods), start = seq_len(nrow(starts)), stringsAsFactors = FALSE)
optimizers <- do.call(rbind, parallel::mclapply(seq_len(nrow(tasks)), function(k)
  fitOne(tasks$method[k], tasks$start[k]), mc.cores = getOption("dMod2.cores"),
  mc.preschedule = FALSE))

## ---- optimizerssave
dir.create("results", showWarnings = FALSE)
# the md5 of the files the result is computed from, which the guide checks on reading
sources <- c("model.R", "sensitivities.R", "optimizers.R", "data/oppelt2018.csv")
attr(optimizers, "stamp") <- list(md5 = tools::md5sum(sources), dMod2 = packageVersion("dMod2"))
saveRDS(optimizers, "results/optimizers.rds")
