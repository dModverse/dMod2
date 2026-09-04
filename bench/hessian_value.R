# How much of a trust step comes from the Hessian: steptype and cos(step, -g)
# per iteration from random starts. A boundary step at small radius tends to
# -g/lambda, so cos near 1 means the Hessian bought nothing at that iteration.
library(dMod2); Sys.setenv(OMP_NUM_THREADS = 1); setwd(tempdir())
so <- list(atol = 1e-9, rtol = 1e-9); cores <- 5
y <- path.expand("~/Documents/Projects/dModverse/cppDE/benchmarks/cache/petab/Benchmark-Models/Bachmann_MSB2011/Bachmann_MSB2011.yaml")
pe <- importPEtab(y, backend = "cppDE", cores = cores, outdir = tempdir(),
                  optionsOde = so, optionsSens = so)
vOpt <- pe$obj(pe$bestfit, cores = cores)$value
npar <- length(pe$bestfit)
cat("optimum", vOpt, " npar", npar, "\n")

# Gradients cached by parameter vector, so reconstructing the step directions
# afterwards costs no extra evaluation.
cache <- new.env(); cache$key <- list(); cache$grad <- list()
obj <- function(p) {
  o <- pe$obj(p, cores = cores)
  cache$key[[length(cache$key) + 1L]]  <- unname(p)
  cache$grad[[length(cache$grad) + 1L]] <- o$gradient
  o
}
lookup <- function(p) {
  hit <- which(vapply(cache$key, function(k) isTRUE(all.equal(k, unname(p), tolerance = 0)), logical(1)))
  if (length(hit)) cache$grad[[hit[1]]] else NULL
}

starts <- msParframe(pe$bestfit, n = 3, seed = 42, sd = 1, keepfirst = FALSE)
inBox <- function(p) pmin(pmax(p, pe$parlower + 1e-3), pe$parupper - 1e-3)

for (k in seq_len(nrow(starts))) {
  p0 <- inBox(unlist(starts[k, ]))[names(pe$bestfit)]
  v0 <- pe$obj(p0, cores = cores)$value
  if (!is.finite(v0)) { cat("\nstart", k, "nicht endlich, uebersprungen\n"); next }
  cache$key <- list(); cache$grad <- list()
  f <- trust(obj, p0, rinit = 1, rmax = 10, iterlim = 300, gtol = 1e-4, blather = TRUE,
             parlower = pe$parlower, parupper = pe$parupper)

  step <- f$argtry - f$argpath
  cosg <- vapply(seq_len(nrow(f$argpath)), function(i) {
    g <- lookup(f$argpath[i, ]); s <- step[i, ]
    if (is.null(g) || sum(s^2) == 0) return(NA_real_)
    g <- g[colnames(f$argpath)]
    -sum(s * g) / sqrt(sum(s^2) * sum(g^2))
  }, numeric(1))

  prog <- 100 * (v0 - f$valpath) / (v0 - vOpt)
  cat(sprintf("\n=== start %d === v0 %.2f -> %.2f (%.1f%%), %d iters, %s\n",
              k, v0, f$value, 100*(v0 - f$value)/(v0 - vOpt), f$iterations, f$stopReason))
  cat("steptype:", paste(names(table(f$steptype)), table(f$steptype), collapse = "  "), "\n")
  phase <- cut(prog, c(-Inf, 50, 90, Inf), labels = c("<50%", "50-90%", ">90%"))
  print(round(do.call(rbind, tapply(cosg, phase, function(z)
    c(n = length(z), median = median(z, na.rm = TRUE),
      min = min(z, na.rm = TRUE)))), 4))
  cat("Anteil cos>0.99:", round(mean(cosg > 0.99, na.rm = TRUE), 3),
      "  cos>0.999:", round(mean(cosg > 0.999, na.rm = TRUE), 3), "\n")
  cat("Newton-Schritte:", sum(f$steptype == "Newton"), "von", length(f$steptype), "\n")
}
