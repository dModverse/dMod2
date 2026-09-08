# -------------------------------------------------------------------------#
# Cost of a gradient relative to a bare value solve
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# `factor` is the price of one objective gradient in value solves. Forward
# sensitivities pay it per evaluation; the adjoint pays a constant instead, so
# `factor_grad` is the number the adjoint has to beat. Measured per model and,
# on one model, as a curve over the number of estimated parameters.
#
# `factor_rev` is the adjoint's own price, from obj(pars, sweep = "reverse").
# It is not expected to be one: a reverse evaluation integrates the states twice
# -- once for the values, once inside the sweep, because a seed only exists
# after the chain above has been walked -- and sweeps a tape one step wide on
# top. What matters is that it stays flat while factor_grad rises, and where the
# two lines cross.
#
# [AUTHOR]
# Simon Beyer
#
# [Date]
# Fri 05 Sep 2026
#
# [Info]
# Needs an idle machine and `OMP_NUM_THREADS=1`: a fit running beside it makes
# every number here meaningless. The two PEtab models ship with the package;
# point `.collection` at the benchmark collection for a wider sweep.
# -------------------------------------------------------------------------#

library(dMod2)
options(dMod.cores = 1)

.outdir <- file.path(tempdir(), "bench_gradientCost")
dir.create(.outdir, recursive = TRUE, showWarnings = FALSE)
.collection <- Sys.getenv("DMOD_PETAB_COLLECTION", "")

# The machine scatters, so report the minimum rather than the mean.
tmin <- function(f, reps = 7L)
  min(vapply(seq_len(reps), function(i) system.time(f())[["elapsed"]], 0.0))

probe <- function(obj, p, label) {
  # A model imported without reverse = TRUE has no reverse object; report NA
  # rather than failing the whole sweep for it.
  t_rev <- tryCatch(tmin(function() obj(p, deriv = TRUE, sweep = "reverse")),
                    error = function(e) NA_real_)
  data.frame(model = label, n_theta = length(p),
             t_value = tmin(function() obj(p, deriv = FALSE)),
             t_grad  = tmin(function() obj(p, deriv = TRUE, hessian = FALSE)),
             t_full  = tmin(function() obj(p, deriv = TRUE, hessian = TRUE)),
             t_rev   = t_rev)
}

importOne <- function(dir, tag) {
  yml <- list.files(dir, pattern = "\\.yaml$", full.names = TRUE)[1]
  od  <- file.path(.outdir, tag); dir.create(od, recursive = TRUE, showWarnings = FALSE)
  importPEtab(yml, backend = "cppDE", cores = 6, modelname = paste0("gc_", tag),
              reverse = TRUE, outdir = od)
}


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Per model
#
# Warnings are surfaced: a model that does not integrate cleanly at its own
# published optimum has a meaningless timing.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
dirs <- c(system.file("extdata", "petab_boehm",    package = "dMod2"),
          system.file("extdata", "petab_bachmann", package = "dMod2"))
if (nzchar(.collection)) dirs <- c(dirs, list.dirs(.collection, recursive = FALSE))

perModel <- NULL
for (d in dirs) {
  tag <- gsub("[^A-Za-z0-9]", "", basename(d))
  pet <- try(importOne(d, tag), silent = TRUE)
  if (inherits(pet, "try-error")) {
    cat("import failed:", basename(d), "\n"); next
  }
  withCallingHandlers(
    perModel <- rbind(perModel, probe(pet$obj, pet$bestfit, basename(d))),
    warning = function(w) {
      cat("  warning [", basename(d), "]:", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    })
  cat("done:", basename(d), "\n"); utils::flush.console()
}

perModel <- perModel[order(perModel$n_theta), ]
perModel$factor_grad <- perModel$t_grad / perModel$t_value
perModel$factor_full <- perModel$t_full / perModel$t_value
perModel$factor_rev  <- perModel$t_rev  / perModel$t_value
perModel$speedup     <- perModel$t_grad / perModel$t_rev
print(perModel, row.names = FALSE, digits = 4)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Over the number of estimated parameters
#
# The same model re-imported with only the first k parameters estimated, the
# rest held at their nominal value. Fixing them at call time does not reduce
# the propagated directions, re-importing does.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
.src <- system.file("extdata", "petab_bachmann", package = "dMod2")
.tab <- list.files(.src, pattern = "^parameters", full.names = TRUE)[1]
estimated <- read.delim(.tab)$parameterId[read.delim(.tab)$estimate == 1]

atK <- function(k) {
  dst <- file.path(.outdir, paste0("pet", k))
  dir.create(dst, recursive = TRUE, showWarnings = FALSE)
  file.copy(list.files(.src, full.names = TRUE), dst, overwrite = TRUE)
  tsv <- file.path(dst, basename(.tab))
  tab <- read.delim(tsv)
  tab$estimate[tab$parameterId %in% estimated[-seq_len(k)]] <- 0L
  write.table(tab, tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  importOne(dst, paste0("k", k))
}

overK <- do.call(rbind, lapply(c(10L, 25L, 50L, 75L, length(estimated)), function(k) {
  pet <- atK(k)
  probe(pet$obj, pet$bestfit, "Bachmann")
}))
overK$factor_grad <- overK$t_grad / overK$t_value
overK$factor_full <- overK$t_full / overK$t_value
overK$factor_rev  <- overK$t_rev  / overK$t_value
overK$speedup     <- overK$t_grad / overK$t_rev
print(overK, row.names = FALSE, digits = 4)

# The whole point in one picture: one line rises with n_theta, the other does
# not, and where they cross is a property of the model.
matplot(overK$n_theta, cbind(overK$factor_grad, overK$factor_rev),
        type = "b", pch = c(1, 4), lty = 1, log = "x", col = c(1, 2),
        xlab = "estimated parameters", ylab = "gradient cost in value solves",
        main = "Bachmann2011")
legend("topleft", c("forward sensitivities", "adjoint"), pch = c(1, 4),
       col = c(1, 2), bty = "n")
