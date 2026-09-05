# Boehm 100-fit multistart: gn vs bfgs vs hybrid (Stage 1 decision experiment).
# Usage: Rscript reverseAD-boehm-experiment.R [nfit]   (default 100)
suppressMessages({ library(dMod2); library(ggplot2) })

nfit  <- { a <- commandArgs(trailingOnly = TRUE); if (length(a)) as.integer(a[1]) else 100L }
# mstrust parallelism forks on Unix (workers inherit the loaded model DLL); on
# Windows the PSOCK workers cannot use the temp-compiled model, so run serial.
ncore <- if (.Platform$OS.type == "windows") 1L else 6L
outdir <- file.path(tempdir(), "boehm"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
petabDir <- system.file("extdata", "petab_boehm", package = "dMod2")

data(boehm); mydataL <- as.datalist(boehm)

epo <- "1.25e-7*exp(-Epo_degradation_BaF3*time)"
reactions <- eqnlist() |>
  assignCompartment(nucpApA = "nuc", nucpApB = "nuc", nucpBpB = "nuc", volume = "0.45") |>
  addReaction("2*STAT5A", "pApA", paste0("k_phos*", epo, "*STAT5A^2"), "p1", compartment = "cyt") |>
  addReaction("STAT5A + STAT5B", "pApB", paste0("k_phos*", epo, "*STAT5A*STAT5B"), "p2", compartment = "cyt") |>
  addReaction("2*STAT5B", "pBpB", paste0("k_phos*", epo, "*STAT5B^2"), "p3", compartment = "cyt") |>
  addReaction("pApA", "nucpApA", "k_imp_homo*pApA", "i1", compartment = "cyt") |>
  addReaction("pApB", "nucpApB", "k_imp_hetero*pApB", "i2", compartment = "cyt") |>
  addReaction("pBpB", "nucpBpB", "k_imp_homo*pBpB", "i3", compartment = "cyt") |>
  addReaction("nucpApA", "2*STAT5A", "k_exp_homo*nucpApA", "e1", compartment = "nuc") |>
  addReaction("nucpApB", "STAT5A + STAT5B", "k_exp_hetero*nucpApB", "e2", compartment = "nuc") |>
  addReaction("nucpBpB", "2*STAT5B", "k_exp_homo*nucpBpB", "e3", compartment = "nuc") |>
  setCompartmentVolume(cyt = "1.4")

opt <- list(atol = 1e-8, rtol = 1e-6, maxattemps = 100L, maxsteps = 1e6)
model <- odemodel(reactions, modelname = "boehm_ode", compile = FALSE, outdir = outdir)
x <- Xs(model, optionsOde = opt, optionsSens = opt)
observables <- eqnvec(
  pSTAT5A_rel = "(100*pApB + 200*pApA*specC17)/(pApB + STAT5A*specC17 + 2*pApA*specC17)",
  pSTAT5B_rel = "-(100*pApB - 200*pBpB*(specC17 - 1))/((STAT5B*(specC17 - 1) - pApB) + 2*pBpB*(specC17 - 1))",
  rSTAT5A_rel = "(100*pApB + 100*STAT5A*specC17 + 200*pApA*specC17)/(2*pApB + STAT5A*specC17 + 2*pApA*specC17 - STAT5B*(specC17 - 1) - 2*pBpB*(specC17 - 1))")
errorModels <- eqnvec(pSTAT5A_rel = "sd_pSTAT5A_rel", pSTAT5B_rel = "sd_pSTAT5B_rel", rSTAT5A_rel = "sd_rSTAT5A_rel")
g <- Y(observables, x, modelname = "boehm_obs", attach.input = FALSE, compile = FALSE, outdir = outdir)
e <- Y(errorModels, g, modelname = "boehm_err", attach.input = FALSE, compile = FALSE, outdir = outdir)

innerpars <- getParameters(model, g, e)
estimated <- c("Epo_degradation_BaF3","k_exp_hetero","k_exp_homo","k_imp_hetero",
               "k_imp_homo","k_phos","sd_pSTAT5A_rel","sd_pSTAT5B_rel","sd_rSTAT5A_rel")
p <- eqnvec() |>
  define("x~x", x = innerpars) |>
  define("x~0", x = c("pApA","pApB","pBpB","nucpApA","nucpApB","nucpBpB")) |>
  insert("STAT5A ~ 207.6*ratio") |> insert("STAT5B ~ 207.6 - 207.6*ratio") |>
  insert("ratio ~ 0.693") |> insert("specC17 ~ 0.107") |>
  insert("x~exp10(x)", x = estimated) |>
  P(modelname = "boehm_trafo", condition = "Boehm2014", compile = FALSE, outdir = outdir)

compile(g, x, p, e, output = "boehm", cores = ncore)
prd <- g * x * p; outerpars <- getParameters(prd)

pars <- read.delim(file.path(petabDir, "parameters_Boehm_JProteomeRes2014.tsv"))
pars <- pars[pars$estimate == 1, ]
lower <- setNames(log10(pars$lowerBound), pars$parameterId)[outerpars]
upper <- setNames(log10(pars$upperBound), pars$parameterId)[outerpars]
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)
dyn <- setdiff(outerpars, grep("^sd_", outerpars, value = TRUE))
obj <- normL2(mydataL, prd, e) + constraintL2(pouter[dyn], sigma = 4)
bestpub <- setNames(log10(pars$nominalValue), pars$parameterId)[outerpars]
cat(sprintf("obj at published optimum: %.4f\n", obj(bestpub, deriv = FALSE)$value))

# center-start sanity fit, one per method
for (hm in c("gn","bfgs","hybrid")) {
  f0 <- try(trust(obj, pouter, rinit = 0.1, rmax = 10, iterlim = 5000,
                  parlower = lower, parupper = upper, hessianMethod = hm), silent = TRUE)
  cat(sprintf("  center %-7s: %s\n", hm,
      if (inherits(f0, "try-error")) "ERROR"
      else sprintf("value=%.4f conv=%s neval=%d qnEval=%d", f0$value, f0$converged, f0$neval, f0$qnEval)))
}

methods <- c("gn","bfgs","hybrid"); seed <- 5555
cat(sprintf("\nRunning %d fits x %d methods (cores=%d) ...\n", nfit, length(methods), ncore))
runs <- lapply(methods, function(hm) {
  set.seed(seed)   # identical start set across methods
  t0 <- proc.time()[["elapsed"]]
  # retry = FALSE: a failed start is not re-sampled, so the RNG stream and thus
  # the start set stay identical across methods, and neval is not inflated.
  fl <- mstrust(obj, center = pouter, sd = 3, fits = nfit, cores = ncore,
                rinit = 0.1, rmax = 10, iterlim = 5000, start1stfromCenter = TRUE,
                retry = FALSE, parlower = lower, parupper = upper, hessianMethod = hm)
  list(frame = as.parframe(fl), secs = proc.time()[["elapsed"]] - t0)
})
names(runs) <- methods

best <- min(vapply(runs, function(r) if (nrow(r$frame)) min(r$frame$value) else Inf, 0.0))
summ <- do.call(rbind, lapply(methods, function(hm) {
  f <- runs[[hm]]$frame
  ok <- if (nrow(f)) f$value <= best + 0.1 else logical(0)   # reached the optimum
  data.frame(method = hm, fits = nfit, completed = nrow(f), n_success = sum(ok),
             success_rate = round(sum(ok) / nfit, 3),
             median_neval_success = if (any(ok)) stats::median(f$neval[ok]) else NA_integer_,
             total_neval = if (nrow(f)) sum(f$neval) else 0L,
             best_value = if (nrow(f)) round(min(f$value), 4) else NA_real_,
             secs = round(runs[[hm]]$secs, 1))
}))
cat("\n==== RESULTS ====\n"); print(summ, row.names = FALSE)
cat(sprintf("\nreference best value: %.4f\n", best))
saveRDS(list(summary = summ, runs = runs, best = best),
        file.path(outdir, sprintf("boehm_experiment_%d.rds", nfit)))
cat("saved:", file.path(outdir, sprintf("boehm_experiment_%d.rds", nfit)), "\n")
