# -------------------------------------------------------------------------#
# Hessian source in a multi-start: gn, hybrid, bfgs
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# Compares the trust-region Hessian sources on Boehm2014 over one shared set of
# starting points. The two metrics disagree, which is the point of the script:
# scored per start a method can look good and still be expensive per hit,
# because it carries hopeless starts a long way before it gives up.
#
# [AUTHOR]
# Simon Beyer
#
# [Date]
# Fri 05 Sep 2026
#
# [Info]
# The model is the one of inst/examples/example_Boehm_JProteomeRes2014.R, built
# here so the script stands alone. The Gauss-Newton-seeded arms dominate the run
# time, because they spend the whole iteration cap on most starts.
# -------------------------------------------------------------------------#

library(dMod2)

.outdir <- file.path(tempdir(), "bench_hessianSource")
dir.create(.outdir, recursive = TRUE, showWarnings = FALSE)
.cores  <- 20
.nstart <- 100

data(boehm)
mydataL <- as.datalist(boehm)

epo <- "1.25e-7*exp(-Epo_degradation_BaF3*time)"
reactions <- eqnlist() |>
  assignCompartment(nucpApA = "nuc", nucpApB = "nuc", nucpBpB = "nuc",
                    volume = "0.45") |>
  addReaction("2*STAT5A", "pApA", paste0("k_phos*", epo, "*STAT5A^2"),
              compartment = "cyt") |>
  addReaction("STAT5A + STAT5B", "pApB", paste0("k_phos*", epo, "*STAT5A*STAT5B"),
              compartment = "cyt") |>
  addReaction("2*STAT5B", "pBpB", paste0("k_phos*", epo, "*STAT5B^2"),
              compartment = "cyt") |>
  addReaction("pApA", "nucpApA", "k_imp_homo*pApA",     compartment = "cyt") |>
  addReaction("pApB", "nucpApB", "k_imp_hetero*pApB",   compartment = "cyt") |>
  addReaction("pBpB", "nucpBpB", "k_imp_homo*pBpB",     compartment = "cyt") |>
  addReaction("nucpApA", "2*STAT5A", "k_exp_homo*nucpApA",   compartment = "nuc") |>
  addReaction("nucpApB", "STAT5A + STAT5B", "k_exp_hetero*nucpApB",
              compartment = "nuc") |>
  addReaction("nucpBpB", "2*STAT5B", "k_exp_homo*nucpBpB",   compartment = "nuc") |>
  setCompartmentVolume(cyt = "1.4")

# The optimiser can only resolve what the integrator delivers, so the tolerances
# are set explicitly rather than left at their defaults.
myOptionsODE <- list(atol = 1e-8, rtol = 1e-6, maxattemps = 100L, maxsteps = 1e6)
model <- odemodel(reactions, modelname = "boehm_ode", compile = FALSE,
                  outdir = .outdir)
x <- Xs(model, optionsOde = myOptionsODE, optionsSens = myOptionsODE)

observables <- eqnvec(
  pSTAT5A_rel = "(100*pApB + 200*pApA*specC17)/(pApB + STAT5A*specC17 + 2*pApA*specC17)",
  pSTAT5B_rel = "-(100*pApB - 200*pBpB*(specC17 - 1))/((STAT5B*(specC17 - 1) - pApB) + 2*pBpB*(specC17 - 1))",
  rSTAT5A_rel = "(100*pApB + 100*STAT5A*specC17 + 200*pApA*specC17)/(2*pApB + STAT5A*specC17 + 2*pApA*specC17 - STAT5B*(specC17 - 1) - 2*pBpB*(specC17 - 1))")
errorModels <- eqnvec(pSTAT5A_rel = "sd_pSTAT5A_rel",
                      pSTAT5B_rel = "sd_pSTAT5B_rel",
                      rSTAT5A_rel = "sd_rSTAT5A_rel")

g <- Y(observables, x, modelname = "boehm_obs", attach.input = FALSE,
       compile = FALSE, outdir = .outdir)
e <- Y(errorModels, g, modelname = "boehm_err", attach.input = FALSE,
       compile = FALSE, outdir = .outdir)

innerpars <- getParameters(model, g, e)
estimated <- c("Epo_degradation_BaF3", "k_exp_hetero", "k_exp_homo",
               "k_imp_hetero", "k_imp_homo", "k_phos",
               "sd_pSTAT5A_rel", "sd_pSTAT5B_rel", "sd_rSTAT5A_rel")
p <- eqnvec() |>
  define("x~x", x = innerpars) |>
  define("x~0", x = c("pApA", "pApB", "pBpB", "nucpApA", "nucpApB", "nucpBpB")) |>
  insert("STAT5A ~ 207.6*ratio") |>
  insert("STAT5B ~ 207.6 - 207.6*ratio") |>
  insert("ratio ~ 0.693") |>
  insert("specC17 ~ 0.107") |>
  insert("x~exp10(x)", x = estimated) |>
  P(modelname = "boehm_trafo", condition = "Boehm2014", compile = FALSE,
    outdir = .outdir)

compile(g, x, p, e, output = "boehm_bench", cores = 2)

prd       <- g*x*p
outerpars <- getParameters(prd)
pouter    <- structure(rep(-1, length(outerpars)), names = outerpars)
dyn       <- setdiff(outerpars, grep("^sd_", outerpars, value = TRUE))
obj       <- normL2(mydataL, prd, e) + constraintL2(pouter[dyn], sigma = 4)

parlower <- structure(rep(-5, length(outerpars)), names = outerpars)
parupper <- structure(rep( 5, length(outerpars)), names = outerpars)

set.seed(20260905)
starts <- msParframe(pouter, n = .nstart, sd = 3)

# One shared iteration cap for every source, set high enough that a run ends on
# a termination test rather than on the cap wherever it is able to.
run <- function(...)
  as.parframe(mstrust(obj, center = starts, fits = nrow(starts), cores = .cores,
                      rinit = 0.1, rmax = 10, iterlim = 5000,
                      parlower = parlower, parupper = parupper, ...))

frames <- list(
  `gn`            = run(hessianMethod = "gn"),
  `hybrid`        = run(hessianMethod = "hybrid"),
  `bfgs, gn`      = run(hessianMethod = "bfgs"),
  `bfgs, id`      = run(hessianMethod = "bfgs", hessianInit = "identity"),
  `sr1, gn`       = run(hessianMethod = "sr1"),
  `sr1, id`       = run(hessianMethod = "sr1", hessianInit = "identity"),
  `bfgs, id, m10` = run(hessianMethod = "bfgs", hessianInit = "identity",
                        qnMemory = 10L))


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Per start
#
# A start "hits" if it lands within 0.1 of the best value seen anywhere.
# `converged` separates a method that misses the optimum from one that never
# stops looking.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
best <- min(vapply(frames, function(f) min(f$value), 0.0))

perStart <- do.call(rbind, lapply(names(frames), function(id) {
  f <- frames[[id]]
  hits <- sum(f$value <= best + 0.1)
  data.frame(source       = id,
             converged    = sum(f$converged),
             hits         = hits,
             median_neval = stats::median(f$neval),
             total_neval  = sum(f$neval),
             median_value = stats::median(f$value),
             best_value   = min(f$value),
             eval_per_hit = round(sum(f$neval) / hits))
}))
print(perStart, row.names = FALSE)

for (id in names(frames)) { cat(id, ": "); print(table(frames[[id]]$stopReason)) }

# The hybrid against the source it starts from, matched per start. A parframe is
# sorted by value, so the starts are lined up on `index` before comparing.
byIndex <- lapply(frames[c("gn", "hybrid")], function(f) f[order(f$index), ])
stopifnot(identical(byIndex$gn$index, byIndex$hybrid$index))
cat("hybrid: switched on", sum(byIndex$hybrid$qnEval > 0), "of", .nstart,
    "starts, lower value than gn on", sum(byIndex$hybrid$value < byIndex$gn$value),
    ", higher on", sum(byIndex$hybrid$value > byIndex$gn$value), "\n")

print(plotValues(frames[["hybrid"]], tol = 0.1, value < 1e4))
