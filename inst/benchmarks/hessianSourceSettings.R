# -------------------------------------------------------------------------#
# Catalogue of trust() Hessian-source variants used by the benchmarks
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# One definition of the optimiser settings compared by bench_hessianSource.R,
# the Optimisation vignette and dev/optimisation/bachmann_hessianSource.R. Each
# selects the subset it can afford to run.
#
# [AUTHOR]
# Simon Beyer
#
# [Date]
# Sat 06 Sep 2026
#
# [Info]
# Keys carry no spaces or commas: they travel as literals into the generated
# SLURM array script. `hessianSourceLabels` holds the display form for tables.
#
# Usage:
#   source(system.file("benchmarks", "hessianSourceSettings.R", package = "dMod2"))
#   settings <- hessianSourceSettings[c("gn", "gn_bfgs", "bfgs_id", "sr1_id")]
# -------------------------------------------------------------------------#

hessianSourceSettings <- list(
  gn          = list(hessianMethod = "gn"),
  gn_bfgs     = list(hessianMethod = "gn", hessianFallback = "bfgs"),
  bfgs_gn     = list(hessianMethod = "bfgs"),
  bfgs_id     = list(hessianMethod = "bfgs",
                     qnControl = list(hessianInit = "identity")),
  sr1_gn      = list(hessianMethod = "sr1"),
  sr1_id      = list(hessianMethod = "sr1",
                     qnControl = list(hessianInit = "identity")),
  bfgs_id_m10 = list(hessianMethod = "bfgs",
                     qnControl = list(hessianInit = "identity",
                                      qnMemory = 10L)),
  gn_sr1      = list(hessianMethod = "gn", hessianFallback = "sr1"),
  gn_bfgs_x3  = list(hessianMethod = "gn", hessianFallback = "bfgs",
                     fallbackLimit = 3L),
  sr1_id_gn   = list(hessianMethod = "sr1", hessianFallback = "gn",
                     qnControl = list(hessianInit = "identity")),
  bfgs_id_gn  = list(hessianMethod = "bfgs", hessianFallback = "gn",
                     qnControl = list(hessianInit = "identity")))

hessianSourceLabels <- c(
  gn          = "gn",
  gn_bfgs     = "gn -> bfgs",
  bfgs_gn     = "bfgs, gn",
  bfgs_id     = "bfgs, id",
  sr1_gn      = "sr1, gn",
  sr1_id      = "sr1, id",
  bfgs_id_m10 = "bfgs, id, m10",
  gn_sr1      = "gn -> sr1",
  gn_bfgs_x3  = "gn <-> bfgs, x3",
  sr1_id_gn   = "sr1, id -> gn",
  bfgs_id_gn  = "bfgs, id -> gn")

stopifnot(setequal(names(hessianSourceSettings), names(hessianSourceLabels)))
