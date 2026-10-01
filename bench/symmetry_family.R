# -------------------------------------------------------------------------#
# A family of signalling networks: symmetryDetection() and symmetryReduction()
# against model size
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# Phosphoproteomic networks of .nMod modules, each a receptor with its ligand,
# a cascade of doubly phosphorylated kinases with their phosphatases, a target
# gene whose protein degrades the receptor, and a factor that the gene product
# makes and the cell both degrades and secretes. Modules cross-talk through the
# last kinase of one module activating the first of the next. Kinase inhibitors
# and the ligands make the conditions, every condition starts at its resting
# state. The secreted factor gives each module one general direction.
#
# [AUTHOR]
# Simon Beyer
#
# [Date]
# Thu 01 Oct 2026
#
# [Info]
# One module has 17 states and 28 parameters. Timings are the minimum over
# .reps runs. Needs dMod2 from devel-symmetry.
# -------------------------------------------------------------------------#

library(dMod2)

.cores <- detectFreeCores()
.depth <- 3L                 # kinases per cascade
.sizes <- c(2L, 4L, 8L)      # modules
.reps  <- 1L


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Generator
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
familyModel <- function(nMod, depth = .depth) {
  r <- eqnlist()
  obs <- character(0)
  pools <- list()
  for (m in seq_len(nMod)) {
    M <- function(x) paste0(x, "_", m)
    K <- function(j, form = "") paste0(form, "K", j, "_", m)
    r <- r |>
      addReaction("", M("R"), paste0("ks_R_", m)) |>
      addReaction(M("R"), "", paste0("(kd_R_", m, " + kfb_", m, " * ", M("Gp"), ") * ", M("R"))) |>
      addReaction(paste(M("L"), "+", M("R")), M("C"), paste0("kon_", m, " * ", M("L"), " * ", M("R"))) |>
      addReaction(M("C"), paste(M("L"), "+", M("R")), paste0("koff_", m, " * ", M("C"))) |>
      addReaction(M("C"), "", paste0("kint_", m, " * ", M("C")))
    # doubly phosphorylated kinases, the first one inhibited
    for (j in seq_len(depth)) {
      act <- if (j == 1L) M("C") else K(j - 1L, "pp")
      if (j == 1L && m > 1L) act <- paste0("(", act, " + kx_", m, " * ", paste0("ppK", depth, "_", m - 1L), ")")
      inh <- if (j == 1L) paste0(" * (1 - bool_inh_", m, ")") else ""
      r <- r |>
        addReaction(K(j), K(j, "p"), paste0("kp1_", j, "_", m, " * ", act, " * ", K(j), inh)) |>
        addReaction(K(j, "p"), K(j, "pp"), paste0("kp2_", j, "_", m, " * ", act, " * ", K(j, "p"), inh)) |>
        addReaction(K(j, "pp"), K(j, "p"), paste0("kdp2_", j, "_", m, " * ", K(j, "pp"))) |>
        addReaction(K(j, "p"), K(j), paste0("kdp1_", j, "_", m, " * ", K(j, "p")))
      pools[[paste0("tK", j, "_", m)]] <- paste(K(j), "+", K(j, "p"), "+", K(j, "pp"))
    }
    # target gene, its protein on the receptor, and the factor it makes: degraded,
    # secreted and cleared as a dimer
    r <- r |>
      addReaction("", M("Gm"), paste0("(kb_G_", m, " + kt_G_", m, " * ", K(depth, "pp"), ")")) |>
      addReaction(M("Gm"), "", paste0("kd_Gm_", m, " * ", M("Gm"))) |>
      addReaction("", M("Gp"), paste0("ktl_G_", m, " * ", M("Gm"))) |>
      addReaction(M("Gp"), "", paste0("kd_Gp_", m, " * ", M("Gp"))) |>
      addReaction("", M("S"), paste0("ks_S_", m, " * ", M("Gp"))) |>
      addReaction(M("S"), "", paste0("(kd_S_", m, " + kq_S_", m, " * ", M("S"), ") * ", M("S"))) |>
      addReaction(M("S"), M("Sx"), paste0("ksec_S_", m, " * ", M("S"))) |>
      addReaction(M("Sx"), "", paste0("kcl_S_", m, " * ", M("Sx"))) |>
      addReaction("", paste0("bool_inh_", m), "0")
    obs <- c(obs,
      setNames(c(paste0("pK1_", m, " + ppK1_", m), K(depth, "pp"),
                 paste0("s_Gm_", m, " * ", M("Gm")), M("Gp"), M("Sx")),
               paste0(c("pK1", paste0("ppK", depth), "Gm", "Gp", "Sx"), "_", m, "_obs")))
  }
  r <- customTotals(r, pools)
  list(f = r, g = as.eqnvec(obs))
}

familyConditions <- function(nMod) {
  cn <- c("ctrl", "stim", paste0("stim_inh", seq_len(nMod)))
  cond <- data.frame(row.names = cn, dose = as.numeric(cn != "ctrl"))
  for (m in seq_len(nMod)) cond[[paste0("var_inh_", m)]] <- as.numeric(cn == paste0("stim_inh", m))
  ev <- eventlist()
  for (m in seq_len(nMod))
    ev <- ev |>
      addEvent(var = paste0("bool_inh_", m), time = -30, value = paste0("var_inh_", m),
               method = "replace") |>
      addEvent(var = paste0("L_", m), time = 0, value = "dose", method = "add")
  list(conditions = cond, events = ev, forcings = paste0("bool_inh_", seq_len(nMod)))
}


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Timings against size
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
timings <- data.frame()
for (n in .sizes) {
  mdl <- familyModel(n)
  cnd <- familyConditions(n)
  tD <- tR <- Inf
  for (rep in seq_len(.reps)) {
    t0 <- Sys.time()
    det <- symmetryDetection(mdl$f, mdl$g, forcings = cnd$forcings, events = cnd$events,
                             conditions = cnd$conditions, equilibrate = TRUE, reduceCQ = TRUE,
                             gaugePreference = NULL, reconstruct = TRUE, cores = .cores,
                             verbose = FALSE)
    tD <- min(tD, as.numeric(Sys.time() - t0, units = "secs"))
    t0 <- Sys.time()
    red <- symmetryReduction(det)
    tR <- min(tR, as.numeric(Sys.time() - t0, units = "secs"))
  }
  gen <- Filter(function(d) d$type == "general", det$symmetries)
  timings <- rbind(timings, data.frame(
    modules = n, states = length(mdl$f$states),
    parameters = length(setdiff(getParameters(mdl$f), mdl$f$states)),
    rank = det$rank, dim = det$dim, scalings = length(det$symmetries) - length(gen),
    general = length(gen), reduced = length(red$removed), remaining = length(red$remaining),
    detect_s = tD, reduce_s = tR))
  print(timings)
}

plot(timings$states, timings$detect_s, log = "xy", type = "b", pch = 19,
     xlab = "states", ylab = "time [s]", ylim = range(c(timings$detect_s, timings$reduce_s)))
lines(timings$states, timings$reduce_s, type = "b", pch = 1)
legend("topleft", c("symmetryDetection", "symmetryReduction"), pch = c(19, 1), bty = "n")
