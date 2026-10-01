# -------------------------------------------------------------------------#
# Multicompartment neuron with voltage-gated channels: symmetries at scale
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# A large conductance-based model for symmetryDetection() and
# symmetryReduction(): soma, axon initial segment and a dendritic cable, each
# compartment with its own conductances, leak and calcium pool. Gates are
# Boltzmann functions of the membrane potential, so states enter through
# exp(). The somatic potential is patched, the dendrites are read by a
# voltage-sensitive dye and by a calcium dye. Channel blockers, two
# neuromodulators and current steps make the conditions. Nothing is fitted.
#
# [AUTHOR]
# Simon Beyer
#
# [Date]
# Wed 30 Sep 2026
#
# [Info]
# Channels after Migliore & Shepherd (2002) and Hu et al. (2009): NaT (Nav1.2
# somatodendritic, Nav1.6 in the AIS), Kdr, KA, KM (Kv7, PIP2-dependent), HCN
# (cAMP shifts its activation), CaL with GHK flux, SK on calcium. Forskolin
# raises cAMP; carbachol activates PLC, which depletes the PIP2 pool of each
# compartment, and KM needs PIP2 (Suh & Hille 2002). Initial values are free:
# every condition starts at the same unknown resting state when the step and
# the drugs are applied at t = 0. .nDend sets the size.
# -------------------------------------------------------------------------#

library(dMod2)

.modelname <- "symmetryIonChannel"
.cores     <- detectFreeCores()

# dendritic compartments of the cable; 12 give about 300 parameters
.nDend <- 12L


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Gates and currents
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Boltzmann steady state and a constant time constant; kinetics are shared by
# every compartment of a region, sign = -1 for an inactivation gate
.ENa <- "50"; .EK <- "-90"; .Eh <- "-30"
.zCa <- "0.0749"   # 2F/RT per mV at 37 C
.Cao <- "2000"     # extracellular calcium, uM

gate <- function(x, V, id, sign = 1, shift = NULL) {
  vh  <- paste0("Vh_", id, if (!is.null(shift)) paste0(" + ", shift))
  arg <- if (sign > 0) paste0("(", vh, " - ", V, ")/k_", id)
         else paste0("(", V, " - (", vh, "))/k_", id)
  paste0("(1/(1 + exp(", arg, ")) - ", x, ")/tau_", id)
}

.regionOf <- function(cp) if (cp == "s") "soma" else if (cp == "a") "ais" else "dend"
.comps <- c("s", "a", paste0("d", seq_len(.nDend)))
# channels per region
.chan <- list(
  soma = c("NaT", "Kdr", "KA", "KM", "HCN", "CaL", "SK"),
  ais  = c("NaT", "Kdr", "KM"),
  dend = c("NaT", "KA", "KM", "HCN", "CaL", "SK"))
# kinetic set of each channel's gates: region-specific for NaT and KA
.kin <- function(ch, reg) switch(ch,
  NaT = if (reg == "ais") "NaTa" else "NaTs",
  KA  = if (reg == "dend") "KAd" else "KAs",
  ch)
# the blocker of each channel
.blocker <- c(NaT = "TTX", Kdr = "TEA", KA = "AP4", KM = "XE991", HCN = "ZD",
              CaL = "Nif", SK = "Apa")

f <- c()
for (cp in .comps) {
  reg <- .regionOf(cp)
  V   <- paste0("V_", cp)
  cur <- character(0)
  blk <- function(ch) paste0("(1 - bool_", .blocker[[ch]], ")")
  for (ch in .chan[[reg]]) {
    ks <- .kin(ch, reg)
    g  <- paste0("g", ch, "_", cp, " * ", blk(ch))
    if (ch == "NaT") {
      f[paste0("m_", cp)] <- gate(paste0("m_", cp), V, paste0(ks, "m"))
      f[paste0("h_", cp)] <- gate(paste0("h_", cp), V, paste0(ks, "h"), sign = -1)
      cur <- c(cur, paste0(g, " * m_", cp, "^3 * h_", cp, " * (", V, " - ", .ENa, ")"))
    }
    if (ch == "Kdr") {
      f[paste0("n_", cp)] <- gate(paste0("n_", cp), V, "Kdrn")
      cur <- c(cur, paste0(g, " * n_", cp, "^4 * (", V, " - (", .EK, "))"))
    }
    if (ch == "KA") {
      f[paste0("a_", cp)] <- gate(paste0("a_", cp), V, paste0(ks, "a"))
      f[paste0("b_", cp)] <- gate(paste0("b_", cp), V, paste0(ks, "b"), sign = -1)
      cur <- c(cur, paste0(g, " * a_", cp, " * b_", cp, " * (", V, " - (", .EK, "))"))
    }
    if (ch == "KM") {
      f[paste0("w_", cp)] <- gate(paste0("w_", cp), V, "KMw")
      pip <- paste0("PIP2_", cp)
      cur <- c(cur, paste0(g, " * w_", cp, " * ", pip, "^2/(K_KM_PIP2^2 + ", pip, "^2) * (",
                           V, " - (", .EK, "))"))
      # the membrane pool of the compartment: resynthesis and PLC hydrolysis
      f[pip] <- paste0("k_PIP2_", cp, " * (PIP2tot_", cp, " - ", pip, ") - (k_PLC_", cp,
                       " + k_PLC_CCh_", cp, " * bool_CCh) * ", pip)
    }
    if (ch == "HCN") {
      f[paste0("y_", cp)] <- gate(paste0("y_", cp), V, "HCNy", sign = -1,
                                  shift = "dVh_HCN_cAMP * cAMP")
      cur <- c(cur, paste0(g, " * y_", cp, " * (", V, " - (", .Eh, "))"))
    }
    if (ch == "CaL") {
      f[paste0("d_", cp)] <- gate(paste0("d_", cp), V, "CaLd")
      f[paste0("f_", cp)] <- gate(paste0("f_", cp), V, "CaLf", sign = -1)
      ghk <- paste0(V, " * (Ca_", cp, " - ", .Cao, " * exp(-", .zCa, " * ", V, "))/(1 - exp(-",
                    .zCa, " * ", V, "))")
      cur <- c(cur, paste0("p", ch, "_", cp, " * ", blk(ch), " * d_", cp, " * f_", cp, " * ", ghk))
      .iCa <- paste0("pCaL_", cp, " * ", blk(ch), " * d_", cp, " * f_", cp, " * ", ghk)
    }
    if (ch == "SK")
      cur <- c(cur, paste0(g, " * Ca_", cp, "^2/(K_SK^2 + Ca_", cp, "^2) * (", V, " - (", .EK, "))"))
  }
  cur <- c(cur, paste0("gL_", cp, " * (", V, " - EL_", cp, ")"))
  # cable: soma - AIS, soma - d1 - d2 - ...
  nb <- switch(cp, s = c("a", if (.nDend) "d1"), a = "s",
               { i <- as.integer(sub("d", "", cp)); c(if (i == 1L) "s" else paste0("d", i - 1L),
                                                         if (i < .nDend) paste0("d", i + 1L)) })
  cpl <- vapply(nb, function(o) {
    gc <- paste0("gc_", paste(sort(c(cp, o)), collapse = ""))
    paste0(gc, " * (V_", o, " - ", V, ")")
  }, "")
  f[V] <- paste0("(", paste(c(if (cp == "s") "I_inj", cpl), collapse = " + "), " - (",
                 paste(cur, collapse = " + "), "))/C_", cp)
  # calcium pool: GHK influx, PMCA, NCX and a buffer with a conserved total
  if ("CaL" %in% .chan[[reg]]) {
    Ca <- paste0("Ca_", cp); CaB <- paste0("CaB_", cp)
    f[Ca] <- paste0("-kF_", cp, " * ", .iCa,
                    " - Vmax_PMCA * ", Ca, "/(Km_PMCA + ", Ca, ") - k_NCX_", cp, " * ", Ca,
                    " + j_leak_", cp,
                    " - kon_B * ", Ca, " * (Btot_", cp, " - ", CaB, ") + koff_B * ", CaB)
    f[CaB] <- paste0("kon_B * ", Ca, " * (Btot_", cp, " - ", CaB, ") - koff_B * ", CaB)
  }
}

# forskolin raises cAMP, which diffuses fast enough to be one pool
f["cAMP"] <- "k_cAMP + k_cAMP_FSK * bool_FSK - k_PDE * cAMP"

# drugs and the injected current are held for the whole record
.drugs <- c(unique(unname(.blocker)), "FSK", "CCh")
for (b in .drugs) f[paste0("bool_", b)] <- "0"
f["I_inj"] <- "0"
f <- as.eqnvec(f)
.states <- names(f)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Observables
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# somatic patch in mV; VSD with one gain and offset for the camera; a calcium
# dye of known affinity (0.345 uM) with one loading scale per region of interest
.vsd <- paste0("d", unique(c(1L, ceiling(.nDend / 2), .nDend)))
.caROI <- c("s", paste0("d", seq_len(.nDend)))
g <- c(Vsoma_obs = "V_s",
       setNames(paste0("s_VSD * V_", .vsd, " + off_VSD"), paste0("VSD_", .vsd, "_obs")),
       setNames(paste0("sF_", .caROI, " * Ca_", .caROI, "/(Ca_", .caROI, " + 0.345)"),
                paste0("Fluo_", .caROI, "_obs")))
g <- as.eqnvec(g)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Conditions
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# two current steps, each drug at the first step
.conds <- c("step1", "step2", paste0("step1_", .drugs))
conditions <- data.frame(row.names = .conds,
                         I_inj = ifelse(.conds == "step2", 0.4, 0.2))
for (b in .drugs) conditions[[paste0("bool_", b)]] <- as.numeric(.conds == paste0("step1_", b))

.pars <- setdiff(cOde::getSymbols(c(as.character(f), as.character(g))), c(.states, "exp"))
cat(sprintf("%d states, %d parameters, %d observables, %d conditions\n",
            length(setdiff(.states, c(paste0("bool_", .drugs), "I_inj"))), length(.pars),
            length(g), nrow(conditions)))

# potentials, half-activations, reversal and dye offsets take either sign
.real <- c(grep("^(V_|Vh_|EL_|dVh_)", c(.states, .pars), value = TRUE), "off_VSD")
.positive <- setdiff(c(.states, .pars), c(.real, paste0("bool_", .drugs), "I_inj"))


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Symmetries
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
t0 <- Sys.time()
idResult <- symmetryDetection(f, g, forcings = NULL, conditions = conditions,
                              positive = .positive, gaugePreference = NULL,
                              reconstruct = TRUE, cores = .cores, verbose = TRUE)
tDetect <- Sys.time() - t0
tDetect
summary(idResult)

t0 <- Sys.time()
idReduction <- symmetryReduction(idResult, positive = .positive)
tReduce <- Sys.time() - t0
tReduce
idReduction
