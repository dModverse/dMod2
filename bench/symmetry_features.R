# -------------------------------------------------------------------------#
# symmetryDetection() and symmetryReduction(): one model per feature
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# Small models that exercise one option each: conserved moieties, steady
# states, events and conditions, free exponents, states inside exponentials,
# the reduction stages and their degree bounds, and a 26-state signalling
# model at the end. Each call prints its result; compare the two engines
# where both apply.
#
# [AUTHOR]
# Simon Beyer
#
# [Date]
# Thu 08 Oct 2026
#
# [Info]
# Needs the Python package symident. Nothing is compiled or written to disk.
# The last section runs for minutes.
# -------------------------------------------------------------------------#

library(dMod2)

.cores <- min(6, detectFreeCores())


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Conserved moiety: reduceCQ, scalingsOnly and the symbolic engine
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
reactions <- eqnlist() |>
  addReaction("A", "B", "k1 * A") |>
  addReaction("B", "A", "k2 * B") |>
  customTotals(list(totC = "A+B"))
g <- eqnvec(Aobs = "alpha * A")

out <- symmetryDetection(reactions, g, reduceCQ = FALSE, reconstruct = TRUE)
out <- symmetryDetection(reactions, g, reduceCQ = TRUE, reconstruct = TRUE)
summary(out)
out <- symmetryDetection(reactions, g, scalingsOnly = TRUE, reduceCQ = FALSE)
out <- symmetryDetection(reactions, g, reduceCQ = TRUE, symEngine = "symbolic")


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Steady state: equilibrate against an explicit steady-state trafo
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
out <- symmetryDetection(reactions, g, equilibrate = TRUE, reconstruct = TRUE)

# The symbolic engine needs the resting state as trafo
mysteadies <- steadyStates(reactions)
out <- symmetryDetection(reactions, g, trafo = as.eqnvec(mysteadies), reduceCQ = FALSE,
                         reconstruct = TRUE, symEngine = "symbolic")

fss <- eqnvec(x = "b - a*x")
gss <- eqnvec(y = "s*x")
out <- symmetryDetection(fss, gss, trafo = eqnvec(x = "b/a"), reconstruct = TRUE)
out <- symmetryDetection(fss, gss, equilibrate = TRUE, reconstruct = TRUE)

# A known dose
dose <- addEvent(eventlist(), var = "x", time = 0, value = "dose",
                 method = "replace")
out <- symmetryDetection(fss, gss, events = dose,
                         conditions = data.frame(dose = 2, row.names = "stim"))
out$identifiable


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Lumped parameters
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
f <- eqnvec(S = "-kcat*Etot*S/(Km + S)")
g.enz <- eqnvec(y = "s*S")
out <- symmetryDetection(f, g.enz, reconstruct = TRUE)
out <- symmetryDetection(f, g.enz, scalingsOnly = TRUE)

gene <- eqnvec(m = "ktx - dm*m", p = "ktl*m - dp*p")
out <- symmetryDetection(gene, eqnvec(y = "p"), symEngine = "symbolic")


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Conditions stacked over a switch
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
fu <- eqnvec(A = "-(k1 + u*k2) * A", u = "0")
events <- addEvent(eventlist(), var = "u", time = -1, value = "var_u",
                   method = "replace")
cond.grid <- data.frame(var_u = c(0, 1), row.names = c("ctrl", "stim"))
out <- symmetryDetection(fu, eqnvec(y = "A"), events = events, conditions = cond.grid)
summary(out)
out <- symmetryDetection(fu, eqnvec(y = "A"), events = events, conditions = cond.grid,
                         symEngine = "symbolic")


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Free Hill and power exponents
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
hill <- eqnlist() |>
  addReaction("0",  "FB", "k_pr_FB")                     |>
  addReaction("FB", "0",  "d_FB * FB")                   |>
  addReaction("0",  "x",  "k_pr_x * K^n / (K^n + FB^n)") |>
  addReaction("x",  "0",  "d_x * x")
out <- symmetryDetection(hill, eqnvec(xobs = "scale * x"), reconstruct = TRUE)
summary(out)

ev <- addEvent(eventlist(), var = "u", time = -1, value = "1", method = "replace")
out <- symmetryDetection(eqnvec(x = "kpr - dp*x^q + kin*u", u = "0"),
                         eqnvec(y = "s*x"), equilibrate = TRUE, forcings = "u", events = ev,
                         conditions = data.frame(var = 1, row.names = "stim"),
                         reconstruct = TRUE)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Partially observed cascade with conserved totals
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
egfr <- eqnlist() |>
  addReaction("EGF + EGFR", "EGF_EGFR", "k_bind * EGF * EGFR") |>
  addReaction("EGF_EGFR", "EGF + EGFR", "k_unbind * EGF_EGFR") |>
  addReaction("MEK", "pMEK", "k_phos_MEK * EGF_EGFR * MEK") |>
  addReaction("pMEK", "MEK", "k_dephos_MEK * pMEK") |>
  addReaction("ERK", "pERK", "k_phos_ERK * pMEK * ERK") |>
  addReaction("pERK", "ERK", "k_dephos_ERK * pERK")
egfr <- customTotals(egfr, list(
  totalEGF  = "EGF + EGF_EGFR", totalEGFR = "EGFR + EGF_EGFR",
  totalMEK  = "MEK + pMEK",     totalERK  = "ERK + pERK"))
out <- symmetryDetection(egfr, eqnvec(pMEK_obs = "scale_pMEK * pMEK",
                                      pERK_obs = "scale_pERK * pERK"),
                         reduceCQ = TRUE, reconstruct = TRUE)

# Steady-state trafo, a stimulus event, and the reduction
egf.event <- addEvent(eventlist(), var = "EGF", time = 0, value = "1", method = "add")
egf <- symmetryDetection(egfr, eqnvec(pMEK_obs = "pMEK", pERK_obs = "pERK"),
                         events = egf.event, trafo = as.eqnvec(steadyStates(egfr)),
                         reduceCQ = FALSE, reconstruct = TRUE)
redEgf <- symmetryReduction(egf, verbose = TRUE)
summary(redEgf)
redEgf$trafo


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Events: segments, gaps and coupling across them
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
f <- eqnvec(R = "kpr - kdg*R + kon*u*R", u = "0")
gR <- eqnvec(y = "scale*R")
ev <- addEvent(eventlist(), var = "u", time = 0, value = "init_u", method = "replace")
cg <- data.frame(init_u = 1, row.names = "Ctrl")
out <- symmetryDetection(f, gR, equilibrate = TRUE,
                         events = ev, conditions = cg, forcings = "u")

ev2 <- ev |> addEvent(var = "u", time = 60, value = "0", method = "replace")
out <- symmetryDetection(f, gR, equilibrate = TRUE,
                         events = ev2, conditions = cg, forcings = "u")
summary(out)

f2 <- eqnvec(x = "kpr/(1 + kinh*inh) - kdeg*x + kstim*stim", inh = "0", stim = "0")
g2 <- eqnvec(y = "s*x")
ev3 <- eventlist() |>
  addEvent(var = "inh",  time = -30, value = "1", method = "replace") |>
  addEvent(var = "stim", time = 0,   value = "1", method = "replace")
out <- symmetryDetection(f2, g2, equilibrate = TRUE,
                         events = ev3, forcings = c("inh", "stim"))

# Both events at the same time
ev.collapsed <- eventlist() |>
  addEvent(var = "inh",  time = 0, value = "1", method = "replace") |>
  addEvent(var = "stim", time = 0, value = "1", method = "replace")
out <- symmetryDetection(f2, g2, equilibrate = TRUE,
                         events = ev.collapsed, forcings = c("inh", "stim"),
                         reconstruct = TRUE)

# Unobserved pre-window against an observed window only
fcpl <- eqnvec(x = "kpr/(1 + kinh*inh) - kdeg*x", inh = "0", obsw = "0")
gcpl <- eqnvec(y = "s * obsw * x")
ev.coupled <- eventlist() |>
  addEvent(var = "inh",  time = -30, value = "1", method = "replace") |>
  addEvent(var = "inh",  time = 0,   value = "0", method = "replace") |>
  addEvent(var = "obsw", time = 0,   value = "1", method = "replace")
out <- symmetryDetection(fcpl, gcpl, equilibrate = TRUE,
                         events = ev.coupled, forcings = "inh", reconstruct = TRUE)
summary(out)

ev.cond <- addEvent(eventlist(), var = "obsw", time = 0, value = "1", method = "replace")
out <- symmetryDetection(fcpl, gcpl, equilibrate = TRUE,
                         events = ev.cond, forcings = "inh", reconstruct = TRUE)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Reduction: fixed coordinates, stages and degree bounds
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
res <- symmetryDetection(gene, eqnvec(y = "p"), reconstruct = TRUE)
symmetryReduction(res, fixed = "ktl")

f2 <- eqnvec(A = "k1 + u*k2 - kdeg*A")
res2 <- symmetryDetection(f2, eqnvec(y = "A"), fixed = "u", reconstruct = TRUE)
red2 <- symmetryReduction(res2)
red2
symmetryReduction(res2, dPoly = 0L, dDarboux = 0L, dExp = 0L)
symmetryReduction(res2, reportZeroCompatibility = TRUE)$zeroCompatibility

f3 <- eqnvec(A = "-k1*A + k2*B", B = "k1*A - k2*B")
g3 <- eqnvec(y = "alpha*A")
res3 <- symmetryDetection(f3, g3, reconstruct = TRUE)
red3 <- symmetryReduction(res3)
symmetryDetection(f3, g3, trafo = red3$trafo)$identifiable

f4 <- eqnvec(x = "-(b - a)/(a*b)*x")
res4 <- symmetryDetection(f4, eqnvec(y = "x"), reconstruct = TRUE)
red4 <- symmetryReduction(res4)
red4$blocks[[1]]$stage
red4$blocks[[1]]$invariants
symmetryReduction(res4, dDarboux = 1L, separable = FALSE)$blocks[[1]]$invariants
symmetryReduction(res4, dDarboux = 0L, separable = FALSE)$blocks[[1]]$invariants


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# States inside exponentials
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
f <- eqnvec(A = "kA - kd*A", B = "kB - exp(A)*B")
g <- eqnvec(y = "s*B")
out <- symmetryDetection(f, g, reconstruct = TRUE)
out <- symmetryDetection(f, g, scalingsOnly = TRUE)

fE <- eqnvec(A = "kA - kd*A", B = "kB - kE*exp(A)*B")
obs <- symmetryDetection(fE, g, reconstruct = TRUE)
red <- symmetryReduction(obs, positive = c("B", "kB", "kd", "kE", "s"))
red
symmetryDetection(fE, g, trafo = red$trafo)$identifiable
out <- symmetryDetection(fE, g, trafo = eqnvec(A = "0"))

out <- symmetryDetection(eqnvec(A = "kA - kd*A", B = "kB - kE*exp10(A)*B"), g,
                         reconstruct = TRUE)

ev <- eventlist() |>
  addEvent(var = "B", time = 0, value = "0", method = "add") |>
  addEvent(var = "A", time = 5, value = "dA", method = "add")
out <- symmetryDetection(fE, g, events = ev, reconstruct = TRUE)

# Hodgkin-Huxley, voltage observed
an <- "0.01*(V + 55)/(1 - exp(-(V + 55)/10))"
bn <- "0.125*exp(-(V + 65)/80)"
am <- "0.1*(V + 40)/(1 - exp(-(V + 40)/10))"
bm <- "4*exp(-(V + 65)/18)"
ah <- "0.07*exp(-(V + 65)/20)"
bh <- "1/(1 + exp(-(V + 35)/10))"
gate <- function(a, b, x) sprintf("(%s)*(1 - %s) - (%s)*%s", a, x, b, x)
hh <- eqnvec(
  V = "(I - gNa*m^3*h*(V - ENa) - gK*n^4*(V - EK) - gL*(V - EL))/C",
  m = gate(am, bm, "m"),
  h = gate(ah, bh, "h"),
  n = gate(an, bn, "n"))
obs <- symmetryDetection(hh, eqnvec(y = "V"), reconstruct = TRUE)
out <- symmetryDetection(hh, eqnvec(y = "V"), scalingsOnly = TRUE)
red <- symmetryReduction(obs, positive = c("C", "gNa", "gK", "gL", "m", "h", "n"))
red
symmetryDetection(hh, eqnvec(y = "V"), trafo = red$trafo)$identifiable

# Morris-Lecar: tanh() and cosh()
ml <- eqnvec(
  V = "(I - gL*(V - VL) - gCa*(1 + tanh((V - V1)/V2))/2*(V - VCa) - gK*w*(V - VK))/C",
  w = "phi*cosh((V - V3)/(2*V4))*((1 + tanh((V - V3)/V4))/2 - w)")
out <- symmetryDetection(ml, eqnvec(y = "V"), reconstruct = TRUE)

# Boltzmann gate with parameters in the exponent
bz <- eqnvec(x = "(1/(1 + exp((Vh - u)/k)) - x)/tau", u = "a - b*u")
out <- symmetryDetection(bz, eqnvec(y = "s*x"), reconstruct = TRUE)
red <- symmetryReduction(out)
red
symmetryDetection(bz, eqnvec(y = "s*x"), trafo = red$trafo)$identifiable


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# TGF-beta/SMAD signalling: 26 states, perturbation conditions
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
addRC <- function(eq, from, to, rate, ...) addReaction(eq, from, to, rate, compartment = "Cell", ...)
addRE <- function(eq, from, to, rate, ...) addReaction(eq, from, to, rate, compartment = "extraCell", ...)

reactions <- eqnlist() |>
  addRC("", "bool_ActD",  "0") |> addRC("", "bool_CHX", "0") |> addRC("", "bool_MG132", "0") |>
  addRE("", "TGFb", "0") |>
  addRC("", "R1mRNA",  "k_pr_R1mRNA * (1 + k_inh_R1mRNA_FB3 * FB3^nhill_R1) * (1 - bool_ActD)") |>
  addRC("", "R2mRNA",  "k_pr_R2mRNA / (1 + k_inh_R2mRNA_FB4 * FB4^nhill_R2) * (1 - bool_ActD)") |>
  addRC("", "FB2mRNA", "k_pr_FB2mRNA * C3^nhill_FB2mRNA / (Km_FB2mRNA^nhill_FB2mRNA + C3^nhill_FB2mRNA) * (1 - bool_ActD)") |>
  addRC("", "FB3mRNA", "k_pr_FB3mRNA * C3^nhill_FB3mRNA / (Km_FB3mRNA^nhill_FB3mRNA + C3^nhill_FB3mRNA) * (1 - bool_ActD)") |>
  addRC("", "FB4mRNA", "k_pr_FB4mRNA * C3^nhill_FB4mRNA / (Km_FB4mRNA^nhill_FB4mRNA + C3^nhill_FB4mRNA) * (1 - bool_ActD)") |>
  addRC("R1mRNA", "", "k_dg_R1mRNA * R1mRNA") |> addRC("R2mRNA", "", "k_dg_R2mRNA * R2mRNA") |>
  addRC("FB2mRNA", "", "k_dg_FB2 * FB2mRNA") |> addRC("FB3mRNA", "", "k_dg_FB3 * FB3mRNA") |>
  addRC("FB4mRNA", "", "k_dg_FB4 * FB4mRNA") |>
  addRC("", "R1",  "k_pr_R1 * R1mRNA * (1 - bool_CHX)") |>
  addRC("", "R2",  "k_pr_R2 * R2mRNA * (1 - bool_CHX)") |>
  addRC("", "FB2", "k_pr_FB2 * FB2mRNA * (1 - bool_CHX)") |>
  addRC("", "FB3", "k_pr_FB3 * FB3mRNA * (1 - bool_CHX)") |>
  addRC("", "FB4", "k_pr_FB4 * FB4mRNA * (1 - bool_CHX)") |>
  addRC("R1", "", "k_dg_R1 * R1 * (1 - bool_MG132)") |> addRC("R2", "", "k_dg_R2 * R2 * (1 - bool_MG132)") |>
  addRC("FB2", "", "k_dg_FB2 * FB2 * (1 - bool_MG132)") |> addRC("FB3", "", "k_dg_FB3 * FB3 * (1 - bool_MG132)") |>
  addRC("FB4", "", "k_dg_FB4 * FB4 * (1 - bool_MG132)") |>
  addRC("R1 + R2", "R1_R2", "k_act_R1_R2 * R1 * R2") |>
  addRC("R1_R2", "R1 + R2", "k_deact_R1_R2 * R1_R2") |>
  addRC("R1_R2", "", "k_dg_R1_R2 * R1_R2 * (1 - bool_MG132)") |>
  addRC("R2 + TGFb", "R2_TGFb", "k_act_R2_TGFb * TGFb * R2 / (km_R2 + R2 + TGFb)", rateCompartment = "Cell") |>
  addRC("R2_TGFb", "R2 + TGFb", "k_deact_R2_TGFb * R2_TGFb") |>
  addRC("R2_TGFb", "R2_TGFb_int", "k_int_R2_TGFb * R2_TGFb") |>
  addRC("R2_TGFb_int", "TGFb", "k_decay_R2_TGFb_int * R2_TGFb_int") |>
  addRC("R2_TGFb_int", "", "k_dg_R2_TGFb_int * R2_TGFb_int * (1 - bool_MG132)") |>
  addRC("R1 + TGFb", "R1_TGFb", "k_act_R1_TGFb * TGFb * R1 / (km_R1 + R1 + TGFb)", rateCompartment = "Cell") |>
  addRC("R1_TGFb", "R1 + TGFb", "k_deact_R1_TGFb * R1_TGFb") |>
  addRC("R1_TGFb", "R1_TGFb_int", "k_int_R1_TGFb * R1_TGFb") |>
  addRC("R1_TGFb_int", "TGFb", "k_decay_R1_TGFb_int * R1_TGFb_int") |>
  addRC("R1_TGFb_int", "", "k_dg_R1_TGFb_int * R1_TGFb_int * (1 - bool_MG132)") |>
  addRC("R1 + R2_TGFb", "R1_R2_TGFb", "k_act_R1_R2_TGFb * R1 * R2_TGFb") |>
  addRC("R2 + R1_TGFb", "R1_R2_TGFb", "k_act_R2_R1_TGFb * R2 * R1_TGFb") |>
  addRC("R1_R2 + TGFb", "R1_R2_TGFb", "k_act_R1_R2_TGFb_direct * TGFb * R1_R2 / (km_R1_R2 + R1_R2 + TGFb)", rateCompartment = "Cell") |>
  addRC("R1_R2_TGFb", "", "(k_dg_R1_R2_TGFb + k_dg_R1_R2_TGFb_FB1 * FB2) * R1_R2_TGFb * (1 - bool_MG132)") |>
  addRC("Smad2", "pSmad2", "(k_phospho_pS2 / (1 + k_inh_pSmad2_FB2 * FB2)) * Smad2 * R1_R2") |>
  addRC("Smad2", "pSmad2", "(k_phospho_pS2 / (1 + k_inh_pSmad2_FB1 * FB2)) * Smad2 * R1_R2_TGFb") |>
  addRC("pSmad2", "Smad2", "k_dephos_S2 * pSmad2") |>
  addRC("Smad3", "pSmad3", "(k_phospho_pS3 / (1 + k_inh_pSmad3_FB2 * FB2)) * Smad3 * R1_R2") |>
  addRC("Smad3", "pSmad3", "(k_phospho_pS3 / (1 + k_inh_pSmad3_FB2 * FB2)) * Smad3 * R1_R2_TGFb") |>
  addRC("pSmad3", "Smad3", "k_dephos_S3 * pSmad3") |>
  addRC("pSmad2 + pSmad3 + Smad4", "C3", "k_form_S4Coip * pSmad2 * pSmad3 * Smad4") |>
  addRC("C3", "Smad2 + pSmad3 + Smad4", "k_dissolve_C3_dp2 * C3") |>
  addRC("C3", "pSmad2 + Smad3 + Smad4", "k_dissolve_C3_dp3 * C3")
reactions$compartments$Cell$volume      <- "1"
reactions$compartments$extraCell$volume <- "volumeEC"
reactions <- customTotals(reactions, list(totalSMAD2 = "Smad2 + pSmad2 + C3",
                                          totalSMAD3 = "Smad3 + pSmad3 + C3",
                                          totalSMAD4 = "Smad4 + C3"))

observables <- eqnvec(
  R1_obs = "scale_R1 * R1", R2_obs = "scale_R2 * R2",
  pSmad2_obs = "scale_pSmad2 * (pSmad2 + C3)", pSmad3_obs = "scale_pSmad3 * (pSmad3 + C3)",
  TSmad2_obs = "scale_TSmad2 * (Smad2 + pSmad2 + C3)", TSmad3_obs = "scale_TSmad3 * (Smad3 + pSmad3 + C3)",
  Smad4_CoIP_obs = "scale_CoIP * C3", TGFBR1_mRNA_obs = "scale_R1mRNA * R1mRNA",
  TGFBR2_mRNA_obs = "scale_R2mRNA * R2mRNA",
  TGFb_obs = "scale_TGFb * TGFb")

events <- eventlist() |>
  addEvent(var = "TGFb",       time = 0,   value = "init_TGFb",      method = "replace") |>
  addEvent(var = "bool_CHX",   time = -30, value = "var_bool_CHX",   method = "replace") |>
  addEvent(var = "bool_MG132", time = -30, value = "var_bool_MG132", method = "replace") |>
  addEvent(var = "bool_ActD",  time = -30, value = "var_bool_ActD",  method = "replace")

# One condition per perturbation; knockdowns rename a synthesis rate
cond.grid <- data.frame(Pertubation = c("Ctrl", "ActD", "CHX", "MG132", "R1Knd", "R2Knd"),
                        init_TGFb = 1, stringsAsFactors = FALSE)
cond.grid$var_bool_ActD  <- ifelse(cond.grid$Pertubation == "ActD",  1, 0)
cond.grid$var_bool_CHX   <- ifelse(cond.grid$Pertubation == "CHX",   1, 0)
cond.grid$var_bool_MG132 <- ifelse(cond.grid$Pertubation == "MG132", 1, 0)
cond.grid$k_pr_R1mRNA <- ifelse(cond.grid$Pertubation == "R1Knd", "k_pr_R1mRNA_R1Knd", "k_pr_R1mRNA")
cond.grid$k_pr_R2mRNA <- ifelse(cond.grid$Pertubation == "R2Knd", "k_pr_R2mRNA_R2Knd", "k_pr_R2mRNA")
rownames(cond.grid) <- cond.grid$Pertubation
cond.grid$Pertubation <- NULL

cond.trafo <- eqnvec() |>
  define("x~x", x = getParameters(reactions, events)) |>
  branch(table = cond.grid, apply = "insert")

out <- symmetryDetection(
  reactions, observables, events = events, trafo = cond.trafo,
  forcings = c("bool_ActD", "bool_CHX", "bool_MG132", "TGFb"),
  equilibrate = TRUE, reduceCQ = TRUE, reconstruct = TRUE,
  cores = .cores)
summary(out)
