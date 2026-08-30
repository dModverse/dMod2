# Worked examples for symmetryDetection(); derivations in vignette("Symmetries").

library(dMod2)
setwd(tempdir())

## 1. canonical scaling, all three engines ------------------------------------
reactions <- eqnlist() |>
  addReaction("A", "B", "k1 * A") |>
  addReaction("B", "A", "k2 * B") |>
  customTotals(list(totC = "A+B"))
g <- eqnvec(Aobs = "alpha * A")

out <- symmetryDetection(reactions, g, method = "observability",
                         reduceCQ = FALSE, reconstruct = TRUE)
out <- symmetryDetection(reactions, g, method = "observability",
                         reduceCQ = TRUE, reconstruct = TRUE)
print(out)     # verdict and generators
summary(out)   # adds the computation block
out <- symmetryDetection(reactions, g, method = "polynomial", reduceCQ = FALSE,
                         polynomial = polynomialControl(ansatz = "multi", pMax = 2L))
out <- symmetryDetection(reactions, g, method = "scaling", reduceCQ = FALSE)
out <- symmetryDetection(reactions, g, method = "observability",
                         reduceCQ = TRUE, symEngine = "symbolic")

## 2. closed-form non-monomial direction (reconstruct = TRUE) -----------------
out <- symmetryDetection(reactions, g, method = "observability", equilibrate = TRUE,
                         reconstruct = TRUE)

# symbolic engine needs an explicit resting state; reduceCQ = FALSE keeps the
# moiety relation the steady state parameterises
mysteadies <- steadyStates(reactions)
out <- symmetryDetection(reactions, g, method = "observability",
                         trafo = as.eqnvec(mysteadies), reduceCQ = FALSE,
                         reconstruct = TRUE, symEngine = "symbolic")

## 3. enzyme assay: kcat*Etot lump and the unit scaling -----------------------
f <- eqnvec(S = "-kcat*Etot*S/(Km + S)")
g.enz <- eqnvec(y = "s*S")
out <- symmetryDetection(f, g.enz, method = "observability", reconstruct = TRUE)
out <- symmetryDetection(f, g.enz, method = "scaling")
out <- symmetryDetection(f, g.enz, method = "polynomial")

## 4. rates confined to ktx*ktl = const ---------------------------------------
gene <- eqnvec(m = "ktx - dm*m", p = "ktl*m - dp*p")
out <- symmetryDetection(gene, eqnvec(y = "p"), method = "observability",
                         reconstruct = TRUE)
out <- symmetryDetection(gene, eqnvec(y = "p"), method = "observability",
                         symEngine = "symbolic")

## 5. stacked conditions separate a switch-gated rate -------------------------
fu <- eqnvec(A = "-(k1 + u*k2) * A", u = "0")
events <- addEvent(eventlist(), var = "u", time = -1, value = "var_u",
                   method = "replace")
cond.grid <- data.frame(var_u = c(0, 1), row.names = c("ctrl", "stim"))
out <- symmetryDetection(fu, eqnvec(y = "A"), method = "observability",
                         events = events, conditions = cond.grid)
summary(out)
out <- symmetryDetection(fu, eqnvec(y = "A"), method = "observability",
                         events = events, conditions = cond.grid,
                         symEngine = "symbolic")

## 6. steady-state initial condition ------------------------------------------
fss <- eqnvec(x = "b - a*x")
gss <- eqnvec(y = "s*x")
out <- symmetryDetection(fss, gss, method = "observability",
                         trafo = eqnvec(x = "b/a"), reconstruct = TRUE)
out <- symmetryDetection(fss, gss, method = "observability",
                         equilibrate = TRUE, reconstruct = TRUE)
# a known dose pins the scale
dose <- addEvent(eventlist(), var = "x", time = 0, value = "dose",
                 method = "replace")
out <- symmetryDetection(fss, gss, method = "observability", events = dose,
                         conditions = data.frame(dose = 2, row.names = "stim"))
out$identifiable

## 7. free Hill and power exponents (recast to rational coordinates) ----------
hill <- eqnlist() |>
  addReaction("0",  "FB", "k_pr_FB")                     |>
  addReaction("FB", "0",  "d_FB * FB")                   |>
  addReaction("0",  "x",  "k_pr_x * K^n / (K^n + FB^n)") |>
  addReaction("x",  "0",  "d_x * x")
out <- symmetryDetection(hill, eqnvec(xobs = "scale * x"),
                         method = "observability", reconstruct = TRUE)
summary(out)   # the n-direction carries a log(base) factor

ev <- addEvent(eventlist(), var = "u", time = -1, value = "1", method = "replace")
out <- symmetryDetection(eqnvec(x = "kpr - dp*x^q + kin*u", u = "0"),
                         eqnvec(y = "s*x"), method = "observability",
                         equilibrate = TRUE, forcings = "u", events = ev,
                         conditions = data.frame(var = 1, row.names = "stim"),
                         reconstruct = TRUE)

## 8. EGF/EGFR -> MEK/ERK cascade, partially observed -------------------------
reactions <- eqnlist() |>
  addReaction("EGF + EGFR", "EGF_EGFR", "k_bind * EGF * EGFR") |>
  addReaction("EGF_EGFR", "EGF + EGFR", "k_unbind * EGF_EGFR") |>
  addReaction("MEK", "pMEK", "k_phos_MEK * EGF_EGFR * MEK") |>
  addReaction("pMEK", "MEK", "k_dephos_MEK * pMEK") |>
  addReaction("ERK", "pERK", "k_phos_ERK * pMEK * ERK") |>
  addReaction("pERK", "ERK", "k_dephos_ERK * pERK")
reactions <- customTotals(reactions, list(
  totalEGF  = "EGF + EGF_EGFR", totalEGFR = "EGFR + EGF_EGFR",
  totalMEK  = "MEK + pMEK",     totalERK  = "ERK + pERK"))
observables <- eqnvec(pMEK_obs = "scale_pMEK * pMEK",
                      pERK_obs = "scale_pERK * pERK")
out <- symmetryDetection(reactions, observables, method = "observability",
                         reduceCQ = TRUE, reconstruct = TRUE)
# rank 11 / 15: three scalings plus one non-monomial direction, in closed form

## 9. events split the timeline into exactly propagated segments --------------
f <- eqnvec(R = "kpr - kdg*R + kon*u*R", u = "0")
gR <- eqnvec(y = "scale*R")
ev <- addEvent(eventlist(), var = "u", time = 0, value = "init_u", method = "replace")
cg <- data.frame(init_u = 1, row.names = "Ctrl")
out <- symmetryDetection(f, gR, method = "observability", equilibrate = TRUE,
                         events = ev, conditions = cg, forcings = "u")

# a washout at t = 60 opens a second segment
ev2 <- ev |> addEvent(var = "u", time = 60, value = "0", method = "replace")
out <- symmetryDetection(f, gR, method = "observability", equilibrate = TRUE,
                         events = ev2, conditions = cg, forcings = "u")
summary(out)   # 2 segments and the gap order

# kinh is seen only through the state propagated across the gap
f2 <- eqnvec(x = "kpr/(1 + kinh*inh) - kdeg*x + kstim*stim", inh = "0", stim = "0")
g2 <- eqnvec(y = "s*x")
ev3 <- eventlist() |>
  addEvent(var = "inh",  time = -30, value = "1", method = "replace") |>
  addEvent(var = "stim", time = 0,   value = "1", method = "replace")
out <- symmetryDetection(f2, g2, method = "observability", equilibrate = TRUE,
                         events = ev3, forcings = c("inh", "stim"))

## 10. coupled segments vs. independent conditions ----------------------------
# kinh acts only in an UNOBSERVED pre-window: only the propagated segment sees it
fcpl <- eqnvec(x = "kpr/(1 + kinh*inh) - kdeg*x", inh = "0", obsw = "0")
gcpl <- eqnvec(y = "s * obsw * x")
ev.coupled <- eventlist() |>
  addEvent(var = "inh",  time = -30, value = "1", method = "replace") |>
  addEvent(var = "inh",  time = 0,   value = "0", method = "replace") |>
  addEvent(var = "obsw", time = 0,   value = "1", method = "replace")
out <- symmetryDetection(fcpl, gcpl, method = "observability", equilibrate = TRUE,
                         events = ev.coupled, forcings = "inh", reconstruct = TRUE)
summary(out)   # rank 4/5, gap order 1

# without the pre-window (the condition view) kinh is gone entirely
ev.cond <- addEvent(eventlist(), var = "obsw", time = 0, value = "1", method = "replace")
out <- symmetryDetection(fcpl, gcpl, method = "observability", equilibrate = TRUE,
                         events = ev.cond, forcings = "inh", reconstruct = TRUE)

# collapsing both events onto t = 0 loses kinh again: a [kinh, kstim] confound
ev.collapsed <- eventlist() |>
  addEvent(var = "inh",  time = 0, value = "1", method = "replace") |>
  addEvent(var = "stim", time = 0, value = "1", method = "replace")
out <- symmetryDetection(f2, g2, method = "observability", equilibrate = TRUE,
                         events = ev.collapsed, forcings = c("inh", "stim"),
                         reconstruct = TRUE)

## 11. TGF-beta / SMAD signalling at scale (26 states, Hill feedback) ---------
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

# switches at t = -30, TGFb dose at t = 0
events <- eventlist() |>
  addEvent(var = "TGFb",       time = 0,   value = "init_TGFb",      method = "replace") |>
  addEvent(var = "bool_CHX",   time = -30, value = "var_bool_CHX",   method = "replace") |>
  addEvent(var = "bool_MG132", time = -30, value = "var_bool_MG132", method = "replace") |>
  addEvent(var = "bool_ActD",  time = -30, value = "var_bool_ActD",  method = "replace")

# one condition per perturbation; knockdowns rename a synthesis rate
cond.grid <- data.frame(Pertubation = c("Ctrl", "ActD", "CHX", "MG132", "R1Knd", "R2Knd"),
                        init_TGFb = 1, stringsAsFactors = FALSE)
cond.grid$var_bool_ActD  <- ifelse(cond.grid$Pertubation == "ActD",  1, 0)
cond.grid$var_bool_CHX   <- ifelse(cond.grid$Pertubation == "CHX",   1, 0)
cond.grid$var_bool_MG132 <- ifelse(cond.grid$Pertubation == "MG132", 1, 0)
cond.grid$k_pr_R1mRNA <- ifelse(cond.grid$Pertubation == "R1Knd", "k_pr_R1mRNA_R1Knd", "k_pr_R1mRNA")
cond.grid$k_pr_R2mRNA <- ifelse(cond.grid$Pertubation == "R2Knd", "k_pr_R2mRNA_R2Knd", "k_pr_R2mRNA")
rownames(cond.grid) <- cond.grid$Pertubation
cond.grid$Pertubation <- NULL

# per-condition trafo list: branch() broadcasts the base trafo over the grid
cond.trafo <- eqnvec() |>
  define("x~x", x = getParameters(reactions, events)) |>
  branch(table = cond.grid, apply = "insert")

out <- symmetryDetection(
  reactions, observables, method = "observability",
  events = events, trafo = cond.trafo,
  forcings = c("bool_ActD","bool_CHX","bool_MG132","TGFb"),
  equilibrate = TRUE, reduceCQ = TRUE, reconstruct = TRUE,
  cores = 6)
summary(out)   # rank 65 / 74: nine scalings, recovered exactly
