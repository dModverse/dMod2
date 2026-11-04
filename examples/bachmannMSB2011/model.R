# JAK2-STAT5 signalling in CFU-E cells with the feedbacks CIS, SOCS3 and SHP1
# (Bachmann et al. 2011, Mol Syst Biol 7:516). Data, model, conditions and objective;
# the other scripts of this project start with source("model.R").

## ---- setup
library(dMod2)

# generated sources and shared objects in build/, reused when nothing changed;
# one evaluation spreads its conditions over the free cores
dir.create("build", showWarnings = FALSE)
options(dMod2.outdir = "build", dMod2.cores = detectFreeCores())


## ---- data
# all observables but the relative pSTAT5B are fitted on log10
data(bachmann)
logFitted <- bachmann$name != "pSTAT5B_rel"
bachmann$value[logFitted] <- log10(bachmann$value[logFitted])

covariates <- c("experiment", "epo_level", "ActD", "CISoe", "SOCS3oe", "SHP1oe")
data <- as.datalist(bachmann, splitBy = "condition", keepCovariates = covariates)
grid <- covariates(data)


## ---- reactions
# rates are written in the compartment of their educts
addCyt <- function(eq, from, to, rate) addReaction(eq, from, to, rate, compartment = "cyt")
addNuc <- function(eq, from, to, rate) addReaction(eq, from, to, rate, compartment = "nuc")

# five nuclear stages delay a transcript before it reaches the cytoplasm
addTranscript <- function(eq, gene) {
  st <- paste0(gene, "nRNA", 1:5)
  eq <- addNuc(eq, "", st[1], sprintf("%sRNAEqc * %sRNATurn * npSTAT5 * (1 - ActD)", gene, gene))
  for (i in 1:4) eq <- addNuc(eq, st[i], st[i + 1], sprintf("%sRNADelay * %s", gene, st[i]))
  eq <- addNuc(eq, st[5], paste0(gene, "RNA"), sprintf("%sRNADelay * %s", gene, st[5]))
  addCyt(eq, paste0(gene, "RNA"), "", sprintf("%sRNATurn * %sRNA", gene, gene))
}

reactions <- eqnlist() |>
  assignCompartment(npSTAT5 = "nuc", CISRNA = "cyt", SOCS3RNA = "cyt") |>
  # receptor: SOCS3 inhibits every phosphorylation, CIS the second site, SHP1 resets all forms
  addCyt("EpoRJAK2", "EpoRpJAK2", "JAK2ActEpo * EpoRJAK2 * epo_level / (1 + SOCS3Inh * SOCS3)") |>
  addCyt("EpoRpJAK2", "p1EpoRpJAK2", "EpoRActJAK2 * EpoRpJAK2 / (1 + SOCS3Inh * SOCS3)") |>
  addCyt("EpoRpJAK2", "p2EpoRpJAK2",
         "3 * EpoRActJAK2 * EpoRpJAK2 / (1 + SOCS3Inh * SOCS3) / (EpoRCISInh * EpoRJAK2_CIS + 1)") |>
  addCyt("p1EpoRpJAK2", "p12EpoRpJAK2",
         "3 * EpoRActJAK2 * p1EpoRpJAK2 / (1 + SOCS3Inh * SOCS3) / (EpoRCISInh * EpoRJAK2_CIS + 1)") |>
  addCyt("p2EpoRpJAK2", "p12EpoRpJAK2", "EpoRActJAK2 * p2EpoRpJAK2 / (1 + SOCS3Inh * SOCS3)") |>
  addCyt("EpoRpJAK2", "EpoRJAK2", "JAK2EpoRDeaSHP1 * EpoRpJAK2 * SHP1Act") |>
  addCyt("p1EpoRpJAK2", "EpoRJAK2", "JAK2EpoRDeaSHP1 * p1EpoRpJAK2 * SHP1Act") |>
  addCyt("p2EpoRpJAK2", "EpoRJAK2", "JAK2EpoRDeaSHP1 * p2EpoRpJAK2 * SHP1Act") |>
  addCyt("p12EpoRpJAK2", "EpoRJAK2", "JAK2EpoRDeaSHP1 * p12EpoRpJAK2 * SHP1Act") |>
  # CIS occupies the receptor; the active receptor clears it
  addCyt("EpoRJAK2_CIS", "", "EpoRCISRemove * EpoRJAK2_CIS * (p1EpoRpJAK2 + p12EpoRpJAK2)") |>
  addCyt("SHP1", "SHP1Act",
         "SHP1ActEpoR * SHP1 * (EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2 + p12EpoRpJAK2)") |>
  addCyt("SHP1Act", "SHP1", "SHP1Dea * SHP1Act") |>
  # STAT5: phosphorylated by JAK2 on any active receptor and by the receptor after site 1
  addCyt("STAT5", "pSTAT5",
         "STAT5ActJAK2 * STAT5 * (EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2 + p12EpoRpJAK2) / (1 + SOCS3Inh * SOCS3)") |>
  addCyt("STAT5", "pSTAT5",
         "STAT5ActEpoR * STAT5 * (p1EpoRpJAK2 + p12EpoRpJAK2)^2 / (1 + SOCS3Inh * SOCS3) / (1 + CISInh * CIS)") |>
  addCyt("pSTAT5", "npSTAT5", "STAT5Imp * pSTAT5") |>
  addNuc("npSTAT5", "STAT5", "STAT5Exp * npSTAT5") |>
  addTranscript("CIS") |>
  addTranscript("SOCS3") |>
  # feedback proteins, each with the overexpression term of its experiment
  addCyt("", "CIS", "CISEqc * CISTurn * CISRNA") |>
  addCyt("CIS", "", "CISTurn * CIS") |>
  addCyt("", "CIS", "CISoe * CISEqcOE * CISTurn") |>
  addCyt("", "SOCS3", "SOCS3Eqc * SOCS3Turn * SOCS3RNA") |>
  addCyt("SOCS3", "", "SOCS3Turn * SOCS3") |>
  addCyt("", "SOCS3", "SOCS3oe * SOCS3EqcOE * SOCS3Turn") |>
  setCompartmentVolume(cyt = "0.4", nuc = "0.275")


## ---- observables
# log10 with a floor that keeps empty species finite
log10Of <- function(x) paste0("log10(", x, " + 1e-15)")
pJAK2 <- "2 * (EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2 + p12EpoRpJAK2)"
pEpoR <- "16 * (p1EpoRpJAK2 + p2EpoRpJAK2 + p12EpoRpJAK2)"

# scale_ and offset_ are resolved per experiment in the condition table
observables <- do.call(eqnvec, c(list(
  pJAK2_au    = log10Of(paste("offset_pJAK2 + scale_pJAK2 / init_EpoRJAK2 *", pJAK2)),
  pEpoR_au    = log10Of(paste("offset_pEpoR + scale_pEpoR / init_EpoRJAK2 *", pEpoR)),
  CIS_au      = log10Of("offset_CIS + scale_CIS / CISEqc / CISRNAEqc / init_STAT5 * CIS"),
  SOCS3_au    = log10Of("offset_SOCS3 + scale_SOCS3 / SOCS3Eqc / SOCS3RNAEqc / init_STAT5 * SOCS3"),
  tSTAT5_au   = log10Of("scale_tSTAT5 / init_STAT5 * (STAT5 + pSTAT5)"),
  pSTAT5_au   = log10Of("offset_pSTAT5 + scale_pSTAT5 / init_STAT5 * pSTAT5"),
  tSHP1_au    = log10Of("scale_tSHP1 / init_SHP1 * (SHP1 + SHP1Act)"),
  CIS_au1     = log10Of("scale1_CIS_dr90 / CISEqc / CISRNAEqc / init_STAT5 * CIS"),
  CIS_au2     = log10Of("scale2_CIS_dr90 / CISEqc / CISRNAEqc / init_STAT5 * CIS"),
  STAT5_abs   = log10Of("STAT5"),
  SHP1_abs    = log10Of("SHP1 + SHP1Act"),
  CIS_abs     = log10Of("CIS"),
  SOCS3_abs   = log10Of("SOCS3"),
  pSTAT5B_rel = "offset_pSTAT5_conc + 100 * pSTAT5 / (pSTAT5 + STAT5)"),
  # three qPCR replicates, each with its own scale
  setNames(as.list(log10Of(paste0("1 + scale_CISRNA_fold", LETTERS[1:3], " / CISRNAEqc / init_STAT5 * CISRNA"))),
           paste0("CISRNA_fold", LETTERS[1:3])),
  setNames(as.list(log10Of(paste0("1 + scale_SOCS3RNA_fold", LETTERS[1:3], " / SOCS3RNAEqc / init_STAT5 * SOCS3RNA"))),
           paste0("SOCS3RNA_fold", LETTERS[1:3]))))


## ---- errors
# one standard deviation per antibody, shared by all experiments
errors <- eqnvec(
  pJAK2_au       = "sd_JAK2EpoR_au",
  pEpoR_au       = "sd_JAK2EpoR_au",
  CIS_au         = "sd_CIS_au",
  CIS_au1        = "sd_CIS_au",
  CIS_au2        = "sd_CIS_au",
  SOCS3_au       = "sd_SOCS3_au",
  tSTAT5_au      = "sd_STAT5_au",
  tSHP1_au       = "sd_SHP1_au",
  pSTAT5_au      = "sd_pSTAT5_au",
  STAT5_abs      = "sd_STAT5_abs",
  SHP1_abs       = "sd_SHP1_abs",
  CIS_abs        = "sd_CIS_abs",
  SOCS3_abs      = "sd_SOCS3_abs",
  pSTAT5B_rel    = "sd_pSTAT5_rel",
  CISRNA_foldA   = "sd_RNA_fold",
  CISRNA_foldB   = "sd_RNA_fold",
  CISRNA_foldC   = "sd_RNA_fold",
  SOCS3RNA_foldA = "sd_RNA_fold",
  SOCS3RNA_foldB = "sd_RNA_fold",
  SOCS3RNA_foldC = "sd_RNA_fold"
)


## ---- conditions
# experiments that recorded each blot; elsewhere its scale is 1 and its offset 0
recordedBy <- list(
  pJAK2  = c("long", "actd", "fine", "cisoe", "socs3oe", "shp1oe", "dr7", "dr30"),
  pEpoR  = c("long", "actd", "fine", "cisoe", "cisoe_pepor", "socs3oe", "shp1oe", "dr7", "dr30"),
  CIS    = c("long", "actd", "cisoe", "socs3oe", "shp1oe"),
  SOCS3  = c("long", "cisoe", "socs3oe"),
  tSTAT5 = c("long", "actd", "shp1oe"),
  pSTAT5 = c("long", "actd", "cisoe", "socs3oe", "shp1oe", "dr10"),
  tSHP1  = "shp1oe")
stemOf <- c(
  pJAK2  = "pJAK2",
  pEpoR  = "pEpoR",
  CIS    = "CIS",
  SOCS3  = "SOCS3",
  tSTAT5 = "tSTAT5",
  pSTAT5 = "pSTAT5",
  tSHP1  = "SHP1"
)
# totals have no background; the pSTAT5 dose response was blotted without one
withOffset <- setdiff(names(stemOf), c("tSTAT5", "tSHP1"))

for (b in names(stemOf)) {
  seen <- grid$experiment %in% recordedBy[[b]]
  grid[[paste0("scale_", b)]] <- ifelse(seen, paste0("scale_", stemOf[[b]], "_", grid$experiment), "1")
  if (b %in% withOffset)
    grid[[paste0("offset_", b)]] <- ifelse(seen & !(b == "pSTAT5" & grid$experiment == "dr10"),
                                           paste0("offset_", stemOf[[b]], "_", grid$experiment), "0")
}
# SOCS3 overexpression widens the pSTAT5 blot
grid$sd_pSTAT5_au <- ifelse(grid$SOCS3oe == 1, "sd_STAT5_au + sd_pSTAT5_socs3oe", "sd_STAT5_au")
grid$perturbation <- with(grid,
  ifelse(ActD == 1, "actinomycin D",
  ifelse(CISoe == 1, "CIS overexpressed",
  ifelse(SOCS3oe == 1, "SOCS3 overexpressed",
  ifelse(SHP1oe == 1, "SHP1 overexpressed", "control")))))
data <- as.datalist(data, condition.grid = grid)

# the columns of the grid that replace symbols; experiment and perturbation are labels
substitutions <- c("epo_level", "ActD", "CISoe", "SOCS3oe", "SHP1oe", "sd_pSTAT5_au",
                   grep("^(scale|offset)_", names(grid), value = TRUE))


## ---- odemodel
solverOptions <- list(atol = 1e-11, rtol = 1e-8, maxsteps = 1e7L, maxattempts = 100L)
model <- odemodel(reactions, modelname = "bachmann_ode", compile = FALSE)
x <- Xs(model, options = solverOptions)
g <- Y(observables, x, modelname = "bachmann_obs", compile = FALSE)
e <- Y(errors, g, modelname = "bachmann_err", compile = FALSE)


## ---- transformation
# only the receptor, SHP1 and STAT5 pools are stocked; overexpressed CIS and SOCS3 start at
# the level their production term holds
zero <- setdiff(reactions$states, c("EpoRJAK2", "SHP1", "STAT5", "EpoRJAK2_CIS", "CIS", "SOCS3"))

trafo <- eqnvec() |>
  define("x <- x", x = getParameters(model, g, e)) |>
  define("x <- 0", x = zero) |>
  insert("EpoRJAK2 <- init_EpoRJAK2") |>
  insert("STAT5 <- init_STAT5") |>
  insert("SHP1 <- init_SHP1 * (1 + SHP1oe * SHP1ProOE)") |>
  # parameters relative to their reference pool; each rule reads the raw parameter
  insert("CISRNAEqc <- CISRNAEqc / init_STAT5") |>
  insert("SOCS3RNAEqc <- SOCS3RNAEqc / init_STAT5") |>
  insert("CISEqc <- CISEqc / CISRNAEqc") |>
  insert("SOCS3Eqc <- SOCS3Eqc / SOCS3RNAEqc") |>
  insert("CISEqcOE <- CISEqcOE * CISEqc") |>
  insert("SOCS3EqcOE <- SOCS3EqcOE * SOCS3Eqc") |>
  insert("CISInh <- CISInh / CISEqc") |>
  insert("SOCS3Inh <- SOCS3Inh / SOCS3Eqc") |>
  insert("JAK2EpoRDeaSHP1 <- JAK2EpoRDeaSHP1 / init_SHP1") |>
  insert("EpoRCISRemove <- EpoRCISRemove / init_EpoRJAK2") |>
  insert("SHP1ActEpoR <- SHP1ActEpoR / init_EpoRJAK2") |>
  insert("STAT5ActJAK2 <- STAT5ActJAK2 / init_EpoRJAK2") |>
  insert("STAT5ActEpoR <- STAT5ActEpoR / init_EpoRJAK2^2") |>
  insert("CIS <- CISoe * CISEqcOE * CISEqc") |>
  insert("EpoRJAK2_CIS <- CISoe") |>
  insert("SOCS3 <- SOCS3oe * SOCS3EqcOE * SOCS3Eqc") |>
  branch(table = grid[substitutions], apply = "insert") |>
  # the RNA equilibrium constants only appear next to a scale
  insert("x <- 1", x = c("CISRNAEqc", "SOCS3RNAEqc")) |>
  insert("x <- exp10(x)", x = .currentSymbols)
p <- P(trafo, modelname = "bachmann_par", compile = FALSE)


## ---- compile
compile(g, x, p, e, output = "bachmann", cores = getOption("dMod2.cores"))
prd <- g * x * p


## ---- objective
# data term and a prior on the receptor pool, on the log10 scale of the parameters
obj <- normL2(data, prd, e) + constraintL2(c(init_EpoRJAK2 = 0.59947), sigma = 0.19346)

# bounds of the fit in log10: from -3 to 3, wider for seven rates
outer <- getParameters(prd)
lower <- setNames(rep(-3, length(outer)), outer)
upper <- setNames(rep(3, length(outer)), outer)
upper[c("CISEqc", "EpoRActJAK2", "JAK2EpoRDeaSHP1", "SOCS3Turn")] <- 4
upper["EpoRCISInh"] <- 6
upper["JAK2ActEpo"] <- 9
upper["CISInh"] <- 12
