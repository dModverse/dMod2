# TNFα-induced NFκB signalling in HepG2 cells (Oppelt et al. 2018, npj Syst Biol Appl 4:23): the
# full model of the publication, before any reduction, with receptors that saturate and a ligand
# the cells secrete, on untreated and TNFα-stimulated cells.
# Data, model, conditions and objective; the other scripts of this project start with
# source("model.R").

## ---- setup
library(dMod2)

# generated sources and shared objects in build/, reused when nothing changed
dir.create("build", showWarnings = FALSE)
options(dMod2.outdir = "build", dMod2.cores = detectFreeCores())


## ---- data
# the data set of the publication, values on the natural log scale; the imaging readouts have
# no standard deviation, the error model supplies one
raw <- read.csv("data/oppelt2018.csv", skip = 15)
readouts <- c(
  "nuclear NFkB"       = "nucNFkB",
  "cytoplasmic NFkB"   = "cytNFkB",
  "A20"                = "A20GFP",
  "IkBa"               = "IkBa",
  "pIkBa"              = "pIkBa",
  "pNFkB"              = "pNFkB",
  "pIKK"               = "pIKK",
  "IP:IkBa | IB:IkBa"  = "ipIkBa_ibIkBa",
  "IP:NFkB | IB:NFkB"  = "ipNFkB_ibNFkB",
  "IP:NFkB | IB:pIkBa" = "ipNFkB_ibpIkBa",
  "IP:NFkB | IB:pNFkB" = "ipNFkB_ibpNFkB",
  "IP:IkBa | IB:pNFkB" = "ipIkBa_ibpNFkB",
  "NFkB/IkBa ratio"    = "ratioNFkBIkBa")
raw$name[grepl("^IP:IkBa \\| IB:NFkB", raw$name)] <- "ipIkBa_ibNFkB"
raw$name <- ifelse(raw$name %in% names(readouts), readouts[raw$name], raw$name)
raw$sigma[raw$name %in% c("nucNFkB", "cytNFkB", "A20GFP")] <- NA

# HepG2 cells without a compound, on the first blot scale
hepg2 <- subset(raw, cell == "HepG2" & scale == 1 & compound == "no")
hepg2$condition <- ifelse(hepg2$treatment == "TNF", "TNF", "control")
data <- as.datalist(hepg2[c("name", "time", "value", "sigma", "condition", "treatment")],
                    splitBy = "condition", keepCovariates = "treatment")


## ---- reactions
reactions <- eqnlist() |>
  assignCompartment(nIkb = "nucleus", pnNfk = "nucleus", nNfk = "nucleus", nNfkIkb = "nucleus",
                    volume = "Vnuc") |>
  # receptors activated by TNFα in the medium, which they do not use up, and by a ligand the
  # cells secrete, saturating in the sum of both
  addReaction("TNF", "TNF + TNFR", "uptake * (TNF + Lig) / (1 + K_R * (TNF + Lig))",
              "receptor activation") |>
  addReaction("TNFR", "", "deact_TNFR * TNFR", "receptor deactivation") |>
  addReaction("", "Lig", "prod_L", "ligand secretion") |>
  addReaction("Lig", "", "degrad_L * Lig", "ligand decay") |>
  # IKK cycle: activation by receptors, inhibited by A20, and basal; inactivation and recovery
  addReaction("Ikk", "pIkk", "act_Ikk_by_TNF / (1 + A20 / inh_Ikk) * TNFR * Ikk", "IKK activation") |>
  addReaction("Ikk", "pIkk", "act_Ikk_Base * Ikk", "basal IKK activation") |>
  addReaction("pIkk", "ppIkk", "act_pIkk * pIkk", "IKK phosphorylation") |>
  addReaction("ppIkk", "iIkk", "deact_ppIkk * ppIkk", "IKK inactivation") |>
  addReaction("iIkk", "Ikk", "trigger_iIkk * iIkk", "IKK recovery") |>
  # IκBα and NFκB phosphorylated in their complex, basal and by pIKK and ppIKK
  addReaction("NfkIkb", "NfkpIkb", "act_Ikb_Base * NfkIkb", "basal IkBa phosphorylation") |>
  addReaction("NfkIkb", "NfkpIkb", "act_Ikb_by_Ikk * pIkk * NfkIkb", "IkBa phosphorylation by pIKK") |>
  addReaction("NfkIkb", "NfkpIkb", "act_Ikb_by_Ikk * eta_s_ppIkk * ppIkk * NfkIkb",
              "IkBa phosphorylation by ppIKK") |>
  addReaction("NfkIkb", "pNfkIkb", "act_Nfk_Base * NfkIkb", "basal NFkB phosphorylation") |>
  addReaction("NfkIkb", "pNfkIkb", "act_Nfk_by_Ikk * pIkk * NfkIkb", "NFkB phosphorylation by pIKK") |>
  addReaction("NfkIkb", "pNfkIkb", "act_Nfk_by_Ikk * eta_s_ppIkk * ppIkk * NfkIkb",
              "NFkB phosphorylation by ppIKK") |>
  addReaction("pNfkIkb", "pNfkpIkb", "act_Ikb_Base * eta_act_Ikb_Base_c_pNfkIkb * pNfkIkb",
              "basal IkBa phosphorylation, pNFkB complex") |>
  addReaction("pNfkIkb", "pNfkpIkb", "act_Ikb_by_Ikk * eta_act_Ikb_c_pNfkIkb * pIkk * pNfkIkb",
              "IkBa phosphorylation by pIKK, pNFkB complex") |>
  addReaction("pNfkIkb", "pNfkpIkb",
              "act_Ikb_by_Ikk * eta_s_ppIkk * eta_act_Ikb_c_pNfkIkb * ppIkk * pNfkIkb",
              "IkBa phosphorylation by ppIKK, pNFkB complex") |>
  addReaction("NfkpIkb", "pNfkpIkb", "act_Nfk_Base * eta_act_Nfk_Base_c_NfkpIkb * NfkpIkb",
              "basal NFkB phosphorylation, pIkBa complex") |>
  addReaction("NfkpIkb", "pNfkpIkb", "act_Nfk_by_Ikk * eta_act_Nfk_c_NfkpIkb * pIkk * NfkpIkb",
              "NFkB phosphorylation by pIKK, pIkBa complex") |>
  addReaction("NfkpIkb", "pNfkpIkb",
              "act_Nfk_by_Ikk * eta_s_ppIkk * eta_act_Nfk_c_NfkpIkb * ppIkk * NfkpIkb",
              "NFkB phosphorylation by ppIKK, pIkBa complex") |>
  # release of NFκB from the complex with phosphorylated IκBα, and complex formation
  addReaction("NfkpIkb", "Nfk + pIkb", "split_NfkpIkb * NfkpIkb", "release of NFkB") |>
  addReaction("pNfkpIkb", "pNfk + pIkb", "split_NfkIkb * eta_split_pNfkpIkb * pNfkpIkb",
              "release of pNFkB") |>
  addReaction("Nfk + Ikb", "NfkIkb", "form_complex * Nfk * Ikb", "complex formation") |>
  # IκBα transcript, induced by nuclear NFκB, its translation, and degradation of pIκBα
  addReaction("", "mIkb", "prod_mIkb_by_nNfk * nNfk", "IkBa transcription") |>
  addReaction("mIkb", "", "degrad_mIkb * mIkb", "IkBa transcript degradation") |>
  addReaction("", "Ikb", "prod_Ikb * mIkb", "IkBa translation") |>
  addReaction("pIkb", "", "degrad_Ikb * pIkb", "pIkBa degradation") |>
  # nuclear import, dephosphorylation of pNFκB, the nuclear complex and its export
  addReaction("Ikb", "nIkb", "int_Nfk * eta_int_Ikb * Ikb", "IkBa import") |>
  addReaction("pNfk", "pnNfk", "int_Nfk * eta_int_pNfk * pNfk", "pNFkB import") |>
  addReaction("Nfk", "nNfk", "int_Nfk * Nfk", "NFkB import") |>
  addReaction("pnNfk", "nNfk", "deact_pnNfk * pnNfk", "nuclear pNFkB dephosphorylation") |>
  addReaction("nIkb + nNfk", "nNfkIkb", "form_complex * eta_form_complex_nuc * nNfk * nIkb",
              "nuclear complex formation") |>
  addReaction("nNfkIkb", "NfkIkb", "ext_nNfkIkb * nNfkIkb", "complex export") |>
  # A20 transcript with two delay stages, and A20 protein
  addReaction("", "RnaA20_1", "build_RnaA20 * nNfk", "A20 transcription") |>
  addReaction("RnaA20_1", "RnaA20_2", "shuttle_RnaA20 * RnaA20_1", "A20 transcript, delay 1") |>
  addReaction("RnaA20_2", "RnaA20", "shuttle_RnaA20 * RnaA20_2", "A20 transcript, delay 2") |>
  addReaction("RnaA20", "", "degrad_RnaA20 * RnaA20", "A20 transcript degradation") |>
  addReaction("", "A20", "build_A20 * RnaA20", "A20 translation") |>
  addReaction("A20", "", "degrad_A20 * A20", "A20 degradation")


## ---- observables
nfkbCyt <- "NfkIkb + NfkpIkb + pNfkIkb + pNfkpIkb + pNfk + Nfk"
ikbCyt <- "NfkIkb + NfkpIkb + pNfkIkb + pNfkpIkb + pIkb + Ikb"
pNfkb <- "pNfkpIkb + pNfkIkb + pNfk"
# log of scale times species plus offset; the imaging signals bleach exponentially, which
# subtracts a term linear in time
logOf <- function(stem, species, bleach = FALSE)
  sprintf(if (bleach) "log(s_%1$s * (%2$s) + off_%1$s) - tau_%1$s * time"
          else "log(s_%1$s * (%2$s) + off_%1$s)", stem, species)

observables <- eqnvec(
  nucNFkB        = logOf("nucNFkB", "pnNfk + nNfk + nNfkIkb", bleach = TRUE),
  cytNFkB        = logOf("cytNFkB", nfkbCyt, bleach = TRUE),
  A20GFP         = logOf("A20GFP", "A20", bleach = TRUE),
  IkBa           = logOf("IkBa", "NfkIkb + pNfkIkb + NfkpIkb + pIkb + Ikb"),
  pNFkB          = logOf("pNFkB", pNfkb),
  pIkBa          = logOf("pIkBa", "NfkpIkb + pNfkpIkb + pIkb"),
  pIKK           = logOf("pIKK", "pIkk + ppIkk"),
  ipIkBa_ibIkBa  = logOf("ipIkBa_ibIkBa", "NfkIkb + pNfkIkb + NfkpIkb + pIkb + Ikb"),
  ipNFkB_ibNFkB  = logOf("ipNFkB_ibNFkB", nfkbCyt),
  ipNFkB_ibpIkBa = logOf("ipNFkB_ibpIkBa", "NfkpIkb + pNfkpIkb"),
  ipNFkB_ibpNFkB = logOf("ipNFkB_ibpNFkB", pNfkb),
  ipIkBa_ibpNFkB = logOf("ipIkBa_ibpNFkB", "pNfkpIkb + pNfkIkb"),
  ipIkBa_ibNFkB  = logOf("ipIkBa_ibNFkB", "NfkIkb + NfkpIkb + pNfkIkb + pNfkpIkb"),
  ratioNFkBIkBa  = sprintf("log((%s) / (%s))", nfkbCyt, ikbCyt))


## ---- errors
# a constant standard deviation on the log scale for each imaging readout
errors <- eqnvec(nucNFkB = "sigma_nucNFkB", cytNFkB = "sigma_cytNFkB", A20GFP = "sigma_A20GFP")


## ---- steadystate
# the cells rest before TNFα is added: a steady state with the amounts of NFκB and IKK conserved;
# one that departs from rest by less than 1e-3 per minute, slower than the six hours of the
# experiment resolve, counts as stable, as just beyond a Hopf bifurcation
pSS <- Pimpl(reactions, modelname = "oppelt_ss", compile = FALSE, controlsPTC = list(stabilityTol = 1e-3))


## ---- odemodel
# TNFα is added at time 0
events <- addEvent(NULL, var = "TNF", time = 0, value = "dose_TNF", method = "replace")
model <- odemodel(reactions, events = events, modelname = "oppelt_ode", compile = FALSE)
x <- Xs(model, options = list(atol = 1e-10, rtol = 1e-8))
g <- Y(observables, x, modelname = "oppelt_obs", compile = FALSE)
e <- Y(errors, g, modelname = "oppelt_err", compile = FALSE)


## ---- transformation
# every parameter on the log scale but the saturation constant, which may be 0; at rest no TNFα,
# the dose is the unit of TNFα
trafo <- eqnvec() |>
  define("x <- x", x = setdiff(getParameters(pSS, x, g, e), reactions$states)) |>
  insert("x <- exp(x)", x = setdiff(.currentSymbols, "K_R")) |>
  define("TNF <- 0") |>
  branch(covariates(data)) |>
  define("dose_TNF <- 0", conditionMatch = "control") |>
  define("dose_TNF <- 1", conditionMatch = "TNF")
p <- P(trafo, modelname = "oppelt_par", compile = FALSE)


## ---- compile
compile(g, x, pSS, p, e, output = "oppelt", cores = getOption("dMod2.cores"))
prd <- g * x * pSS * p


## ---- objective
# a weak prior on the log scale: the offsets below the signals, the eta factors around 1
outer <- getParameters(prd)
mu <- structure(rep(-1, length(outer)), names = outer)
sigma <- structure(rep(4, length(outer)), names = outer)
mu[grepl("^off_", outer)] <- -2
sigma[grepl("^off_", outer)] <- 2
mu[grepl("^eta_", outer)] <- 0
obj <- normL2(data, prd, e) + constraintL2(mu, sigma = sigma)


## ---- tolerances
# the stopping tolerances of trust(), an order of magnitude above the error that the solver
# settings leave in the value (below 1e-6) and the gradient (below 1e-4)
tolerances <- list(ftol = 1e-5, mtol = 1e-5, gtol = 1e-3)


## ---- gauge
# one coordinate of each symmetry (symmetries.R) fixed: a scaling at 1, at 0 on the log scale,
# the saturation of the receptors at 0
gauge <- c(Vnuc = 0, totalNfk = 0, totalIkk = 0, uptake = 0, prod_Ikb = 0, build_RnaA20 = 0,
           split_NfkIkb = 0, s_A20GFP = 0, prod_L = 0, K_R = 0)
