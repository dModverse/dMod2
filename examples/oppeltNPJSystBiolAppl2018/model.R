# TNFα-induced NFκB signaling in HepG2 cells (Oppelt et al. 2018, npj Syst Biol Appl 4:23): the
# full model of the publication, before any reduction, on untreated and TNFα-stimulated cells.
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
  # receptors activated by TNFα in the medium, which they do not use up
  addReaction("TNF", "TNF + TNFR", "uptake * TNF", "receptor activation") |>
  addReaction("TNFR", "", "deact_TNFR * TNFR", "receptor deactivation") |>
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
  # IκBα transcript, induced by nuclear NFκB, the translation, and degradation of pIκBα
  addReaction("", "mIkb", "prod_mIkb_by_nNfk * nNfk", "IkBa transcription") |>
  addReaction("mIkb", "", "degrad_mIkb * mIkb", "IkBa transcript degradation") |>
  addReaction("", "Ikb", "prod_Ikb * mIkb", "IkBa translation") |>
  addReaction("pIkb", "", "degrad_Ikb * pIkb", "pIkBa degradation") |>
  # nuclear import, dephosphorylation of pNFκB, the nuclear complex and the complex export
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
# an unstable one is returned rather than refused, the objective below keeps it stable
pSS <- Pimpl(reactions, modelname = "oppelt_ss", compile = FALSE, deriv2 = TRUE,
             controlsPTC = list(stability = FALSE))


## ---- odemodel
# TNFα is added at time 0
events <- addEvent(NULL, var = "TNF", time = 0, value = "dose_TNF", method = "replace")
model <- odemodel(reactions, events = events, modelname = "oppelt_ode", compile = FALSE)
x <- Xs(model, options = list(atol = 1e-10, rtol = 1e-8))
g <- Y(observables, x, modelname = "oppelt_obs", compile = FALSE)
e <- Y(errors, g, modelname = "oppelt_err", compile = FALSE)


## ---- transformation
# every parameter on the log scale; at rest no TNFα, the dose is the unit of TNFα
trafo <- eqnvec() |>
  define("x <- x", x = setdiff(getParameters(pSS, x, g, e), reactions$states)) |>
  insert("x <- exp(x)", x = .currentSymbols) |>
  define("TNF <- 0") |>
  branch(covariates(data)) |>
  define("dose_TNF <- 0", conditionMatch = "control") |>
  define("dose_TNF <- 1", conditionMatch = "TNF")
p <- P(trafo, modelname = "oppelt_par", compile = FALSE)


## ---- compile
compile(g, x, pSS, p, e, output = "oppelt", cores = 6)
prd <- g * x * pSS * p


## ---- objective
# a weak prior on the log scale: mean 0, standard deviation 4, a factor of about 55, narrower
# where the order of magnitude is known before the fit
outer <- getParameters(prd)
mu <- structure(rep(0, length(outer)), names = outer)
sigma <- structure(rep(4, length(outer)), names = outer)
# bleaching removes at most a few percent of the signal over the time of the data
bleach <- grep("^tau_", outer, value = TRUE)
mu[bleach] <- log(0.1 / max(as.data.frame(data)$time))
sigma[bleach] <- 2
# the imaging readouts scatter by about 5%
noise <- grep("^sigma_", outer, value = TRUE)
mu[noise] <- log(0.05)
sigma[noise] <- 1
# without TNFα, IKK, IκBα and NFκB are activated far more slowly than after it
basal <- grep("_Base$", outer, value = TRUE)
mu[basal] <- -8
sigma[basal] <- 2
# transcripts and A20 are made and degraded within about an hour
turnover <- c("shuttle_RnaA20", "degrad_RnaA20", "build_A20", "degrad_A20", "degrad_mIkb")
mu[turnover] <- log(1 / 60)
sigma[turnover] <- 2
# IKK recovers from its inactive form over hours
mu["trigger_iIkk"] <- log(1 / 240)
sigma["trigger_iIkk"] <- 2
# the cells rest: a term that is zero on a stable resting state and reaches 100 where it loses
# stability, on a ramp of a tenth of the inverse time span of the data, wide enough for the model
# of a fit step to see the boundary before the step crosses it
rampWidth <- 0.1 / max(as.data.frame(data)$time)
stability <- constraintStability(pSS, p, width = rampWidth)
obj <- normL2(data, prd, e) + constraintL2(mu, sigma = sigma) + stability


## ---- tolerances
# the stopping tolerances of trust(), above the error that the solver settings leave: ftol and
# mtol an order of magnitude above that of the value (about 2e-6), gtol two above that of the
# gradient (about 8e-5)
tolerances <- list(ftol = 1e-4, mtol = 1e-4, gtol = 1e-2)


## ---- gauge
# one coordinate of each scaling (symmetries.R) fixed at 1, at 0 on the log scale
gauge <- c(Vnuc = 0, totalNfk = 0, totalIkk = 0, uptake = 0, prod_Ikb = 0, build_RnaA20 = 0,
           split_NfkIkb = 0, s_A20GFP = 0)

