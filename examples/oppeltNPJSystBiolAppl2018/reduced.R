# The core model of the publication, reduced by the limits the profiles of profiles.R point to. Data, observables and
# errors are those of model.R; the scripts of the reduced model start with source("reduced.R").

# model.R, unless it ran in this session
if (!exists("obj")) source("model.R")

## ---- reducedreactions
# the core model of the publication (Table S5 of the supplement), with the two delay steps of the
# A20 transcript of the full model and a separate factor for the release of pNFκB
reactionsReduced <- eqnlist() |>
  assignCompartment(nIkb = "nucleus", pnNfk = "nucleus", nNfk = "nucleus", nNfkIkb = "nucleus",
                    volume = "Vnuc") |>
  addReaction("TNF", "TNF + TNFR", "uptake * TNF", "receptor activation") |>
  addReaction("TNFR", "", "deact_TNFR * TNFR", "receptor deactivation") |>
  # IKK activated by the receptors alone, without basal activation and inhibition by A20
  addReaction("Ikk", "pIkk", "act_Ikk_by_TNF * TNFR * Ikk", "IKK activation") |>
  addReaction("pIkk", "ppIkk", "act_pIkk * pIkk", "IKK phosphorylation") |>
  addReaction("ppIkk", "iIkk", "deact_ppIkk * ppIkk", "IKK inactivation") |>
  addReaction("iIkk", "Ikk", "trigger_iIkk * iIkk", "IKK recovery") |>
  # IκBα phosphorylated by pIKK in both complexes and basally only in the complex with pNFκB;
  # NFκB phosphorylated by pIKK only; ppIKK phosphorylates neither
  addReaction("NfkIkb", "NfkpIkb", "act_Ikb_by_Ikk * pIkk * NfkIkb", "IkBa phosphorylation by pIKK") |>
  addReaction("pNfkIkb", "pNfkpIkb", "act_Ikb_complex * pNfkIkb",
              "basal IkBa phosphorylation, pNFkB complex") |>
  addReaction("pNfkIkb", "pNfkpIkb", "act_Ikb_by_Ikk * pIkk * pNfkIkb",
              "IkBa phosphorylation by pIKK, pNFkB complex") |>
  addReaction("NfkIkb", "pNfkIkb", "act_Nfk_by_Ikk * pIkk * NfkIkb", "NFkB phosphorylation by pIKK") |>
  addReaction("NfkpIkb", "pNfkpIkb", "act_Nfk_by_Ikk * eta_act_Nfk_c_NfkpIkb * pIkk * NfkpIkb",
              "NFkB phosphorylation by pIKK, pIkBa complex") |>
  addReaction("NfkpIkb", "Nfk + pIkb", "split_NfkpIkb * NfkpIkb", "release of NFkB") |>
  addReaction("pNfkpIkb", "pNfk + pIkb", "split_NfkIkb * eta_split_pNfkpIkb * pNfkpIkb",
              "release of pNFkB") |>
  addReaction("Nfk + Ikb", "NfkIkb", "form_complex * Nfk * Ikb", "complex formation") |>
  addReaction("", "mIkb", "prod_mIkb_by_nNfk * nNfk", "IkBa transcription") |>
  addReaction("mIkb", "", "degrad_mIkb * mIkb", "IkBa transcript degradation") |>
  addReaction("", "Ikb", "prod_Ikb * mIkb", "IkBa translation") |>
  addReaction("pIkb", "", "degrad_Ikb * pIkb", "pIkBa degradation") |>
  addReaction("Ikb", "nIkb", "int_Nfk * eta_int_Ikb * Ikb", "IkBa import") |>
  addReaction("pNfk", "pnNfk", "int_Nfk * eta_int_pNfk * pNfk", "pNFkB import") |>
  addReaction("Nfk", "nNfk", "int_Nfk * Nfk", "NFkB import") |>
  addReaction("pnNfk", "nNfk", "deact_pnNfk * pnNfk", "nuclear pNFkB dephosphorylation") |>
  addReaction("nIkb + nNfk", "nNfkIkb", "form_complex_nuc * nNfk * nIkb",
              "nuclear complex formation") |>
  addReaction("nNfkIkb", "NfkIkb", "ext_nNfkIkb * nNfkIkb", "complex export") |>
  addReaction("", "RnaA20_1", "build_RnaA20 * nNfk", "A20 transcription") |>
  addReaction("RnaA20_1", "RnaA20_2", "shuttle_RnaA20 * RnaA20_1", "A20 transcript, delay 1") |>
  addReaction("RnaA20_2", "RnaA20", "shuttle_RnaA20 * RnaA20_2", "A20 transcript, delay 2") |>
  addReaction("RnaA20", "", "degrad_RnaA20 * RnaA20", "A20 transcript degradation") |>
  addReaction("", "A20", "build_A20 * RnaA20", "A20 translation") |>
  addReaction("A20", "", "degrad_A20 * A20", "A20 degradation")


## ---- reducedobservables
# the readout whose scale goes to 0 is the offset
observablesReduced <- observables
observablesReduced["ipNFkB_ibNFkB"] <- "log(off_ipNFkB_ibNFkB)"


## ---- reducedmodel
# at rest all NFκB is bound to IκBα and all IKK inactive, the other states are 0: without basal
# activation this is a steady state for every parameter value, as in the publication; TNFα added
# at time 0 as in model.R
modelReduced <- odemodel(reactionsReduced, events = events, modelname = "oppelt_red_ode", compile = FALSE)
xReduced <- Xs(modelReduced, options = list(atol = 1e-10, rtol = 1e-8))
gReduced <- Y(observablesReduced, xReduced, modelname = "oppelt_red_obs", compile = FALSE)
eReduced <- Y(errors, gReduced, modelname = "oppelt_red_err", compile = FALSE)
resting <- c(NfkIkb = "totalNfk", Ikk = "totalIkk")
resting <- c(resting, structure(rep("0", length(setdiff(reactionsReduced$states, names(resting)))),
                                names = setdiff(reactionsReduced$states, names(resting))))

# the fixed values of the publication (Tables S8 and S9): fast processes at 1000 per minute, slow
# receptor deactivation and NFκB import, offsets and bleaching that vanish
fixedReduced <- c(deact_TNFR = 1e-3, int_Nfk = 1e-2, deact_pnNfk = 1000, form_complex_nuc = 1000,
                  ext_nNfkIkb = 1000, off_IkBa = 0, off_ipIkBa_ibIkBa = 0, tau_A20GFP = 0)
trafoReduced <- eqnvec() |>
  define("x <- x", x = setdiff(getParameters(xReduced, gReduced, eReduced),
                               reactionsReduced$states)) |>
  define("x <- y", x = names(fixedReduced), y = fixedReduced) |>
  define("x <- y", x = names(resting), y = resting) |>
  insert("x <- exp(x)", x = .currentSymbols) |>
  define("TNF <- 0") |>
  branch(covariates(data)) |>
  define("dose_TNF <- 0", conditionMatch = "control") |>
  define("dose_TNF <- 1", conditionMatch = "TNF")
pReduced <- P(trafoReduced, modelname = "oppelt_red_par", compile = FALSE)


## ---- reducedcompile
compile(gReduced, xReduced, pReduced, eReduced, output = "oppelt_red",
        cores = getOption("dMod2.cores"))
prdReduced <- gReduced * xReduced * pReduced


## ---- reducedobjective
# the prior of model.R on the parameters that remain
outerReduced <- getParameters(prdReduced)
# the basal phosphorylation in the pNFκB complex at the product of the two rates it replaces
muReduced <- c(mu, act_Ikb_complex = mu[["act_Ikb_Base"]] + mu[["eta_act_Ikb_Base_c_pNfkIkb"]])[outerReduced]
sigmaReduced <- c(sigma, act_Ikb_complex = sqrt(sigma[["act_Ikb_Base"]]^2 +
                                                  sigma[["eta_act_Ikb_Base_c_pNfkIkb"]]^2))[outerReduced]
gaugeReduced <- gauge[intersect(names(gauge), outerReduced)]
objReduced <- normL2(data, prdReduced, eReduced) + constraintL2(muReduced, sigma = sigmaReduced)
