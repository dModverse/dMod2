# Fits of the core model of the publication and of variants that release one of its settings,
# each against the full model by the likelihood ratio. Saves results/variants.rds, which the guide
# reads.

# reduced.R, unless it ran in this session
if (!exists("objReduced")) source("reduced.R")

## ---- variants
# the core model: one delay step of the A20 transcript and the factor of the release of NFκB also
# for pNFκB; the variants add the second delay step, the separate factor, or estimate one of the
# fixed values
reactionsCore <- subset(reactionsReduced, !Description %in% c("A20 transcript, delay 1",
                                                              "A20 transcript, delay 2",
                                                              "release of pNFkB")) |>
  addReaction("RnaA20_1", "RnaA20", "shuttle_RnaA20 * RnaA20_1", "A20 transcript, delay") |>
  addReaction("pNfkpIkb", "pNfk + pIkb", "split_NfkIkb * pNfkpIkb", "release of pNFkB")
reactionsDelays <- subset(reactionsReduced, Description != "release of pNFkB") |>
  addReaction("pNfkpIkb", "pNfk + pIkb", "split_NfkIkb * pNfkpIkb", "release of pNFkB")
reactionsSplit <- subset(reactionsReduced, !Description %in% c("A20 transcript, delay 1",
                                                               "A20 transcript, delay 2")) |>
  addReaction("RnaA20_1", "RnaA20", "shuttle_RnaA20 * RnaA20_1", "A20 transcript, delay")

# trajectory, observables and errors of a set of reactions
odeOf <- function(name, reactions) {
  model <- odemodel(reactions, events = events, modelname = paste0("oppelt_var_", name, "_ode"),
                    compile = FALSE)
  x <- Xs(model, options = list(atol = 1e-10, rtol = 1e-8))
  g <- Y(observablesReduced, x, modelname = paste0("oppelt_var_", name, "_obs"), compile = FALSE)
  e <- Y(errors, g, modelname = paste0("oppelt_var_", name, "_err"), compile = FALSE)
  list(reactions = reactions, x = x, g = g, e = e)
}
# the parameter transformation of reduced.R with the fixed values `fixed`
withFixed <- function(ode, name, fixed) {
  rest <- c(NfkIkb = "totalNfk", Ikk = "totalIkk")
  zero <- setdiff(ode$reactions$states, names(rest))
  rest <- c(rest, structure(rep("0", length(zero)), names = zero))
  trafo <- eqnvec() |>
    define("x <- x", x = setdiff(getParameters(ode$x, ode$g, ode$e), ode$reactions$states)) |>
    define("x <- y", x = names(fixed), y = fixed) |>
    define("x <- y", x = names(rest), y = rest) |>
    insert("x <- exp(x)", x = .currentSymbols) |>
    define("TNF <- 0") |>
    branch(covariates(data)) |>
    define("dose_TNF <- 0", conditionMatch = "control") |>
    define("dose_TNF <- 1", conditionMatch = "TNF")
  c(ode, list(p = P(trafo, modelname = paste0("oppelt_var_", name, "_par"), compile = FALSE)))
}
odes <- list(core = odeOf("core", reactionsCore), delays = odeOf("delays", reactionsDelays),
             split = odeOf("split", reactionsSplit))
variants <- c(list(core = withFixed(odes$core, "core", fixedReduced),
                   delays = withFixed(odes$delays, "delays", fixedReduced),
                   split = withFixed(odes$split, "split", fixedReduced)),
              sapply(names(fixedReduced), function(r)
                withFixed(odes$core, r, fixedReduced[names(fixedReduced) != r]), simplify = FALSE))
do.call(compile, c(unlist(lapply(odes, `[`, c("x", "g", "e")), recursive = FALSE),
                   lapply(variants, `[[`, "p"),
                   list(output = "oppelt_var", cores = getOption("dMod2.cores"))))

## ---- variantfits
# each variant from the refitted full model, as refit.R starts the reduced model, then 20 fits
# around that start; the prior of model.R
# a rate that replaces a product of two rates of the full model, on the log scale their sum
full <- readRDS("results/profiles.rds")$refit
products <- list(act_Ikb_complex = c("act_Ikb_Base", "eta_act_Ikb_Base_c_pNfkIkb"),
                 form_complex_nuc = c("form_complex", "eta_form_complex_nuc"))
startFull <- c(full$argument, gauge,
               vapply(products, function(f) sum(full$argument[f]), 0))
muAll <- c(mu, vapply(products, function(f) sum(mu[f]), 0))
sigmaAll <- c(sigma, vapply(products, function(f) sqrt(sum(sigma[f]^2)), 0))
dataFull <- normL2(data, prd, e)(full$argument, fixed = gauge)$value
fitVariant <- function(v) {
  prdV <- v$g * v$x * v$p
  outer <- getParameters(prdV)
  gaugeV <- gauge[intersect(names(gauge), outer)]
  free <- setdiff(outer, names(gaugeV))
  objV <- normL2(data, prdV, v$e) + constraintL2(muAll[outer], sigma = sigmaAll[outer])
  start <- ifelse(free %in% names(startFull), startFull[free], muAll[free])
  names(start) <- free
  set.seed(1)
  fits <- mstrust(objV, start, fixed = gaugeV, sd = 0.5, fits = 20, startFromCenter = TRUE,
                  iterlim = 1000, tolControl = tolerances, cores = getOption("dMod2.cores")) |>
    as.parframe()
  best <- as.parvec(fits)
  dataV <- normL2(data, prdV, v$e)(best, fixed = gaugeV)$value
  removed <- length(full$argument) - length(free)
  c(data = dataV, difference = dataV - dataFull, removed = removed,
    threshold = qchisq(0.95, removed))
}
variantTable <- do.call(rbind, lapply(variants, fitVariant))

## ---- variantsave
sources <- c("model.R", "reduced.R", "variants.R", "data/oppelt2018.csv", "results/profiles.rds")
attr(variantTable, "stamp") <- list(md5 = tools::md5sum(sources), dMod2 = packageVersion("dMod2"))
saveRDS(variantTable, "results/variants.rds")
