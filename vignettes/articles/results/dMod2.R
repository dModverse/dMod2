# The model, multistart and profiles of the Get started page, dMod2.Rmd, which shows the labeled
# parts and reads the fits and profiles from dMod2.rds. Run in this folder.

## ---- setup
library(dMod2)
options(dMod2.outdir = tempdir())

## ---- data
data(badata)
data <- as.datalist(badata)

## ---- reactions
reactions <- eqnlist() |>
  addReaction("TCA_buffer", "TCA_cell", "k_import*TCA_buffer", "uptake") |>
  addReaction("TCA_cell", "TCA_buffer", "k_export_sinus*TCA_cell", "sinusoidal export") |>
  addReaction("TCA_cell", "TCA_cana", "k_export_cana*TCA_cell", "canalicular export") |>
  addReaction("TCA_cana", "TCA_buffer", "k_reflux*TCA_cana", "reflux into the buffer")
reactions

## ---- model
x <- Xs(odemodel(reactions, modelname = "start_x", compile = FALSE))
observables <- eqnvec(buffer = "s*TCA_buffer", cellular = "s*(TCA_cana + TCA_cell)")
g <- Y(observables, x, modelname = "start_g", compile = FALSE)

## ---- trafo
trafo <- eqnvec() |>
  define("x <- x", x = getParameters(x, g)) |>
  define("TCA_buffer <- 0") |>
  insert("x <- 10^y", x = .currentSymbols, y = toupper(.currentSymbols)) |>
  branch(conditions = c("closed", "open")) |>
  insert("K_REFLUX <- K_REFLUX_OPEN", conditionMatch = "open")
trafo
p <- P(trafo, modelname = "start_p", compile = FALSE)
compile(x, g, p, output = "start", cores = 6)

## ---- objective
outerpars <- getParameters(p)
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)
obj <- normL2(data, g * x * p, attrName = "data") +
  constraintL2(pouter, sigma = 4, attrName = "prior")

## ---- fit
set.seed(5555)
fits <- mstrust(obj, pouter, sd = 4, fits = 50, cores = 6)
fits <- as.parframe(fits)

## ---- best
best <- as.parvec(fits)

## ---- profile
profiles <- profile(obj, best, whichPar = names(best), method = "integrate", cores = 6,
                    limits = c(lower = -5, upper = 5), stepControl = list(stop = "data"))

## ---- save
result <- list(fits = fits, profiles = profiles)
attr(result, "stamp") <- list(md5 = tools::md5sum("dMod2.R"), dMod2 = as.character(packageVersion("dMod2")))
saveRDS(result, "dMod2.rds")
