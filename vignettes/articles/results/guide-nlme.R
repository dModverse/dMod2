# The model and fits of guide chapter 7, guide-nlme.Rmd, which shows the labeled parts and reads
# the fits from guide-nlme.rds. Run in this folder.

## ---- setup
library(dMod2)
options(dMod2.outdir = file.path(tempdir(), "guide-mixed-effects"))
dir.create(getOption("dMod2.outdir"), showWarnings = FALSE)

## ---- reactions
reactions <- eqnlist() |>
  assignCompartment(Ag = "depot") |>
  assignCompartment(Cc = "central", volume = "V") |>
  setCompartmentVolume(unit = "1") |>
  addReaction("Ag", "Cc", "Ka * Ag", "absorption") |>
  addReaction("Cc", "", "Cl * Cc", "elimination", rateCompartment = "unit")
as.eqnvec(reactions)

## ---- model
dm <- c("forward", "reverse", "forward-reverse")
x <- Xs(odemodel(reactions, modelname = "theo", compile = FALSE, deriv2 = TRUE, derivMode = dm))
g <- Y(c(cp = "Cc"), x, modelname = "theo_g", compile = FALSE, deriv2 = TRUE, derivMode = dm)
e <- Y(c(cp = "sigma"), g, attachInput = FALSE, modelname = "theo_e", compile = FALSE, deriv2 = TRUE,
       derivMode = dm)

## ---- data
data(Theoph, package = "datasets")
Theoph$Subject <- as.character(Theoph$Subject)
data <- as.datalist(data.frame(name = "cp", time = Theoph$Time, value = Theoph$conc, sigma = NA_real_,
                               condition = Theoph$Subject))
subjects <- sort(unique(Theoph$Subject))
table <- data.frame(ID = subjects, dose = tapply(Theoph$Dose * Theoph$Wt, Theoph$Subject, `[`, 1)[subjects],
                    weight = tapply(Theoph$Wt, Theoph$Subject, `[`, 1)[subjects], row.names = subjects)

## ---- trafo
trafo <- eqnvec() |>
  define("x <- x", x = c("Ka", "V", "Cl", "Ag", "Cc", "sigma")) |>
  branch(table = table) |>
  insert("Ag <- dose", dose = dose) |>
  insert("Cc <- 0") |>
  insert("x <- exp(x)", x = c("Ka", "V", "Cl", "sigma")) |>
  insert(Ka + V + Cl ~ 1 + (1 | ID))
p <- P(trafo, method = "explicit", modelname = "theo_p", compile = FALSE, deriv2 = TRUE, derivMode = dm)
compile(x, g, e, p, output = "theo", cores = 6)
random <- randomEffects(p)
random

## ---- focei
init <- c(Ka = 0, V = 3, Cl = 0.7, sigma = log(0.5))
fit <- em(obj = normL2(data, g * x * p, errmodel = e), random = random, init = init, verbose = FALSE)
summary(fit)

## ---- compare
obj <- normL2(data, g * x * p, errmodel = e)
set.seed(1)
fits <- list(focei = fit,
             quadrature = em(obj, random, init, method = "quadrature", verbose = FALSE),
             saem = em(obj, random, init, method = "saem", verbose = FALSE))
round(sapply(fits, function(f) c(f$argument[names(fit$argument)], minus2logL = f$value)), 3)

## ---- laplaceem
laplaceEM <- em(obj, random, init, method = "laplaceEM", verbose = FALSE)
round(c(laplaceEM$argument[names(fit$argument)], minus2logL = laplaceEM$value), 3)

## ---- covariates
table$logWeight <- log(table$weight / 70)
trafo2 <- eqnvec() |>
  define("x <- x", x = c("Ka", "V", "Cl", "Ag", "Cc", "sigma")) |>
  branch(table = table) |>
  insert("Ag <- dose", dose = dose) |>
  insert("Cc <- 0") |>
  insert("x <- exp(x)", x = c("Ka", "V", "Cl", "sigma")) |>
  insert(Ka + V ~ 1 + (1 | ID)) |>
  insert(Cl ~ 1 + logWeight + (1 | ID))
p2 <- P(trafo2, method = "explicit", modelname = "theo_p2", compile = FALSE, deriv2 = TRUE, derivMode = dm)
compile(p2, output = "theo2", cores = 6)
fit2 <- em(normL2(data, g * x * p2, errmodel = e), randomEffects(p2), init, verbose = FALSE)
comparison <- anova(fit, fit2)
comparison

## ---- save
# the fits without the model, which an evaluation at the estimates restores
model <- c("prdfn", "data", "errfn", "design")
result <- lapply(list(fit = fit, quadrature = fits$quadrature, saem = fits$saem, laplaceEM = laplaceEM,
                      fit2 = fit2), function(f) { f[model] <- NULL; f })
attr(result, "stamp") <- list(md5 = tools::md5sum("guide-nlme.R"), dMod2 = packageVersion("dMod2"))
saveRDS(result, "guide-nlme.rds")
