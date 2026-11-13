# Multistarts of the theophylline model of guide chapter 7 by FOCEI, quadrature and SAEM, and
# the two optima they reach, each evaluated by every method. Saves results/multistart.rds, which
# the guide reads.

library(dMod2)

## ---- theophylline
reactions <- eqnlist() |>
  assignCompartment(Ag = "depot") |>
  assignCompartment(Cc = "central", volume = "V") |>
  setCompartmentVolume(unit = "1") |>
  addReaction("Ag", "Cc", "Ka * Ag", "absorption") |>
  addReaction("Cc", "", "Cl * Cc", "elimination", rateCompartment = "unit")
dm <- c("forward", "reverse", "forward-reverse")
x <- Xs(odemodel(reactions, modelname = "theo", compile = FALSE, deriv2 = TRUE, derivMode = dm))
g <- Y(c(cp = "Cc"), x, modelname = "theo_g", compile = FALSE, deriv2 = TRUE, derivMode = dm)
e <- Y(c(cp = "sigma"), g, attachInput = FALSE, modelname = "theo_e", compile = FALSE, deriv2 = TRUE,
       derivMode = dm)
data(Theoph, package = "datasets")
Theoph$Subject <- as.character(Theoph$Subject)
data <- as.datalist(data.frame(name = "cp", time = Theoph$Time, value = Theoph$conc, sigma = NA_real_,
                               condition = Theoph$Subject))
subjects <- sort(unique(Theoph$Subject))
table <- data.frame(ID = subjects, dose = tapply(Theoph$Dose * Theoph$Wt, Theoph$Subject, `[`, 1)[subjects],
                    weight = tapply(Theoph$Wt, Theoph$Subject, `[`, 1)[subjects], row.names = subjects)
trafo <- branch(define(NULL, "x <- x", x = c("Ka", "V", "Cl", "Ag", "Cc", "sigma")), table = table) |>
  insert("Ag <- dose", dose = dose) |>
  insert("Cc <- 0") |>
  insert("x <- exp(x)", x = c("Ka", "V", "Cl", "sigma")) |>
  insert(Ka + V + Cl ~ 1 + (1 | ID))
p <- P(trafo, method = "explicit", modelname = "theo_p", compile = FALSE, deriv2 = TRUE, derivMode = dm)
compile(x, g, e, p, output = "theo")
random <- randomEffects(p)
obj <- normL2(data, g * x * p, errmodel = e)
init <- c(Ka = 0, V = 3, Cl = 0.7, sigma = log(0.5))

## ---- multistart
# 200 starts per method, from the same draws, each run to the end; SAEM with up to 500 burn-in and
# 1000 averaging iterations
methods <- c(FOCEI = "focei", quadrature = "quadrature", SAEM = "saem")
controls <- list(FOCEI = list(), quadrature = list(), SAEM = list(saem = list(nBurnin = 500, nEM = 1000)))
ms <- lapply(names(methods), function(m) {
  set.seed(1)
  msem(obj, random, init, method = methods[[m]], fits = 200, sd = 1,
       control = controls[[m]], cores = getOption("dMod2.cores", 1))
})
names(ms) <- names(methods)

## ---- crosscheck
# the best fit of every method refitted by every other method from the estimates of that fit, so that the
# values of the optima compare on the scale of each method
best <- lapply(ms, function(f) f$argument)
cross <- do.call(rbind, lapply(names(best), function(from) do.call(rbind, lapply(names(methods), function(by) {
  fit <- tryCatch(em(obj, random, best[[from]], method = methods[[by]], control = controls[[by]],
                     verbose = FALSE),
                  error = function(err) NULL)
  data.frame(from = from, by = by, value = if (is.null(fit)) NA else fit$value,
             converged = !is.null(fit) && isTRUE(fit$converged),
             Ka = if (is.null(fit)) NA else exp(fit$argument[["Ka"]]),
             V = if (is.null(fit)) NA else exp(fit$argument[["V"]]))
}))))

## ---- save
dir.create("results", showWarnings = FALSE)
saved <- list(multistart = lapply(ms, `[[`, "multistart"), best = best, cross = cross)
# the md5 of the script, which the guide checks on reading
attr(saved, "stamp") <- list(md5 = tools::md5sum("theophyllineMultistart.R"), dMod2 = packageVersion("dMod2"))
saveRDS(saved, "results/multistart.rds")
