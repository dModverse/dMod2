# -------------------------------------------------------------------------#
# JAK2-STAT5 signalling in erythroid progenitor cells
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# Bachmann et al. (2011) in dMod2: Epo induced JAK2-STAT5 signalling with the
# three negative feedbacks CIS, SOCS3 and SHP1. 541 measurements from thirteen
# experiments, 113 estimated parameters. The reaction network follows the
# Data2Dynamics model `jak2_stat5_feedbacks`; the objective is evaluated at the
# optimum the benchmark collection publishes.
#
# [AUTHOR]
# Simon Beyer
#
# [Date]
# Wed 03 Sep 2026
#
# [Info]
# The data ship with the package as `bachmann`, the PEtab form of the same
# problem as inst/extdata/petab_bachmann. The last section imports that one and
# checks the two against each other.
# -------------------------------------------------------------------------#

library(dMod2)
library(ggplot2)

.modelname <- "bachmann"
# every generated source, object and shared library goes here, never into the
# working directory
.outdir    <- file.path(tempdir(), .modelname)
.fit       <- FALSE   # multi-start search, hours on this model

if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)
.petabDir <- system.file("extdata", "petab_bachmann", package = "dMod2")


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Load data
#
# Values are the recorded signal divided by the maximum of its own column, the
# gauge the scale parameters are estimated in. All but one observable are fitted
# on log10, which is where their estimated standard deviations live as well.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
data(bachmann)
.logFitted <- bachmann$name != "pSTAT5B_rel"
bachmann$value[.logFitted] <- log10(bachmann$value[.logFitted])

.covariates <- c("experiment", "epo_level", "ActD", "CISoe", "SOCS3oe", "SHP1oe")
mydataL   <- as.datalist(bachmann, split.by = "condition",
                         keep.covariates = .covariates)
cond.grid <- covariates(mydataL)

cat(sprintf("%d points, %d conditions, %d observables\n", nrow(bachmann),
            length(mydataL), length(unique(bachmann$name))))


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Define model
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Rate expressions are written in the compartment of their educts. npSTAT5 and
# the nuclear transcripts are placed first, so that the reactions producing them
# do not claim them for the cytoplasm.

addCyt <- function(eq, from, to, rate) addReaction(eq, from, to, rate, compartment = "cyt")
addNuc <- function(eq, from, to, rate) addReaction(eq, from, to, rate, compartment = "nuc")

.nRNA <- function(eq, gene) {
  # five nuclear stages give the transcript the measured delay before it reaches
  # the cytoplasm, without a delay differential equation
  st <- paste0(gene, "nRNA", 1:5)
  eq <- addNuc(eq, "", st[1],
               sprintf("%sRNAEqc * %sRNATurn * npSTAT5 * (1 - ActD)", gene, gene))
  for (i in 1:4)
    eq <- addNuc(eq, st[i], st[i + 1], sprintf("%sRNADelay * %s", gene, st[i]))
  eq <- addNuc(eq, st[5], paste0(gene, "RNA"),
               sprintf("%sRNADelay * %s", gene, st[5]))
  addCyt(eq, paste0(gene, "RNA"), "", sprintf("%sRNATurn * %sRNA", gene, gene))
}

reactions <- eqnlist() |>
  assignCompartment(npSTAT5 = "nuc", CISRNA = "cyt", SOCS3RNA = "cyt") |>

  # receptor activation. SOCS3 blocks every phosphorylation step, CIS only the
  # second site, and SHP1 returns all four forms to the unphosphorylated one
  addCyt("EpoRJAK2", "EpoRpJAK2",
         "JAK2ActEpo * EpoRJAK2 * epo_level / (1 + SOCS3Inh * SOCS3)") |>
  addCyt("EpoRpJAK2", "p1EpoRpJAK2",
         "EpoRActJAK2 * EpoRpJAK2 / (1 + SOCS3Inh * SOCS3)") |>
  addCyt("EpoRpJAK2", "p2EpoRpJAK2",
         "3 * EpoRActJAK2 * EpoRpJAK2 / (1 + SOCS3Inh * SOCS3) / (EpoRCISInh * EpoRJAK2_CIS + 1)") |>
  addCyt("p1EpoRpJAK2", "p12EpoRpJAK2",
         "3 * EpoRActJAK2 * p1EpoRpJAK2 / (1 + SOCS3Inh * SOCS3) / (EpoRCISInh * EpoRJAK2_CIS + 1)") |>
  addCyt("p2EpoRpJAK2", "p12EpoRpJAK2",
         "EpoRActJAK2 * p2EpoRpJAK2 / (1 + SOCS3Inh * SOCS3)") |>
  addCyt("EpoRpJAK2",    "EpoRJAK2", "JAK2EpoRDeaSHP1 * EpoRpJAK2 * SHP1Act") |>
  addCyt("p1EpoRpJAK2",  "EpoRJAK2", "JAK2EpoRDeaSHP1 * p1EpoRpJAK2 * SHP1Act") |>
  addCyt("p2EpoRpJAK2",  "EpoRJAK2", "JAK2EpoRDeaSHP1 * p2EpoRpJAK2 * SHP1Act") |>
  addCyt("p12EpoRpJAK2", "EpoRJAK2", "JAK2EpoRDeaSHP1 * p12EpoRpJAK2 * SHP1Act") |>

  # CIS occupies the receptor rather than being consumed by it, so the pool is
  # counted as an occupied fraction that the active receptor clears
  addCyt("EpoRJAK2_CIS", "",
         "EpoRCISRemove * EpoRJAK2_CIS * (p1EpoRpJAK2 + p12EpoRpJAK2)") |>

  addCyt("SHP1", "SHP1Act",
         "SHP1ActEpoR * SHP1 * (EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2 + p12EpoRpJAK2)") |>
  addCyt("SHP1Act", "SHP1", "SHP1Dea * SHP1Act") |>

  # two routes to pSTAT5: JAK2 phosphorylates on any active receptor, the
  # receptor itself only on the doubly phosphorylated one and cooperatively
  addCyt("STAT5", "pSTAT5",
         "STAT5ActJAK2 * STAT5 * (EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2 + p12EpoRpJAK2) / (1 + SOCS3Inh * SOCS3)") |>
  addCyt("STAT5", "pSTAT5",
         "STAT5ActEpoR * STAT5 * (p1EpoRpJAK2 + p12EpoRpJAK2)^2 / (1 + SOCS3Inh * SOCS3) / (1 + CISInh * CIS)") |>
  addCyt("pSTAT5", "npSTAT5", "STAT5Imp * pSTAT5") |>
  addNuc("npSTAT5", "STAT5", "STAT5Exp * npSTAT5") |>

  .nRNA("CIS") |>
  .nRNA("SOCS3") |>

  # feedback proteins, each with an overexpression term the OE experiments switch on
  addCyt("", "CIS",   "CISEqc * CISTurn * CISRNA") |>
  addCyt("CIS", "",   "CISTurn * CIS") |>
  addCyt("", "CIS",   "CISoe * CISEqcOE * CISTurn") |>
  addCyt("", "SOCS3", "SOCS3Eqc * SOCS3Turn * SOCS3RNA") |>
  addCyt("SOCS3", "", "SOCS3Turn * SOCS3") |>
  addCyt("", "SOCS3", "SOCS3oe * SOCS3EqcOE * SOCS3Turn") |>

  setCompartmentVolume(cyt = "0.4", nuc = "0.275")

reactions

# the floor keeps log10 finite where a species is still empty. It matches what
# importPEtab() applies, so the two objectives further down stay comparable.
.log10 <- function(x) paste0("log10(", x, " + 1e-10)")
.pJAK2 <- "2 * (EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2 + p12EpoRpJAK2)"
.pEpoR <- "16 * (p1EpoRpJAK2 + p2EpoRpJAK2 + p12EpoRpJAK2)"

# scale_ and offset_ carry no experiment here; the transformation below resolves
# them per condition, which is what keeps one observation function for all of them
observables <- do.call(eqnvec, c(list(
  pJAK2_au    = .log10(paste("offset_pJAK2 + scale_pJAK2 / init_EpoRJAK2 *", .pJAK2)),
  pEpoR_au    = .log10(paste("offset_pEpoR + scale_pEpoR / init_EpoRJAK2 *", .pEpoR)),
  CIS_au      = .log10("offset_CIS + scale_CIS / CISEqc / CISRNAEqc / init_STAT5 * CIS"),
  SOCS3_au    = .log10("offset_SOCS3 + scale_SOCS3 / SOCS3Eqc / SOCS3RNAEqc / init_STAT5 * SOCS3"),
  tSTAT5_au   = .log10("scale_tSTAT5 / init_STAT5 * (STAT5 + pSTAT5)"),
  pSTAT5_au   = .log10("offset_pSTAT5 + scale_pSTAT5 / init_STAT5 * pSTAT5"),
  tSHP1_au    = .log10("scale_tSHP1 / init_SHP1 * (SHP1 + SHP1Act)"),
  CIS_au1     = .log10("scale1_CIS_dr90 / CISEqc / CISRNAEqc / init_STAT5 * CIS"),
  CIS_au2     = .log10("scale2_CIS_dr90 / CISEqc / CISRNAEqc / init_STAT5 * CIS"),
  STAT5_abs   = .log10("STAT5"),
  SHP1_abs    = .log10("SHP1 + SHP1Act"),
  CIS_abs     = .log10("CIS"),
  SOCS3_abs   = .log10("SOCS3"),
  pSTAT5B_rel = "offset_pSTAT5_conc + 100 * pSTAT5 / (pSTAT5 + STAT5)"),
  # the qPCR panels were run three times, each replicate with its own scale
  setNames(as.list(.log10(paste0("1 + scale_CISRNA_fold", LETTERS[1:3],
                                 " / CISRNAEqc / init_STAT5 * CISRNA"))),
           paste0("CISRNA_fold", LETTERS[1:3])),
  setNames(as.list(.log10(paste0("1 + scale_SOCS3RNA_fold", LETTERS[1:3],
                                 " / SOCS3RNAEqc / init_STAT5 * SOCS3RNA"))),
           paste0("SOCS3RNA_fold", LETTERS[1:3]))))

stopifnot(setequal(names(observables), unique(bachmann$name)))

# one standard deviation per antibody, shared by every experiment using it
errorModels <- eqnvec(
  pJAK2_au = "sd_JAK2EpoR_au", pEpoR_au = "sd_JAK2EpoR_au",
  CIS_au   = "sd_CIS_au", CIS_au1 = "sd_CIS_au", CIS_au2 = "sd_CIS_au",
  SOCS3_au = "sd_SOCS3_au", tSTAT5_au = "sd_STAT5_au", tSHP1_au = "sd_SHP1_au",
  pSTAT5_au = "sd_pSTAT5_au",
  STAT5_abs = "sd_STAT5_abs", SHP1_abs = "sd_SHP1_abs",
  CIS_abs = "sd_CIS_abs", SOCS3_abs = "sd_SOCS3_abs",
  pSTAT5B_rel = "sd_pSTAT5_rel",
  CISRNA_foldA = "sd_RNA_fold", CISRNA_foldB = "sd_RNA_fold",
  CISRNA_foldC = "sd_RNA_fold", SOCS3RNA_foldA = "sd_RNA_fold",
  SOCS3RNA_foldB = "sd_RNA_fold", SOCS3RNA_foldC = "sd_RNA_fold")


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Condition table
#
# Every generic scale and offset is resolved here. An observable an experiment
# did not record is pinned, so it contributes no free direction to the fit.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
.recordedBy <- list(
  pJAK2  = c("long", "actd", "fine", "cisoe", "socs3oe", "shp1oe", "dr7", "dr30"),
  pEpoR  = c("long", "actd", "fine", "cisoe", "cisoe_pepor", "socs3oe", "shp1oe",
             "dr7", "dr30"),
  CIS    = c("long", "actd", "cisoe", "socs3oe", "shp1oe"),
  SOCS3  = c("long", "cisoe", "socs3oe"),
  tSTAT5 = c("long", "actd", "shp1oe"),
  pSTAT5 = c("long", "actd", "cisoe", "socs3oe", "shp1oe", "dr10"),
  tSHP1  = "shp1oe")
# the parameter stem is the observable name except for total SHP1
.stemOf   <- c(pJAK2 = "pJAK2", pEpoR = "pEpoR", CIS = "CIS", SOCS3 = "SOCS3",
               tSTAT5 = "tSTAT5", pSTAT5 = "pSTAT5", tSHP1 = "SHP1")
# a total has no background, and the pSTAT5 dose response was blotted without one
.offsetIn <- setdiff(names(.stemOf), c("tSTAT5", "tSHP1"))

.exp <- cond.grid$experiment
for (b in names(.stemOf)) {
  seen <- .exp %in% .recordedBy[[b]]
  cond.grid[[paste0("scale_", b)]] <-
    ifelse(seen, paste0("scale_", .stemOf[[b]], "_", .exp), "1")
  if (b %in% .offsetIn)
    cond.grid[[paste0("offset_", b)]] <-
      ifelse(seen & !(b == "pSTAT5" & .exp == "dr10"),
             paste0("offset_", .stemOf[[b]], "_", .exp), "0")
}
# SOCS3 overexpression widens the pSTAT5 blot beyond its usual spread
cond.grid$sd_pSTAT5_au <- ifelse(cond.grid$SOCS3oe == 1,
                                 "sd_STAT5_au + sd_pSTAT5_socs3oe", "sd_STAT5_au")

# Each perturbation experiment holds a control arm and a perturbed one; the
# label separates the two and is what the plots colour by. Writing the grid
# back makes it the covariate table the plotting functions read.
cond.grid$perturbation <- with(cond.grid,
  ifelse(ActD == 1, "actinomycin D",
  ifelse(CISoe == 1, "CIS overexpressed",
  ifelse(SOCS3oe == 1, "SOCS3 overexpressed",
  ifelse(SHP1oe == 1, "SHP1 overexpressed", "control")))))
attr(mydataL, "condition.grid") <- cond.grid

# Which columns of the grid are parameter substitutions. The rest, `experiment`
# and `perturbation`, are labels the plots read and must not reach `branch()`.
.substitutions <- c("epo_level", "ActD", "CISoe", "SOCS3oe", "SHP1oe",
                    "sd_pSTAT5_au",
                    grep("^(scale|offset)_", names(cond.grid), value = TRUE))


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Build
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
myOptionsODE  <- list(atol = 1e-11, rtol = 1e-8, maxsteps = 1e7L, maxattemps = 100L)
myOptionsSens <- list(atol = 1e-11, rtol = 1e-8, maxsteps = 1e7L, maxattemps = 100L)

model <- odemodel(reactions, modelname = "bachmann_ode", compile = FALSE, outdir = .outdir)
x <- Xs(model, optionsOde = myOptionsODE, optionsSens = myOptionsSens)
g <- Y(observables, x, modelname = "bachmann_obs", attach.input = FALSE,
       compile = FALSE, outdir = .outdir)
e <- Y(errorModels, g, modelname = "bachmann_err", attach.input = FALSE,
       compile = FALSE, outdir = .outdir)

# only the receptor, SHP1 and STAT5 pools are stocked; the overexpression
# experiments start with their protein already present
.zero <- setdiff(reactions$states,
                 c("EpoRJAK2", "SHP1", "STAT5", "EpoRJAK2_CIS", "CIS", "SOCS3"))

trafo <- eqnvec() |>
  define("x~x", x = getParameters(model, g, e)) |>
  define("x~0", x = .zero) |>
  insert("EpoRJAK2 ~ init_EpoRJAK2") |>
  insert("STAT5 ~ init_STAT5") |>
  insert("SHP1 ~ init_SHP1 * (1 + SHP1oe * SHP1ProOE)") |>

  # The model was written in parameters that carry their own reference scale, so
  # the estimated ones are dimensionless. Each rule uses the raw parameter on its
  # right hand side, which is what fixes the order: a parameter is rewritten
  # before it is used, never after.
  insert("CISRNAEqc ~ CISRNAEqc / init_STAT5") |>
  insert("SOCS3RNAEqc ~ SOCS3RNAEqc / init_STAT5") |>
  insert("CISEqc ~ CISEqc / CISRNAEqc") |>
  insert("SOCS3Eqc ~ SOCS3Eqc / SOCS3RNAEqc") |>
  insert("CISEqcOE ~ CISEqcOE * CISEqc") |>
  insert("SOCS3EqcOE ~ SOCS3EqcOE * SOCS3Eqc") |>
  insert("CISInh ~ CISInh / CISEqc") |>
  insert("SOCS3Inh ~ SOCS3Inh / SOCS3Eqc") |>
  insert("JAK2EpoRDeaSHP1 ~ JAK2EpoRDeaSHP1 / init_SHP1") |>
  insert("EpoRCISRemove ~ EpoRCISRemove / init_EpoRJAK2") |>
  insert("SHP1ActEpoR ~ SHP1ActEpoR / init_EpoRJAK2") |>
  insert("STAT5ActJAK2 ~ STAT5ActJAK2 / init_EpoRJAK2") |>
  insert("STAT5ActEpoR ~ STAT5ActEpoR / init_EpoRJAK2^2") |>

  # the two overexpressed pools start at the level their production term holds
  insert("CIS ~ CISoe * CISEqcOE * CISEqc") |>
  insert("EpoRJAK2_CIS ~ CISoe") |>
  insert("SOCS3 ~ SOCS3oe * SOCS3EqcOE * SOCS3Eqc") |>

  branch(table = cond.grid[.substitutions], apply = "insert") |>
  # both RNA equilibrium constants only ever appear next to a scale parameter
  insert("x~1", x = c("CISRNAEqc", "SOCS3RNAEqc")) |>
  insert("x~exp10(x)", x = .currentSymbols)

p <- P(trafo, modelname = "bachmann_trafo", compile = FALSE, outdir = .outdir)

compile(g, x, p, e, output = .modelname, cores = 12)

prd       <- g*x*p
outerpars <- getParameters(prd)
cat(length(outerpars), "estimated parameters\n")


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Objective
#
# Bounds, the published optimum and the one prior in the problem all come from
# the PEtab parameter table, so this script and the import below start from the
# same numbers. The prior is on the receptor pool, which the paper measured at
# 4.15 nM with 30 percent spread; everything else the bounds hold in.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
.pars  <- read.delim(file.path(.petabDir, "parameters_Bachmann_MSB2011.tsv"))
.prior <- as.numeric(strsplit(
  .pars$objectivePriorParameters[.pars$parameterId == "init_EpoRJAK2"], ";")[[1]])
.pars  <- .pars[.pars$estimate == 1, ]

obj <- normL2(mydataL, prd, e) +
  constraintL2(setNames(.prior[1], "init_EpoRJAK2"), sigma = .prior[2],
               attr.name = "prior")

# the published optimum, on the log10 scale the table declares
bestfit <- setNames(log10(.pars$nominalValue), .pars$parameterId)
stopifnot(setequal(names(bestfit), outerpars))

obj(bestfit, cores = 10) |> system.time()


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Multi-start fit
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
if (.fit) {
  pouter <- structure(rep(-1, length(outerpars)), names = outerpars)
  fits <- mstrust(obj, pouter, sd = 4, fits = 500, cores = 10, rinit = 0.1, rmax = 10,  iterlim = 5000,
                  parlower = setNames(log10(.pars$lowerBound), .pars$parameterId),
                  parupper = setNames(log10(.pars$upperBound), .pars$parameterId))
  bestfit <- as.parvec(as.parframe(fits))
}


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Predictions over data
#
# One observation function serves all 36 conditions, so a prediction also holds
# observables an experiment never recorded. Those are dropped everywhere below,
# which is what keeps the axes on the measured range. Time courses and dose
# responses do not share an x axis and are drawn apart.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
recorded <- unique(paste(bachmann$condition, bachmann$name))

# prediction as a data.frame, with the error model as a `sigma` column and the
# condition grid attached, which is what the ribbons and the colours read from
.frame <- function(conditions, times)
  as.data.frame(prd(times, bestfit, conditions = conditions, deriv = FALSE),
                data = mydataL[conditions], errfn = e) |>
  subset(paste(condition, name) %in% recorded)

# One standard deviation of the error model, in the colour of its own line and
# dashed at the edge. plotCombined() already carries the package colour scale,
# so only the fill has to be set to the same palette. `alpha` reaches the fill
# only, which keeps the edge crisp.
.band <- function(pl, df, colour = "condition", x = "time")
  pl + geom_ribbon(data = df,
                   aes(x = .data[[x]], ymin = value - sigma, ymax = value + sigma,
                       group = condition, colour = .data[[colour]],
                       fill = .data[[colour]]),
                   inherit.aes = FALSE, alpha = 0.15, linetype = "dashed",
                   linewidth = 0.3) +
    scale_fill_dMod()

times <- seq(0, 250, length.out = 500)

## --- the reference time course, all six observables ---
plotCombined(prd(times, bestfit, conditions = "long", deriv = FALSE),
             mydataL["long"], paste(condition, name) %in% recorded) |>
  .band(.frame("long", times)) +
  labs(x = "time [min]", y = "log10 signal [a.u.]", colour = NULL, fill = NULL,
       title = "CFU-E, 5 units/ml Epo")

## --- what each perturbation does, control against perturbed ---
.pert <- rownames(cond.grid)[cond.grid$experiment %in%
                             c("actd", "cisoe", "socs3oe", "shp1oe")]
plotCombined(prd(times, bestfit, conditions = .pert, deriv = FALSE),
             mydataL[.pert],
             paste(condition, name) %in% recorded,
             name %in% c("pJAK2_au", "pEpoR_au", "pSTAT5_au", "CIS_au"),
             aesthetics = list(color = "perturbation", group = "condition")) |>
  .band(subset(.frame(.pert, times),
               name %in% c("pJAK2_au", "pEpoR_au", "pSTAT5_au", "CIS_au")),
        colour = "perturbation") +
  facet_grid(name ~ experiment, scales = "free_y") +
  labs(x = "time [min]", y = "log10 signal [a.u.]", colour = NULL, fill = NULL)

## --- dose response: the dose swept through the model, data on top ---
# Each dose is a condition of its own, so the recorded points alone would give a
# four-point curve. `epo_level` is an inner parameter, so the dose can instead be
# swept through one condition of the experiment, which keeps that experiment's
# scales and its measurement time and draws the response the model actually
# predicts between the doses.
.drData <- subset(bachmann, grepl("^dr", experiment))
.dose   <- 10^seq(log10(min(.drData$epo_level)), log10(max(.drData$epo_level)),
                  length.out = 300)
.inner  <- p(bestfit)

.sweep <- function(ex) {
  cn   <- rownames(cond.grid)[cond.grid$experiment == ex][1]
  tt   <- unique(.drData$time[.drData$experiment == ex])
  seen <- unique(.drData$name[.drData$experiment == ex])
  pin  <- setNames(as.numeric(.inner[[cn]]), names(.inner[[cn]]))
  do.call(rbind, lapply(.dose, function(d) {
    pin["epo_level"] <- d
    out <- as.data.frame((g * x)(sort(unique(c(0, tt))), pin, deriv = FALSE),
                         errfn = e)
    cbind(out[out$time == tt & out$name %in% seen, c("name", "value", "sigma")],
          epo_level = d, experiment = ex)
  }))
}
.drPred <- do.call(rbind, lapply(unique(.drData$experiment), .sweep))

ggplot(.drPred, aes(epo_level, value, colour = experiment, fill = experiment)) +
  geom_ribbon(aes(ymin = value - sigma, ymax = value + sigma),
              alpha = 0.15, linetype = "dashed", linewidth = 0.3) +
  geom_line() +
  geom_point(data = .drData) +
  facet_wrap(~name, scales = "free_y") +
  scale_x_log10() + scale_color_dMod() + scale_fill_dMod() +
  theme_dMod() + 
  labs(x = "Epo [units/cell]", y = "log10 signal [a.u.]", colour = NULL, fill = NULL)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Cross-check against the PEtab form of the same problem
#
# The benchmark collection distributes Bachmann as PEtab, and the package ships
# that copy. Importing it builds the model a second time, from SBML and the four
# tables instead of from the reactions above, so the two are independent up to
# the data they share.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
petab <- importPEtab(file.path(.petabDir, "Bachmann_MSB2011.yaml"),
                     backend = "cppDE", cores = 4, outdir = .outdir,
                     optionsOde = myOptionsODE, optionsSens = myOptionsSens)

stopifnot(setequal(names(petab$bestfit), outerpars),
          isTRUE(all.equal(petab$bestfit[outerpars], bestfit[outerpars])))

# The two values differ by a constant: PEtab defines the likelihood on the
# linear measurement, so importPEtab() adds the Jacobian of the log10 transform,
# and it normalises the prior density. The sum of squares carries neither, so
# that is what the two models have to agree on.
chi2 <- c(hand  = attr(obj(bestfit, deriv = FALSE),       "chi2"),
          PEtab = attr(petab$obj(bestfit, deriv = FALSE), "chi2"))
chi2
