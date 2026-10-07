## L1 selection of cell-type-specific parameters with scanL1(), after Hauber,
## Rosenblatt and Timmer (2023). Becker Epo receptor model, six cell lines,
## line A is the reference. Ground truth: kon shared by all lines, ke split
## {A,B,C} vs {D,E,F}, kt in three pairs {A,B} {C,D} {E,F}.
library(dMod2)
outdir <- tempdir()
set.seed(2)

## Model ----------------------------------------------------------------------
reactions <- NULL
reactions <- addReaction(reactions, "Epo + EpoR", "Epo_EpoR",
                         "Epo * EpoR * kon / init_Epo", "binding")
reactions <- addReaction(reactions, "Epo_EpoR", "Epo + EpoR", "Epo_EpoR * koff", "unbinding")
reactions <- addReaction(reactions, "", "EpoR", "4 * init_Epo * init_EpoR_rel * kt", "synthesis")
reactions <- addReaction(reactions, "EpoR", "", "EpoR * kt", "turnover")
reactions <- addReaction(reactions, "Epo_EpoR", "Epo_EpoR_i", "Epo_EpoR * ke", "internalisation")
reactions <- addReaction(reactions, "Epo_EpoR_i", "Epo + EpoR", "Epo_EpoR_i * kex", "recycling")
reactions <- addReaction(reactions, "Epo_EpoR_i", "dEpo_i", "Epo_EpoR_i * kdi", "intra degr.")
reactions <- addReaction(reactions, "Epo_EpoR_i", "dEpo_e", "Epo_EpoR_i * kde", "extra degr.")

x <- Xs(odemodel(reactions, modelname = "beckerL1", compile = FALSE, outdir = outdir),
        compile = FALSE)
g <- Y(eqnvec(y_ext = "log(Epo + dEpo_e + 1)", y_mem = "log(Epo_EpoR + 1)",
              y_int = "log(Epo_EpoR_i + dEpo_i + 1)"),
       x, modelname = "beckerL1_obs", compile = FALSE, attach.input = FALSE,
       outdir = outdir)
e <- Y(eqnvec(y_ext = "sigma", y_mem = "sigma", y_int = "sigma"), g,
       modelname = "beckerL1_err", compile = FALSE, attach.input = FALSE,
       outdir = outdir)

## Fold changes r_<par>_<line> on the log scale, zero for the reference line.
dose  <- 1347.49
lines <- c("A", "B", "C", "D", "E", "F")
cand  <- c("kon", "ke", "kt")
fc    <- function(par, line) if (line == "A") "0" else paste0("r_", par, "_", line)
trafo <- setNames(lapply(lines, function(l) eqnvec(
  kon  = paste0("exp(log_kon + ", fc("kon", l), ")"),
  ke   = paste0("exp(log_ke + ", fc("ke", l), ")"),
  kt   = paste0("exp(log_kt + ", fc("kt", l), ")"),
  koff = "exp(log_koff)", kex = "exp(log_kex)", kdi = "exp(log_kdi)",
  kde = "exp(log_kde)", init_Epo = as.character(dose),
  init_EpoR_rel = "exp(log_EpoR_rel)", Epo = as.character(dose),
  EpoR = paste0("4 * ", dose, " * exp(log_EpoR_rel)"),
  Epo_EpoR = "0", Epo_EpoR_i = "0", dEpo_i = "0", dEpo_e = "0",
  sigma = "exp(log_sigma)")), lines)
p <- P(trafo, modelname = "beckerL1_p", compile = FALSE, outdir = outdir)
compile(g, x, p, e, output = "beckerL1")
prd <- g * x * p

## Data -----------------------------------------------------------------------
st <- c(log_kon = log(0.1513), log_koff = log(0.08069), log_kt = log(0.01601),
        log_ke = log(0.05546), log_kex = log(0.000578), log_kdi = log(0.001273),
        log_kde = log(0.01199), log_EpoR_rel = log(0.09205), log_sigma = log(0.03))
truth <- rbind(kon = c(A = 0, B = 0,   C = 0,   D = 0,    E = 0,    F = 0),
               ke  = c(A = 0, B = 0,   C = 0,   D = -1.6, E = -1.6, F = -1.6),
               kt  = c(A = 0, B = 0,   C = -0.7, D = -0.7, E = -1.4, F = -1.4))
r_true <- unlist(lapply(cand, function(pp)
  setNames(truth[pp, -1], paste0("r_", pp, "_", lines[-1]))))
times <- c(3, 6, 12, 20, 30, 45, 60, 90, 120, 180, 240, 300)
sim <- wide2long(prd(times, c(st, r_true), deriv = FALSE))
sim$value <- sim$value + rnorm(nrow(sim), 0, 0.03)
sim$sigma <- NA
data <- as.datalist(sim)
plot(prd(times, c(st, r_true)), data)

obj    <- normL2(data, prd, errmodel = e)
fixed  <- st[c("log_koff", "log_kex", "log_kdi", "log_kde", "log_EpoR_rel")]
center <- c(st[c(paste0("log_", cand), "log_sigma")], r_true * 0)

## Cluster penalty: one block per parameter, anchored at the reference line.
blocks <- lapply(cand, function(pp)
  list(pars = paste0("r_", pp, "_", lines[-1]), anchor = 0))
grid <- 10^seq(-1, 3, length.out = 41)

t_cl <- system.time(
  scanCl <- scanL1(obj, center, groups = blocks, fixed = fixed, lambda = grid,
                   fits = 10, cores = 4, sd = 0.3))
t_cl
scanCl
plot(scanCl)
plot(scanCl, type = "test")
scanCl$structure[[scanCl$selected]]$groups

## The same with the L0.8 penalty of the paper.
t_q <- system.time(
  scanQ <- scanL1(obj, center, groups = blocks, fixed = fixed, lambda = grid,
                  fits = 10, cores = 4, sd = 0.3, q = 0.8))
t_q
scanQ
plot(scanQ)

## Reference penalty only: which fold changes differ from line A at all.
t_ref <- system.time(
  scanRef <- scanL1(obj, center, reference = names(r_true), fixed = fixed,
                    lambda = grid, fits = 10, cores = 4, sd = 0.3))
t_ref
scanRef
plot(scanRef)
setdiff(names(r_true), scanRef$structure[[scanRef$selected]]$removed)
