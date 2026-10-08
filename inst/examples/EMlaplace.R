\dontrun{

## ---------------------------------------------------------------------------
## Marginal-likelihood L1/Laplace model selection: which parameters are
## cell-line specific? (Becker et al. 2010, Epo receptor; 2 cell lines.)
##
## Classical workflow: L1-penalise the per-line parameter deviations, SCAN the
## penalty strength lambda over a grid, run a multistart fit at every grid
## point, and pick lambda by BIC. Cost: O(grid x multistart) fits.
##
## EM reframes the penalised fit as a nonlinear mixed-effects model with a
## Laplace random-effect density: the deviations are random effects and the
## single lambda is ESTIMATED by marginal maximum likelihood (one fit, no grid).
## sparsify then walks a short nested chain of supports and keeps the one with
## the smallest marginal -2 log L. It recovers the same parsimonious model the
## BIC-scan would, at O(K) fits instead of O(grid x multistart). Theory:
## notes/laplace_nlme_theory.Rmd.
##
## Reference encoding (Hauber et al. 2023 style): line A is the baseline
## (theta = exp(mu)); line B carries the deviation d (theta = exp(mu + d)) and
## |d| is penalised. d is a single coordinate, so trustL1 drives the deviation
## of a shared parameter cleanly to zero.
##
## Ground truth here: kon and ke differ between the two lines (individual),
## koff and kt are shared. Candidates = {kon, koff, kt, ke}.
## ---------------------------------------------------------------------------

library(dMod2)
outdir <- tempdir()          # keep generated C/C++ + shared objects out of the tree
set.seed(7)

## 1. Model: the Becker Epo receptor network, built inline (no SBML / Python). -
reactions <- NULL
reactions <- addReaction(reactions, "Epo + EpoR", "Epo_EpoR",
                         "Epo * EpoR * kon / init_Epo", "Epo-receptor binding")
reactions <- addReaction(reactions, "Epo_EpoR", "Epo + EpoR",
                         "Epo_EpoR * koff", "unbinding")
reactions <- addReaction(reactions, "", "EpoR",
                         "4 * init_Epo * init_EpoR_rel * kt", "receptor synthesis")
reactions <- addReaction(reactions, "EpoR", "", "EpoR * kt", "receptor turnover")
reactions <- addReaction(reactions, "Epo_EpoR", "Epo_EpoR_i",
                         "Epo_EpoR * ke", "internalisation")
reactions <- addReaction(reactions, "Epo_EpoR_i", "Epo + EpoR",
                         "Epo_EpoR_i * kex", "recycling")
reactions <- addReaction(reactions, "Epo_EpoR_i", "dEpo_i",
                         "Epo_EpoR_i * kdi", "intracellular degradation")
reactions <- addReaction(reactions, "Epo_EpoR_i", "dEpo_e",
                         "Epo_EpoR_i * kde", "extracellular degradation")

m <- odemodel(reactions, modelname = "becker", backend = "cppDE",
              compile = FALSE, deriv2 = TRUE, outdir = outdir)
x <- Xs(m, compile = FALSE)
g <- Y(eqnvec(y_ext = "log(Epo + dEpo_e + 1)",
              y_mem = "log(Epo_EpoR + 1)",
              y_int = "log(Epo_EpoR_i + dEpo_i + 1)"),
       x, modelname = "becker_obs", compile = FALSE, deriv2 = TRUE,
       attachInput = FALSE, outdir = outdir)
e <- Y(eqnvec(y_ext = "sigma", y_mem = "sigma", y_int = "sigma"), g,
       modelname = "becker_err", compile = FALSE, deriv2 = TRUE,
       attachInput = FALSE, outdir = outdir)

## 2. Reference encoding: candidates kon,koff,kt,ke carry a line-B deviation.
dose  <- 1347.49                                       # init_Epo
lines <- c("A", "B")
trafo <- eqnvec(
  kon  = "exp(log_kon  + eta_kon)",  koff = "exp(log_koff + eta_koff)",
  kt   = "exp(log_kt   + eta_kt)",   ke   = "exp(log_ke   + eta_ke)",
  kex  = "exp(log_kex)", kdi = "exp(log_kdi)", kde = "exp(log_kde)",
  init_Epo = "dose", init_EpoR_rel = "exp(log_EpoR_rel)",
  Epo = "dose", EpoR = "4 * dose * exp(log_EpoR_rel)",
  Epo_EpoR = "0", Epo_EpoR_i = "0", dEpo_i = "0", dEpo_e = "0",
  sigma = "exp(log_sigma)")
## line A: eta = "0" (baseline); line B: eta = the estimated difference.
cl_table <- data.frame(
  eta_kon = c("0", "eta_kon_B"), eta_koff = c("0", "eta_koff_B"),
  eta_kt  = c("0", "eta_kt_B"),  eta_ke   = c("0", "eta_ke_B"),
  dose = c(dose, dose), row.names = lines, stringsAsFactors = FALSE)
p <- P(branch(trafo, table = cl_table, apply = "insert"), method = "explicit",
       modelname = "becker_p", compile = FALSE, deriv2 = TRUE, outdir = outdir)
prd <- g * x * p
compile(prd, e, cores = 4)

## 3. Simulate two cell lines: kon, ke individual; koff, kt shared. -----------
st <- c(log_kon = log(0.1513), log_koff = log(0.08069), log_kt = log(0.01601),
        log_ke = log(0.05546), log_kex = log(0.000578), log_kdi = log(0.001273),
        log_kde = log(0.01199), log_EpoR_rel = log(0.09205),
        log_sigma = log(0.05))
d_true <- c(eta_kon_B = 0.8, eta_koff_B = 0, eta_kt_B = 0, eta_ke_B = 1.0)
times  <- c(3, 6, 12, 20, 30, 45, 60, 90, 120, 180, 240, 300)
sim <- prd(times, c(st, d_true), deriv = FALSE)
dl_rows <- do.call(rbind, lapply(lines, function(s) {
  pr <- sim[[s]]
  do.call(rbind, lapply(c("y_ext", "y_mem", "y_int"), function(nm)
    data.frame(name = nm, time = pr[, "time"],
               value = pr[, nm] + rnorm(nrow(pr), 0, 0.05),
               sigma = NA_real_, condition = s, stringsAsFactors = FALSE)))
}))
dlist <- as.datalist(dl_rows)

## 4. Penalised objective + start (emInit adds the single lambda). -----------
cand   <- c("kon", "koff", "kt", "ke")
pen    <- penaltyL1(cand, subjects = "B")
obj    <- normL2(dlist, prd, errmodel = e) + constraintL1(pen)
fixed  <- st[c("log_kex", "log_kdi", "log_kde", "log_EpoR_rel")]
center <- emInit(c(st[paste0("log_", cand)], log_sigma = log(0.1)), pen,
                  lambda = 1)

## 5a. EM: ONE mixed-effects fit estimates lambda and ranks |d_hat|. ------
fit <- msEM(obj, center, fixed = fixed, method = "focei",
                control = list(cm1 = list(iterlim = 50L)),
                fits = 16, cores = 4, sd = 0.5, verbose = TRUE)
print(fit)
cat("estimated lambda:", fit$lambda, "\n")
print(round(fit$etaModes, 4))              # line-B deviations d_hat

## 5b. sparsify: nested-support chain, pick min marginal -2 log L. -----------
##     The marginal integrates the deviations out, so a spurious individual
##     parameter LOWERS the data -2logL but RAISES the marginal (Occam factor);
##     the minimum of the chain is the parsimonious model.
sel <- sparsify(obj, center, fixed = fixed, method = "focei",
                 control = list(cm1 = list(iterlim = 50L)),
                 fits = 12, cores = 4, sd = 0.5, verbose = TRUE)
print(sel)
stopifnot(setequal(sel$support, c("kon", "ke")))     # the parsimonious model

## 6. Baseline: the classical L1 scan (Hauber et al. 2023) with scanL1().
##    A penalised multistart fit per lambda, an unpenalised refit per distinct
##    support, the largest lambda the likelihood ratio test does not reject.
eta_names <- as.vector(pen$subjectEtas)              # line-B deviations
scan <- scanL1(normL2(dlist, prd, errmodel = e),
               c(st[paste0("log_", cand)], log_sigma = log(0.1),
                 setNames(rep(0, length(eta_names)), eta_names)),
               reference = eta_names, fixed = fixed,
               lambda = 10^seq(-1, 4, length.out = 16), fits = 10, cores = 4,
               sd = 0.3)
print(scan)
print(plot(scan))
indiv <- setdiff(eta_names, scan$structure[[scan$selected]]$removed)
cat(sprintf("\nsparsify  S = {%s}\nscanL1    S = {%s}\n",
            paste(sort(sel$support), collapse = ","),
            paste(sort(gsub("^eta_|_B$", "", indiv)), collapse = ",")))
}
