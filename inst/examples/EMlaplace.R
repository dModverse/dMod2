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
       attach.input = FALSE, outdir = outdir)
e <- Y(eqnvec(y_ext = "sigma", y_mem = "sigma", y_int = "sigma"), g,
       modelname = "becker_err", compile = FALSE, deriv2 = TRUE,
       attach.input = FALSE, outdir = outdir)

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

## 6. Baseline: the classical trustL1 lambda-scan + BIC (what sparsify saves).
##    20 lambdas x 50 starts = 1000 fits, versus sparsify's ~K+1 marginal fits.
obj_data  <- normL2(dlist, prd, errmodel = e)
eta_names <- as.vector(pen$subjectEtas)              # line-B deviations
mu_pen    <- setNames(rep(0, length(eta_names)), eta_names)
N         <- nrow(dl_rows)
p_struct  <- length(cand) + 1L                       # candidate mu's + sigma
lam_grid  <- 10^seq(-5, 5, length.out = 20)
scan_center <- c(st[paste0("log_", cand)], log_sigma = log(0.1),
                 setNames(rep(0, length(eta_names)), eta_names))

## 6a. regularisation path: which deviations survive at each lambda. ----------
path <- do.call(rbind, lapply(lam_grid, function(lam) {
  ms <- mstrust(obj_data, center = scan_center, studyname = "scan",
                optmethod = "trustL1", lambda = lam, mu = mu_pen, fixed = fixed,
                fits = 50, cores = 4, sd = 0.3, rinit = 1, rmax = 10,
                iterlim = 400, resultPath = tempdir(), output = FALSE)
  keep <- abs(as.parvec(as.parframe(ms), index = 1)[eta_names]) > 1e-4
  data.frame(lambda = lam, nnz = sum(keep),
             support = paste(sort(gsub("^eta_|_B$", "", eta_names[keep])),
                             collapse = ","), stringsAsFactors = FALSE)
}))
print(path, row.names = FALSE)

## 6b. BIC over the DISTINCT supports the path visits. Compute BIC from an
##     UNPENALISED (relaxed) refit of each support: at a shrinking lambda the
##     deviations are biased toward 0, so scoring supports at the penalised fit
##     is apples-to-oranges across lambda; the post-selection refit is the
##     standard way to BIC-score a lasso support (relaxed lasso, Meinshausen 2007).
relaxRefit <- function(S) {                          # clean data -2 log L for S
  eta_in  <- if (length(S)) as.vector(pen$subjectEtas[, paste0("eta_", S),
                                                      drop = FALSE]) else character(0)
  eta_out <- setdiff(eta_names, eta_in)
  fixedS  <- c(fixed, setNames(rep(0, length(eta_out)), eta_out))
  start   <- c(st[paste0("log_", cand)], log_sigma = log(0.1),
               setNames(rep(0, length(eta_in)), eta_in))
  ms <- mstrust(obj_data, center = start, studyname = "relax", optmethod = "trust",
                fixed = fixedS, fits = 12, cores = 4, sd = 0.3, rinit = 1,
                rmax = 10, iterlim = 400, resultPath = tempdir(), output = FALSE)
  obj_data(as.parvec(as.parframe(ms), index = 1), fixed = fixedS,
           deriv = FALSE)$value
}
supports <- unique(lapply(strsplit(path$support, ","), function(s) s[nzchar(s)]))
bic_vals <- vapply(supports,
                   function(S) relaxRefit(S) + log(N) * (p_struct + length(S)), 0.0)
bic_tab  <- data.frame(
  support = vapply(supports, function(S)
    if (length(S)) paste(sort(S), collapse = ",") else "(none)", ""),
  size = lengths(supports), BIC = bic_vals, stringsAsFactors = FALSE)
bic_tab  <- bic_tab[order(bic_vals), ]
bic_tab$sel <- ifelse(bic_tab$BIC == min(bic_tab$BIC), "<= BIC min", "")
print(bic_tab, row.names = FALSE)

## sparsify finds the SAME parsimonious model the BIC-scan selects. ----------
bic_support <- supports[[which.min(bic_vals)]]
cat(sprintf(paste0("\nregSelect  S_hat = {%s}\nBIC-scan   S     = {%s}\n",
                   "(sparsify: ~%d marginal fits;  scan: %d x %d = %d fits)\n"),
            paste(sort(sel$support), collapse = ","),
            paste(sort(bic_support), collapse = ","),
            length(cand) + 1L, length(lam_grid), 50L, length(lam_grid) * 50L))
}
