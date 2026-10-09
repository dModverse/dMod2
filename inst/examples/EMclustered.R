\dontrun{

## Which of six cell lines share a parameter value (Becker et al. 2010, Epo
## receptor)? The clustered penalty fuses all pairwise differences, and sparsify
## selects each parameter's grouping by exact-marginal model selection.

library(dMod2)
outdir <- tempdir()          # keep generated C/C++ + shared objects out of the tree
set.seed(2)

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

m <- odemodel(reactions, modelname = "beckerc", backend = "cppDE",
              compile = FALSE, deriv2 = TRUE, outdir = outdir)
x <- Xs(m, compile = FALSE)
g <- Y(eqnvec(y_ext = "log(Epo + dEpo_e + 1)",
              y_mem = "log(Epo_EpoR + 1)",
              y_int = "log(Epo_EpoR_i + dEpo_i + 1)"),
       x, modelname = "beckerc_obs", compile = FALSE, deriv2 = TRUE,
       attachInput = FALSE, outdir = outdir)
e <- Y(eqnvec(y_ext = "sigma", y_mem = "sigma", y_int = "sigma"), g,
       modelname = "beckerc_err", compile = FALSE, deriv2 = TRUE,
       attachInput = FALSE, outdir = outdir)

## 2. Clustered encoding: every line has a deviation for each candidate.
##    Unlike the 2-line reference encoding (EMlaplace.R), the complete-graph fusion
##    needs an eta for ALL subjects; the mean level is anchored into log_<par>.
dose  <- 1347.49                                       # init_Epo
lines <- c("A", "B", "C", "D", "E", "F")
cand  <- c("kon", "ke", "kt")
trafo <- eqnvec(
  kon  = "exp(log_kon + eta_kon)", ke = "exp(log_ke + eta_ke)",
  kt   = "exp(log_kt  + eta_kt)",  koff = "exp(log_koff)",
  kex  = "exp(log_kex)", kdi = "exp(log_kdi)", kde = "exp(log_kde)",
  init_Epo = "dose", init_EpoR_rel = "exp(log_EpoR_rel)",
  Epo = "dose", EpoR = "4 * dose * exp(log_EpoR_rel)",
  Epo_EpoR = "0", Epo_EpoR_i = "0", dEpo_i = "0", dEpo_e = "0",
  sigma = "exp(log_sigma)")
tab <- data.frame(
  eta_kon = paste0("eta_kon_", lines), eta_ke = paste0("eta_ke_", lines),
  eta_kt  = paste0("eta_kt_",  lines), dose = dose,
  row.names = lines, stringsAsFactors = FALSE)
p <- P(branch(trafo, table = tab, apply = "insert"), method = "explicit",
       modelname = "beckerc_p", compile = FALSE, deriv2 = TRUE, outdir = outdir)
prd <- g * x * p
compile(prd, e, cores = 4)

## 3. Simulate six cell lines with the three ground-truth groupings. -----------
st <- c(log_kon = log(0.1513), log_koff = log(0.08069), log_kt = log(0.01601),
        log_ke = log(0.05546), log_kex = log(0.000578), log_kdi = log(0.001273),
        log_kde = log(0.01199), log_EpoR_rel = log(0.09205), log_sigma = log(0.05))
truth <- rbind(                                        # anchored: rows sum to 0
  eta_kon = c(A = 0,    B = 0,    C = 0,    D = 0,    E = 0,    F = 0),
  eta_ke  = c(A = 0.8,  B = 0.8,  C = 0.8,  D = -0.8, E = -0.8, F = -0.8),
  eta_kt  = c(A = 0.7,  B = 0.7,  C = 0,    D = 0,    E = -0.7, F = -0.7))
eta_true <- unlist(lapply(cand, function(pp)
  stats::setNames(truth[paste0("eta_", pp), ], paste0("eta_", pp, "_", lines))))
times <- c(3, 6, 12, 20, 30, 45, 60, 90, 120, 180, 240, 300)
sim <- prd(times, c(st, eta_true), deriv = FALSE)
dl_rows <- do.call(rbind, lapply(lines, function(s) {
  pr <- sim[[s]]
  do.call(rbind, lapply(c("y_ext", "y_mem", "y_int"), function(nm)
    data.frame(name = nm, time = pr[, "time"],
               value = pr[, nm] + rnorm(nrow(pr), 0, 0.03),
               sigma = NA_real_, condition = s, stringsAsFactors = FALSE)))
}))
dlist <- as.datalist(dl_rows)

## 4. Clustered penalty + objective + start (emInit adds the single lambda).
pen   <- penaltyL1(cand, method = "clustered", subjects = lines)
obj   <- normL2(dlist, prd, errmodel = e) + constraintL1(pen)
fixed <- st[c("log_koff", "log_kex", "log_kdi", "log_kde", "log_EpoR_rel")]
init  <- emInit(c(st[paste0("log_", cand)], log_sigma = log(0.03)), pen,
                 lambda = 5)

## 5. Recover the groupings. One all-separate fit gives the FOCE linearisation;
##    the clustered marginal scores the cluster path and the joint selector fixes
##    each parameter's grouping from the full cross-parameter covariance.
sel <- sparsify(obj, init, fixed = fixed,
                        fits = 8, cores = 4, sd = 0.3, verbose = TRUE)
print(sel)

## 6. Diagnostic plots: the recovered grouping and the score chain. ------------
print(plot(sel, type = "grouping"))   # eta_hat per line, clusters coloured
print(plot(sel, type = "chain"))      # marginal -2 log L vs G, selected marked

## 7. Recovery check against the ground truth. --------------------------------
truth_grp <- list(kon = list(lines),
                  ke  = list(c("A", "B", "C"), c("D", "E", "F")),
                  kt  = list(c("A", "B"), c("C", "D"), c("E", "F")))
norm_grp <- function(P) sort(vapply(lapply(P, sort), paste, "", collapse = ","))
cat("\n== recovered vs truth ==\n")
for (pp in cand) {
  rec <- sel$perParam[[paste0("eta_", pp)]]$clusters
  cat(sprintf("  %-4s recovered %-24s truth %-24s MATCH = %s\n", pp,
              paste(vapply(rec, paste, "", collapse = ","), collapse = " | "),
              paste(vapply(truth_grp[[pp]], paste, "", collapse = ","), collapse = " | "),
              identical(norm_grp(rec), norm_grp(truth_grp[[pp]]))))
}
}
