\dontrun{

## ---------------------------------------------------------------------------
## Nonlinear mixed-effects walkthrough: Theophylline (12 subjects, oral PK).
##
## Structural model, on the log scale so positivity is automatic:
##   Ka_i = exp(tka + eta_Ka_i), V_i = exp(tv + eta_V_i), Cl_i = exp(tcl + eta_Cl_i)
## Residual model, the usual lognormal / constant-CV PK convention:
##   log(DV_ij) = log(Cc_ij) + eps_ij,   eps_ij ~ N(0, sigma_add^2)
##
## Reference values (Pinheiro & Bates 2000): Ka_pop ~ 1.5 1/h, V_pop ~ 32 L,
## Cl_pop ~ 3.0 L/h; sd(eta) ~ 0.5 / 0.1 / 0.3.
## ---------------------------------------------------------------------------

library(dMod2)

## Keep generated C/C++ sources and shared objects out of the project tree.
## Every constructor below takes an `outdir`; point it at a temporary folder.
outdir <- tempdir()

## 1. Data. Pre-dose rows are dropped: the model predicts Cc(0) = 0 by
## construction, and log(0) is not a usable observation.
data(Theoph, package = "datasets")
Theoph$Subject <- as.character(Theoph$Subject)
Theoph_pos     <- Theoph[Theoph$Time > 0, , drop = FALSE]
subjects       <- sort(unique(Theoph_pos$Subject))

doses <- vapply(subjects, function(s) {
  rec <- Theoph[Theoph$Subject == s, ][1, ]
  rec$Dose * rec$Wt
}, 0.0)

## sigma = NA leaves the residual SD to be estimated jointly with the
## structural parameters and the omegas (the NONMEM / nlmixr2 convention).
dlist <- as.datalist(data.frame(
  name      = "y",
  time      = Theoph_pos$Time,
  value     = log(Theoph_pos$conc),
  sigma     = NA_real_,
  condition = Theoph_pos$Subject,
  stringsAsFactors = FALSE))

## 2. Prediction chain: one-compartment oral PK.
reactions <- eqnlist()
reactions <- addReaction(reactions, "Ag", "",   "Ka * Ag",     "absorption")
reactions <- addReaction(reactions, "",   "Cc", "Ka * Ag / V", "appearance")
reactions <- addReaction(reactions, "Cc", "",   "Cl/V * Cc",   "elimination")

## deriv2 = TRUE is required throughout: the NLME kernels need second-order
## sensitivities. The pieces are built with compile = FALSE and compiled
## together below, once the whole chain exists.
m <- odemodel(reactions, modelname = "theoph_ode", backend = "cppDE",
              compile = FALSE, deriv2 = TRUE, outdir = outdir)
x <- Xs(m)
g <- Y(c(y = "log(Cc + 1e-9)"), x, modelname = "theoph_obs",
       compile = FALSE, deriv2 = TRUE, outdir = outdir)
e <- Y(eqnvec(y = "sigma_add"), g, attach.input = FALSE,
       modelname = "theoph_err", compile = FALSE, deriv2 = TRUE, outdir = outdir)

## 3. Per-subject parameter transformation. branch(apply = "insert") expands
## the shared trafo into one condition per subject, routing each subject's
## own eta and dose into the shared structural form.
trafo <- eqnvec(Ka        = "exp(tka + eta_Ka)",
                V         = "exp(tv  + eta_V)",
                Cl        = "exp(tcl + eta_Cl)",
                Ag        = "Ag_init",
                Cc        = "0",
                sigma_add = "exp(log_sigma_add)")
subj_table <- data.frame(
  eta_Ka  = paste0("eta_Ka_", subjects),
  eta_V   = paste0("eta_V_",  subjects),
  eta_Cl  = paste0("eta_Cl_", subjects),
  Ag_init = doses,
  row.names = subjects, stringsAsFactors = FALSE)
p   <- P(branch(trafo, table = subj_table, apply = "insert"),
         method = "explicit", modelname = "theoph_p",
         compile = FALSE, deriv2 = TRUE, outdir = outdir)
prd <- g * x * p

## Compile and load the generated C/C++ for the whole chain in one pass.
compile(prd, e, cores = 4)

## 4. Marginal likelihood. normL2() and constraintL2() stamp the model pieces
## onto the objective, so EM() recovers prd, dlist, e and om from `obj`
## and they are never passed again.
om  <- omega(eta = c("eta_Ka", "eta_V", "eta_Cl"), subjects = subjects)
obj <- normL2(dlist, prd, errmodel = e) +
         constraintL2(mu = 0, Omega = om)

## emInit() fills in the omega Cholesky entries (here sd = 0.3 on each
## diagonal) around the named structural starting values.
init <- emInit(c(tka = 0, tv = 3, tcl = 0.7, log_sigma_add = log(0.2)),
                 om, sd = 0.3)

## 5. Fit. FOCEI is the default and the right starting point for exploratory
## work: Laplace E-step plus an exact d log|H_i| / d theta correction.
fit <- EM(obj, init, method = "focei")

## Population estimates with SE and RSE%, Omega SDs and correlations, and eta
## shrinkage. coef(), vcov() and confint() (Wald) are also available.
summary(fit)

## fit$value is the plain-ML -2 log L (== nlmixr2's -2LL, == NONMEM's
## "OFV with constant"). For a raw NONMEM .lst OBJ compare fit$value_nonmem.
fit$value
fit$value_nonmem

## 6. Diagnostics. All return ggplot objects.
plot          (fit)   # DV vs IPRED and DV vs PRED
plotIndivs    (fit)   # per-subject IPRED / PRED curves over the data
plotResiduals (fit)   # IWRES vs IPRED and IWRES vs TIME
plotHistIndivs(fit)   # eta histograms with QQ line vs N(0, Omega_kk)

## ---------------------------------------------------------------------------
## Other estimators. All share this objective, this init, and the fit layout;
## only the E-step and the outer loop differ.
## ---------------------------------------------------------------------------

## Sparse-grid Gauss-Hermite refinement, warm-started from FOCEI. Use it as a
## confirmatory step once the model structure is settled: it removes the
## Laplace bias in Omega, at a cost that grows with the number of etas.
fit_quad <- EM(obj, init, method = "foceiQuadrature")
plotTrace(fit_quad)   # per-stage ECM convergence

## Deterministic EM at the Laplace level: the closed-form Omega update keeps
## the H_i^-1 posterior-variance term that ITS drops.
fit_em <- EM(obj, init, method = "laplaceEM",
                  control = list(quadrature = list(maxEcm = 200L)))

## Stochastic approximation EM. MCMC E-step, so no Gaussian assumption on the
## posterior and a cost per iteration that does not grow with the number of
## etas. Preferred for high-dimensional or strongly nonlinear problems.
fit_saem <- EM(obj, init, method = "saem",
                    control = list(saem = list(nBurnin = 200L, nEM = 200L)))

## Multi-start from a perturbed center, to check for local optima.
fit_ms <- msEM(obj, center = init, method = "focei",
                    fits = 8, sd = 0.5, cores = 4)

}
