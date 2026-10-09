# Reverse-mode derivatives through the whole chain, checked on toy models.
# Needs `derivMode = c("forward", "reverse")` on the odemodel, `compile = TRUE` on
# g and p, and then `obj(pars, sweep = "reverse")`.

library(dMod2)

.outdir <- file.path(tempdir(), "reverseAD")
if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)
setwd(.outdir)

# Tight, so the discretisation gap between the two modes stays below the
# comparisons below. Section 5 is the one that varies it on purpose.
tight <- list(atol = 1e-11, rtol = 1e-11)

reactions <- eqnlist() |>
  addReaction("A", "B", "k1*A", "conversion") |>
  addReaction("B", "", "k2*B", "decay")


# 1. The solver alone ---------------------------------------------------------
# One value solve with checkpoints, then one backward sweep from a cotangent on
# the prediction to one on the model parameters; the oracle is w' S from forward.
model <- odemodel(reactions, modelname = "rev_ode", deriv = TRUE, deriv2 = TRUE,
                  derivMode = c("forward", "reverse", "forward-forward",
                                "forward-reverse"),
                  outdir = .outdir, compile = TRUE)

x     <- Xs(model, options = tight)
inner <- c(A = 2, B = 0, k1 = 0.6, k2 = 0.3)
times <- seq(0, 8, length.out = 41)

fwd <- x(times, inner)[[1]]
set.seed(1)
w <- matrix(rnorm(nrow(fwd) * 2), nrow(fwd), 2, dimnames = list(NULL, c("A", "B")))

ref <- apply(attr(fwd, "deriv") * as.vector(w), 3, sum)
# A vjp answers a cotangent with a direction axis: column 1 is the gradient,
# and second order would put its derivatives beside it.
got <- attr(x, "mappings")[[1]] |> attr("vjpfn") |> (\(f) f(times, inner, NULL, w))()
got <- got[, 1L]

print(rbind(forward = ref, reverse = got[names(ref)]))
cat("1. solver          max |difference| =",
    format(max(abs(ref - got[names(ref)])), digits = 3), "\n\n")


# 2. The whole chain, g * x * p -----------------------------------------------
# normL2 seeds the residual cotangent and pushes it back through g, x and p in
# one pass.
g <- Y(c(obsA = "s*A", obsB = "s*B"), reactions, compile = TRUE, deriv2 = TRUE,
       modelname = "rev_obs", outdir = .outdir)
p <- P(c(A = "exp(logA)", B = "0", k1 = "exp(logk1)", k2 = "exp(logk2)",
         s = "exp(logs)"),
       condition = "C1", compile = TRUE, deriv2 = TRUE,
       modelname = "rev_p", outdir = .outdir)

pars <- c(logA = log(2), logk1 = log(0.6), logk2 = log(0.3), logs = log(1.5))
truth <- (g * x * p)(times, pars)[[1]]

set.seed(4)
idx <- c(4L, 9L, 16L, 26L, 36L)
dat <- do.call(rbind, lapply(c("obsA", "obsB"), function(nm)
  data.frame(name = nm, time = truth[idx, "time"],
             value = truth[idx, nm] * exp(rnorm(length(idx), 0, 0.05)),
             sigma = 0.1, condition = "C1", stringsAsFactors = FALSE)))
obj <- normL2(as.datalist(dat), g * x * p)

f <- obj(pars, deriv = TRUE)
r <- obj(pars, deriv = TRUE, sweep = "reverse")

print(rbind(forward = f$gradient, reverse = r$gradient[names(f$gradient)]))
cat("2. g * x * p       max |difference| =",
    format(max(abs(f$gradient - r$gradient[names(f$gradient)])), digits = 3),
    "\n   the reverse objective has no Hessian:", is.null(r$hessian), "\n\n")


# 3. Several conditions, and a parameter only one of them has -----------------
# Cotangents of the `+` branches meet on shared parameters; a condition-specific
# parameter receives a contribution from its own condition only.
conds <- c("C1", "C2")
trafo2 <- c(A = "exp(logA)", B = "0", k1 = "exp(logk1)", k2 = "exp(logk2)",
            s = "exp(logs)")
p2 <- Reduce("+", lapply(conds, function(cn)
  P(repar(trafo2, paste0("logk1 ~ logk1 + dk_", cn)), condition = cn,
    compile = TRUE, modelname = paste0("rev_p2_", cn), outdir = .outdir)))

pars2 <- c(pars, dk_C1 = 0.1, dk_C2 = -0.15)
pred2 <- (g * x * p2)(times, pars2)
set.seed(7)
dat2 <- do.call(rbind, lapply(conds, function(cn)
  do.call(rbind, lapply(c("obsA", "obsB"), function(nm)
    data.frame(name = nm, time = pred2[[cn]][idx, "time"],
               value = pred2[[cn]][idx, nm] * exp(rnorm(length(idx), 0, 0.05)),
               sigma = 0.1, condition = cn, stringsAsFactors = FALSE)))))
obj2 <- normL2(as.datalist(dat2), g * x * p2)

f2 <- obj2(pars2, deriv = TRUE)
r2 <- obj2(pars2, deriv = TRUE, sweep = "reverse")
print(rbind(forward = f2$gradient, reverse = r2$gradient[names(f2$gradient)]))
cat("3. two conditions  max |difference| =",
    format(max(abs(f2$gradient - r2$gradient[names(f2$gradient)])), digits = 3),
    "\n\n")


# 4. An estimated error model -------------------------------------------------
# A sigma that depends on theta seeds the error model's output as well, which
# reaches the prediction a second time and the error parameters.
e <- Y(c(obsA = "sd_rel*obsA + sd_abs", obsB = "sd_rel*obsB + sd_abs"), g,
       states = c("obsA", "obsB"), parameters = c("sd_rel", "sd_abs"),
       compile = TRUE, modelname = "rev_err", outdir = .outdir)

pe <- P(c(A = "exp(logA)", B = "0", k1 = "exp(logk1)", k2 = "exp(logk2)",
          s = "exp(logs)", sd_rel = "exp(logsdrel)", sd_abs = "exp(logsdabs)"),
        condition = "C1", compile = TRUE, modelname = "rev_pe", outdir = .outdir)

pars_e <- c(pars, logsdrel = log(0.08), logsdabs = log(0.02))
# The same data, but with sigma unknown so the error model has to supply it.
dat_e <- dat; dat_e$sigma <- NULL
obj_e <- normL2(as.datalist(dat_e), g * x * pe, e)

fe <- obj_e(pars_e, deriv = TRUE)
re <- obj_e(pars_e, deriv = TRUE, sweep = "reverse")
print(rbind(forward = fe$gradient, reverse = re$gradient[names(fe$gradient)]))
cat("4. error model     max |difference| =",
    format(max(abs(fe$gradient - re$gradient[names(fe$gradient)])), digits = 3),
    "\n\n")


# 5. Why the two do not agree to machine precision ----------------------------
# Each mode differentiates its own discretisation exactly, so the gap falls with
# the tolerance; a missing channel in the reverse pass would leave a floor.
gap <- t(vapply(10^-c(4, 6, 8, 10, 12), function(tt) {
  o  <- list(atol = tt, rtol = tt)
  xx <- Xs(model, options = o)
  oo <- normL2(as.datalist(dat), g * xx * p)
  a  <- oo(pars, deriv = TRUE)
  b  <- oo(pars, deriv = TRUE, sweep = "reverse")
  c(tol      = tt,
    value    = abs(a$value - b$value) / abs(a$value),
    gradient = max(abs(a$gradient - b$gradient[names(a$gradient)])) /
               max(abs(a$gradient)))
}, numeric(3)))
print(format(as.data.frame(gap), digits = 3), row.names = FALSE)
cat("\n5. the gap tracks the tolerance over eight decades, so it is the\n",
    "  discretisation and not the adjoint.\n\n")


# 5b. A log observable with an epsilon below the solver's noise floor amplifies
# integration noise in both modes; a forward/reverse mismatch exposes it. The fix
# is in the model: an epsilon above the noise floor, a tolerance below it.


# 6. Events -------------------------------------------------------------------
# A jump is replayed backwards through the same scalar-generic function that
# applied it; the reverse pass only reads back which events fired and when.
ev <- eventlist(var = "A", time = "t_dose", value = "d_amt", method = "add")
model_ev <- odemodel(reactions, events = ev, modelname = "rev_ev",
                     deriv = TRUE, derivMode = c("forward", "reverse"), outdir = .outdir,
                     compile = TRUE)
x_ev <- Xs(model_ev, options = tight)
p_ev <- P(c(A = "exp(logA)", B = "0", k1 = "exp(logk1)", k2 = "exp(logk2)",
            s = "exp(logs)", t_dose = "3", d_amt = "exp(logdose)"),
          condition = "C1", compile = TRUE, modelname = "rev_pev",
          outdir = .outdir)

pars_ev <- c(pars, logdose = log(0.8))
obj_ev  <- normL2(as.datalist(dat), g * x_ev * p_ev)
fv <- obj_ev(pars_ev, deriv = TRUE)
rv <- obj_ev(pars_ev, deriv = TRUE, sweep = "reverse")
print(rbind(forward = fv$gradient, reverse = rv$gradient[names(fv$gradient)]))
cat("6. dosing event    max |difference| =",
    format(max(abs(fv$gradient - rv$gradient[names(fv$gradient)])), digits = 3),
    "\n   the dose itself is estimated, so the jump needs a derivative.\n\n")


# 7. What it costs ------------------------------------------------------------
# Forward cost grows with n_theta, reverse cost does not; four parameters are
# too few for reverse to win, the ratio is what falls as they grow.
reps <- 20
tf <- system.time(for (i in seq_len(reps)) obj(pars, deriv = TRUE))[["elapsed"]]
tr <- system.time(for (i in seq_len(reps))
                    obj(pars, deriv = TRUE, sweep = "reverse"))[["elapsed"]]
cat(sprintf("7. %d gradients    forward %.3fs   reverse %.3fs   ratio %.2f\n",
            reps, tf, tr, tr / tf))
cat("   with", length(pars), "parameters. The ratio is what falls as they grow.\n")

# 8. The exact Hessian --------------------------------------------------------
# `sweep` picks forward over forward or forward over reverse; the model needs the
# matching `derivMode` entries, and g and p need `deriv2 = TRUE`.
hf <- obj(pars, deriv2 = TRUE)
hr <- obj(pars, sweep = "reverse", deriv2 = TRUE)
cat(sprintf("8. exact Hessian   forward vs reverse: max|difference| %.2e\n",
            max(abs(hf$hessian - hr$hessian))))
cat("   attr(., \"sweep\") says which answered:",
    attr(hr, "sweep"), "\n")

# `hessian` asks for a Hessian at all, `deriv2` for the exact one rather than
# Gauss-Newton; a contradiction resolves to the cheaper answer with a warning.
gn <- obj(pars, hessian = TRUE)
cat(sprintf("   Gauss-Newton vs exact: max|difference| %.2e  (they are not the ",
            max(abs(gn$hessian - hf$hessian))))
cat("same matrix)\n")

# 9. An optimiser that chooses per evaluation ---------------------------------
# One exact Hessian at the start, then quasi-Newton on reverse gradients. It pays
# only past the crossover of section 7, at tens of parameters.
start <- pars + 0.3
fit_gn <- trust(obj, start, iterlim = 100L)
fit_ex <- trust(obj, start, iterlim = 100L, hessianMethod = "sr1",
                sweep = "reverse",
                qnControl = list(hessianInit = "exact",
                                 hessianReseed = "stall"))
cat(sprintf("9. gn: value %.6f in %d iterations\n", fit_gn$value,
            fit_gn$iterations))
cat(sprintf("   sr1 from an exact seed, reverse: value %.6f in %d iterations\n",
            fit_ex$value, fit_ex$iterations))
