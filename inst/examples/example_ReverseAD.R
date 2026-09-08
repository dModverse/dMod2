# -------------------------------------------------------------------------#
# Reverse-mode derivatives through the whole chain
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# The forward mode threads a Jacobian: every object in a chain receives
# dP/dtheta on its input and returns dX/dtheta on its output, so everything in
# flight is n_theta wide and the ODE integrates 1 + n_theta copies of itself.
# The reverse mode threads a cotangent the other way and is as wide as the
# number of seeds, which for a gradient is one. Its cost does not grow with the
# number of parameters.
#
# This file walks the chain piece by piece on toy models small enough to check
# by hand, and each section states what its own oracle is. Every section stands
# alone.
#
# [WHAT IT TAKES]
# One extra compilation on the ODE, `odemodel(..., reverse = TRUE)`, and
# `compile = TRUE` on the observation and transformation functions, because the
# reverse path has no interpreted fallback. Then:
#
#     obj(pars, sweep = "reverse")
#
# [WHY THE TWO DO NOT AGREE TO MACHINE PRECISION]
# They differ by O(tol), and not because of the adjoint. A forward-sensitivity
# solve carries n_theta tangent columns, and the step-size controller's error
# norm takes the maximum over the state norm AND every one of those columns. A
# maximum over a larger set is larger, so that run takes a finer step sequence
# than a value-only run does. The two modes therefore differentiate two
# different discretisations, each exactly. The gap falls with the tolerance,
# which section 5 measures.
#
# The corollary is in reverse's favour: its gradient belongs to the trajectory a
# value-only prediction produces, so value and gradient are consistent with each
# other. Under forward sensitivities they are not.
#
# [AUTHOR]
# Simon Beyer
#
# [Date]
# Mon 08 Sep 2026
# -------------------------------------------------------------------------#

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


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 1. The solver alone
#
# Xs() backwards is one value solve with a checkpoint per step, then one
# backward sweep. What goes in is a cotangent on the prediction; what comes out
# is a cotangent on the model's own parameters, which is the row set a forward
# sens1ini seeds.
#
# The oracle is the forward sensitivity contracted with the same cotangent: the
# two compute w' S from opposite ends.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
model <- odemodel(reactions, modelname = "rev_ode", deriv = TRUE,
                  reverse = TRUE, outdir = .outdir, compile = TRUE)

x     <- Xs(model, optionsOde = tight, optionsSens = tight)
inner <- c(A = 2, B = 0, k1 = 0.6, k2 = 0.3)
times <- seq(0, 8, length.out = 41)

fwd <- x(times, inner)[[1]]
set.seed(1)
w <- matrix(rnorm(nrow(fwd) * 2), nrow(fwd), 2, dimnames = list(NULL, c("A", "B")))

ref <- apply(attr(fwd, "deriv") * as.vector(w), 3, sum)
got <- attr(x, "mappings")[[1]] |> attr("vjpfn") |> (\(f) f(times, inner, NULL, w))()

print(rbind(forward = ref, reverse = got[names(ref)]))
cat("1. solver          max |difference| =",
    format(max(abs(ref - got[names(ref)])), digits = 3), "\n\n")


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 2. The whole chain, g * x * p
#
# normL2 turns the residuals into a seed -- not the residual vector itself, see
# section 4 -- and pushes it back through the observation function, the solver
# and the transformation in one pass. `sweep = "reverse"` is the only thing the
# caller says.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
g <- Y(c(obsA = "s*A", obsB = "s*B"), reactions, compile = TRUE,
       modelname = "rev_obs", outdir = .outdir)
p <- P(c(A = "exp(logA)", B = "0", k1 = "exp(logk1)", k2 = "exp(logk2)",
         s = "exp(logs)"),
       condition = "C1", compile = TRUE, modelname = "rev_p", outdir = .outdir)

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
    "\n   the reverse objective carries no Hessian:", is.null(r$hessian), "\n\n")


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 3. Several conditions, and a parameter that only one of them has
#
# Conditions are independent, so the backward pass walks each branch of the `+`
# and the cotangents meet again on the shared parameters. dk_C1 and dk_C2 reach
# one condition each and must come back with a contribution from that one only.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
conds <- c("C1", "C2")
trafo2 <- c(A = "exp(logA)", B = "0", k1 = "exp(logk1)", k2 = "exp(logk2)",
            s = "exp(logs)")
p2 <- Reduce("+", lapply(conds, function(cn)
  P(repar(paste0("logk1 ~ logk1 + dk_", cn), trafo2), condition = cn,
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


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 4. An estimated error model
#
# This is where the seed stops being the residual vector. sigma depends on theta
# as well, so a data row seeds two things: the prediction, and the error model's
# own output. The second one then travels back through the error model, which
# reads the prediction, so it lands on the prediction a second time as well as
# on sd_rel and sd_abs.
#
# Leaving that second channel out would still produce a gradient, and it would
# be wrong only in the sd_* directions -- which is exactly the kind of error a
# finite-difference check at 1e-3 does not catch.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
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


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 5. Why the two do not agree to machine precision
#
# The gap is the difference between two discretisations, so it falls with the
# solver tolerance and not with anything about the adjoint. A missing channel in
# the reverse pass would leave a floor here instead.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
gap <- t(vapply(10^-c(4, 6, 8, 10, 12), function(tt) {
  o  <- list(atol = tt, rtol = tt)
  xx <- Xs(model, optionsOde = o, optionsSens = o)
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


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 6. Events, and why they need no second derivation
#
# A jump is replayed backwards through the very function that applied it
# forwards -- the saltation transport is written over the scalar type, so the
# reverse type takes the same branch the forward one does. What the reverse pass
# reads back is which events fired and where a root sat, because those are
# control decisions and not arithmetic.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
ev <- eventlist(var = "A", time = "t_dose", value = "d_amt", method = "add")
model_ev <- odemodel(reactions, events = ev, modelname = "rev_ev",
                     deriv = TRUE, reverse = TRUE, outdir = .outdir,
                     compile = TRUE)
x_ev <- Xs(model_ev, optionsOde = tight, optionsSens = tight)
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
    "\n   the dose itself is estimated, so the jump has to carry a derivative.\n\n")


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 7. What it costs
#
# The forward gradient integrates 1 + n_theta copies of the system; the reverse
# one integrates the states twice -- once for the values, once inside the sweep
# -- and sweeps a tape one step wide. So the forward line rises with n_theta and
# the reverse line does not, and where they cross is a property of the model.
#
# Four parameters is far too few for reverse to win. The point of the number
# below is the shape, not the winner.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
reps <- 20
tf <- system.time(for (i in seq_len(reps)) obj(pars, deriv = TRUE))[["elapsed"]]
tr <- system.time(for (i in seq_len(reps))
                    obj(pars, deriv = TRUE, sweep = "reverse"))[["elapsed"]]
cat(sprintf("7. %d gradients    forward %.3fs   reverse %.3fs   ratio %.2f\n",
            reps, tf, tr, tr / tf))
cat("   with", length(pars), "parameters. The ratio is what falls as they grow.\n")
