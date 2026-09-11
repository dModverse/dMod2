# -------------------------------------------------------------------------#
# The discrete adjoint against CVODES ASA on Bachmann
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# Both reverse routes answer the same interface, `obj(pars, sweep = "reverse")`,
# and differ only in the ODE object underneath: our own discrete adjoint, which
# replays every accepted forward step, and CVODES adjoint sensitivity analysis,
# which integrates the adjoint as its own ODE over checkpointed states. This
# times them against each other on the whole chain up to normL2, over all
# conditions at once and then over one condition at a time.
#
# The two granularities answer different questions. The whole objective is what
# an optimiser pays per iteration. The per-condition table says where that time
# sits, and whether the whole is the sum of its parts: a ratio near one means
# the chain is the sum of its conditions, well under one means the whole shares
# work that the parts each repeat.
#
# The model is the hand-built one from
# inst/examples/example_BachmannMSB2011.R, sourced up to its objective rather
# than transcribed. Its observables floor log10 at 1e-10, which sits above
# what the solver resolves at these tolerances; the PEtab form of the same
# problem floors at 1e-15, and a floor below the tolerance multiplies whatever
# noise is left by several decades and makes every gradient comparison on it
# meaningless.
#
# Only the ODE object differs between the two chains. The observation, error
# and transformation functions are the example's own, so the likelihood, the
# per-condition data groups and the parameter scale are shared by construction.
#
# Every route is called once before it is timed. A first call builds dMod2's
# prepared-batch handle and fills the object caches, so timing it cold measures
# the marshalling rather than the solve.
#
# [WHAT IT TAKES]
# SUNDIALS, for the ASA column; without it that column is NA and the rest still
# runs. An idle machine: every number here is a wall time. Two full Bachmann
# compilations on top of the example's own.
#
# [AUTHOR]
# Simon Beyer
#
# [Date]
# Wed 10 Sep 2026
# -------------------------------------------------------------------------#

library(dMod2)
library(microbenchmark)

Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", GOTO_NUM_THREADS = "1")
options(dMod.cores = 1, cppDE.cores = 1)

# atol under rtol under the observables' own floor: the three have to hold in
# that order, or the gradient is measuring the floor.
TOL <- list(atol = 1e-11, rtol = 1e-9, maxsteps = 1e7L, maxattemps = 100L)

# The machine scatters badly, so report the minimum over repetitions rather
# than the mean or the median. One condition is two orders of magnitude cheaper
# than all of them together and affords more of them.
REPS_WHOLE <- 7L
REPS_COND  <- 15L

tmin <- function(f, reps) {
  f()
  min(microbenchmark(f(), times = reps, unit = "ms")$time) / 1e6
}


# -----------------------------------------------------------------------------
# The example's own chain, up to but not including its first objective call
#
# `reactions`, `observables`, `errorModels`, `trafo`, `mydataL` and `bestfit`
# come from there, along with the compiled g, e and p this benchmark reuses.
# -----------------------------------------------------------------------------
.example <- system.file("examples", "example_BachmannMSB2011.R", package = "dMod2")
.src <- readLines(.example)
.upto <- grep("^obj\\(bestfit", .src)[1] - 1L
eval(parse(text = .src[seq_len(.upto)]), envir = globalenv())

.bdir <- file.path(tempdir(), "bench_adjoint")
dir.create(.bdir, recursive = TRUE, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# One ODE object per backend, both directions in each
# -----------------------------------------------------------------------------
.cfg <- tryCatch(get("cvodeConfig", envir = asNamespace("cppDE")),
                 error = function(e) NULL)
hasASA <- isTRUE(.cfg$available) &&
  "derivMode" %in% names(formals(cppDE::cvode)) &&
  "reverse" %in% eval(formals(cppDE::cvode)$derivMode)

mC <- odemodel(reactions, modelname = "bench_cpp", backend = "cppDE",
               derivMode = c("forward", "reverse"), compile = FALSE,
               outdir = .bdir)
xC <- Xs(mC, optionsOde = TOL, optionsSens = TOL)

# The same model on the Rosenbrock stepper. Its adjoint is a different piece of
# arithmetic: six direct solves transposed, against a corrector's implicit
# function theorem, so the two are worth seeing side by side.
mR <- odemodel(reactions, modelname = "bench_rb4", backend = "cppDE",
               method = "rb4", derivMode = c("forward", "reverse"),
               compile = FALSE, outdir = .bdir)
xR <- Xs(mR, optionsOde = TOL, optionsSens = TOL)

xS <- NULL
if (hasASA) {
  mS <- odemodel(reactions, modelname = "bench_sun", backend = "Sundials",
                 derivMode = c("forward", "reverse"), compile = FALSE,
                 outdir = .bdir)
  xS <- Xs(mS, optionsOde = TOL, optionsSens = TOL)
  compile(xC, xR, xS, output = "bench_adjoint", cores = 12)
} else {
  cat("SUNDIALS absent or cvode() has no reverse direction: ASA column is NA.\n")
  compile(xC, xR, output = "bench_adjoint", cores = 12)
}

# normL2 alone, not the example's objective: its prior term is one scalar over
# the whole problem and would be counted once per condition below.
mkobj  <- function(dat) normL2(dat, g * xC * p, e)
mkobjR <- function(dat) normL2(dat, g * xR * p, e)
mkobjS <- function(dat) if (is.null(xS)) NULL else normL2(dat, g * xS * p, e)

pars <- bestfit


# -----------------------------------------------------------------------------
# All conditions at once, the whole chain up to normL2
#
# The value run is the unit the other three are quoted in: how many value runs
# a gradient costs is the number that decides which direction is worth taking.
# -----------------------------------------------------------------------------
objC_all <- mkobj(mydataL)
objR_all <- mkobjR(mydataL)
objS_all <- mkobjS(mydataL)

# A reverse evaluation returns no Hessian: the one invariant that says the
# direction arrived rather than being swallowed by a wrapper in the chain.
chk <- objC_all(pars, deriv = TRUE, sweep = "reverse")
if (!is.null(chk$hessian)) stop("sweep = \"reverse\" returned a Hessian")

# The two chains must be the same model before their times mean anything.
if (!is.null(objS_all)) {
  v1 <- objC_all(pars, deriv = FALSE)$value
  v2 <- objS_all(pars, deriv = FALSE)$value
  cat("\nlikelihood, cppDE against Sundials:", format(v1, digits = 12), "/",
      format(v2, digits = 12), " relative gap",
      format(abs(v2 - v1) / abs(v1), digits = 3), "\n")
}

w_val <- tmin(function() objC_all(pars, deriv = FALSE), REPS_WHOLE)
w_fwd <- tmin(function() objC_all(pars, deriv = TRUE, hessian = FALSE), REPS_WHOLE)
w_rev <- tmin(function() objC_all(pars, deriv = TRUE, sweep = "reverse"), REPS_WHOLE)
w_rrev <- tmin(function() objR_all(pars, deriv = TRUE, sweep = "reverse"), REPS_WHOLE)
w_rval <- tmin(function() objR_all(pars, deriv = FALSE), REPS_WHOLE)
w_asa <- if (is.null(objS_all)) NA_real_ else
  tmin(function() objS_all(pars, deriv = TRUE, sweep = "reverse"), REPS_WHOLE)

whole <- data.frame(
  route     = c("value", "forward", "reverse", "value rb4", "reverse rb4", "ASA"),
  ms        = c(w_val, w_fwd, w_rev, w_rval, w_rrev, w_asa),
  in_values = c(w_val, w_fwd, w_rev, w_rval, w_rrev, w_asa) / w_val,
  stringsAsFactors = FALSE)

cat("\nAll conditions, whole chain to normL2, one core\n\n")
print(format(whole, digits = 3), row.names = FALSE)


# -----------------------------------------------------------------------------
# Do the two adjoints agree?
#
# Speed without this says nothing. Both are held against the forward gradient,
# which is the third discretisation of the same derivative and the only oracle
# either of them has here.
# -----------------------------------------------------------------------------
g_fwd <- objC_all(pars, deriv = TRUE, hessian = FALSE)$gradient
g_rev <- objC_all(pars, deriv = TRUE, sweep = "reverse")$gradient
g_rrev <- objR_all(pars, deriv = TRUE, sweep = "reverse")$gradient
g_asa <- if (is.null(objS_all)) NULL else
  objS_all(pars, deriv = TRUE, sweep = "reverse")$gradient

# Against the largest component rather than each component's own size: a
# forward component near zero would otherwise set the whole scale.
relgap <- function(a, b) {
  nm <- names(a)
  max(abs(b[nm] - a[nm])) / max(abs(a[nm]))
}
# The angle, which is what a line search feels.
cosgap <- function(a, b) {
  nm <- names(a)
  1 - sum(a[nm] * b[nm]) / sqrt(sum(a[nm]^2) * sum(b[nm]^2))
}

agree <- data.frame(
  against_forward    = c("reverse", "reverse rb4", "ASA"),
  worst_over_largest = c(relgap(g_fwd, g_rev), relgap(g_fwd, g_rrev),
                         if (is.null(g_asa)) NA_real_ else relgap(g_fwd, g_asa)),
  one_minus_cos      = c(cosgap(g_fwd, g_rev), cosgap(g_fwd, g_rrev),
                         if (is.null(g_asa)) NA_real_ else cosgap(g_fwd, g_asa)),
  stringsAsFactors = FALSE)

cat("\nThe two adjoints against the forward gradient\n\n")
print(format(agree, digits = 3), row.names = FALSE)


# -----------------------------------------------------------------------------
# One condition at a time
#
# `steps` is the accepted step count of the reverse solve for that condition,
# read off the sweep's own grid. It is the thing our cost is proportional to and
# ASA's is not: we replay every accepted forward step, ASA integrates the
# adjoint on a grid of its own choosing.
# -----------------------------------------------------------------------------
conds <- names(mydataL)
rows <- lapply(conds, function(cn) {
  objC <- mkobj(mydataL[cn])
  objS <- mkobjS(mydataL[cn])

  t_val <- tmin(function() objC(pars, deriv = FALSE), REPS_COND)
  t_fwd <- tmin(function() objC(pars, deriv = TRUE, hessian = FALSE), REPS_COND)
  t_rev <- tmin(function() objC(pars, deriv = TRUE, sweep = "reverse"), REPS_COND)
  t_asa <- if (is.null(objS)) NA_real_ else
    tmin(function() objS(pars, deriv = TRUE, sweep = "reverse"), REPS_COND)

  steps <- tryCatch({
    inner <- p(pars)[[cn]]
    tt <- sort(unique(c(0, mydataL[[cn]]$time)))
    gr <- cppDE::solveODE(mC$reversed, tt, inner,
                          seed = array(1, c(length(tt),
                                            length(attr(mC$func, "variables")), 1L)),
                          adjointGrid = TRUE, abstol = TOL$atol, reltol = TOL$rtol,
                          maxsteps = TOL$maxsteps)
    length(gr$adjointGrid$h)
  }, error = function(e) NA_integer_)

  data.frame(condition = cn, n_times = nrow(mydataL[[cn]]), steps = steps,
             t_value = t_val, t_fwd = t_fwd, t_rev = t_rev, t_asa = t_asa,
             stringsAsFactors = FALSE)
})
per <- do.call(rbind, rows)
per$rev_over_asa <- per$t_rev / per$t_asa
per <- per[order(-per$t_rev), ]

cat("\nPer condition, one objective each, milliseconds\n\n")
print(format(per, digits = 3), row.names = FALSE)


# -----------------------------------------------------------------------------
# The parts against the whole
# -----------------------------------------------------------------------------
cmp <- data.frame(
  route        = c("value", "forward", "reverse", "ASA"),
  sum_of_parts = c(sum(per$t_value), sum(per$t_fwd), sum(per$t_rev), sum(per$t_asa)),
  whole        = c(w_val, w_fwd, w_rev, w_asa),
  stringsAsFactors = FALSE)
cmp$whole_over_parts <- cmp$whole / cmp$sum_of_parts

cat("\n\nThe parts against the whole, milliseconds\n\n")
print(format(cmp, digits = 3), row.names = FALSE)

cat("\n  conditions:", length(conds), "  n_theta:", length(pars),
    "  atol", TOL$atol, " rtol", TOL$rtol, "\n")
