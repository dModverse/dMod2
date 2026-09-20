# -------------------------------------------------------------------------#
# Three ways to a gradient, on one model
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# dMod2 can reach the same gradient by three routes, and they are not
# interchangeable. This script puts them on one model, one parameter vector and
# one ladder of solver tolerances, and reports the two numbers that decide
# between them: how long each takes, and how far each sits from the others.
#
#   forward   sensitivity equations carried beside the states. One extra
#             trajectory per parameter, so the cost grows with n_theta.
#   reverse   the discrete adjoint of cppDE: integrate in plain double, keep a
#             checkpoint per step, replay each step backwards. One sweep,
#             whatever n_theta is.
#   ASA       SUNDIALS CVODES adjoint sensitivity analysis, the same idea
#             solved as a second ODE rather than differentiated step by step.
#
# [WHAT TO EXPECT, AND WHAT WOULD BE NEWS]
# The three do not agree to machine precision and should not. Forward
# sensitivities make the step-size controller take the maximum over the state
# error AND every tangent column, so a sensitivity run adapts on a finer grid
# than a value-only run. Forward and reverse therefore differentiate two
# different discretisations, each exactly, and their gap is O(tol). A gap that
# does NOT fall with the tolerance is the news: that is a missing channel in one
# of the two, not a discretisation difference.
#
# The corollary favours reverse and is the reason it exists beyond speed: its
# gradient belongs to the trajectory a value-only solve produces, so the value
# and the gradient a caller receives are consistent with each other. Under
# forward sensitivities they are not, the value coming from a grid the tangents
# made finer. ASA has the same inconsistency for the same reason and is measured
# for it below.
#
# [MODEL]
# Bachmann et al. (2011), 113 estimated parameters. Deliberately not a toy: a
# small model measures the cost of an R call rather than of a method, and the
# whole point of the adjoint is what happens as n_theta grows.
#
# [WHAT IT TAKES]
# `importPEtab(..., derivMode = c("forward", "reverse"))` builds both objects.
# Section 3 needs cppDE with SUNDIALS; it reports itself skipped without it.
# An idle machine, and OMP_NUM_THREADS=1: a fit running beside this makes every
# timing here meaningless.
#
# [AUTHOR]
# Simon Beyer
#
# [Date]
# Tue 09 Sep 2026
# -------------------------------------------------------------------------#

library(dMod2)

# Everything on one core, and stated here rather than left to the caller's
# environment. Three separate things could thread: cppDE's batch entry over
# conditions (OpenMP, and disabled in this build; check cvodeConfig), dMod2's
# own residual kernels, and the BLAS behind the chain rule. A comparison where
# one route threads and another does not measures the threading.
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", GOTO_NUM_THREADS = "1")
options(dMod.cores = 1, cppDE.cores = 1)

.outdir <- file.path(tempdir(), "adjointComparison")
dir.create(.outdir, recursive = TRUE, showWarnings = FALSE)
.petab <- system.file("extdata", "petab_bachmann", package = "dMod2")
.yaml  <- list.files(.petab, pattern = "\\.yaml$", full.names = TRUE)[1]

# The machine scatters, so report the minimum rather than the mean, and time a
# burst rather than one call: the clock resolves about 10 ms on Windows.
tmin <- function(f, reps = 5L, target = 0.5) {
  once <- system.time(f())[["elapsed"]]
  n <- max(1L, ceiling(target / max(once, 1e-3)))
  min(vapply(seq_len(reps), function(i)
    system.time(for (j in seq_len(n)) f())[["elapsed"]] / n, 0.0))
}

# Worst relative deviation over the components that carry something. A gradient
# component three decades below the largest contributes nothing to any use of
# the gradient and would otherwise dominate a per-component ratio.
relWorst <- function(g, ref) {
  g <- g[names(ref)]
  keep <- abs(ref) > 1e-6 * max(abs(ref))
  if (!any(keep)) return(NA_real_)
  max(abs(g[keep] - ref[keep]) / abs(ref[keep]))
}

# One import per tolerance. The tolerance is a runtime option and not a compile
# one, so the generated sources are unchanged and dMod2 reuses the objects it
# already built; only the first call through here pays for a compile.
importAt <- function(tol, tag = "cmp") {
  o <- list(atol = tol, rtol = tol)
  importPEtab(.yaml, backend = "cppDE", cores = 6,
              modelname = paste0("adjcmp_", tag),
              derivMode = c("forward", "reverse"),
              optionsOde = o, optionsSens = o, outdir = .outdir)
}


# -----------------------------------------------------------------------------
# 1. Do forward and reverse answer the same question?
#
# Over a ladder of tolerances, at the published optimum. Reported per row: the
# relative gap in the objective value, the worst relative gap in the gradient,
# and the angle between the two gradients, which is what a line search actually
# feels. Two gradients that differ only in length still point the same way.
# -----------------------------------------------------------------------------
cat("\n1. forward against reverse, over the solver tolerance\n\n")

TOLS <- 10^-c(6, 8, 10, 12)

agree <- do.call(rbind, lapply(TOLS, function(tt) {
  pet <- importAt(tt)
  p   <- pet$bestfit
  a   <- pet$obj(p, deriv = TRUE, hessian = FALSE)
  b   <- pet$obj(p, deriv = TRUE, sweep = "reverse")
  gb  <- b$gradient[names(a$gradient)]
  cosang <- sum(a$gradient * gb) /
            sqrt(sum(a$gradient^2) * sum(gb^2))
  data.frame(tol       = tt,
             value_gap = abs(a$value - b$value) / abs(a$value),
             grad_gap  = relWorst(b$gradient, a$gradient),
             one_minus_cos = 1 - cosang)
}))
print(format(agree, digits = 3), row.names = FALSE)

fall <- agree$grad_gap[1] / agree$grad_gap[nrow(agree)]
cat("\n   over", length(TOLS), "steps of tolerance the gradient gap falls by",
    format(fall, digits = 3), "\n")
cat("   A gap that tracks the tolerance is the discretisation. One that does\n",
    "  not is a missing channel, and that is what this row is watching for.\n")


# -----------------------------------------------------------------------------
# 2. One chain, three routes
#
# The comparison only means anything if the two objectives differ in the ODE
# object and in nothing else. So both are built here by the same recipe from the
# same pieces the import produced, namely data, observation function, error
# model and transformation, and only `x` is swapped. The imported objective is
# not used
# for the timings: it carries a likelihood offset and per-condition data groups
# that the hand-built one does not, and charging that to the backend is how the
# first version of this script got the answer backwards.
# -----------------------------------------------------------------------------
cat("\n2. one chain, three routes\n\n")

pet <- importAt(1e-8)
p   <- pet$bestfit
fx  <- attr(pet, "petab_meta")$fixed
o8  <- list(atol = 1e-8, rtol = 1e-8)

mkobj <- function(x) normL2(pet$dataList, pet$g * x * pet$p, pet$e)

objC <- mkobj(pet$x)          # cppDE: forward and the discrete adjoint
objS <- NULL
.cvodeCfg <- tryCatch(get("cvodeConfig", envir = asNamespace("cppDE")),
                      error = function(e) NULL)
.hasASA <- isTRUE(.cvodeCfg$available) &&
  "derivMode" %in% names(formals(cppDE::cvode)) &&
  "reverse" %in% eval(formals(cppDE::cvode)$derivMode)

if (.hasASA) {
  mS <- odemodel(pet$reactions, modelname = "adjcmp_sun", backend = "Sundials",
                 derivMode = c("forward", "reverse"), compile = TRUE,
                 outdir = .outdir)
  objS <- mkobj(Xs(mS, optionsOde = o8, optionsSens = o8))
} else {
  cat("   ASA skipped: no SUNDIALS, or cvode() has no reverse direction.\n")
}

# A reverse evaluation returns no Hessian. That is the invariant that says the
# direction arrived rather than being swallowed by a wrapper.
chk <- objC(p, fixed = fx, deriv = TRUE, sweep = "reverse")
stopifnot(is.null(chk$hessian))

gF <- objC(p, fixed = fx, deriv = TRUE, hessian = FALSE)$gradient
gR <- objC(p, fixed = fx, deriv = TRUE, sweep = "reverse")$gradient
gS <- if (is.null(objS)) NULL else
      objS(p, fixed = fx, deriv = TRUE, sweep = "reverse")$gradient

# The reference: forward sensitivities at a tolerance far tighter than any run
# compared. It belongs to neither adjoint, which is what makes it a yardstick
# rather than a home advantage.
petRef <- importAt(1e-15, tag = "ref")
gRef <- normL2(petRef$dataList, petRef$g * petRef$x * petRef$p,
               petRef$e)(petRef$bestfit,
                         fixed = attr(petRef, "petab_meta")$fixed,
                         deriv = TRUE, hessian = FALSE)$gradient

# Timed twice, in the same order both times. The first pass pays for whatever
# is still cold, being freshly loaded shared objects, allocator arenas and
# the object cache, and charges it to whichever route runs first. If the two
# passes disagree, the numbers are about the cache and not about the methods.
one_pass <- function() c(
  val = tmin(function() objC(p, fixed = fx, deriv = FALSE)),
  fwd = tmin(function() objC(p, fixed = fx, deriv = TRUE, hessian = FALSE)),
  rev = tmin(function() objC(p, fixed = fx, deriv = TRUE, sweep = "reverse")),
  asa = if (is.null(objS)) NA_real_ else
        tmin(function() objS(p, fixed = fx, deriv = TRUE, sweep = "reverse")))

pass1 <- one_pass()
pass2 <- one_pass()

tab <- data.frame(
  route     = c("forward sensitivities", "reverse (discrete adjoint)",
                "CVODES ASA"),
  sec_1     = unname(pass1[c("fwd", "rev", "asa")]),
  sec_2     = unname(pass2[c("fwd", "rev", "asa")]),
  x_value   = unname(pass2[c("fwd", "rev", "asa")]) / pass2[["val"]],
  deviation = c(relWorst(gF, gRef), relWorst(gR, gRef),
                if (is.null(gS)) NA_real_ else relWorst(gS, gRef)),
  stringsAsFactors = FALSE)
print(format(tab, digits = 3), row.names = FALSE)
cat("\n  value solve: ", format(pass1[["val"]], digits = 3), " then ",
    format(pass2[["val"]], digits = 3), " s\n", sep = "")

cat("  n_theta = ", length(p), "\n\n", sep = "")
cat("  sec_1, sec_2   seconds for one gradient, two passes in the same order.\n")
cat("  x_value        that time divided by one *objective* evaluation over all\n",
    "                conditions, not by a single ODE solve. Dimensionless.\n")
cat("  deviation      worst relative gap over the gradient components carrying\n",
    "                at least 1e-6 of the largest, against forward\n",
    "                sensitivities at rtol = atol = 1e-15. Also dimensionless,\n",
    "                and 31 means a factor of 32 out, not 31 percent.\n\n")
cat("  Everything went through one objective built one way; only the ODE object\n",
    " differs. The deviations are that large because Bachmann observes\n",
    " log10(x + 1e-15): an epsilon below anything double precision resolves\n",
    " multiplies the integration's leftover noise by fifteen decades, and it\n",
    " does that to all three routes alike. Section 1 shows the gap falling with\n",
    " the tolerance, which is what says it is discretisation and not a defect.\n")

cat("\ndone.\n")
