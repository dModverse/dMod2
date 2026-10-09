# One gradient four ways on Bachmann et al. (2011): forward sensitivities, the
# cppDE adjoint on the bdf and rb4 steppers, and CVODES ASA (if SUNDIALS is
# available). Timings need an idle machine.

library(dMod2)

# Single-threaded everywhere, so the comparison does not measure threading.
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", GOTO_NUM_THREADS = "1")
options(dMod.cores = 1, cppDE.cores = 1)

.outdir <- file.path(tempdir(), "bachmannReverse")
dir.create(.outdir, recursive = TRUE, showWarnings = FALSE)

# atol below rtol below the observables' floor, or the gradient measures the floor.
TOL <- list(atol = 1e-11, rtol = 1e-9, maxsteps = 1e7L, maxattempts = 100L)


# The model, from example_BachmannMSB2011.R --------------------------------------
# It supplies `reactions`, `mydataL`, `bestfit` and the compiled g, e and p;
# only the ODE object differs between routes.
.example <- system.file("examples", "example_BachmannMSB2011.R", package = "dMod2")
.src     <- readLines(.example)
.upto    <- grep("^obj\\(bestfit", .src)[1] - 1L
eval(parse(text = .src[seq_len(.upto)]), envir = globalenv())

cat(length(mydataL), "conditions,", length(bestfit), "estimated parameters\n")


# One ODE object per route ------------------------------------------------------
# `derivMode` selects the directions and their entry points; the stepper
# `method` is independent of it.
mBDF <- odemodel(reactions, modelname = "bachRev_bdf", backend = "cppDE",
                 method = "bdf", derivMode = c("forward", "reverse"),
                 compile = FALSE, outdir = .outdir)

mRB4 <- odemodel(reactions, modelname = "bachRev_rb4", backend = "cppDE",
                 method = "rb4", derivMode = c("forward", "reverse"),
                 compile = FALSE, outdir = .outdir)

# cvodeConfig records what cppDE's configure found; it is read, not called.
.cfg   <- tryCatch(get("cvodeConfig", envir = asNamespace("cppDE")),
                   error = function(e) NULL)
hasASA <- isTRUE(.cfg$available)
mASA <- if (hasASA)
  odemodel(reactions, modelname = "bachRev_asa", backend = "Sundials",
           derivMode = c("forward", "reverse"), compile = FALSE,
           outdir = .outdir) else NULL

xBDF <- Xs(mBDF, options = TOL)
xRB4 <- Xs(mRB4, options = TOL)
xASA <- if (hasASA) Xs(mASA, options = TOL) else NULL

# One compile call for all of them, so the sources share a build.
if (hasASA) {
  compile(xBDF, xRB4, xASA, output = "bachmannReverse", cores = 4)
} else {
  cat("SUNDIALS absent: the ASA route is skipped.\n")
  compile(xBDF, xRB4, output = "bachmannReverse", cores = 4)
}


# The chain, three times over the same g, e and p -------------------------------
# normL2 without the example's prior, which is the same scalar on every route.
objBDF <- normL2(mydataL, g * xBDF * p, e)
objRB4 <- normL2(mydataL, g * xRB4 * p, e)
objASA <- if (hasASA) normL2(mydataL, g * xASA * p, e) else NULL

pars <- bestfit


# Calling it: the same objective, three answers ---------------------------------
# `deriv = FALSE` is the value, `deriv = TRUE` the forward direction, and
# `sweep = "reverse"` the adjoint, which returns a gradient and no Hessian.
val <- objBDF(pars, deriv = FALSE)
fwd <- objBDF(pars, deriv = TRUE, hessian = FALSE)
rev <- objBDF(pars, deriv = TRUE, sweep = "reverse")

cat("\nvalue      ", format(val$value, digits = 12), "\n")
cat("forward    ", format(fwd$value, digits = 12),
    " gradient of length", length(fwd$gradient), "\n")
cat("reverse    ", format(rev$value, digits = 12),
    " gradient of length", length(rev$gradient),
    if (is.null(rev$hessian)) " and no hessian, by construction\n" else "\n")

# The gradient is named and not in parameter order: index it by name.
cat("\nfirst five gradient entries, by name\n")
print(round(rev$gradient[names(pars)[1:5]], 6))


# What the four routes cost, and how far apart they sit -------------------------
tmin <- function(f, n = 5L) {
  f()                                     # warm the caches and the handles
  min(replicate(n, system.time(f())[["elapsed"]])) * 1000
}

t_val <- tmin(function() objBDF(pars, deriv = FALSE))
t_fwd <- tmin(function() objBDF(pars, deriv = TRUE, hessian = FALSE))
t_bdf <- tmin(function() objBDF(pars, deriv = TRUE, sweep = "reverse"))
t_rb4 <- tmin(function() objRB4(pars, deriv = TRUE, sweep = "reverse"))
t_asa <- if (hasASA) tmin(function() objASA(pars, deriv = TRUE, sweep = "reverse")) else NA

timing <- data.frame(
  route     = c("value", "forward", "reverse bdf", "reverse rb4", "ASA"),
  ms        = c(t_val, t_fwd, t_bdf, t_rb4, t_asa),
  in_values = c(t_val, t_fwd, t_bdf, t_rb4, t_asa) / t_val)

cat("\nOne objective call, milliseconds\n\n")
print(format(timing, digits = 3), row.names = FALSE)

# The forward gradient is the oracle the adjoints are held against.
g_fwd <- fwd$gradient
g_bdf <- rev$gradient
g_rb4 <- objRB4(pars, deriv = TRUE, sweep = "reverse")$gradient
g_asa <- if (hasASA) objASA(pars, deriv = TRUE, sweep = "reverse")$gradient else NULL

# Relative to the largest component, so a near-zero component does not set the
# scale; the angle is what a line search feels.
relgap <- function(a, b) { nm <- names(a); max(abs(b[nm] - a[nm])) / max(abs(a[nm])) }
cosgap <- function(a, b) { nm <- names(a)
  1 - sum(a[nm] * b[nm]) / sqrt(sum(a[nm]^2) * sum(b[nm]^2)) }

agree <- data.frame(
  against_forward    = c("reverse bdf", "reverse rb4", "ASA"),
  worst_over_largest = c(relgap(g_fwd, g_bdf), relgap(g_fwd, g_rb4),
                         if (hasASA) relgap(g_fwd, g_asa) else NA),
  one_minus_cos      = c(cosgap(g_fwd, g_bdf), cosgap(g_fwd, g_rb4),
                         if (hasASA) cosgap(g_fwd, g_asa) else NA))

cat("\nThe adjoints against the forward gradient\n\n")
print(format(agree, digits = 3), row.names = FALSE)


# The prediction underneath, for one condition ----------------------------------
# The chain the objective holds is callable on its own.
cond  <- names(mydataL)[1]
times <- sort(unique(c(0, mydataL[[cond]]$time)))
pred  <- (g * xBDF * p)(times, pars, conditions = cond)

cat("\nprediction for condition", cond, ":", nrow(pred[[cond]]), "rows,",
    ncol(pred[[cond]]) - 1L, "observables\n")
print(head(pred[[cond]][, 1:4], 3))


# The backward pass can be told how accurate to be ------------------------------
# `optionsReverse$refine` puts each backward step under an error test; `gradtol`
# adds the step's share of the gradient to it.
xR   <- Xs(mBDF, options = TOL,
           optionsReverse = list(refine = TRUE, gradtol = 1e-6))
objR <- normL2(mydataL, g * xR * p, e)
g_r  <- objR(pars, deriv = TRUE, sweep = "reverse")$gradient

cat("\nchecked sweep:",
    format(tmin(function() objR(pars, deriv = TRUE, sweep = "reverse")), digits = 3),
    "ms, 1 - cos against forward", format(cosgap(g_fwd, g_r), digits = 3), "\n")
