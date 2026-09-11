# -------------------------------------------------------------------------#
# One gradient, four ways: forward, two written adjoints, and CVODES ASA
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# dMod2 reaches the same gradient by several routes, and the choice is made in
# one place: the ODE object. Everything above it, the observation function, the
# error model, the parameter transformation and normL2, is the same chain. This
# script builds the routes side by side on one model and one parameter vector,
# so the API that selects them is visible in isolation.
#
#   forward     sensitivity equations carried beside the states. One extra
#               trajectory per parameter, so the cost grows with n_theta.
#   reverse     cppDE's written adjoint. The forward run keeps a checkpoint per
#               accepted step; the backward walk applies each step's adjoint in
#               closed form. One sweep, whatever n_theta is.
#   ASA         SUNDIALS CVODES adjoint sensitivity analysis, the same idea
#               solved as a second ODE rather than differentiated step by step.
#
# The reverse route is shown on both stepper families, because their adjoints
# are different pieces of arithmetic:
#
#   bdf   a corrector method. Its step solves an equation, so the adjoint is the
#         implicit function theorem on that equation, one transposed solve.
#   rb4   a Rosenbrock method. Its stages are direct solves, so the adjoint is
#         their transposed recursion. Nothing is truncated, which is why its
#         gradient sits closer to the forward one.
#
# [WHAT TO EXPECT]
# On this model the gradient costs about three value runs whichever adjoint
# takes it, and the forward route costs about thirty. rb4 is slower here than
# bdf, but its adjoint is not: Bachmann is stiff and rb4 pays for that in the
# forward direction, which both routes share.
#
# The four do not agree to machine precision and should not. A sensitivity run
# lets the step-size controller take the maximum over the state error AND every
# tangent column, so it adapts on a finer grid than a value-only run. Forward
# and reverse therefore differentiate two different discretisations, each
# exactly, and their gap is O(tol). A gap that does not fall with the tolerance
# would be news; a gap that does is arithmetic.
#
# The corollary favours reverse beyond speed: its gradient belongs to the
# trajectory a value-only solve produces, so the value and the gradient a caller
# receives are consistent with each other.
#
# [WHAT IT TAKES]
# The model, the data and the fit come from example_BachmannMSB2011.R, which is
# sourced up to its first objective call rather than transcribed. The ASA route
# needs cppDE built against SUNDIALS; without it that column says so and the
# rest still runs. An idle machine, and one thread: a fit running beside this
# makes every timing here meaningless.
#
# [AUTHOR]
# Simon Beyer
#
# [Date]
# Thu 11 Sep 2026
# -------------------------------------------------------------------------#

library(dMod2)

# Three separate things could thread: cppDE's batch entry over conditions,
# dMod2's own residual kernels, and the BLAS behind the chain rule. A comparison
# where one route threads and another does not measures the threading.
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", GOTO_NUM_THREADS = "1")
options(dMod.cores = 1, cppDE.cores = 1)

.outdir <- file.path(tempdir(), "bachmannReverse")
dir.create(.outdir, recursive = TRUE, showWarnings = FALSE)

# atol under rtol under the observables' own floor: the three have to hold in
# that order, or the gradient is measuring the floor rather than the model.
TOL <- list(atol = 1e-11, rtol = 1e-9, maxsteps = 1e7L, maxattemps = 100L)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# The model, from the example next door
#
# `reactions`, `mydataL` and `bestfit` come from there, along with the compiled
# observation function g, error model e and parameter transformation p that this
# script reuses unchanged. Only the ODE object differs between the routes.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
.example <- system.file("examples", "example_BachmannMSB2011.R", package = "dMod2")
.src     <- readLines(.example)
.upto    <- grep("^obj\\(bestfit", .src)[1] - 1L
eval(parse(text = .src[seq_len(.upto)]), envir = globalenv())

cat(length(mydataL), "conditions,", length(bestfit), "estimated parameters\n")


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# One ODE object per route
#
# `derivMode` is the argument that decides. It is matched, so both directions in
# one object is `c("forward", "reverse")`; asking for only one omits the other's
# entry points and its compile time with them. The stepper is `method`, and it
# is separate: every method carries both directions.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
mBDF <- odemodel(reactions, modelname = "bachRev_bdf", backend = "cppDE",
                 method = "bdf", derivMode = c("forward", "reverse"),
                 compile = FALSE, outdir = .outdir)

mRB4 <- odemodel(reactions, modelname = "bachRev_rb4", backend = "cppDE",
                 method = "rb4", derivMode = c("forward", "reverse"),
                 compile = FALSE, outdir = .outdir)

# cvodeConfig is the package's own record of what ./configure found. It is not
# exported, so this reads it rather than calling it.
.cfg   <- tryCatch(get("cvodeConfig", envir = asNamespace("cppDE")),
                   error = function(e) NULL)
hasASA <- isTRUE(.cfg$available)
mASA <- if (hasASA)
  odemodel(reactions, modelname = "bachRev_asa", backend = "Sundials",
           derivMode = c("forward", "reverse"), compile = FALSE,
           outdir = .outdir) else NULL

xBDF <- Xs(mBDF, optionsOde = TOL, optionsSens = TOL)
xRB4 <- Xs(mRB4, optionsOde = TOL, optionsSens = TOL)
xASA <- if (hasASA) Xs(mASA, optionsOde = TOL, optionsSens = TOL) else NULL

# One compile call for all of them, so the sources share a build.
if (hasASA) {
  compile(xBDF, xRB4, xASA, output = "bachmannReverse", cores = 4)
} else {
  cat("SUNDIALS absent: the ASA route is skipped.\n")
  compile(xBDF, xRB4, output = "bachmannReverse", cores = 4)
}


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# The chain, three times over the same g, e and p
#
# normL2 rather than the example's own objective: that one adds a prior term,
# which is the same scalar on every route and would only pad the comparison.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
objBDF <- normL2(mydataL, g * xBDF * p, e)
objRB4 <- normL2(mydataL, g * xRB4 * p, e)
objASA <- if (hasASA) normL2(mydataL, g * xASA * p, e) else NULL

pars <- bestfit


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Calling it: the same objective, three answers
#
# `deriv = FALSE` is the value alone. `deriv = TRUE` is the forward direction.
# `sweep = "reverse"` asks for the adjoint instead, and returns no Hessian:
# a reverse pass gives the gradient and nothing above it.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
val <- objBDF(pars, deriv = FALSE)
fwd <- objBDF(pars, deriv = TRUE, hessian = FALSE)
rev <- objBDF(pars, deriv = TRUE, sweep = "reverse")

cat("\nvalue      ", format(val$value, digits = 12), "\n")
cat("forward    ", format(fwd$value, digits = 12),
    " gradient of length", length(fwd$gradient), "\n")
cat("reverse    ", format(rev$value, digits = 12),
    " gradient of length", length(rev$gradient),
    if (is.null(rev$hessian)) " and no hessian, by construction\n" else "\n")

# The gradient comes back named and in its own order, which is not the order of
# the parameter vector. Index it by name, never by position.
cat("\nfirst five gradient entries, by name\n")
print(round(rev$gradient[names(pars)[1:5]], 6))


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# What the four routes cost, and how far apart they sit
#
# The machine scatters, so the minimum over repetitions rather than the mean.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
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

# The forward gradient is the oracle: it is the third discretisation of the same
# derivative and the only one the adjoints can be held against here.
g_fwd <- fwd$gradient
g_bdf <- rev$gradient
g_rb4 <- objRB4(pars, deriv = TRUE, sweep = "reverse")$gradient
g_asa <- if (hasASA) objASA(pars, deriv = TRUE, sweep = "reverse")$gradient else NULL

# Against the largest component rather than each component's own size: a forward
# component near zero would otherwise set the whole scale. The angle is what a
# line search actually feels.
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


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# The prediction underneath, for one condition
#
# The chain is callable on its own, which is where to look when a gradient
# surprises: the same g * x * p that the objective holds.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
cond  <- names(mydataL)[1]
times <- sort(unique(c(0, mydataL[[cond]]$time)))
pred  <- (g * xBDF * p)(times, pars, conditions = cond)

cat("\nprediction for condition", cond, ":", nrow(pred[[cond]]), "rows,",
    ncol(pred[[cond]]) - 1L, "observables\n")
print(head(pred[[cond]][, 1:4], 3))


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# The backward pass can be told how accurate to be
#
# `optionsReverse` weights the backward step size by the adjoint of the previous
# evaluation: where lambda is large the grid becomes finer, and it can only
# become finer, so a weight left over from a parameter the optimiser has since
# moved away from costs steps and never accuracy. `gradtol` turns it on.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
xW   <- Xs(mBDF, optionsOde = TOL, optionsSens = TOL,
           optionsReverse = list(gradtol = 1e-6, floor = 1e-3))
objW <- normL2(mydataL, g * xW * p, e)
g_w  <- objW(pars, deriv = TRUE, sweep = "reverse")$gradient

cat("\nweighted backward grid:",
    format(tmin(function() objW(pars, deriv = TRUE, sweep = "reverse")), digits = 3),
    "ms, 1 - cos against forward", format(cosgap(g_fwd, g_w), digits = 3), "\n")
