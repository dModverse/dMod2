# -------------------------------------------------------------------------#
# Forward, reverse and ASA, one condition at a time
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# Two measurements of the same two methods disagree by a factor of five, and
# this script exists to find out why.
#
#   one condition, one bare solve   our reverse is 2.6x faster than ASA
#   36 conditions, whole objective  ASA is 1.9x faster than our reverse
#
# Threading is already excluded: this cppDE has no OpenMP, the BLAS is R's own
# reference build, and dMod2's kernels take their thread count from
# `dMod.cores`. All three are pinned below anyway, because a comparison in
# which one route threads and another does not measures the threading.
#
# What is left is the conditions. The bare solve used the first one, which need
# not be representative: they differ in time points, in stiffness, and in how
# many steps the integration takes. So build one objective per condition, run
# all three routes through each, and see whether the per-condition times add up
# to the whole-objective time. If they do not, the chain is doing something
# that is not the sum of its conditions, and that is the thing to look at.
#
# [WHAT IT TAKES]
# SUNDIALS, for the ASA column; without it that column is NA and the rest still
# runs. An idle machine: every number here is a wall time.
#
# [AUTHOR]
# Simon Beyer
#
# [Date]
# Tue 09 Sep 2026
# -------------------------------------------------------------------------#

library(dMod2)

Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", GOTO_NUM_THREADS = "1")
options(dMod.cores = 1, cppDE.cores = 1)

.outdir <- file.path(tempdir(), "bench_perCondition")
dir.create(.outdir, recursive = TRUE, showWarnings = FALSE)
.yaml <- list.files(system.file("extdata", "petab_bachmann", package = "dMod2"),
                    pattern = "[.]yaml$", full.names = TRUE)[1]
TOL <- list(atol = 1e-8, rtol = 1e-8)

# The machine scatters, so report the minimum rather than the mean, and time a
# burst rather than one call: the clock resolves about 10 ms on Windows.
tmin <- function(f, reps = 3L, target = 0.15) {
  once <- system.time(f())[["elapsed"]]
  n <- max(1L, ceiling(target / max(once, 1e-3)))
  min(vapply(seq_len(reps),
             function(i) system.time(for (j in seq_len(n)) f())[["elapsed"]] / n, 0.0))
}


# -----------------------------------------------------------------------------
# One import, two prediction chains
#
# Both objectives are built from the same pieces by the same recipe, and only
# the ODE object differs. Building one of them from the import and the other by
# hand is how an earlier version of this comparison got the sign wrong: the
# imported objective carries a likelihood offset and per-condition data groups
# that a hand-built one does not.
# -----------------------------------------------------------------------------
pet <- importPEtab(.yaml, backend = "cppDE", cores = 6, modelname = "bpc",
                   derivMode = c("forward", "reverse"),
                   optionsOde = TOL, optionsSens = TOL, outdir = .outdir)
pars  <- pet$bestfit
fixed <- attr(pet, "petab_meta")$fixed

.cfg <- tryCatch(get("cvodeConfig", envir = asNamespace("cppDE")),
                 error = function(e) NULL)
hasASA <- isTRUE(.cfg$available) &&
  "sweep" %in% names(formals(cppDE::cvode)) &&
  "reverse" %in% eval(formals(cppDE::cvode)$sweep)

xS <- NULL
if (hasASA) {
  mS <- odemodel(pet$reactions, modelname = "bpc_sun", backend = "Sundials",
                 derivMode = c("forward", "reverse"), compile = TRUE,
                 outdir = .outdir)
  xS <- Xs(mS, optionsOde = TOL, optionsSens = TOL)
} else {
  cat("SUNDIALS absent or cvode() has no reverse direction: ASA column is NA.\n")
}

mkobj <- function(x, dat) normL2(dat, pet$g * x * pet$p, pet$e)


# -----------------------------------------------------------------------------
# Per condition
#
# `steps` is the accepted step count of the reverse solve for that condition,
# read off the sweep's own grid. It is the thing our cost is proportional to and
# ASA's is not: we replay every accepted forward step, ASA integrates the
# adjoint on a grid of its own choosing.
# -----------------------------------------------------------------------------
conds <- names(pet$dataList)
rows <- lapply(conds, function(cn) {
  dat  <- pet$dataList[cn]
  objC <- mkobj(pet$x, dat)
  objS <- if (is.null(xS)) NULL else mkobj(xS, dat)

  # A reverse evaluation returns no Hessian: the one invariant that says the
  # direction arrived rather than being swallowed by a wrapper in the chain.
  chk <- objC(pars, fixed = fixed, deriv = TRUE, sweep = "reverse")
  if (!is.null(chk$hessian))
    stop("condition ", cn, ": sweep = \"reverse\" returned a Hessian")

  t_val <- tmin(function() objC(pars, fixed = fixed, deriv = FALSE))
  t_fwd <- tmin(function() objC(pars, fixed = fixed, deriv = TRUE, hessian = FALSE))
  t_rev <- tmin(function() objC(pars, fixed = fixed, deriv = TRUE, sweep = "reverse"))
  t_asa <- if (is.null(objS)) NA_real_ else
           tmin(function() objS(pars, fixed = fixed, deriv = TRUE, sweep = "reverse"))

  steps <- tryCatch({
    inner <- pet$p(pars, fixed = fixed)[[cn]]
    tt <- sort(unique(c(0, dat[[cn]]$time)))
    g <- cppDE::solveODE(pet$odemodel$reversed, tt, inner,
                         seed = array(1, c(length(tt),
                                           length(attr(pet$odemodel$func, "variables")), 1L)),
                         adjointGrid = TRUE, abstol = TOL$atol, reltol = TOL$rtol)
    length(g$adjointGrid$h)
  }, error = function(e) NA_integer_)

  data.frame(condition = cn, n_times = nrow(dat[[cn]]), steps = steps,
             t_value = t_val, t_fwd = t_fwd, t_rev = t_rev, t_asa = t_asa,
             stringsAsFactors = FALSE)
})
per <- do.call(rbind, rows)
per$rev_over_asa <- per$t_rev / per$t_asa
per <- per[order(-per$t_rev), ]

cat("\nPer condition, one objective each, everything on one core\n\n")
print(format(per, digits = 3), row.names = FALSE)


# -----------------------------------------------------------------------------
# Do the parts add up to the whole?
#
# If the sum over conditions matches the objective built over all of them, the
# chain is the sum of its conditions and the per-condition table explains the
# whole. If it does not, the difference is the thing to chase.
# -----------------------------------------------------------------------------
objC_all <- mkobj(pet$x, pet$dataList)
objS_all <- if (is.null(xS)) NULL else mkobj(xS, pet$dataList)

w_val <- tmin(function() objC_all(pars, fixed = fixed, deriv = FALSE))
w_fwd <- tmin(function() objC_all(pars, fixed = fixed, deriv = TRUE, hessian = FALSE))
w_rev <- tmin(function() objC_all(pars, fixed = fixed, deriv = TRUE, sweep = "reverse"))
w_asa <- if (is.null(objS_all)) NA_real_ else
         tmin(function() objS_all(pars, fixed = fixed, deriv = TRUE, sweep = "reverse"))

cmp <- data.frame(
  route      = c("value", "forward", "reverse", "ASA"),
  sum_of_parts = c(sum(per$t_value), sum(per$t_fwd), sum(per$t_rev), sum(per$t_asa)),
  whole        = c(w_val, w_fwd, w_rev, w_asa),
  stringsAsFactors = FALSE)
cmp$whole_over_parts <- cmp$whole / cmp$sum_of_parts

cat("\n\nThe parts against the whole\n\n")
print(format(cmp, digits = 3), row.names = FALSE)

cat("\n  conditions:", length(conds), "  n_theta:", length(pars), "\n")
cat("  A ratio near one means the objective is the sum of its conditions and\n",
    " the table above accounts for it. Well under one means the whole shares\n",
    " work the parts each repeat, one compiled prediction call for all\n",
    " conditions rather than 36, and then per-condition timings overstate\n",
    " every route, though not necessarily by the same factor.\n")
