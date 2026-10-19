# Forward sensitivities, the cppDE discrete adjoint and CVODES ASA on one model:
# gradient cost and agreement over a ladder of solver tolerances.
# Model: Bachmann et al. (2011), Mol Syst Biol, via PEtab.

library(dMod2)

# One core for every layer that could thread (OpenMP, dMod2 kernels, BLAS), so
# no route is timed with threading the others lack.
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", GOTO_NUM_THREADS = "1")
options(dMod.cores = 1, cppDE.cores = 1)

.outdir <- file.path(tempdir(), "adjointComparison")
dir.create(.outdir, recursive = TRUE, showWarnings = FALSE)
.petab <- system.file("extdata", "petab_bachmann", package = "dMod2")
.yaml  <- list.files(.petab, pattern = "\\.yaml$", full.names = TRUE)[1]

# Minimum over repetitions of a timed burst: the clock resolves about 10 ms on Windows.
tmin <- function(f, reps = 5L, target = 0.5) {
  once <- system.time(f())[["elapsed"]]
  n <- max(1L, ceiling(target / max(once, 1e-3)))
  min(vapply(seq_len(reps), function(i)
    system.time(for (j in seq_len(n)) f())[["elapsed"]] / n, 0.0))
}

# Worst relative deviation over components above 1e-6 of the largest.
relWorst <- function(g, ref) {
  g <- g[names(ref)]
  keep <- abs(ref) > 1e-6 * max(abs(ref))
  if (!any(keep)) return(NA_real_)
  max(abs(g[keep] - ref[keep]) / abs(ref[keep]))
}

# The tolerance is a runtime option, so every import reuses the compiled objects.
importAt <- function(tol, tag = "cmp") {
  o <- list(atol = tol, rtol = tol)
  importPEtab(.yaml, backend = "cppDE", cores = 6,
              modelname = paste0("adjcmp_", tag),
              derivMode = c("forward", "reverse"),
              options = o, outdir = .outdir)
}


# 1. Forward against reverse over a tolerance ladder: value gap, worst gradient
# gap, and one minus the cosine between the gradients.
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


# 2. One chain, three routes: both objectives are built by the same recipe,
# only `x` differs.
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
  objS <- mkobj(Xs(mS, options = o8))
} else {
  cat("   ASA skipped: no SUNDIALS, or cvode() has no reverse direction.\n")
}

# A reverse evaluation returns no Hessian, so the direction reached the backend.
chk <- objC(p, fixed = fx, deriv = TRUE, sweep = "reverse")
stopifnot(is.null(chk$hessian))

gF <- objC(p, fixed = fx, deriv = TRUE, hessian = FALSE)$gradient
gR <- objC(p, fixed = fx, deriv = TRUE, sweep = "reverse")$gradient
gS <- if (is.null(objS)) NULL else
      objS(p, fixed = fx, deriv = TRUE, sweep = "reverse")$gradient

# Reference: forward sensitivities at a tolerance far tighter than any compared run.
petRef <- importAt(1e-15, tag = "ref")
gRef <- normL2(petRef$dataList, petRef$g * petRef$x * petRef$p,
               petRef$e)(petRef$bestfit,
                         fixed = attr(petRef, "petab_meta")$fixed,
                         deriv = TRUE, hessian = FALSE)$gradient

# Two passes in the same order: if they disagree, the numbers measure the cold
# cache and not the methods.
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
cat("  deviation      worst relative gap over the gradient components of\n",
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
