# Counts Steihaug-Toint CG iterations on the trust-region subproblem along a
# real optimisation path. A matrix-free step costs 2k linear passes for k CG
# iterations, against n_theta forward sensitivity directions today, so k(n_theta)
# decides whether the adjoint has a caller at all.
#
# normL2's Hessian is Gauss-Newton only when no noise parameter is estimated.
# With dsigma present it also carries -2*wr/s^2 (dpred x dsigma + dsigma x dpred)
# and (4*wr^2 - 2)/s^2 dsigma x dsigma, which are indefinite. Moré-Sorensen
# handles that natively, Steihaug only by stopping at the first negative
# curvature direction, so the exit reason is part of the measurement.
#
# The count is taken on the plain subproblem min g'p + p'Hp/2 s.t. ||p|| <= r.
# trust()'s reflective boundary solves a Coleman-Li scaled variant; that changes
# the step, not the conditioning this measures.

library(dMod2)

# OMP_NUM_THREADS above 1 fails internally. Parallelism goes through `cores`.
Sys.setenv(OMP_NUM_THREADS = 1)
cores <- 6
setwd(tempdir())

# The backend defaults let trust() stagnate: the radius collapses while the
# value stands still, so a path has almost no accepted iterates to measure.
solverOpts <- list(atol = 1e-9, rtol = 1e-9)

# Displacement of the start from the published optimum, on the PEtab scale.
# Large enough for a path worth measuring, small enough to converge back.
startSd <- 0.1

## --- Steihaug-Toint truncated CG -------------------------------------------

# Positive root of ||z + tau d|| = radius.
boundaryStep <- function(z, d, radius) {
  dd <- sum(d * d)
  zd <- sum(z * d)
  zz <- sum(z * z)
  (-zd + sqrt(zd * zd - dd * (zz - radius * radius))) / dd
}

steihaug <- function(H, g, radius, tol = 1e-8, maxit = 10 * length(g)) {
  z  <- numeric(length(g))
  r  <- -g
  d  <- r
  rr <- sum(r * r)
  gnorm <- sqrt(rr)
  if (gnorm == 0) return(list(step = z, iter = 0L, exit = "zerograd"))
  for (j in seq_len(maxit)) {
    Hd  <- as.vector(H %*% d)
    dHd <- sum(d * Hd)
    if (dHd <= 0)
      return(list(step = z + boundaryStep(z, d, radius) * d,
                  iter = j, exit = "curvature"))
    alpha <- rr / dHd
    znext <- z + alpha * d
    if (sqrt(sum(znext * znext)) >= radius)
      return(list(step = z + boundaryStep(z, d, radius) * d,
                  iter = j, exit = "boundary"))
    z   <- znext
    r   <- r - alpha * Hd
    rrn <- sum(r * r)
    if (sqrt(rrn) <= tol * gnorm)
      return(list(step = z, iter = j, exit = "tolerance"))
    d  <- r + (rrn / rr) * d
    rr <- rrn
  }
  list(step = z, iter = maxit, exit = "maxit")
}

# Per-iterate row: k inside the true radius, k unconstrained as a conditioning
# proxy, and the spectrum that decides whether Steihaug can use the curvature.
modelDecrease <- function(H, g, p) -(sum(g * p) + 0.5 * sum(p * (H %*% p)))

cgRow <- function(model, obj, fit, i, ...) {
  path <- fit$argpath
  o  <- obj(structure(path[i, ], names = colnames(path)), ...)
  ev <- eigen(o$hessian, symmetric = TRUE, only.values = TRUE)$values
  tr <- steihaug(o$hessian, o$gradient, radius = fit$r[i])
  fu <- steihaug(o$hessian, o$gradient, radius = Inf, tol = 1e-6)
  pos <- ev[ev > 0]
  # trust()'s own trial step, scored on the same model, is the quality baseline.
  stepMS <- fit$argtry[i, ] - path[i, ]
  data.frame(model = model, npar = ncol(path), iterate = i,
             radius = fit$r[i], value = o$value,
             kTrust = tr$iter, exitTrust = tr$exit,
             kFull = fu$iter, exitFull = fu$exit,
             evMin = min(ev), nNeg = sum(ev <= 0),
             condPos = max(pos) / min(pos),
             mdCG = modelDecrease(o$hessian, o$gradient, tr$step),
             mdMS = modelDecrease(o$hessian, o$gradient, stepMS))
}

## --- Correctness of the CG loop ---------------------------------------------

# On a positive definite system without an active radius, CG must reproduce the
# Newton step; with a radius it must stop exactly on the boundary.
set.seed(7)
Arand <- matrix(rnorm(40 * 40), 40)
Hpd   <- crossprod(Arand) + diag(40)
gpd   <- rnorm(40)
cgPD  <- steihaug(Hpd, gpd, radius = Inf, tol = 1e-14)
stopifnot(max(abs(cgPD$step - solve(Hpd, -gpd))) < 1e-8)
cgTR  <- steihaug(Hpd, gpd, radius = 0.01)
stopifnot(identical(cgTR$exit, "boundary"),
          abs(sqrt(sum(cgTR$step^2)) - 0.01) < 1e-12)

## --- Boehm, 9 parameters ----------------------------------------------------

petabBoehm <- importPEtab(
  system.file("extdata/petab_boehm/Boehm.yaml", package = "dMod2"),
  backend = "cppDE", cores = cores, outdir = tempdir(),
  optionsOde = solverOpts, optionsSens = solverOpts)

objBoehm <- petabBoehm$obj
set.seed(1)
startBoehm <- pmin(pmax(petabBoehm$bestfit +
                          rnorm(length(petabBoehm$bestfit), sd = startSd),
                        petabBoehm$parlower + 1e-3),
                   petabBoehm$parupper - 1e-3)

fitBoehm <- trust(objBoehm, startBoehm, rinit = 1, rmax = 10,
                  iterlim = 200, gtol = 1e-4, blather = TRUE,
                  parlower = petabBoehm$parlower,
                  parupper = petabBoehm$parupper)

cgBoehm <- do.call(rbind, lapply(which(fitBoehm$accept), function(i)
  cgRow("Boehm", objBoehm, fitBoehm, i)))

cgBoehm

## --- Bachmann, 113 parameters -----------------------------------------------

yamlBachmann <- path.expand(file.path(
  "~/Documents/Projects/dModverse/cppDE/benchmarks/cache/petab/Benchmark-Models",
  "Bachmann_MSB2011/Bachmann_MSB2011.yaml"))

petabBachmann <- importPEtab(yamlBachmann, backend = "cppDE",
                             cores = cores, outdir = tempdir(),
                             optionsOde = solverOpts,
                             optionsSens = solverOpts)

objBachmann <- petabBachmann$obj
set.seed(1)
startBachmann <- pmin(pmax(petabBachmann$bestfit +
                             rnorm(length(petabBachmann$bestfit), sd = startSd),
                           petabBachmann$parlower + 1e-3),
                      petabBachmann$parupper - 1e-3)

fitBachmann <- trust(objBachmann, startBachmann, rinit = 1, rmax = 10,
                     iterlim = 50, gtol = 1e-4, blather = TRUE,
                     parlower = petabBachmann$parlower,
                     parupper = petabBachmann$parupper,
                     cores = cores)

cgBachmann <- do.call(rbind, lapply(which(fitBachmann$accept), function(i)
  cgRow("Bachmann", objBachmann, fitBachmann, i, cores = cores)))

cgBachmann

## --- Lang, 294 parameters ---------------------------------------------------
## Heavy: 124 states, generating and compiling the sensitivity model takes long.

yamlLang <- path.expand(file.path(
  "~/Documents/Projects/dModverse/cppDE/benchmarks/cache/petab/Benchmark-Models",
  "Lang_PLOSComputBiol2024/Lang_PLOSComputBiol2024.yaml"))

petabLang <- importPEtab(yamlLang, backend = "cppDE",
                         cores = cores, outdir = tempdir(),
                         optionsOde = solverOpts, optionsSens = solverOpts)

objLang <- petabLang$obj
set.seed(1)
startLang <- pmin(pmax(petabLang$bestfit +
                         rnorm(length(petabLang$bestfit), sd = startSd),
                       petabLang$parlower + 1e-3),
                  petabLang$parupper - 1e-3)

fitLang <- trust(objLang, startLang, rinit = 1, rmax = 10,
                 iterlim = 20, gtol = 1e-4, blather = TRUE,
                 parlower = petabLang$parlower,
                 parupper = petabLang$parupper,
                 cores = cores)

cgLang <- do.call(rbind, lapply(which(fitLang$accept), function(i)
  cgRow("Lang", objLang, fitLang, i, cores = cores)))

cgLang

## --- Verdict ----------------------------------------------------------------
## Matrix-free wins where 2k stays below npar. kFull is the conditioning-bound
## count, kTrust what a Steihaug step would actually spend.

cg <- rbind(cgBoehm, cgBachmann, cgLang)
cg$ratioTrust <- 2 * cg$kTrust / cg$npar
cg$ratioFull  <- 2 * cg$kFull  / cg$npar
cg$quality    <- cg$mdCG / cg$mdMS

with(cg, boxplot(kFull ~ model, ylab = "CG iterations to tolerance", xlab = ""))
with(cg, plot(npar, ratioFull, log = "xy", pch = 19,
              xlab = "n_theta", ylab = "2k / n_theta"))
abline(h = 1, lty = 2)

aggregate(cbind(kTrust, kFull, ratioFull, nNeg, quality) ~ model + npar, cg,
          median)
table(cg$model, cg$exitFull)
