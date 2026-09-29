# numSteady bench: Pimpl at random parameters.
#
#   Rscript bench/numSteady/battery.R [ndraw] [nwalk] [outdir]
#
# ndraw parameter vectors, log10 uniform in [-2, 1], each solved cold in the
# five arms; nwalk random walks of 25 steps (sd 0.15 in log10) solved warm, as
# in a fit. Reports failures, time, the largest residual relative to the
# turnover (states above 1e-12 of the largest), conservation error, the
# largest eigenvalue on ker C, and for the first draws the IFT sensitivities
# against central differences.

args   <- commandArgs(TRUE)
ndraw  <- if (length(args) >= 1L) as.integer(args[1]) else 60L
nwalk  <- if (length(args) >= 2L) as.integer(args[2]) else 10L
outdir <- if (length(args) >= 3L) args[3] else file.path(tempdir(), "numSteady_battery")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressMessages({library(dMod2); library(data.table)})
here <- if (file.exists("bench/numSteady/model.R")) "bench/numSteady" else "."
source(file.path(here, "model.R"))
m <- numSteadyModel(here)
setwd(outdir)

pss <- Pimpl(m$reactions, forcings = m$forcings, compile = TRUE, modelname = "nsb_ss")
E   <- environment(attr(pss, "mappings")[[1]])
dep <- E$dependent
C   <- E$C_dep

kin  <- setdiff(getParameters(pss), c(m$perturbation, "volumeC", "volumeEC", "tSmad2", "tSmad3", "tSmad4"))
arms <- lapply(list(Ctrl = c(0, 0, 0, 0), R1Knd = c(0, 0, 1, 0), R2Knd = c(0, 0, 0, 1),
                    OE_R1 = c(1, 0, 0, 0), OE_R2 = c(0, 1, 0, 0)), setNames, m$perturbation)
draw <- function() c(setNames(10^runif(length(kin), -2, 1), kin), volumeC = 1, volumeEC = 10,
                     tSmad2 = 10^runif(1, -1, 0), tSmad3 = 10^runif(1, -1, 0), tSmad4 = 10^runif(1, -1, 0))

check <- function(out, p) {
  x  <- unclass(out)[dep]
  pv <- setNames(as.numeric(p[E$parms_all]), E$parms_all)
  J  <- E$eval_J(x, pv)[dep, dep]
  f  <- E$eval_G(x, pv)[1, dep]
  act <- dep[x > 0]
  big <- dep[x > 1e-12 * max(x)]
  Ca <- C[, act, drop = FALSE]
  Nb <- MASS::Null(t(Ca))
  list(relres = max(abs(f[big]) / pmax(abs(J[big, ]) %*% abs(x), 1e-300)),
       cq = max(abs(drop(C %*% x) - pv[rownames(C)]) / pv[rownames(C)]),
       maxRe = max(Re(eigen(t(Nb) %*% J[act, act] %*% Nb, only.values = TRUE)$values)))
}
solve1 <- function(p, cond, deriv = TRUE) {
  t <- system.time(o <- tryCatch(pss(p, deriv = deriv, condition = cond)[[1]],
                                 error = function(e) conditionMessage(e)))[3]
  list(out = o, t = t)
}

set.seed(1); cold <- list()
for (i in seq_len(ndraw)) {
  p0 <- draw()
  for (a in names(arms)) {
    p <- c(p0, arms[[a]])
    resetWarmStarts(pss, verbose = FALSE)
    s <- solve1(p, a)
    ok <- !is.character(s$out)
    ck <- if (ok) check(s$out, p) else list(relres = NA, cq = NA, maxRe = NA)
    fd <- NA_real_
    if (ok && i <= 5L) {
      # IFT against central differences in log10 parameters, states above 1e-8
      x  <- unclass(s$out)[dep]; big <- dep[x > 1e-8]
      Jd <- attr(s$out, "deriv")[big, kin] %*% diag(p[kin] * log(10))
      Jf <- sapply(kin, function(k) {
        h <- 1e-5; up <- p; dn <- p; up[k] <- p[k] * 10^h; dn[k] <- p[k] * 10^-h
        (unclass(pss(up, deriv = FALSE, condition = a)[[1]])[big] -
         unclass(pss(dn, deriv = FALSE, condition = a)[[1]])[big]) / (2 * h) })
      fd <- max(abs(Jd - Jf) / abs(x[big]))
    }
    cold[[length(cold) + 1L]] <- data.table(draw = i, arm = a, ok = ok, time = s$t,
                                            relres = ck$relres, cq = ck$cq, maxRe = ck$maxRe, fd = fd)
  }
}
cold <- rbindlist(cold)

set.seed(2); walk <- list()
for (w in seq_len(nwalk)) {
  p0 <- draw(); a <- sample(names(arms), 1); lp <- log10(p0[kin])
  resetWarmStarts(pss, verbose = FALSE)
  for (s in 1:25) {
    if (s > 1) lp <- lp + rnorm(length(lp), 0, 0.15)
    p <- c(p0, arms[[a]]); p[kin] <- 10^lp
    r <- solve1(p, a)
    walk[[length(walk) + 1L]] <- data.table(walk = w, step = s, ok = !is.character(r$out), time = r$t)
  }
}
walk <- rbindlist(walk)
saveRDS(list(cold = cold, walk = walk, stats = as.list(E$stats)), file.path(outdir, "battery.rds"))

cat(sprintf("cold: %d/%d solved, median %.3f s, max relres %.1e, max cq %.1e, max Re %.1e, max fd %.1e\n",
            sum(cold$ok), nrow(cold), median(cold$time), max(cold$relres, na.rm = TRUE),
            max(cold$cq, na.rm = TRUE), max(cold$maxRe, na.rm = TRUE), max(cold$fd, na.rm = TRUE)))
cat(sprintf("warm: %d/%d solved, median %.3f s\n", sum(walk$ok), nrow(walk), median(walk$time)))
s <- E$stats
cat(sprintf("solves %d, iterations per solve %.1f, starts: %s\n", s$solves, s$iter / max(1, s$solves),
            paste(names(s$how), s$how, sep = "=", collapse = " ")))
