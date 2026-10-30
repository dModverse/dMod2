# Sparse against dense Hessians on simulated alignments. trust() over the number of
# experiments sets the thresholds of sparse = "auto"; em() by FOCEI over the number of groups;
# one alignment of the size of a large blotIt preprocessing, its runtime and memory.

library(dMod2)

# Data ----
#
# Targets measured on gels at the times of a time course, with replicates. Every gel measures
# a few targets; a target is measured on several gels. The scale is per gel and target, the
# experiment, and the values per target and time are shared over the gels.
simulate <- function(nGel, targetsPerGel = 2, gelsPerTarget = 10, nTime = 8, reps = 2, seed = 1) {
  set.seed(seed)
  nTarget <- max(2, round(nGel * targetsPerGel / gelsPerTarget))
  targets <- sprintf("T%04d", seq_len(nTarget))
  times <- sprintf("t%02d", seq_len(nTime))
  truth <- matrix(rlnorm(nTarget * nTime, 0, 1), nTarget, dimnames = list(targets, times))
  gels <- sprintf("G%04d", seq_len(nGel))
  pairs <- do.call(rbind, lapply(seq_len(nGel), function(g)
    data.frame(gel = gels[g], target = targets[(g * targetsPerGel + seq_len(targetsPerGel) - 1) %% nTarget + 1])))
  pairs <- unique(pairs)
  pairs$scale <- rlnorm(nrow(pairs), 3, 0.5)
  d <- pairs[rep(seq_len(nrow(pairs)), each = nTime * reps), ]
  d$time <- rep(rep(times, each = reps), nrow(pairs))
  d$value <- truth[cbind(d$target, d$time)] * d$scale * rlnorm(nrow(d), 0, 0.1)
  d$experiment <- paste(d$gel, d$target, sep = "_")
  d$scale <- NULL
  rownames(d) <- NULL
  d
}
model <- function(d, scaling = sj ~ experiment)
  alignmentModel(d, model = "yi * sj", errmodel = "sigmaR * prediction", biological = yi ~ target * time,
                 scaling = scaling, error = sigmaR ~ 1, trafo = c("yi <- exp(yi)", "sj <- exp(sj)",
                                                                  "sigmaR <- exp(sigmaR)"))
dense <- list(dMod2.sparseAuto = list(parameters = Inf))
sparse <- list(dMod2.sparseAuto = list(parameters = 0, fill = 1))

# trust(): dense against sparse over the number of experiments ----
#
# The minimum over three fits of each; the fill of the Hessian at the start.
sizes <- c(5, 10, 20, 50, 100, 200, 500, 1000)
timing <- do.call(rbind, lapply(sizes, function(nGel) {
  m <- model(simulate(nGel))
  P <- length(dMod2:::.alignmentPars(m)$fixed)
  obj <- dMod2:::.alignmentObjective(m) + dMod2:::.alignmentNormalization(m)
  st <- dMod2:::.alignmentStart(m, dMod2:::.alignmentPars(m)$fixed)
  fill <- Matrix::nnzero(obj(st, hessian = "sparse")$hessian) / P^2
  run <- function(o) min(replicate(3, system.time(withr::with_options(o, align(m, fits = 1)))[["elapsed"]]))
  data.frame(gels = nGel, parameters = P, fill = fill, readings = nrow(m$data),
             dense = if (P <= 2500) run(dense) else NA, sparse = run(sparse))
}))
timing
plot(sparse ~ parameters, timing, log = "xy", type = "b", ylab = "seconds", ylim = range(timing[c("dense", "sparse")], na.rm = TRUE))
lines(dense ~ parameters, timing, type = "b", col = 2)
legend("topleft", c("sparse", "dense"), col = 1:2, lty = 1)

# em(): FOCEI dense against sparse over the number of groups ----
#
# A random scale per gel on top of the scale per experiment, so every gel is a group.
focei <- do.call(rbind, lapply(c(50, 200, 1000), function(nGel) {
  m <- model(simulate(nGel), scaling = sj ~ experiment + (1 | gel))
  P <- length(dMod2:::.alignmentPars(m)$fixed)
  run <- function(o) min(replicate(2, system.time(withr::with_options(o, align(m, fits = 1)))[["elapsed"]]))
  data.frame(groups = nGel, parameters = P, readings = nrow(m$data),
             dense = if (P <= 2500) run(dense) else NA, sparse = run(sparse))
}))
focei

# An alignment of blotIt size ----
#
# 106 gels with 100 targets each, 52 combinations of time and treatment in triplicate: 1.65e6
# readings, 5200 values and 10600 scales. Run it alone; gc() reports the peak of the R heap.
big <- simulate(106, targetsPerGel = 100, gelsPerTarget = 106, nTime = 52, reps = 3)
nrow(big)
invisible(gc(reset = TRUE))
t0 <- Sys.time()
mBig <- model(big)
Sys.time() - t0
fitBig <- align(mBig, fits = 1)
Sys.time() - t0
fitBig
fitBig$fit$iterations
gc()
