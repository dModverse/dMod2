# Plot functions on profiles, fits and predictions, without compiled code.

pdf(NULL)

# a is coupled to b and c, d is independent
.coupled <- function(pars, ...) {
  a <- pars[["a"]]; b <- pars[["b"]]; c <- pars[["c"]]; d <- pars[["d"]]
  v <- (a - b)^2 + (a - 2*c)^2 + 0.1*(a - 1)^2 + (d - 1)^2
  g <- c(a = 2*(a - b) + 2*(a - 2*c) + 0.2*(a - 1), b = -2*(a - b), c = -4*(a - 2*c), d = 2*(d - 1))
  H <- matrix(0, 4, 4, dimnames = list(names(g), names(g)))
  H["a", "a"] <- 4.2; H["a", "b"] <- H["b", "a"] <- -2; H["a", "c"] <- H["c", "a"] <- -4
  H["b", "b"] <- 2; H["c", "c"] <- 8; H["d", "d"] <- 2
  objlist(value = v, gradient = g, hessian = H)
}
.coupledProfile <- function()
  profile(.coupled, c(a = 1, b = 1, c = 0.5, d = 1), whichPar = "a", limits = c(-1, 1), cores = 1)

.fitFrame <- function(values) {
  parframe(data.frame(value = values, converged = TRUE, iterations = 1, index = seq_along(values),
                      a = seq_along(values)),
           parameters = "a", metanames = c("value", "converged", "iterations", "index"))
}

test_that("plotPathsMulti normalizes every path to the range -1 to 1", {
  skip_if_not_installed("cowplot")
  d <- as.data.frame(plotPathsMulti(.coupledProfile(), "a", normalizePaths = TRUE)$data)
  for (pt in c("b", "c")) {
    y <- d$y[d$partner == pt]
    expect_equal(max(abs(y)), 1)
    expect_gt(length(unique(y)), 2)
  }
})

test_that("plotPars detects steps also on a frame not sorted by value", {
  d <- as.data.frame(plotPars(.fitFrame(c(40, 10, 15, 10.2)))$data)
  steps <- unique(d[, c("value", "step")])
  expect_equal(as.integer(as.character(steps$step[order(steps$value)])), c(1, 1, 2, 3))
})

test_that("a layer added to plotCombined inherits the plotCombined data", {
  x <- Xd(data.frame(name = "A", time = 0:2, row.names = paste0("A", 0:2)), condition = "C1")
  pars <- c(A0 = 1, A1 = 0.5, A2 = 0.2)
  data <- as.datalist(data.frame(name = "A", time = 0:2, value = pars + 0.01, sigma = 0.1,
                                 condition = "C1"))
  p <- plotCombined(x(seq(0, 2, 0.5), pars), data) + ggplot2::geom_point(ggplot2::aes(y = value))
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plotResiduals labels each fit with the residuals of the fit row", {
  x <- Xd(data.frame(name = "A", time = 0:1, row.names = c("a0", "a1")), condition = "C1")
  data <- as.datalist(data.frame(name = "A", time = 0:1, value = c(1, 2), sigma = 1, condition = "C1"))
  pf <- parframe(data.frame(value = c(8, 0), index = 1:2, a0 = c(3, 1), a1 = c(4, 2)),
                 parameters = c("a0", "a1"), metanames = c("value", "index"))
  d <- as.data.frame(plotResiduals(pf, x, data, split = c("index", "condition"))$data)
  expect_equal(d$res[order(d$index)], c(8, 0))
})

test_that("plotPars and plotValues take a positional subset", {
  fs <- .fitFrame(c(1, 2, 3))
  fs$converged <- c(TRUE, FALSE, TRUE)
  expect_equal(nrow(plotPars(fs, converged)$data), 2L)
  expect_no_error(ggplot2::ggplot_build(plotValues(fs, converged)))
})

test_that("plotValues takes the fits of mstrust() and a named list, each collection ranked separately", {
  fits <- mstrust(.coupled, c(a = 1, b = 1, c = 0.5, d = 1), fits = 3, sd = 1, cores = 1, iterlim = 50)
  expect_equal(attr(plotValues(fits), "data")$value, sort(as.parframe(fits)$value))
  expect_equal(names(as.parvec(fits)), c("a", "b", "c", "d"))
  lim <- 2.5
  p <- plotValues(list(first = .fitFrame(c(3, 1)), second = .fitFrame(c(5, 2, 4))), value < lim + 2,
                  relative = TRUE)
  d <- attr(p, "data")
  expect_equal(levels(d$set), c("first", "second"))
  expect_equal(d$delta[d$set == "second"], c(0, 2))
  expect_equal(d$index[d$set == "second"], 1:2)
  expect_no_error(ggplot2::ggplot_build(p))
  expect_error(plotValues(list(.fitFrame(1))), "name")
})

test_that("plotArray takes a positional subset and needs no covtable", {
  times <- 0:5
  x <- Xd(data.frame(name = "A", time = times, row.names = paste0("A", times)))
  pars <- structure(exp(-times / 2), names = getParameters(x))
  profs <- profile(constraintL2(mu = pars, sigma = 0.1), pars, whichPar = "A2",
                   limits = c(-1, 1), cores = 1)
  p <- plotArray("A2", profs, x, seq(0, 5, 0.1), time < 2)
  expect_lt(max(p$data$time), 2)
})

test_that("plotArray draws nSim trajectories from the optimum, also for small values", {
  times <- 0:5
  x <- Xd(data.frame(name = "A", time = times, row.names = paste0("A", times)))
  pars <- structure(exp(-times / 2) / 1000, names = getParameters(x))
  profs <- suppressWarnings(profile(constraintL2(mu = pars, sigma = 1e-4), pars,
                                    whichPar = "A2", limits = c(-2e-4, 2e-4), cores = 1))
  nUp <- nrow(profs) - which.min(abs(profs$constraint)) + 1L
  expect_gt(nUp, 2L)
  for (n in c(1, 2, 3, 100))
    expect_equal(length(unique(plotArray("A2", profs, x, times, nSim = n)$data$ParValue)),
                 min(n, nUp))
  d <- plotArray("A2", profs, x, times, direction = "down")$data
  expect_true(any(d$ParValue == pars[["A2"]]))
})

test_that("plotPaths draws relative paths and names an unprofiled parameter", {
  pr <- .coupledProfile()
  d <- as.data.frame(plotPaths(pr)$data)
  expect_setequal(unique(d$combination), paste("a", c("b", "c", "d"), sep = " - \n"))
  origin <- which.min(abs(pr$constraint))
  expect_true(all(abs(d$y[d$x == 0]) < 1e-12))
  expect_equal(nrow(plotPaths(pr, sort = TRUE)$data), nrow(d))
  expect_error(plotPaths(pr, whichPar = "b"), "not profiled: b")
})

test_that("coordTransform transforms values and sigmas, whatever the other columns", {
  d <- data.frame(name = c("A", "A", "B"), time = 0:2, value = c(1, 4, 9), sigma = c(1, 2, NA),
                  n = 99)
  out <- coordTransform(d, list(A = "log(x)", B = "abs(x)"))
  expect_equal(nrow(out), 3L)
  expect_equal(out$value, c(log(1), log(4), 9))
  expect_equal(out$sigma[1:2], c(1, 2) / c(1, 4))
})

test_that("plotProfilesAndPaths draws the profiles of several parameters", {
  skip_if_not_installed("cowplot")
  pars <- c(a = 1, b = 0.5, c = 2)
  profs <- profile(constraintL2(mu = pars, sigma = 0.1), pars, whichPar = c("a", "b"),
                   limits = c(-1, 1), cores = 1)
  expect_no_error(plotProfilesAndPaths(profs, c("a", "b"), nPars = 2, ncol = 2))
})

test_that("coordTransform keeps the order of the rows", {
  d <- data.frame(name = c("A", "B", "A", "B"), time = 0:3, value = 1:4)
  out <- coordTransform(d, "log(x)")
  expect_equal(out$time, 0:3)
  expect_equal(out$value, log(1:4))
})

test_that("plotProfilesAndPaths draws the total besides the contributions", {
  skip_if_not_installed("cowplot")
  pars <- c(a = 1, b = 0.5)
  profs <- profile(constraintL2(mu = pars, sigma = 0.1), pars, whichPar = "a",
                   limits = c(-1, 1), cores = 1)
  panels <- list()
  # the first call aligns the profile panel with the paths
  local_mocked_bindings(align_plots = function(...) {
    if (!length(panels)) panels <<- list(...)
    list(...)
  }, .package = "cowplot")
  plotProfilesAndPaths(profs, "a", nPars = 2)
  expect_true("total" %in% panels[[1]]$data$mode)
})

test_that("plotProfilesAndPaths draws the thresholds and the cut-off it is given", {
  skip_if_not_installed("cowplot")
  pars <- c(a = 1, b = 0.5)
  profs <- profile(constraintL2(mu = pars, sigma = 0.1), pars, whichPar = "a",
                   limits = c(-1, 1), cores = 1)
  panels <- list()
  local_mocked_bindings(align_plots = function(...) {
    if (!length(panels)) panels <<- list(...)
    list(...)
  }, .package = "cowplot")
  plotProfilesAndPaths(profs, "a", nPars = 2, threshold = c(F = 4.5), maxvalue = 8)
  y <- panels[[1]]$scales$get_scales("y")
  expect_equal(unname(y$breaks), c(0, 4.5))
  expect_identical(y$labels, c("optimum", "F"))
  expect_equal(y$limits, c(NA, 8))
})

test_that("the parlist overlay of plotProfile draws parameters only", {
  pars <- c(a = 1, b = 0.5)
  profs <- profile(constraintL2(mu = pars, sigma = 0.1), pars, whichPar = c("a", "b"),
                   limits = c(-1, 1), cores = 1)
  pf <- parframe(data.frame(value = 0, converged = TRUE, iterations = 3L,
                            stopReason = "converged", neval = 4L, a = 1, b = 0.5),
                 parameters = c("a", "b"),
                 metanames = c("value", "converged", "iterations", "stopReason", "neval"))
  p <- plotProfile(profs, parlist = pf)
  pts <- p$layers[[length(p$layers)]]$data
  expect_setequal(unique(pts$name), c("a", "b"))
  lim <- "a"
  f <- function() { lim <- "b"; plotProfile(profs, name == lim, mode == "total") }
  expect_equal(unique(as.character(attr(f(), "data")$name)), "b")
})

test_that("plot() of a prdframe takes a data frame", {
  times <- 0:5
  grid <- data.frame(name = "A", time = times, row.names = paste0("A", times))
  x <- Xd(grid, condition = "C1")
  pred <- x(times, structure(exp(-times / 2), names = getParameters(x)))
  data <- data.frame(name = "A", time = times, value = exp(-times / 2), sigma = 0.1)
  expect_s3_class(plot(pred[[1]], data), "ggplot")
})

# A fit of spike-and-slab random effects k and A over three groups, as em() returns it, without
# the model: A of group G3 reaches no data.
.slabFit <- function(method = "focei", shift = 0) {
  groups <- paste0("G", 1:3)
  random <- data.frame(group = rep(groups, 2), eta = rep(c("k", "A"), each = 3),
                       parameter = paste0("eta_", rep(c("k", "A"), each = 3), "_", groups),
                       mode = c(0, 0.8, -0.5, 0, 0, 0), se = c(0.05, 0.1, 0.2, 0.3, 0.02, 0.5),
                       data = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE), cluster = c(0, 1, 2, 0, 0, 0),
                       p0 = c(0.9, 0.01, 0.2, 0.6, 0.95, 0.7), stringsAsFactors = FALSE)
  modes <- matrix(random$mode, 3, dimnames = list(groups, c("eta_k", "eta_A")))
  inner <- modes
  inner[, "eta_A"] <- c(0.2, -0.1, 0)
  argument <- c(k = -0.5 + shift, A = 2.3, prob_ID = 0.3, tau_ID = 0.7)
  information <- diag(c(50, 40, 30, 20)) + 1
  dimnames(information) <- list(names(argument), names(argument))
  structure(list(argument = argument, value = -40, gradient = 0 * argument, method = method,
                 information = information, etaModes = modes, innerModes = inner,
                 random = random,
                 slab = data.frame(name = c("prob_ID", "tau_ID"), distribution = "spikeSlab", symbols = "k, A",
                                   estimate = c(0.3, 0.7), fixed = FALSE, stringsAsFactors = FALSE),
                 stageTrace = data.frame(iter = 1:3, phase = "trust", OFV = c(-30, -39, -40),
                                         deltaPar = c(NA, 0.2, 0.01), prob_ID = c(0.5, 0.35, 0.3),
                                         tau_ID = c(0.3, 0.6, 0.7)),
                 multistart = data.frame(start = 1:4, value = c(-40, -40, -35, NA), converged = TRUE,
                                         iterations = 10L)),
            class = c("em", "list"))
}

test_that("the inclusion heatmap shows 1 - p0 with the reported values of the deviations", {
  fit <- .slabFit()
  d <- plotRandomEffects(fit, type = "inclusion")$data
  expect_equal(d$inclusion, 1 - fit$random$p0)
  expect_equal(d$label, c("", "0.80", "-0.50", "", "", "no data"))
  expect_equal(levels(d$eta), c("A", "k"))
})

test_that("conditional draws of spike-and-slab effects: spike mass p0, the conditional variance, the prior", {
  fit <- .slabFit()
  set.seed(1)
  v <- dMod2:::.spikeSlabValues(fit, "conditional", 20000L)
  # the effect without data stays out
  expect_equal(nrow(v), 5 * 20000)
  r <- fit$random[fit$random$data, ]
  for (i in seq_len(nrow(r))) {
    x <- v$value[(i - 1) * 20000 + 1:20000]
    expect_equal(mean(x == 0), r$p0[i], tolerance = 0.02)
    expect_equal(sd(x), r$se[i], tolerance = 0.05)
  }
  expect_equal(unique(v$spike), 0.7)
  expect_equal(unique(v$sd), 0.7)
  # histogram bars have area 1, the exact zeros the mean p0; the prior spike has area 1 - prob
  p <- plotRandomEffects(fit, draws = 2000L)
  bars <- p$data
  for (e in c("k", "A")) {
    b <- bars[bars$eta == e, ]
    area <- (b$xmax - b$xmin) * (b$ymax - b$ymin)
    expect_equal(sum(area), 1)
    expect_equal(sum(area[b$part == "spike"]), mean(r$p0[r$eta == e]), tolerance = 0.03)
  }
  spike <- p$layers[[3]]$data
  expect_equal(spike$width * spike$density, c(0.7, 0.7))
  # with only sparse effects there is no QQ plot
  expect_s3_class(p, "ggplot")
  expect_s3_class(plotRandomEffects(fit, type = "modes"), "ggplot")
})

test_that("plotTrace draws the hyperparameters of a sparse fit", {
  d <- plotTrace(.slabFit())$data
  expect_equal(levels(d$panel), c("OFV", "|delta par|", "prob_ID", "tau_ID"))
  expect_true(all(is.finite(d$value)))
})

test_that("plotTrace leaves out iterations without a value", {
  fit <- .slabFit()
  fit$stageTrace$deltaPar[1] <- NA
  d <- plotTrace(fit)$data
  expect_equal(sum(d$panel == "|delta par|"), nrow(fit$stageTrace) - 1L)
  expect_true(all(is.finite(d$value)))
})

test_that("plotPars draws the estimates of several em fits side by side with Wald intervals", {
  fits <- list(focei = .slabFit(), saem = .slabFit("saem", shift = 0.1))
  fits$saem$information <- NULL
  p <- plotPars(fits, which = c("k", "A", "prob_ID"))
  d <- p$data
  expect_equal(nrow(d), 6)
  expect_equal(levels(d$fit), c("focei", "saem"))
  expect_equal(d$estimate[d$fit == "saem" & d$parameter == "k"], -0.4)
  # Wald intervals where the fit has standard errors only
  se <- sqrt(diag(solve(fits$focei$information)))[1:3]
  expect_equal(d$upper[d$fit == "focei"] - d$estimate[d$fit == "focei"], unname(qnorm(0.975) * se))
  expect_true(all(is.na(d$lower[d$fit == "saem"])))
  # the first fit of a row on top
  expect_gt(d$position[1], d$position[4])
  # unnamed fits are labeled by their method
  expect_equal(levels(plotPars(unname(fits))$data$fit), c("focei", "saem"))
  expect_error(plotPars(list(a = 1)), "em fits")
})

test_that("plotValues draws the multistart of an em fit and of a list of them", {
  fit <- .slabFit()
  expect_equal(plotValues(fit)$data$value, c(-40, -40, -35))
  d <- plotValues(list(a = fit, b = fit))$data
  expect_equal(as.vector(table(d$set)), c(3, 3))
})
