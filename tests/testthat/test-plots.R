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

test_that("a layer added to plotCombined inherits its data", {
  x <- Xd(data.frame(name = "A", time = 0:2, row.names = paste0("A", 0:2)), condition = "C1")
  pars <- c(A0 = 1, A1 = 0.5, A2 = 0.2)
  data <- as.datalist(data.frame(name = "A", time = 0:2, value = pars + 0.01, sigma = 0.1,
                                 condition = "C1"))
  p <- plotCombined(x(seq(0, 2, 0.5), pars), data) + ggplot2::geom_point(ggplot2::aes(y = value))
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plotResiduals labels each fit with the residuals of its own row", {
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

test_that("plotPathsMulti colors the nPars largest paths and puts the rest in Others", {
  skip_if_not_installed("cowplot")
  pars <- c(A0 = 1, A1 = 0.5, A2 = 0.3, A3 = 0.2, A4 = 0.1)
  profs <- profile(constraintL2(mu = pars, sigma = 0.1), pars, whichPar = "A1",
                   limits = c(-1, 1), cores = 1)
  d <- as.data.frame(plotPathsMulti(profs, "A1", nPars = 2)$data)
  expect_lte(length(unique(d$label)), 3L)
  expect_true("Others" %in% d$label)
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
  # the first call aligns the profile panel with its paths
  local_mocked_bindings(align_plots = function(...) {
    if (!length(panels)) panels <<- list(...)
    list(...)
  }, .package = "cowplot")
  plotProfilesAndPaths(profs, "a", nPars = 2)
  expect_true("total" %in% panels[[1]]$data$mode)
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
