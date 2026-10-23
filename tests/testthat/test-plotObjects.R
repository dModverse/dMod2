# The ggplot objects of the prediction, data, flux, residual and path plots: layers, mappings,
# facets, scales and the plotted data. No compiled code.

pdf(NULL)

.grid2 <- function(times = 0:3) {
  g <- data.frame(name = rep(c("A", "B"), each = length(times)), time = rep(times, 2))
  rownames(g) <- paste0(g$name, g$time)
  g
}
.x2 <- function(conditions = "C1") Reduce(`+`, lapply(conditions, function(C) Xd(.grid2(), C)))
.pars2 <- c(A0 = 1, A1 = 0.5, A2 = 0.25, A3 = 0.125, B0 = 0, B1 = 0.3, B2 = 0.6, B3 = 0.9)
.data2 <- function(conditions = "C1") as.datalist(do.call(rbind, lapply(conditions, function(C)
  data.frame(name = "A", time = 0:3, value = c(1, 0.6, 0.2, 0.1), sigma = 0.1, condition = C,
             dose = match(C, c("C1", "C2"))))),
  splitBy = "condition", keepCovariates = "dose")

.geoms <- function(p) unname(vapply(p$layers, function(l) class(l$geom)[1], ""))
.map <- function(m) vapply(m, rlang::as_label, "")
.facetVars <- function(p) {
  f <- p$facet$params
  list(rows = names(c(f$facets, f$rows)), cols = names(f$cols))
}
.scale <- function(p, aes) p$scales$get_scales(aes)

# sigma from a constant error model, one prdframe per condition
.errfn <- function(out, pars, conditions = NULL, ...) {
  s <- out
  s[, -1] <- 0.2
  setNames(list(prdframe(s, parameters = pars)), conditions)
}

test_that("plotCombined draws prediction lines and data points with error bars", {
  pred <- .x2(c("C1", "C2"))(0:3, .pars2)
  p <- plotCombined(pred, .data2(c("C1", "C2")))
  expect_s3_class(p, "ggplot")
  expect_equal(.geoms(p), c("GeomLine", "GeomPoint", "GeomErrorbar"))
  expect_equal(.map(p$mapping), c(x = "time", y = "value", ymin = "value - sigma",
                                  ymax = "value + sigma", group = "condition",
                                  colour = "condition"))
  expect_equal(.map(p$layers[[2]]$mapping), c(shape = "bloq"))
  expect_equal(p$layers[[3]]$geom_params$width, 0)
  expect_s3_class(p$facet, "FacetWrap")
  expect_equal(.facetVars(p)$rows, "name")
  expect_equal(p$facet$params$free, list(x = TRUE, y = TRUE))
  expect_equal(.scale(p, "shape")$name, "BLoQ")
  expect_equal(p$guides$guides$shape, "none")

  d <- p$data
  expect_equal(names(d), c("time", "name", "value", "sigma", "condition", "dose"))
  expect_equal(nrow(d), 2 * 8 + 2 * 4)
  expect_true(all(is.na(d$sigma[1:16])))
  expect_equal(d$value[d$name == "B" & d$condition == "C2"], c(0, 0.3, 0.6, 0.9))
  expect_equal(d$sigma[17:24], rep(0.1, 8))
  expect_equal(d$dose[17:24], rep(1:2, each = 4))
  expect_equal(p$layers[[1]]$data$value, d$value[1:16])
  expect_equal(p$layers[[2]]$data$bloq, rep("no", 8))
})

test_that("plotCombined facets by grid and by name and condition", {
  pred <- .x2()(0:3, .pars2)
  p <- plotCombined(pred, .data2(), facet = "grid", scales = "free_y")
  expect_s3_class(p$facet, "FacetGrid")
  expect_equal(.facetVars(p), list(rows = "name", cols = "condition"))
  expect_equal(p$facet$params$free, list(x = FALSE, y = TRUE))
  expect_equal(.map(p$mapping), c(x = "time", y = "value", ymin = "value - sigma",
                                  ymax = "value + sigma"))
  p <- plotCombined(pred, .data2(), facet = "wrap_plain")
  expect_s3_class(p$facet, "FacetWrap")
  expect_equal(.facetVars(p)$rows, c("name", "condition"))
  expect_error(plotCombined(pred, facet = "none"), "should be one of")
})

test_that("plotCombined subsets, transforms and takes extra aesthetics", {
  pred <- .x2()(0:3, .pars2)
  p <- plotCombined(pred, .data2(), name == "A", time > 0,
                    transform = list(A = "log(x)"),
                    aesthetics = list(linetype = "name", color = "dose"))
  d <- p$data
  expect_equal(unique(as.character(d$name)), "A")
  expect_equal(d$value, log(c(0.5, 0.25, 0.125, 0.6, 0.2, 0.1)))
  expect_equal(d$sigma, c(NA, NA, NA, 0.1 / c(0.6, 0.2, 0.1)))
  expect_mapequal(.map(p$mapping), c(x = "time", y = "value", ymin = "value - sigma",
                                     ymax = "value + sigma", group = "condition",
                                     linetype = "name", colour = "dose"))
  # an aesthetic replaces the default of the same name, in either spelling
  for (a in list(list(colour = "dose"), list(color = "dose"))) {
    m <- .map(plotCombined(pred, aesthetics = a)$mapping)
    expect_equal(unname(m[names(m) == "colour"]), "dose")
  }
})

test_that("plotCombined marks data below the limit of quantification", {
  pred <- .x2()(0:3, .pars2)
  data <- .data2()
  data[[1]]$lloq <- 0.15
  p <- plotCombined(pred, data)
  expect_equal(p$layers[[2]]$data$bloq, c("no", "no", "no", "yes"))
  expect_null(p$guides$guides$shape)
  expect_equal(.scale(p, "shape")$palette(2), c(yes = 4, no = 19))
})

test_that("plotCombined without data draws the prediction only", {
  p <- plotCombined(.x2()(0:3, .pars2))
  expect_equal(.geoms(p), "GeomLine")
  expect_equal(names(p$data), c("time", "name", "value", "sigma", "condition"))
})

test_that("plot() of a prdlist names unnamed predictions and data", {
  pred <- .x2(c("C1", "C2"))(0:3, .pars2)
  data <- .data2(c("C1", "C2"))
  attr(data, "names") <- NULL
  p <- plot(unname(pred), data)
  expect_equal(unique(p$data$condition), c("C1", "C2"))
  expect_equal(p$layers[[2]]$data$condition, rep(c("C1", "C2"), each = 4))
  p <- plot(pred, .data2(c("C1", "C2")), facet = "grid", scales = "fixed")
  expect_s3_class(p$facet, "FacetGrid")
  expect_equal(p$data$dose, c(rep(1:2, each = 8), rep(1:2, each = 4)))
})

test_that("plot() of a prdframe names it C1, also its data frame", {
  pred <- .x2()(0:3, .pars2)[[1]]
  d <- data.frame(name = "A", time = 0:3, value = 1, sigma = 0.1)
  p <- plot(pred, d, time < 3, facet = "wrap_plain")
  expect_equal(unique(p$data$condition), "C1")
  expect_equal(max(p$data$time), 2)
  expect_equal(.facetVars(p)$rows, c("name", "condition"))
  expect_equal(.geoms(plot(pred)), "GeomLine")
})

test_that("plotPrediction draws one line per condition and keeps its data", {
  pred <- .x2(c("C1", "C2"))(0:3, .pars2)
  p <- plotPrediction(pred)
  expect_equal(.geoms(p), "GeomLine")
  expect_mapequal(.map(p$mapping), c(x = "time", y = "value", group = "condition",
                                     colour = "condition"))
  expect_s3_class(p$facet, "FacetWrap")
  expect_equal(.facetVars(p)$rows, "name")
  d <- attr(p, "data")
  expect_equal(names(d), c("time", "name", "value", "condition", "sigma"))
  expect_true(all(is.nan(d$sigma)))
  expect_equal(nrow(d), 16L)

  p <- plotPrediction(pred, condition == "C2", name == "A", facet = "grid", scales = "free_x",
                      transform = list(A = "2*x"))
  expect_s3_class(p$facet, "FacetGrid")
  expect_equal(.facetVars(p), list(rows = "name", cols = "condition"))
  expect_equal(p$facet$params$free, list(x = TRUE, y = FALSE))
  expect_equal(.map(p$mapping), c(x = "time", y = "value"))
  expect_equal(attr(p, "data")$value, 2 * c(1, 0.5, 0.25, 0.125))
})

test_that("plotPrediction draws the sigma of an error model as a ribbon", {
  pred <- .x2()(0:3, .pars2)
  p <- plotPrediction(pred, errfn = .errfn)
  expect_equal(.geoms(p), c("GeomRibbon", "GeomLine"))
  expect_equal(.map(p$layers[[1]]$mapping), c(ymin = "value - sigma", ymax = "value + sigma",
                                              fill = "condition"))
  expect_equal(p$layers[[1]]$aes_params, list(linetype = 0, alpha = 0.3))
  expect_equal(attr(p, "data")$sigma, rep(0.2, 8))
})

test_that("plotData of a data frame splits it into conditions", {
  d <- data.frame(name = "A", time = 0:3, value = 1:4, sigma = 0.1,
                  condition = rep(c("u", "v"), each = 2))
  p <- plotData(d, time > 0)
  expect_equal(.geoms(p), c("GeomPoint", "GeomErrorbar"))
  expect_equal(p$data$condition, c("u", "v", "v"))
  expect_equal(p$data$value, 2:4)
})

test_that("plotFluxes stacks the fluxes of every condition", {
  x <- .x2(c("C1", "C2"))
  p <- plotFluxes(.pars2, x, 0:3, c(prod = "A0", deg = "0.5*A*B", "Heaviside(time - 1)"))
  expect_equal(.geoms(p), "GeomArea")
  expect_s3_class(p$layers[[1]]$position, "PositionStack")
  expect_equal(p$layers[[1]]$aes_params, list(alpha = 0.3, colour = "darkgrey", linewidth = 0.4))
  expect_equal(.map(p$mapping), c(x = "time", y = "value", group = "name", fill = "name"))
  expect_equal(.facetVars(p)$rows, "condition")
  expect_equal(.scale(p, "fill")$name, "Fluxes:")
  expect_equal(p$labels[c("x", "y")], list(x = "time", y = "flux contribution"))

  out <- attr(p, "out")
  expect_equal(names(out), c("time", "name", "value", "condition"))
  expect_equal(levels(out$name), c("prod", "deg", "Heaviside(time - 1)"))
  expect_equal(out$condition, rep(c("C1", "C2"), each = 12))
  A <- c(1, 0.5, 0.25, 0.125); B <- c(0, 0.3, 0.6, 0.9)
  expect_equal(out$value[1:12], c(rep(1, 4), 0.5 * A * B, c(0, 0.5, 1, 1)))
  expect_identical(p$data, out)
})

test_that("plotFluxes names unnamed predictions and passes its arguments to x", {
  x <- .x2(c("C1", "C2"))
  p <- plotFluxes(.pars2, function(...) unname(x(...)), 0:1,
                  list("piecewise(1, time < 1, 2)", "exp10(B1)"), legendTitle = "F")
  out <- attr(p, "out")
  expect_equal(unique(out$condition), c("C1", "C2"))
  expect_equal(levels(out$name), c("piecewise(1, time < 1, 2)", "exp10(B1)"))
  expect_equal(out$value[1:4], c(1, 2, 10^0.3, 10^0.3))
  expect_equal(.scale(p, "fill")$name, "F")
  p <- plotFluxes(.pars2, x, 0:1, c(f = "A"), conditions = "C2")
  expect_equal(unique(attr(p, "out")$condition), "C2")
})

test_that("plotResiduals sums squared weighted residuals per split", {
  x <- .x2(c("C1", "C2"))
  data <- .data2(c("C1", "C2"))
  pf <- parframe(rbind(data.frame(value = 1, as.list(.pars2)),
                       data.frame(value = 2, as.list(.pars2 * 2))),
                 parameters = names(.pars2), metanames = "value")
  r0 <- ((c(1, 0.5, 0.25, 0.125) - c(1, 0.6, 0.2, 0.1)) / 0.1)^2
  r1 <- ((2 * c(1, 0.5, 0.25, 0.125) - c(1, 0.6, 0.2, 0.1)) / 0.1)^2

  p <- plotResiduals(pf, x, data)
  out <- attr(p, "out")
  expect_equal(out, data.frame(condition = names(data), res = sum(r0) + sum(r1)))
  expect_equal(.geoms(p), c("GeomPoint", "GeomLine"))
  expect_equal(.map(p$mapping), c(x = "condition", y = "res", colour = "condition",
                                  group = "condition"))
  expect_s3_class(p$facet, "FacetNull")

  p <- plotResiduals(pf, x, data, split = c("time", "index", "condition"))
  out <- attr(p, "out")
  expect_equal(names(out), c("time", "index", "condition", "res"))
  expect_equal(out$index, rep(rep(c("1", "2"), each = 2), 4))
  expect_equal(out$res[out$condition == "C1"], c(rbind(r0, r1)))
  expect_equal(.map(p$mapping), c(x = "time", y = "res", colour = "index", group = "index"))
  expect_s3_class(p$facet, "FacetWrap")
  expect_equal(.facetVars(p)$rows, "condition")
})

test_that("plotResiduals adds log(sigma^2) of an error model and uses the index column", {
  x <- .x2()
  data <- .data2()
  data[[1]]$sigma[2:4] <- NA
  pf <- parframe(data.frame(value = 1, index = 7L, as.list(.pars2)),
                 parameters = names(.pars2), metanames = c("value", "index"))
  s <- c(0.1, 0.2, 0.2, 0.2)
  r <- ((c(1, 0.5, 0.25, 0.125) - c(1, 0.6, 0.2, 0.1)) / s)^2 + log(s^2)
  out <- attr(plotResiduals(pf, x, data, split = c("name", "index"), errmodel = .errfn), "out")
  expect_equal(out$name, "A")
  expect_equal(out$index, "7")
  expect_equal(out$res, sum(r))
})

.profA2 <- function(scale = 1) {
  times <- 0:5
  x <- Xd(data.frame(name = "A", time = times, row.names = paste0("A", times)))
  pars <- structure(exp(-times / 2), names = getParameters(x))
  list(x = x, pars = pars,
       profs = profile(constraintL2(mu = pars, sigma = 0.1), pars, whichPar = "A2",
                       limits = c(-1, 1), cores = 1))
}

test_that("plotArray colors trajectories along the profile by the parameter value", {
  pr <- .profA2()
  p <- plotArray("A2", pr$profs, pr$x, c(0, 2, 4), nSim = 3)
  expect_equal(.geoms(p), "GeomLine")
  expect_equal(.map(p$mapping), c(x = "time", y = "value", group = "set", colour = "ParValue"))
  expect_s3_class(p$facet, "FacetGrid")
  expect_equal(.facetVars(p), list(rows = "name", cols = "condition"))
  expect_equal(p$facet$params$free, list(x = FALSE, y = TRUE))
  expect_equal(.scale(p, "colour")$name, "A2")
  expect_equal(p$labels[c("x", "y")], list(x = "time", y = "value"))
  d <- p$data
  expect_equal(names(d), c("time", "name", "value", "condition", "set", "ParValue"))
  expect_equal(d$set, rep(1:3, each = 3))
  expect_equal(d$condition, rep("1", 9))
  expect_equal(d$ParValue[1], pr$pars[["A2"]])
  expect_true(all(diff(unique(d$ParValue)) > 0))
  expect_equal(d$value[d$time == 2], unique(d$ParValue))
  down <- plotArray("A2", pr$profs, pr$x, 4, direction = "down", nSim = 2)$data
  expect_true(all(diff(down$ParValue) < 0))
  expect_error(plotArray("A9", pr$profs, pr$x, 0:5), "'A9' is not profiled")
})

test_that("plotArray merges a covariate table and warns about a failing set", {
  pr <- .profA2()
  cov <- data.frame(dose = 3, row.names = "1")
  p <- plotArray("A2", pr$profs, pr$x, 0:1, dose == 3, covtable = cov, nSim = 2)
  expect_equal(names(p$data), c("condition", "time", "name", "value", "set", "ParValue", "dose"))
  expect_equal(nrow(p$data), 4L)
  calls <- 0
  prd <- function(...) {
    calls <<- calls + 1
    if (calls == 2) stop("boom")
    pr$x(...)
  }
  expect_warning(p <- plotArray("A2", pr$profs, prd, 0:1, nSim = 3),
                 "parameter set 2 failed: boom")
  expect_equal(unique(p$data$set), c(1, 3))
})

test_that("plotPaths draws the paths of each pair, relative or absolute", {
  pars <- c(a = 1, b = 0.5)
  profs <- profile(constraintL2(mu = pars, sigma = 0.1), pars, whichPar = c("a", "b"),
                   limits = c(-1, 1), cores = 1)
  p <- plotPaths(profs)
  expect_equal(.geoms(p), "GeomPath")
  expect_equal(.map(p$mapping), c(x = "x", y = "y", group = "interaction(name, proflist)",
                                  colour = "name", linetype = "proflist"))
  expect_s3_class(p$facet, "FacetWrap")
  expect_equal(.facetVars(p)$rows, "combination")
  expect_equal(p$facet$params$free, list(x = FALSE, y = FALSE))
  expect_equal(p$labels$x, expression(paste(Delta, "parameter 1")))
  expect_equal(p$labels$y, expression(paste(Delta, "parameter 2")))
  expect_equal(.scale(p, "linetype")$name, "profile\nlist")
  expect_equal(.scale(p, "colour")$name, "profiled\nparameter")
  d <- attr(p, "data")
  expect_identical(d, p$data)
  expect_equal(levels(factor(d$combination)), c("a - \nb", "b - \na"))

  p <- plotPaths(profs, name == "a", relative = FALSE, scales = "free")
  expect_equal(p$labels[c("x", "y")], list(x = "parameter 1", y = "parameter 2"))
  expect_equal(p$facet$params$free, list(x = TRUE, y = TRUE))
  d <- attr(p, "data")
  expect_equal(unique(d$name), "a")
  expect_setequal(unique(d$combination), c("a - \nb", "b - \na"))
  expect_equal(range(d$x[d$combination == "a - \nb"]), range(profs$a[profs$whichPar == "a"]))
})

test_that("predict() returns one block per parameter row with the data as attribute", {
  x <- .x2(c("C1", "C2"))
  pf <- parframe(rbind(data.frame(value = 3, as.list(.pars2)),
                       data.frame(value = 4, as.list(.pars2 * 2))),
                 parameters = names(.pars2), metanames = "value")
  r <- predict(x, times = c(0, 3), pars = pf, data = .data2(c("C1", "C2")))
  expect_equal(names(r), c("time", "name", "value", "condition", "sigma", "dose", ".value"))
  expect_equal(nrow(r), 2 * 2 * 2 * 2)
  expect_equal(r$value[r$.value == 4 & r$name == "B" & r$condition == "C2"], c(0, 1.8))
  expect_equal(r$dose, ifelse(r$condition == "C1", 1, 2))
  d <- attr(r, "data")
  expect_equal(nrow(d), 8L)
  expect_true(all(c(".value", "lloq", "dose") %in% names(d)))
  expect_true(all(is.na(d$.value)))
  expect_equal(rownames(d), as.character(17:24))

  r <- predict(x, times = 0, pars = pf, conditions = "C2", data = .data2(c("C1", "C2")))
  expect_equal(unique(r$condition), "C2")
  expect_equal(unique(attr(r, "data")$condition), "C2")

  r <- predict(x, times = 0, pars = pf[1, ])
  expect_equal(attr(r, "data"), data.frame())
  expect_equal(r$.value, rep(3, 4))
})
