# Observable behaviour of mstrust(), reduceReplicates(), fitErrorModel() and
# the result parsing of steadyStates(): printed output, files, routing, values.

quadObj <- function(pars, fixed = NULL, ...) {
  p <- c(pars, fixed)
  if (p[["a"]] > 1.5) stop("outside")
  r <- c(a = p[["a"]] - 1, b = 2 * (p[["b"]] + 0.5))
  list(value = sum(r^2), gradient = c(a = 2 * r[["a"]], b = 4 * r[["b"]]),
       hessian = matrix(c(2, 0, 0, 8), 2, 2, dimnames = list(c("a", "b"), c("a", "b"))))
}

stripRunDate <- function(x) gsub("trial-([0-9]+)-[0-9-]+", "trial-\\1-DATE", x)


# ---- mstrust ------------------------------------------------------------

test_that("mstrust prints the multi-start summary with stats = TRUE", {
  set.seed(4)
  expect_snapshot(fits <- mstrust(quadObj, c(a = 0, b = 0), fits = 6, sd = 2,
                                  retry = FALSE, stats = TRUE))
  expect_identical(sum(dMod2:::.statParlist(fits) == "error"), 3L)
  set.seed(4)
  expect_snapshot(invisible(mstrust(quadObj, c(a = 0, b = 0), fits = 2, sd = 0.5,
                                    iterlim = 1, stats = TRUE)))
  expect_silent(mstrust(quadObj, c(a = 0, b = 0), fits = 2, sd = 0.5))
})


test_that("mstrust writes the run folder, log and parameter list with output = TRUE", {
  d <- withr::local_tempdir()
  set.seed(4)
  fits <- mstrust(quadObj, c(a = 0, b = 0), fits = 4, sd = 2, retry = FALSE,
                  resultPath = d, name = "run", output = TRUE)
  files <- stripRunDate(list.files(d, recursive = TRUE, include.dirs = TRUE))
  expect_identical(sort(files), c("run", "run/trial-1-DATE", "run/trial-1-DATE/mstrust.log",
                                  "run/trial-1-DATE/parameterList.Rda"))
  log <- list.files(d, "mstrust.log", recursive = TRUE, full.names = TRUE)
  expect_snapshot(writeLines(readLines(log)))
  saved <- readRDS(list.files(d, "parameterList.Rda", recursive = TRUE, full.names = TRUE))
  expect_identical(unclass(saved)[1:4], unclass(fits)[1:4])

  # A second run in the same folder is the next trial; cautiousMode keeps the
  # deparsed fits in interRes/.
  fits2 <- mstrust(quadObj, c(a = 0, b = 0), fits = 2, sd = 0.1,
                   resultPath = d, name = "run", output = TRUE, cautiousMode = TRUE)
  files <- stripRunDate(list.files(d, recursive = TRUE))
  expect_identical(sort(files[grepl("trial-2", files)]),
                   c("run/trial-2-DATE/interRes/fit-1.R", "run/trial-2-DATE/interRes/fit-2.R",
                     "run/trial-2-DATE/mstrust.log", "run/trial-2-DATE/parameterList.Rda"))
  fit1 <- dget(list.files(d, "fit-1.R$", recursive = TRUE, full.names = TRUE))
  expect_identical(names(fit1), c("value", "argument", "iterations", "converged", "elapsed"))
  expect_identical(fit1$argument, fits2[[1]]$argument)

  # cautiousMode alone keeps an empty interRes/ and a log without fits
  d2 <- withr::local_tempdir()
  invisible(mstrust(quadObj, c(a = 0, b = 0), fits = 1, sd = 0.1, resultPath = d2,
                    cautiousMode = TRUE))
  expect_identical(sort(stripRunDate(list.files(d2, recursive = TRUE, include.dirs = TRUE))),
                   c("mstrust", "mstrust/trial-1-DATE", "mstrust/trial-1-DATE/interRes",
                     "mstrust/trial-1-DATE/mstrust.log"))
})


test_that("mstrust routes its arguments to the optimiser, the sampler and the objective", {
  seen <- new.env()
  myopt <- function(objfun, parinit, iterlim = 1, ...) {
    seen$iterlim <- iterlim
    seen$dots <- names(list(...))
    list(argument = parinit, value = objfun(parinit, ...)$value, converged = TRUE,
         iterations = 0L)
  }
  mysample <- function(n, spread) seq_len(n) * spread
  obj <- function(pars, scale = 1, ...) list(value = scale * sum(pars^2))
  fits <- mstrust(obj, c(a = 0, b = 0), fits = 2, optmethod = myopt,
                  samplefun = mysample, spread = 0.5, iterlim = 7, scale = 10)
  expect_identical(seen$iterlim, 7)
  expect_identical(sort(seen$dots), c("rinit", "rmax", "scale"))
  expect_equal(fits[[1]]$parinit, c(a = 0.5, b = 1))
  expect_equal(fits[[2]]$value, 12.5)
  run <- attr(fits, "run")
  expect_identical(run$optimiser, "myopt")
  expect_identical(run$sampler, "mysample")
  expect_identical(run$optimiserArgs, list(iterlim = 7))
  expect_identical(run$samplerArgs, list(n = 2L, spread = 0.5))
  expect_identical(run$fits, 2)
  expect_equal(run$cores, c(fits = 1))

  # startFromCenter puts the first fit at the centre
  fits <- mstrust(obj, c(a = 0, b = 0), fits = 2, optmethod = myopt,
                  samplefun = mysample, spread = 0.5, startFromCenter = TRUE)
  expect_equal(fits[[1]]$parinit, c(a = 0, b = 0))
  expect_equal(fits[[2]]$parinit, c(a = 0.5, b = 1))

  # a parframe centre gives one fit per row
  starts <- msParframe(c(a = 0, b = 0), n = 3, sd = 1)
  fits <- mstrust(obj, starts, fits = 10, optmethod = myopt)
  expect_length(fits, 3L)
  expect_equal(fits[[3]]$parinit, as.parvec(starts, 3))
  expect_error(mstrust(obj, c(a = 0, b = 0), fits = 0, optmethod = myopt),
               "fits must be a positive number")
})


# ---- reduceReplicates ---------------------------------------------------

replicateData <- function() {
  data.frame(name = rep(c("A", "B"), each = 6),
             time = rep(c(0, 0, 0, 1, 1, 2), 2),
             value = c(1, 1.2, 0.9, 2, 2.4, 3, 5, 5.5, 4.5, 6, 7, 8),
             sigma = c(0.1, 0.2, 0.1, 0.2, 0.3, 0.1, 0.5, 0.4, 0.5, 0.3, 0.2, 0.6),
             condition = "ctrl",
             dose = rep(c(1, 2), each = 6),
             batch = c(1, 2, 3, 1, 2, 1, 1, 2, 3, 1, 2, 1),
             stringsAsFactors = FALSE)
}

# Columns, row names and class, independent of the order of the attributes.
pinFrame <- function(x) dput(list(columns = as.list(x), rows = attr(x, "row.names"), class = class(x)))

test_that("reduceReplicates gives mean, standard error and the stable columns", {
  d <- replicateData()
  expect_snapshot(r <- reduceReplicates(d))
  expect_snapshot(pinFrame(r))
  expect_snapshot(pinFrame(suppressMessages(reduceReplicates(d, select = c("condition", "dose")))))
  expect_snapshot(pinFrame(suppressMessages(reduceReplicates(d, keep = c("batch", "nothere")))))
  expect_snapshot(pinFrame(suppressMessages(reduceReplicates(d, datatrans = "log(x)"))))
  expect_snapshot(pinFrame(suppressMessages(reduceReplicates(d, weighted = TRUE))))
  expect_snapshot(pinFrame(suppressMessages(reduceReplicates(d[, -7]))))
})

test_that("reduceReplicates checks its input", {
  d <- replicateData()
  expect_error(reduceReplicates(d[, -1]),
               "Mandatory column names are: name, time, value, condition")
  expect_error(reduceReplicates(d[, -4], weighted = TRUE), "sigma")
})


# ---- fitErrorModel ------------------------------------------------------

errorModelData <- function() {
  set.seed(11)
  mu <- rep(c(0.1, 1, 10, 100), each = 5)
  d <- data.frame(name = rep(c("y", "z"), each = 20), time = rep(rep(1:4, each = 5), 2),
                  condition = "ctrl",
                  value = c(mu, 2 * mu) + rnorm(40, sd = sqrt(0.01 + 0.04 * c(mu, 2 * mu)^2)))
  suppressMessages(reduceReplicates(d))
}

roundNumeric <- function(df) {
  df[] <- lapply(df, function(x) if (is.numeric(x)) signif(x, 4) else x)
  df
}

test_that("fitErrorModel fits per group and returns the standard error", {
  skip_if_not_installed("optimx")
  red <- errorModelData()
  out <- fitErrorModel(red, factors = "name", par = c(s0 = -4, srel = -3))
  expect_identical(names(out), names(red))
  expect_identical(out[setdiff(names(red), "sigma")], red[setdiff(names(red), "sigma")])
  expect_snapshot(print(roundNumeric(out)))

  full <- fitErrorModel(red, factors = "name", par = c(s0 = -4, srel = -3), blather = TRUE)
  expect_identical(names(full), c(names(red), "s0", "srel", "cbLower68", "cbUpper68",
                                  "cbLower95", "cbUpper95", "condidnt", "sigmaLS"))
  expect_identical(attr(full, "errorModel"), "exp(s0)+exp(srel)*x^2")
  expect_identical(full$sigma, out$sigma)
  expect_identical(full$sigmaLS, red$sigma)
  expect_type(full$condidnt, "character")
  expect_snapshot(print(roundNumeric(full)))

  # bounds, named or positional; two factors give one group per combination
  bnd <- fitErrorModel(red, factors = "name", par = c(s0 = -4, srel = -1.5),
                       lower = c(srel = -2), blather = TRUE)
  expect_equal(unique(bnd$srel), -2)
  pos <- fitErrorModel(red, factors = "name", par = c(s0 = -4, srel = -1.5),
                       lower = c(-Inf, -2), blather = TRUE)
  expect_equal(pos$srel, bnd$srel)
  two <- fitErrorModel(red, factors = c("name", "condition"), par = c(s0 = -4, srel = -3),
                       blather = TRUE)
  expect_identical(unique(as.character(two$condidnt)), c("y ctrl", "z ctrl"))
})

test_that("fitErrorModel plots variances and bands per group", {
  skip_if_not_installed("optimx")
  red <- errorModelData()
  pdf(NULL); on.exit(dev.off())
  invisible(fitErrorModel(red, factors = "name", par = c(s0 = -4, srel = -3), plotting = TRUE))
  p <- ggplot2::last_plot()
  expect_identical(unname(vapply(p$layers, function(l) class(l$geom)[1], "")),
                   c("GeomPoint", "GeomLine", "GeomRibbon", "GeomRibbon"))
  expect_identical(p$labels$y, "variance")
  expect_identical(names(p$facet$params$facets), "condidnt")
  b <- ggplot2::ggplot_build(p)
  expect_equal(10^b$data[[1]]$y, red$sigma^2 * red$n)
  expect_identical(as.character(b$layout$layout$condidnt), c("y", "z"))
})


# ---- steadyStates -------------------------------------------------------

test_that("steadyStates refuses named givenCQs and returns lhs = rhs as a named vector", {
  skip_if_not(reticulate::py_module_available("sympy"))
  f <- addReaction(eqnlist(), "A", "B", "k1*A")
  f <- addReaction(f, "B", "A", "k2*B")
  expect_error(steadyStates(f, givenCQs = c(tot = "A + B = tot")), "unname")
  ss <- steadyStates(f)
  expect_type(ss, "character")
  expect_identical(names(ss), trimws(names(ss)))
  expect_snapshot(ss)
})
