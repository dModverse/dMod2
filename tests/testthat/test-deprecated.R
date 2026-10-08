# Renamed arguments and functions: the old name works with one warning, both
# names together are an error.

test_that("renamed arguments of the objectives warn and keep working", {
  mu <- c(a = 0, b = 1)
  expect_warning(o <- constraintL2(mu, attr.name = "p"), "'attr.name' is deprecated")
  expect_true(!is.null(attr(o(c(a = 1, b = 1)), "p")))
  expect_error(constraintL2(mu, attrName = "p", attr.name = "q"),
               "give 'attrName' only")
  expect_warning(constraintCauchy(mu, attr.name = "p"), "'attr.name' is deprecated")
  expect_error(constraintGamma(c(a = 1), bogus = 1), "unused argument")

  expect_warning(v <- datapointL2("A", 1, value = "d", attrName = "val",
                                  condition = "C1"),
                 "'value' is deprecated")
  expect_identical(attr(v, "parameters"), "d")
  expect_identical(controls(v, "mu"), c(d = "A"))
  expect_identical(suppressWarnings(controls(v, "attr.name")), "val")
  expect_error(datapointL2("A", 1, parameter = "d", value = "e", condition = "C1"),
               "give 'parameter' only")
  expect_error(normL2(NULL, NULL, attrName = "x", attr.name = "y"),
               "give 'attrName' only")
  expect_error(normL2(NULL, NULL, optBLOQ = "M1", opt.BLOQ = "M3"),
               "give 'optBLOQ' only")
})

test_that("renamed arguments of data and parameter frames warn", {
  df <- data.frame(name = "A", time = 0:1, value = 1, sigma = 1,
                   condition = c("x", "y"), dose = 1:2)
  expect_warning(d <- as.datalist(df, split.by = "condition",
                                  keepCovariates = "dose"),
                 "'split.by' is deprecated")
  expect_warning(as.datalist(df, splitBy = "condition", keep.covariates = "dose"),
                 "'keep.covariates' is deprecated")
  expect_identical(names(d), c("x", "y"))
  expect_true("dose" %in% names(covariates(d)))

  expect_warning(pf <- parframe(data.frame(value = 1, a = 2), parameters = "a",
                                metanames = "value", obj.attributes = "value"),
                 "'obj.attributes' is deprecated")
  expect_identical(attr(pf, "obj.attributes"), "value")
  fits <- list(list(value = 2, argument = c(a = 1), converged = TRUE, iterations = 1L),
               list(value = 1, argument = c(a = 2), converged = TRUE, iterations = 1L))
  class(fits) <- c("parlist", "list")
  expect_warning(pf <- as.parframe(fits, sort.by = "value"), "'sort.by' is deprecated")
  expect_identical(pf$a, c(2, 1))
})

test_that("renamed arguments of utilities warn", {
  x <- structure(1:3, names = c("a", "b", "c"), note = "n")
  expect_warning(y <- attrs(x, atr = "note", keep = FALSE), "'atr' is deprecated")
  expect_null(attr(y, "note"))
  expect_warning(r <- profileParsPerNode(1:10, fits_per_node = 4),
                 "'fits_per_node' is deprecated")
  expect_identical(r, profileParsPerNode(1:10, 4))
  expect_warning(ctl <- reconstControl(perprimeCap = 50L), "'perprimeCap' is deprecated")
  expect_identical(ctl$perPrimeCap, 50L)
  expect_identical(dMod2:::.symControlSnake(c("perPrimeCap", "minSupportCandCap",
                                              "degreeCap")),
                   c("perprime_cap", "minsupport_cand_cap", "degree_cap"))
})

test_that("rref names its result and ignores fractions with a warning", {
  A <- rbind(c(1, 2, 3), c(2, 4, 7))
  r <- rref(A)
  expect_named(r, c("rref", "pivots"))
  expect_identical(r[[1]], r$rref)
  expect_warning(rref(A, fractions = TRUE), "'fractions' is deprecated")
})

test_that("repar takes the transformation first and recognises the old order", {
  tr <- c(a = "a", b = "b")
  expect_identical(repar(tr, "x ~ exp(x)", x = "a"), c(a = "exp(a)", b = "b"))
  expect_warning(old <- repar("x ~ exp(x)", tr, x = "a"), "repar\\(expr, trafo\\)")
  expect_identical(old, repar(tr, "x ~ exp(x)", x = "a"))
  expect_warning(id <- repar("x ~ x", x = c("a", "b")), "deprecated")
  expect_equal(unclass(id), c(a = "a", b = "b"))
  expect_equal(unclass(repar(NULL, "x ~ x", x = "a")), c(a = "a"))
})

test_that("trust returns a trustfit; vcov and profile are stats methods", {
  obj <- constraintL2(mu = c(a = 0.5, b = -0.3), sigma = 1)
  fit <- trust(obj, c(a = 0, b = 0), rinit = 1, rmax = 10)
  expect_s3_class(fit, "trustfit")
  expect_true(is.list(fit))
  expect_identical(vcov(fit), stats::vcov(fit))
  expect_warning(V <- vcov(fit = fit), "'fit' is deprecated")
  expect_equal(V, vcov(fit))
  expect_false(exists("profile", envir = asNamespace("dMod2"), inherits = FALSE))

  expect_warning(
    p1 <- profile(objfun = obj, pars = c(a = 0.5, b = -0.3), whichPar = "a",
                  limits = c(-1, 1), cores = 1),
    "'objfun' is deprecated")
  p2 <- profile(obj, c(a = 0.5, b = -0.3), "a", limits = c(-1, 1), cores = 1)
  expect_equal(p1, p2)
  plain <- function(p, ...) obj(p, ...)
  expect_s3_class(profile(plain, c(a = 0.5, b = -0.3), "a", limits = c(-1, 1),
                          cores = 1), "parframe")

  expect_warning(trust(obj, c(a = 0, b = 0), parscale = c(1, 1)),
                 "'parscale' is deprecated")
  expect_warning(f2 <- trust(obj, c(a = 0, b = 0), stepControl = list(theta.max = 0.9)),
                 "theta.max is deprecated")
  expect_equal(f2$argument, fit$argument, tolerance = 1e-6)
  expect_error(trust(obj, c(a = 0, b = 0),
                     stepControl = list(theta.max = 0.9, thetaMax = 0.9)),
               "thetaMax only")

  fits <- mstrust(obj, c(a = 0, b = 0), fits = 2, cores = 1, sd = 1,
                  samplefun = stats::rnorm)
  expect_s3_class(fits[[1]], "trustfit")
  expect_warning(mstrust(obj, c(a = 0, b = 0), fits = 1, cores = 1,
                         start1stfromCenter = TRUE),
                 "'start1stfromCenter' is deprecated")
})

test_that("msParframe leaves the global RNG alone", {
  set.seed(7)
  before <- .Random.seed
  msParframe(c(a = 0, b = 1), n = 3, seed = 1)
  expect_identical(.Random.seed, before)
})

test_that("readPetab* are deprecated aliases", {
  skip_if_not_installed("yaml")
  yaml <- system.file("extdata/petab_boehm/Boehm.yaml", package = "dMod2")
  skip_if(!nzchar(yaml))
  expect_warning(t1 <- readPetabTables(yaml), "'readPetabTables' is deprecated")
  expect_identical(t1, readPEtabTables(yaml))
  expect_warning(readPetabYaml(yaml), "'readPetabYaml' is deprecated")
})

test_that("plot functions and remote helpers take the renamed arguments", {
  expect_error(plotPathsMulti(NULL, whichPar = "a", whichpars = "a"),
               "give 'whichPar' only")
  expect_error(plotProfilesAndPaths(NULL, whichPar = "a", ncol = 1, ncols = 2),
               "give 'ncol' only")
  times <- 0:5
  grid <- data.frame(name = "A", time = times, row.names = paste0("A", times))
  x <- Xd(grid)
  pars <- structure(exp(-times / 2), names = getParameters(x))
  expect_warning(p <- plotFluxes(pars, x, seq(0, 5, 0.5), c(prod = "0.2"),
                                 nameFlux = "F"),
                 "'nameFlux' is deprecated")
  expect_s3_class(p, "ggplot")
  expect_error(distributedComputing(1, jobname = "x", nRep = 1, no_rep = 2),
               "give 'nRep' only")
  expect_error(symmetryDetection(eqnvec(x = "-k*x")), "give the observables")
  expect_false("rates" %in% names(formals(steadyStates)))
  expect_warning(dMod2:::.droppedArgs(list(rates = 1), "rates", "steadyStates"),
                 "'rates' is deprecated and ignored")
  expect_error(dMod2:::.droppedArgs(list(rates = 1, bogus = 2), "rates", "f") |>
                 suppressWarnings(), "unused argument")
})
