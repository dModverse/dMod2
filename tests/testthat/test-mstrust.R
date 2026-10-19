# mstrust(), profile(), vcov() and confint.parframe() on quadratic objectives
# (constraintL2), so the expected behaviour has a closed form.


# ---- mstrust ------------------------------------------------------------

test_that("mstrust returns the global minimum across random starts on a convex quadratic", {
  target <- c(a = 1.0, b = -0.5)
  obj <- constraintL2(mu = target, sigma = 1)

  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd), add = TRUE)
  fits <- mstrust(objfun = obj, center = c(a = 0, b = 0),
                  name = "test_mstrust_convex",
                  rinit = 1, rmax = 10, iterlim = 100,
                  fits = 5, sd = 1, cores = 1, output = FALSE)

  pf <- as.parframe(fits)
  best <- pf[which.min(pf$value), ]
  expect_lt(abs(best$a - target[["a"]]), 1e-6)
  expect_lt(abs(best$b - target[["b"]]), 1e-6)
  expect_equal(min(pf$value), 0, tolerance = 1e-8)
})


test_that("a multi-start reports the Hessian source handovers per fit", {
  # nSwitch is the gate quantity for a switching Hessian source, so it has to
  # survive the trip from trust() through the parlist into the parframe.
  target <- c(a = 1.0, b = -0.5)
  obj <- constraintL2(mu = target, sigma = 1)

  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd), add = TRUE)
  fits <- mstrust(objfun = obj, center = c(a = 0, b = 0),
                  name = "test_mstrust_switch",
                  rinit = 1, rmax = 10, iterlim = 100,
                  phases = list(list(hessianMethod = "gn"), list(hessianMethod = "bfgs")),
                  fits = 3, sd = 1, cores = 1, output = FALSE)

  pf <- as.parframe(fits)
  expect_true("nSwitch" %in% names(pf))
  expect_true(is.integer(pf$nSwitch))
  expect_false(anyNA(pf$nSwitch))
  expect_identical(pf$iterations.gn + pf$iterations.bfgs, pf$iterations)
})


test_that("a multi-start records wall times and run settings, not Hessians", {
  obj <- constraintL2(mu = c(a = 1.0, b = -0.5), sigma = 1)
  fits <- mstrust(objfun = obj, center = c(a = 0, b = 0),
                  rinit = 1, rmax = 10, iterlim = 100,
                  fits = 3, sd = 1, cores = 1)

  expect_false(any(vapply(fits, function(f) "hessian" %in% names(f), TRUE)))
  expect_true(all(vapply(fits, function(f) f$elapsed >= 0, TRUE)))

  run <- attr(fits, "run")
  expect_equal(run$fits, 3)
  expect_equal(run$optimiser, "trust")
  expect_equal(run$optimiserArgs$iterlim, 100)
  expect_equal(run$samplerArgs$sd, 1)
  expect_gte(run$elapsed, 0)
  expect_named(run$versions, c("R", "dMod2", "cppDE"))

  pf <- as.parframe(fits)
  expect_true("elapsed" %in% attr(pf, "metanames"))
  expect_false(anyNA(pf$elapsed))
  expect_identical(attr(pf[1:2, ], "run"), run)
})


test_that("mstrust refuses an argument both the optimiser and the sampler take", {
  # An argument is routed by name, so one that both could take has no single
  # destination. The clash is caught before any fit starts.
  obj <- constraintL2(mu = c(a = 1.0, b = -0.5), sigma = 1)
  seen <- new.env()
  myopt <- function(objfun, parinit, sd = 1, ...) {
    seen$sd <- sd
    list(argument = parinit, value = 0, converged = TRUE, iterations = 0L)
  }

  expect_error(mstrust(objfun = obj, center = c(a = 0, b = 0), fits = 1,
                       optmethod = myopt, sd = 2, cores = 1),
               "Argument names of the optimiser and rnorm\\(\\) clash")
  expect_null(seen$sd)

  # The sample size is always an argument, so an optimiser that takes `n`
  # clashes with every sampler even when the caller passes nothing more.
  nopt <- function(objfun, parinit, n = 1, ...)
    list(argument = parinit, value = 0, converged = TRUE, iterations = 0L)
  expect_error(mstrust(objfun = obj, center = c(a = 0, b = 0), fits = 1,
                       optmethod = nopt, samplefun = "runif", cores = 1),
               "Argument names of the optimiser and runif\\(\\) clash")

  # A sampler without `sd` leaves it to the optimiser alone.
  fits <- mstrust(objfun = obj, center = c(a = 0, b = 0), fits = 1,
                  optmethod = myopt, samplefun = "runif", sd = 2, cores = 1)
  expect_length(fits, 1L)
  expect_equal(seen$sd, 2)
})


# ---- profile ------------------------------------------------------------

test_that("profile on a 1D quadratic increases monotonically on both sides", {
  obj <- constraintL2(mu = c(theta = 0.0, nuisance = 0.0), sigma = 1)

  prof <- profile(obj, pars = c(theta = 0, nuisance = 0),
                  whichPar = "theta",
                  limits = c(lower = -2, upper = 2),
                  method = "integrate",
                  verbose = FALSE, cores = 1)
  pf <- as.data.frame(prof)
  centered <- pf[order(pf$theta), ]
  left  <- centered[centered$theta <= 0, ]
  right <- centered[centered$theta >= 0, ]
  expect_true(all(diff(left$value)  <= 1e-6))
  expect_true(all(diff(right$value) >= -1e-6))
})


test_that("profile on a 1D Gaussian crosses chi^2 = 3.84 at +/- z(0.95)*sigma", {
  sigma <- 0.4
  obj <- constraintL2(mu = c(theta = 0, nuisance = 0), sigma = sigma)
  prof <- profile(obj, pars = c(theta = 0, nuisance = 0),
                  whichPar = "theta",
                  limits = c(lower = -3, upper = 3),
                  method = "integrate", verbose = FALSE, cores = 1)
  ci <- confint(prof, level = 0.95, val.column = "value")
  expected_half_width <- qnorm(0.975) * sigma
  half_width <- (ci$upper - ci$lower) / 2
  expect_lt(abs(half_width - expected_half_width) / expected_half_width,
            0.10)
})


test_that("profile respects limits on a flat direction", {
  # Along an exactly flat direction the stepsize doubles every step; limits
  # clamp the last step rather than being overshot.
  nm <- c("a", "b", "c")
  obj <- function(pars, fixed = NULL, deriv = TRUE, ...) {
    pp <- c(pars, fixed)[nm]
    v <- as.numeric((pp[["a"]] + pp[["b"]])^2 + pp[["c"]]^2)
    g <- c(a = 2 * (pp[["a"]] + pp[["b"]]),
           b = 2 * (pp[["a"]] + pp[["b"]]),
           c = 2 * pp[["c"]])
    H <- matrix(c(2, 2, 0, 2, 2, 0, 0, 0, 2), 3, 3, dimnames = list(nm, nm))
    free <- names(pars)
    structure(list(value = v, gradient = g[free],
                   hessian = H[free, free, drop = FALSE]),
              class = c("objlist", "list"))
  }
  class(obj) <- c("objfn", "fn")
  attr(obj, "modelname") <- character(0)

  prof <- suppressWarnings(profile(
    obj, pars = c(a = 0.3, b = -0.3, c = 0), whichPar = "a",
    method = "integrate", stepControl = list(limit = 20),
    limits = c(lower = -1.5, upper = 1.5), cores = 1))

  expect_equal(range(prof$constraint), c(-1.5, 1.5), tolerance = 1e-6)
})


# ---- vcov ---------------------------------------------------------------

test_that("vcov(fit) equals solve(0.5 * H) for a quadratic with known Hessian", {
  mu <- c(a = 0.5, b = -0.3)
  obj <- constraintL2(mu = mu, sigma = 1)
  fit <- trust(obj, parinit = c(a = 0, b = 0),
               rinit = 1, rmax = 10, iterlim = 50, printIter = FALSE)
  V <- vcov(fit)
  expect_equal(unname(V), diag(2), tolerance = 1e-6)

  mu2 <- c(a = 0.0, b = 0.0)
  obj2 <- constraintL2(mu = mu2, sigma = 0.5)
  fit2 <- trust(obj2, parinit = c(a = 0.2, b = -0.2),
                rinit = 1, rmax = 10, iterlim = 50, printIter = FALSE)
  V2 <- vcov(fit2)
  expect_equal(unname(V2), 0.5^2 * diag(2), tolerance = 1e-6)
})


test_that("vcov of an objective evaluation matches the fit and zeroes fixed parameters", {
  obj <- constraintL2(mu = c(a = 0.5, b = -0.3), sigma = 0.5)
  fit <- trust(obj, parinit = c(a = 0, b = 0),
               rinit = 1, rmax = 10, iterlim = 50, printIter = FALSE)
  expect_equal(vcov(obj(fit$argument)), vcov(fit))

  V <- vcov(obj(fit$argument), fixed = "b")
  expect_equal(V["a", "a"], 0.25, tolerance = 1e-6)
  expect_equal(unname(c(V["b", ], V[, "b"])), rep(0, 4))
})


test_that("vcov of an mstrust fit asks for an objective evaluation", {
  obj <- constraintL2(mu = c(a = 0.5, b = -0.3), sigma = 1)
  fits <- mstrust(obj, center = c(a = 0, b = 0), fits = 1, sd = 1,
                  rinit = 1, rmax = 10, iterlim = 50)
  expect_error(vcov(fits[[1]]), "vcov\\(obj\\(fit\\$argument\\)\\)")
})


# ---- confint ------------------------------------------------------------

# ---- retry / nTries -----------------------------------------------------

make_flaky_obj <- function(fail_first = 3L) {
  counter <- 0L
  fn <- function(pars, fixed = NULL, deriv = TRUE, ...) {
    counter <<- counter + 1L
    if (counter <= fail_first)
      stop("flaky objfn: deliberate failure ", counter, "/", fail_first)
    v <- sum(pars * pars)
    g <- 2 * pars
    H <- 2 * diag(length(pars))
    dimnames(H) <- list(names(pars), names(pars))
    structure(list(value = v, gradient = g, hessian = H),
              class = c("objlist", "list"))
  }
  class(fn) <- c("objfn", "fn")
  attr(fn, "parameters") <- c("a", "b")
  fn
}


test_that("mstrust retries a failing fit up to nTries times", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd), add = TRUE)

  obj <- make_flaky_obj(fail_first = 3L)
  out <- mstrust(objfun = obj, center = c(a = 0, b = 0),
                 name = "test_mstrust_retry_ok",
                 rinit = 1, rmax = 10, iterlim = 50,
                 fits = 1L, sd = 1, cores = 1L,
                 retry = TRUE, nTries = 10L, output = FALSE)

  expect_length(out, 1L)
  expect_true(is.null(out[[1]]$error))
  expect_true(isTRUE(out[[1]]$converged))
})


test_that("mstrust gives up after nTries and records an error", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd), add = TRUE)

  obj <- make_flaky_obj(fail_first = 100L)
  out <- mstrust(objfun = obj, center = c(a = 0, b = 0),
                 name = "test_mstrust_retry_fail",
                 rinit = 1, rmax = 10, iterlim = 50,
                 fits = 1L, sd = 1, cores = 1L,
                 retry = TRUE, nTries = 3L, output = FALSE)

  expect_length(out, 1L)
  expect_true(!is.null(out[[1]]$error))
})


test_that("mstrust with retry = FALSE does not re-attempt", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd), add = TRUE)

  obj <- make_flaky_obj(fail_first = 1L)
  out <- mstrust(objfun = obj, center = c(a = 0, b = 0),
                 name = "test_mstrust_no_retry",
                 rinit = 1, rmax = 10, iterlim = 50,
                 fits = 1L, sd = 1, cores = 1L,
                 retry = FALSE, output = FALSE)

  expect_length(out, 1L)
  expect_true(!is.null(out[[1]]$error))
})


# ---- on-disk results are opt-in -----------------------------------------

test_that("mstrust writes nothing unless asked", {
  target <- c(a = 1.0, b = -0.5)
  obj <- constraintL2(mu = target, sigma = 1)
  d <- file.path(tempdir(), paste0("ms_quiet_", as.integer(runif(1, 1e6, 9e6))))
  dir.create(d)

  invisible(mstrust(objfun = obj, center = c(a = 0, b = 0), resultPath = d,
                    rinit = 1, rmax = 10, iterlim = 100, fits = 3, sd = 1,
                    cores = 1))
  expect_equal(list.files(d, recursive = TRUE, include.dirs = TRUE),
               character(0))

  invisible(mstrust(objfun = obj, center = c(a = 0, b = 0), resultPath = d,
                    name = "asked", rinit = 1, rmax = 10, iterlim = 100,
                    fits = 3, sd = 1, cores = 1, output = TRUE))
  written <- list.files(d, recursive = TRUE)
  expect_true(any(grepl("^asked/", written)))
  expect_true(any(grepl("mstrust\\.log$", written)))
  expect_true(any(grepl("parameterList\\.Rda$", written)))
})


test_that("cautiousMode removes the .Rda fits of its own run only", {
  obj <- constraintL2(mu = c(a = 1.0, b = -0.5), sigma = 1)
  d <- file.path(tempdir(), paste0("ms_cautious_", as.integer(runif(1, 1e6, 9e6))))
  dir.create(d)
  oldwd <- setwd(d); on.exit(setwd(oldwd), add = TRUE)
  saveRDS(1, "fit-1.Rda")

  fits <- mstrust(objfun = obj, center = c(a = 0, b = 0), resultPath = d,
                  name = "cautious", rinit = 1, rmax = 10, iterlim = 100,
                  fits = 2, sd = 1, cores = 1, output = TRUE, cautiousMode = TRUE)
  inter <- list.files(d, recursive = TRUE)
  inter <- inter[grepl("^cautious/.*/interRes/", inter)]
  expect_true(file.exists(file.path(d, "fit-1.Rda")))
  expect_false(any(grepl("\\.Rda$", inter)))
  expect_true(all(c("fit-1.R", "fit-2.R") %in% basename(inter)))

  # The deparsed fit restores the parlist entry bit for bit.
  restored <- dget(file.path(d, inter[basename(inter) == "fit-1.R"]))
  expect_identical(restored$argument, fits[[1]]$argument)
  expect_identical(restored$value, fits[[1]]$value)
})


test_that("subset() on a parframe sees the calling frame", {
  pf <- parframe(data.frame(value = c(1, 5, 9), a = c(0.1, 0.2, 0.3)),
                 parameters = "a", metanames = "value")

  # The condition is written where the caller stands, not where subset() runs.
  pick <- function(cutoff) subset(pf, value < cutoff)

  expect_equal(nrow(pick(6)), 2L)
  expect_equal(nrow(subset(pf, value < 6)), 2L)
  expect_s3_class(pick(6), "parframe")
})


# ---- PSOCK workers --------------------------------------------------------

# Two conditions behind a Pimpl steady state, each source compiled into its
# own shared object, so a worker has to find all of them by itself.
.psock_env <- new.env(parent = emptyenv())

psock_model <- function() {
  if (!is.null(.psock_env$obj)) return(.psock_env)
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)
  mn <- function(x) paste0("test_psock_", x, "_", as.integer(Sys.time()))

  el <- eqnlist() |>
    addReaction("", "A", "k_in") |> addReaction("A", "B", "k1*A") |>
    addReaction("B", "", "k2*B")
  x  <- Xs(odemodel(el, modelname = mn("x"), compile = FALSE), compile = FALSE)
  g  <- Y(c(yB = "s*B"), x, parameters = "s", modelname = mn("g"), compile = FALSE)
  ss <- Pimpl(el, modelname = mn("ss"), compile = FALSE, verbose = FALSE)
  trafo <- define(NULL, "x~exp(x)", x = getParameters(x, g))
  p1 <- P(trafo, condition = "c1", modelname = mn("p1"), compile = FALSE)
  p2 <- P(insert(trafo, "x ~ 2*x", x = "k_in"), condition = "c2",
          modelname = mn("p2"), compile = FALSE)
  compile(g, x, ss, p1, p2, cores = test_cores())

  p <- ss * (p1 + p2)
  pars <- structure(rep(0, length(getParameters(p))), names = getParameters(p))
  pred <- (g * x * p)(seq(0, 10, 1), pars)
  data <- as.datalist(lapply(pred, function(m)
    data.frame(name = "yB", time = m[, "time"], value = m[, "yB"], sigma = 0.1)))

  .psock_env$pars <- pars
  .psock_env$obj  <- normL2(data, g * x * p) + constraintL2(pars, sigma = 10)
  .psock_env
}


test_that("mstrust fits a Pimpl model on PSOCK workers with a condition axis", {
  skip_on_cran()
  skip_if(dMod2:::.availableCores() < 4L, "needs four cores")
  m <- psock_model()

  fits <- mstrust(m$obj, m$pars, fits = 2, sd = 0.5, iterlim = 30,
                  cores = c(fits = 2, conditions = 2))

  errs <- unlist(lapply(fits, `[[`, "error"))
  expect_null(errs)
  expect_true(all(is.finite(as.parframe(fits)$value)))
})


test_that("profile runs on PSOCK workers with a condition axis", {
  skip_on_cran()
  skip_if(dMod2:::.availableCores() < 4L, "needs four cores")
  m <- psock_model()

  prof <- profile(m$obj, m$pars, whichPar = c("k1", "k2"), verbose = FALSE,
                  stepControl = list(limit = 5),
                  cores = c(pars = 2, conditions = 2))

  expect_setequal(unique(prof$whichPar), c("k1", "k2"))
  expect_true(all(is.finite(prof$value)))
})


test_that("trust reports why the objective failed at parinit", {
  obj <- function(pars, ...) stop("no shared object for this model")
  expect_error(trust(obj, c(a = 0), rinit = 1, rmax = 10),
               "objfun failed: no shared object for this model")
})
