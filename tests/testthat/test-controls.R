# Behavioral tests for controls() / controls<-.
#
# Every construction option a function reads when it runs lives in the
# `controls` list of its kernel, is read from there at every call, and is
# reachable through controls() on the function itself and on anything summed
# or composed from it. The tests change controls, so every object they touch
# is built here rather than taken from the shared fixtures.

skip_on_cran()
testthat::skip_if_not_installed("cppDE")


## ---- Shared models -----------------------------------------------------

.ctl_env <- new.env(parent = emptyenv())

# Every model the file evaluates, generated first and linked into one shared
# object. The prediction functions are built from the compiled odemodels by
# each test that needs one, so a test starts from the constructor's controls.
ctl_models <- function() {
  if (!is.null(.ctl_env$models)) return(.ctl_env$models)
  d <- file.path(tempdir(), "dmod_ctl")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  owd <- setwd(d); on.exit(setwd(owd), add = TRUE)

  re <- addReaction(eqnlist(), "A", "", "k*A", "decay")
  m  <- odemodel(re, modelname = "ctl_ode", derivMode = c("forward", "reverse"),
                 outdir = d, compile = FALSE)
  mf <- odemodel(eqnvec(A = "u - k*A"), forcings = "u", modelname = "ctl_forc",
                 outdir = d, compile = FALSE)
  g  <- Y(c(y = "s*A"), f = re, derivMode = c("forward", "reverse"),
          modelname = "ctl_obs", outdir = d)
  p  <- P(list(C1 = c(A = "exp(A_log)", k = "exp(k_log)", s = "1"),
               C2 = c(A = "exp(A_log)", k = "exp(k_log)", s = "2")),
          modelname = "ctl_p", outdir = d)
  # Leaves k and s to be passed through, which only attachInput does.
  pa <- P(list(C1 = c(A = "exp(A_log)"), C2 = c(A = "2*exp(A_log)")),
          derivMode = c("forward", "reverse"), modelname = "ctl_pa", outdir = d)
  pim <- Pimpl(c(x = "x^2 - a"), parameters = "a",
               modelname = "ctl_pim", outdir = d, verbose = FALSE)

  compile(m, mf, g, p, pa, pim, output = "ctl_all", cores = test_cores())

  .ctl_env$models <- list(m = m, mf = mf, g = g, p = p, pa = pa, pim = pim)
}

# The kernel of an unspecific leaf.
.ctl_kernel <- function(f) attr(f, "mappings")[[1L]]

.ctl_times <- c(0, 2, 5)


## ---- Xs (cppDE): options and forcings -------------------------------------

test_that("Xs solves with the options set through controls, defaults kept", {
  fx <- ctl_models()
  x <- Xs(fx$m)
  inner <- c(A = 1, k = 0.5)
  exact <- exp(-0.5 * .ctl_times)
  err <- function() max(abs(x(.ctl_times, inner, deriv = FALSE)[[1]][, "A"] - exact))

  controls(x, NULL, "options") <- list(atol = 1e-12, rtol = 1e-12)
  expect_lt(err(), 1e-9)
  controls(x, NULL, "options") <- list(atol = 1e-2, rtol = 1e-2)
  expect_gt(err(), 1e-5)
  # The control holds what was set, not the merged list.
  expect_identical(controls(x, NULL, "options"), list(atol = 1e-2, rtol = 1e-2))

  # A replacement names the entries it changes. The others are the defaults,
  # among them onFailure = "stop".
  controls(x, NULL, "options") <- list(maxsteps = 2L)
  expect_error(x(.ctl_times, inner, deriv = FALSE), "Maximum number of steps")
  controls(x, NULL, "options") <- list(maxsteps = 2L, onFailure = "warn")
  expect_warning(x(.ctl_times, inner, deriv = FALSE), "Maximum number of steps")

  # options reaches the sensitivity solve, optionsSens overrides it there.
  controls(x, NULL, "options") <- list(maxsteps = 2L)
  expect_error(x(.ctl_times, inner, deriv = TRUE), "Maximum number of steps")
  controls(x, NULL, "optionsSens") <- list(maxsteps = 1e6L)
  expect_silent(x(.ctl_times, inner, deriv = TRUE))
  controls(x, NULL, "options") <- list()
  controls(x, NULL, "optionsSens") <- list(maxsteps = 2L)
  expect_error(x(.ctl_times, inner, deriv = TRUE), "Maximum number of steps")
  expect_silent(x(.ctl_times, inner, deriv = FALSE))
})

test_that("optionsOde is a deprecated alias of options", {
  fx <- ctl_models()
  inner <- c(A = 1, k = 0.5)
  expect_warning(x <- Xs(fx$m, optionsOde = list(maxsteps = 2L)),
                 "'optionsOde' is deprecated")
  expect_identical(controls(x, NULL, "options"), list(maxsteps = 2L))
  expect_error(x(.ctl_times, inner, deriv = FALSE), "Maximum number of steps")
  expect_warning(Xs(fx$m, options = list(), optionsOde = list(rtol = 1e-8)),
                 "'optionsOde' is deprecated")
  expect_error(Xs(fx$m, options = list(rtol = 1e-8), optionsOde = list()),
               "give 'options' only")

  expect_warning(controls(x, NULL, "optionsOde") <- list(),
                 "'optionsOde' is deprecated")
  expect_identical(controls(x, NULL, "options"), list())
  expect_warning(expect_identical(controls(x, NULL, "optionsOde"), list()),
                 "'optionsOde' is deprecated")

  expect_warning(xf <- Xf(fx$m, optionsOde = list(maxsteps = 2L)),
                 "'optionsOde' is deprecated")
  expect_error(xf(.ctl_times, inner), "Maximum number of steps")
})

test_that("Xs warns about an unknown option set through controls, once", {
  fx <- ctl_models()
  x <- Xs(fx$m)
  inner <- c(A = 1, k = 0.5)
  controls(x, NULL, "options") <- list(rtol = 1e-8, bogus = 1)
  expect_warning(x(.ctl_times, inner, deriv = FALSE),
                 "options: Ignoring unknown option\\(s\\): bogus")
  expect_no_warning(x(.ctl_times, inner, deriv = FALSE))
})

test_that("Xs keeps the forcings as given and integrates the ones set later", {
  fx <- ctl_models()
  u1 <- data.frame(name = "u", time = c(0, 100), value = c(1, 1))
  xf <- Xs(fx$mf, forcings = u1)
  inner <- c(A = 0, k = 1)
  expect_identical(controls(xf, NULL, "forcings"), u1)

  a1 <- xf(.ctl_times, inner, deriv = FALSE)[[1]][, "A"]
  controls(xf, NULL, "forcings") <- data.frame(name = "u", time = c(0, 100),
                                               value = c(2, 2))
  a2 <- xf(.ctl_times, inner, deriv = FALSE)[[1]][, "A"]
  expect_equal(a2, 2 * a1, tolerance = 1e-5)

  # The batched path keys its prepared handle on the forcings, so a change
  # has to reach it as well.
  many <- function() dMod2:::.predictMany(xf, .ctl_times,
                                          list(inner, c(A = 0, k = 2)),
                                          c("a", "b"), deriv = FALSE)
  b2 <- many()$a[, "A"]
  controls(xf, NULL, "forcings") <- data.frame(name = "u", time = c(0, 100),
                                               value = c(3, 3))
  b3 <- many()$a[, "A"]
  expect_equal(b2, a2, tolerance = 1e-8)
  expect_equal(b3, 3 * a1, tolerance = 1e-5)

  # Readers of the forcings see the ones the leaf runs with.
  s <- NULL
  capture.output(s <- summary(xf))
  expect_equal(s[[1]]$forcings$value, c(3, 3))

  # Checked like the argument, when the leaf runs.
  controls(xf, NULL, "forcings") <- data.frame(name = "u", time = c(0, 1),
                                               value = c(NA, 1))
  expect_error(xf(.ctl_times, inner, deriv = FALSE), "NA values")
})

test_that("Xf keeps its forcings and options as given, too", {
  fx <- ctl_models()
  u1 <- data.frame(name = "u", time = c(0, 100), value = c(1, 1))
  xf <- Xf(fx$mf, forcings = u1)
  inner <- c(A = 0, k = 1)
  a1 <- xf(.ctl_times, inner)[[1]][, "A"]
  controls(xf, NULL, "forcings") <- transform(u1, value = 2)
  expect_equal(xf(.ctl_times, inner)[[1]][, "A"], 2 * a1, tolerance = 1e-5)

  # A partial replacement keeps the default onFailure = "stop".
  controls(xf, NULL, "options") <- list(maxsteps = 2L)
  expect_error(xf(.ctl_times, inner), "Maximum number of steps")
})

test_that("optionsReverse set later meets the check the constructor makes", {
  fx <- ctl_models()
  expect_error(Xs(fx$mf, optionsReverse = list(refine = TRUE)), "optionsReverse")
  xf <- Xs(fx$mf, forcings = data.frame(name = "u", time = c(0, 100), value = 1))
  controls(xf, NULL, "optionsReverse") <- list(refine = TRUE)
  expect_error(xf(.ctl_times, c(A = 0, k = 1)), "optionsReverse")
  controls(xf, NULL, "optionsReverse") <- NULL
  expect_silent(xf(.ctl_times, c(A = 0, k = 1)))
})

test_that("optionsReverse warns about an unknown entry, once", {
  fx <- ctl_models()
  expect_warning(x <- Xs(fx$m, optionsReverse = list(refine = TRUE, gradtoll = 1)),
                 "optionsReverse: Ignoring unknown option\\(s\\): gradtoll")
  k <- .ctl_kernel(x)
  expect_true(environment(k)$sweepCtl()$refine)

  controls(x, NULL, "optionsReverse") <- list(refin = TRUE)
  expect_warning(ctl <- environment(k)$sweepCtl(),
                 "optionsReverse: Ignoring unknown option\\(s\\): refin")
  expect_null(ctl)
  expect_no_warning(environment(k)$sweepCtl())
})

test_that("sensErrCon reaches the forward solver from optionsSens", {
  fx <- ctl_models()
  x <- Xs(fx$m)
  k <- .ctl_kernel(x)
  p <- c(A = 1, k = 0.5)
  seen <- list()
  real <- cppDE::solveODE
  local_mocked_bindings(solveODE = function(...) {
    seen <<- c(seen, list(list(...)$sensErrCon))
    real(...)
  }, .package = "cppDE")

  k(.ctl_times, p, deriv = FALSE)
  k(.ctl_times, p, deriv = TRUE)
  k(.ctl_times, p, deriv = FALSE, keepStore = TRUE)
  w <- matrix(1, length(.ctl_times), 1L, dimnames = list(NULL, "A"))
  attr(k, "vjpfn")(.ctl_times, p, NULL, w)
  expect_identical(seen, list(NULL, TRUE, NULL, NULL))

  seen <- list()
  controls(x, NULL, "optionsSens") <- list(sensErrCon = FALSE)
  k(.ctl_times, p, deriv = FALSE)
  k(.ctl_times, p, deriv = TRUE)
  k(.ctl_times, p, deriv = FALSE, keepStore = TRUE)
  attr(k, "vjpfn")(.ctl_times, p, NULL, w)
  expect_identical(seen, list(NULL, FALSE, NULL, NULL))

  # sensErrCon belongs to optionsSens alone, so options does not know it.
  expect_warning(controls(x, NULL, "options") <- list(sensErrCon = FALSE),
                 NA)
  expect_warning(k(.ctl_times, p, deriv = FALSE), "sensErrCon")
})

test_that("the reverse control uses refine and gradtol as they are when it is used", {
  fx <- ctl_models()
  x <- Xs(fx$m, optionsReverse = list(refine = TRUE, gradtol = 1e-4))
  k <- .ctl_kernel(x)
  ctl <- environment(k)$sweepCtl()
  expect_true(ctl$refine)
  expect_equal(ctl$gradtol, 1e-4)

  controls(x, NULL, "optionsReverse") <- list(refine = TRUE, gradtol = 1e-6)
  expect_equal(environment(k)$sweepCtl()$gradtol, 1e-6)
  controls(x, NULL, "optionsReverse") <- NULL
  expect_null(environment(k)$sweepCtl())
})

test_that("Xs keeps only user options in controls", {
  fx <- ctl_models()
  x <- Xs(fx$m)
  capture.output(nm <- controls(x))
  expect_setequal(nm[[1]], c("forcings", "names", "options", "optionsSens",
                             "optionsReverse"))
})


## ---- Pexpl --------------------------------------------------------------

test_that("Pexpl reads attachInput from controls in every entry", {
  fx <- ctl_models()
  pa <- fx$pa
  on.exit(controls(pa, NULL, "attachInput") <- FALSE, add = TRUE)
  pars <- c(A_log = 0, k = 0.5, s = 2)

  expect_false(controls(pa, "C1", "attachInput"))
  out <- pa(pars)
  expect_false(any(c("k", "s") %in% names(out$C1)))

  controls(pa, "C1", "attachInput") <- TRUE
  out <- pa(pars)
  expect_true(all(c("k", "s") %in% names(out$C1)))
  expect_false(any(c("k", "s") %in% names(out$C2)))

  # The batched entry of the same kernel.
  bf <- attr(attr(pa, "mappings")$C1, "batchfn")
  res <- bf(list(pars, pars * 2), list(NULL, NULL), deriv = TRUE, deriv2 = FALSE,
            conditions = list("C1", "C1"), cores = 1L)
  expect_true(all(vapply(res, function(r) all(c("k", "s") %in% names(r)), TRUE)))
  expect_true(all(vapply(res, function(r)
    all(c("k", "s") %in% rownames(attr(r, "deriv"))), TRUE)))

  # The reverse entry: a pass-through parameter keeps the cotangent put on it.
  vjp <- attr(attr(pa, "mappings")$C1, "vjpfn")
  ct <- vjp(pars, NULL, c(A = 1, k = 3, s = 0))
  expect_equal(unname(ct["k", 1L]), 3)
  controls(pa, "C1", "attachInput") <- FALSE
  ct <- vjp(pars, NULL, c(A = 1, k = 3, s = 0))
  expect_equal(unname(ct["k", 1L]), 0)
})

test_that("attachInput set through controls reaches both sweeps", {
  fx <- ctl_models()
  pa <- fx$pa
  on.exit(controls(pa, NULL, "attachInput") <- FALSE, add = TRUE)
  controls(pa, NULL, "attachInput") <- TRUE

  tight <- list(atol = 1e-11, rtol = 1e-11)
  x <- Xs(fx$m, options = tight)
  prd <- fx$g * x * pa
  pars <- c(A_log = 0, k = 0.5, s = 2)
  data <- as.datalist(data.frame(
    name = "y", time = c(1, 2, 4, 1, 2, 4),
    value = c(1.1, 0.8, 0.3, 2.3, 1.4, 0.6), sigma = 0.1,
    condition = rep(c("C1", "C2"), each = 3)), splitBy = "condition")
  obj <- normL2(data, prd)

  # k and s reach the prediction only through attachInput, in the value pass
  # as in either derivative pass. The oracle is the value itself.
  f <- obj(pars, deriv = TRUE)
  r <- obj(pars, deriv = TRUE, sweep = "reverse")
  fd <- vapply(names(pars), function(nm) {
    h <- 1e-5
    up <- dn <- pars
    up[nm] <- up[nm] + h
    dn[nm] <- dn[nm] - h
    (obj(up, deriv = FALSE)$value - obj(dn, deriv = FALSE)$value) / (2 * h)
  }, 0)
  expect_true(all(abs(fd[c("k", "s")]) > 1e-3))
  expect_equal(unname(f$gradient[names(pars)]), unname(fd), tolerance = 1e-5)
  expect_equal(unname(r$gradient[names(pars)]), unname(fd), tolerance = 1e-5)
})


## ---- Pimpl --------------------------------------------------------------

test_that("Pimpl reads keepRoot and its solver options from controls", {
  fx <- ctl_models()
  pim <- fx$pim
  on.exit({
    controls(pim, NULL, "keepRoot") <- TRUE
    controls(pim, NULL, "controlsPTC") <- list()
  }, add = TRUE)
  expect_identical(controls(pim, NULL, "controlsPTC"), list())

  resetWarmStarts(pim, verbose = FALSE)
  reg <- environment(.ctl_kernel(pim))$reg
  controls(pim, NULL, "keepRoot") <- FALSE
  expect_equal(as.numeric(pim(c(a = 4, x = 1))[[1]]["x"]), 2, tolerance = 1e-8)
  expect_null(reg$get(NULL)$arch)
  controls(pim, NULL, "keepRoot") <- TRUE
  pim(c(a = 4, x = 1))
  expect_false(is.null(reg$get(NULL)$arch))

  # One damped Newton step from far away does not reach the root, and every
  # attempt says so.
  resetWarmStarts(pim, verbose = FALSE)
  controls(pim, NULL, "controlsPTC") <- list(maxit = 1L)
  expect_error(pim(c(a = 4, x = 1e3)), "no convergence in 1 iteration")

  # A partial replacement keeps the other defaults.
  controls(pim, NULL, "controlsPTC") <- list(reltol = 1e-13)
  expect_equal(as.numeric(pim(c(a = 4, x = 1e3))[[1]]["x"]), 2, tolerance = 1e-12)
})


## ---- The accessor ---------------------------------------------------------

test_that("condition = NULL sets every condition, and gets the first", {
  fx <- ctl_models()
  p <- fx$p
  on.exit(controls(p, NULL, "attachInput") <- FALSE, add = TRUE)

  controls(p, "C2", "attachInput") <- TRUE
  expect_false(controls(p, NULL, "attachInput"))
  expect_false(controls(p, "C1", "attachInput"))
  expect_true(controls(p, "C2", "attachInput"))
  expect_true(controls(p, 2, "attachInput"))

  controls(p, NULL, "attachInput") <- TRUE
  expect_true(controls(p, "C1", "attachInput"))
  expect_true(controls(p, "C2", "attachInput"))
})

test_that("the setter refuses unknown controls and conditions", {
  fx <- ctl_models()
  p <- fx$p
  expect_error(controls(p, NULL, "bogus") <- 1, "no function .* control 'bogus'")
  expect_error(controls(p, "C1", "bogus") <- 1, "Available: attachInput")
  expect_error(controls(p, "C3", "attachInput") <- TRUE, "unknown condition 'C3'")
  expect_error(controls(p, "C3", "attachInput"), "unknown condition 'C3'")
  expect_error(controls(p, 3, "attachInput") <- TRUE, "out of range")
  expect_null(controls(p, "C1", "bogus"))
  capture.output(nm <- controls(p))
  expect_identical(names(nm), c("C1", "C2"))

  # A kernel without controls does not get any.
  xt <- Xt()
  expect_error(controls(xt, NULL, "foo") <- 1, "no function")
  expect_false(exists("controls", envir = environment(.ctl_kernel(xt)),
                      inherits = FALSE))
})

test_that("controls on g * x * p reach the leaves of the composition", {
  fx <- ctl_models()
  x <- Xs(fx$m)
  g <- fx$g; p <- fx$p
  on.exit({
    controls(g, NULL, "attachInput") <- FALSE
    controls(p, NULL, "attachInput") <- FALSE
  }, add = TRUE)
  gxp <- g * x * p
  pars <- c(A_log = 0, k_log = log(0.5))

  capture.output(nm <- controls(gxp))
  expect_identical(names(nm), c("obsfn", "prdfn", "parfn C1", "parfn C2"))
  expect_true("options" %in% nm$prdfn)
  capture.output(nm <- controls(gxp, "C2"))
  expect_identical(names(nm), c("obsfn", "prdfn", "parfn C2"))

  # Set on the composition, seen on the leaf and used by the composition.
  controls(gxp, NULL, "options") <- list(maxsteps = 2L, onFailure = "stop")
  expect_identical(controls(x, NULL, "options"),
                   list(maxsteps = 2L, onFailure = "stop"))
  expect_error(gxp(.ctl_times, pars, deriv = FALSE), "did not complete")
  controls(x, NULL, "options") <- list()
  expect_identical(controls(gxp, "C1", "options"), list())

  # A name several factors share is set on every one that answers the
  # condition: the observation function, which answers all of them, and the
  # C2 transformation, not the C1 one.
  controls(gxp, "C2", "attachInput") <- TRUE
  expect_true(controls(g, NULL, "attachInput"))
  expect_true(controls(p, "C2", "attachInput"))
  expect_false(controls(p, "C1", "attachInput"))
  out <- gxp(.ctl_times, pars, deriv = FALSE)
  expect_true("A" %in% colnames(out$C1))

  expect_error(controls(gxp, "C3", "attachInput") <- TRUE, "unknown condition")
  expect_error(controls(gxp, NULL, "bogus") <- TRUE, "no function")
})

test_that("controls reach through %.*% and objfn * parfn, which stay one term", {
  fx <- ctl_models()
  x <- Xs(fx$m)
  data <- as.datalist(data.frame(
    name = "y", time = c(1, 2, 1, 2), value = c(0.6, 0.4, 1.2, 0.7),
    sigma = 0.1, condition = c("C1", "C1", "C2", "C2")), splitBy = "condition")
  obj <- normL2(data, fx$g * x * fx$p)

  sobj <- 2 %.*% obj
  expect_true("multipleShootingControl" %in% controls(sobj))
  controls(sobj, "multipleShootingControl") <- list(nodes = c(1, 3))
  expect_identical(controls(obj, "multipleShootingControl"), list(nodes = c(1, 3)))
  expect_identical(controls(sobj, "multipleShootingControl"), list(nodes = c(1, 3)))
  expect_error(controls(sobj, "bogus") <- 1, "no control 'bogus'")

  op <- obj * fx$p
  controls(op, "multipleShootingControl") <- list(nodes = 2)
  expect_identical(controls(obj, "multipleShootingControl"), list(nodes = 2))

  # Multiple shooting reads the terms of an objective as its summands; a
  # scaled or composed objective is one of them, not the objective inside.
  terms <- dMod2:::.objTerms(sobj)
  expect_length(terms, 1L)
  expect_identical(terms[[1]], sobj)
  expect_identical(dMod2:::.objTerms(op)[[1]], op)
  expect_length(dMod2:::.objTerms(sobj + obj), 2L)
  expect_length(dMod2:::.objTerms(sobj + op + obj), 3L)

  # A sum of wrappers reaches every objective inside.
  vali <- datapointL2("A", 2, "newpoint", condition = "C1")
  both <- sobj + 3 %.*% vali
  expect_setequal(controls(both), c("multipleShootingControl", "mu", "time",
                                    "sigma", "attrName"))
  controls(both, "sigma") <- 2
  expect_identical(controls(vali, "sigma"), 2)
})
