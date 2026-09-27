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
  # Leaves k and s to be passed through, which only attach.input does.
  pa <- P(list(C1 = c(A = "exp(A_log)"), C2 = c(A = "2*exp(A_log)")),
          derivMode = c("forward", "reverse"), modelname = "ctl_pa", outdir = d)
  pd <- Pequil(c(A = "k_in - k_out * A"), parameters = c("k_in", "k_out"),
               modelname = "ctl_pd", outdir = d, verbose = FALSE)
  el_ERK <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")
  erk <- Pequil(el_ERK, parameters = c("k1", "k2"), expressInTotals = TRUE,
                attach.input = TRUE, controlsMS = list(nStarts = 1L),
                modelname = "ctl_erk", outdir = d, verbose = FALSE)
  pim <- Pimpl(c(x = "x^2 - a"), parameters = "a", controlsMS = list(nStarts = 1L),
               modelname = "ctl_pim", outdir = d, verbose = FALSE)

  compile(m, mf, g, p, pa, pd, erk, pim, output = "ctl_all", cores = 4L)

  .ctl_env$models <- list(m = m, mf = mf, g = g, p = p, pa = pa, pd = pd,
                          erk = erk, pim = pim)
}

# The kernel of an unspecific leaf.
.ctl_kernel <- function(f) attr(f, "mappings")[[1L]]

.ctl_times <- c(0, 2, 5)


## ---- Xs (cppDE): options and forcings -------------------------------------

test_that("Xs solves with the optionsOde set through controls, defaults kept", {
  fx <- ctl_models()
  x <- Xs(fx$m)
  inner <- c(A = 1, k = 0.5)
  exact <- exp(-0.5 * .ctl_times)
  err <- function() max(abs(x(.ctl_times, inner, deriv = FALSE)[[1]][, "A"] - exact))

  controls(x, NULL, "optionsOde") <- list(atol = 1e-12, rtol = 1e-12)
  expect_lt(err(), 1e-9)
  controls(x, NULL, "optionsOde") <- list(atol = 1e-2, rtol = 1e-2)
  expect_gt(err(), 1e-5)
  # The control holds what was set, not the merged list.
  expect_identical(controls(x, NULL, "optionsOde"), list(atol = 1e-2, rtol = 1e-2))

  # A replacement names the entries it changes. The others are the defaults,
  # among them onFailure = "warn", which turns a failed solve into a warning.
  controls(x, NULL, "optionsOde") <- list(maxsteps = 2L, onFailure = "stop")
  expect_error(x(.ctl_times, inner, deriv = FALSE), "Maximum number of steps")
  controls(x, NULL, "optionsOde") <- list(maxsteps = 2L)
  expect_warning(x(.ctl_times, inner, deriv = FALSE), "Maximum number of steps")

  # The sensitivity solve reads its own entry.
  controls(x, NULL, "optionsOde") <- list()
  controls(x, NULL, "optionsSens") <- list(maxsteps = 2L, onFailure = "stop")
  expect_error(x(.ctl_times, inner, deriv = TRUE), "Maximum number of steps")
  expect_silent(x(.ctl_times, inner, deriv = FALSE))
})

test_that("Xs warns about an unknown option set through controls, once", {
  fx <- ctl_models()
  x <- Xs(fx$m)
  inner <- c(A = 1, k = 0.5)
  controls(x, NULL, "optionsOde") <- list(rtol = 1e-8, bogus = 1)
  expect_warning(x(.ctl_times, inner, deriv = FALSE),
                 "optionsOde: Ignoring unknown option\\(s\\): bogus")
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

  # Xf defaults to onFailure = "stop", which a partial replacement keeps.
  controls(xf, NULL, "optionsOde") <- list(maxsteps = 2L)
  expect_error(xf(.ctl_times, inner), "Maximum number of steps")
})

test_that("optionsReverse set later meets the check the constructor makes", {
  fx <- ctl_models()
  expect_error(Xs(fx$mf, optionsReverse = list(gradtol = 1e-3)), "optionsReverse")
  xf <- Xs(fx$mf, forcings = data.frame(name = "u", time = c(0, 100), value = 1))
  controls(xf, NULL, "optionsReverse") <- list(gradtol = 1e-3)
  expect_error(xf(.ctl_times, c(A = 0, k = 1)), "optionsReverse")
  controls(xf, NULL, "optionsReverse") <- NULL
  expect_silent(xf(.ctl_times, c(A = 0, k = 1)))
})

test_that("the reverse weight uses gradtol and floor as they are when it is used", {
  fx <- ctl_models()
  x <- Xs(fx$m, optionsReverse = list(gradtol = 1e-4))
  k <- .ctl_kernel(x)
  w <- matrix(1, length(.ctl_times), 1L, dimnames = list(NULL, "A"))
  attr(k, "vjpfn")(.ctl_times, c(A = 1, k = 0.5), NULL, w)

  wt <- environment(k)$weightGet(NULL, .ctl_times)
  expect_false(is.null(wt$lambda))
  expect_equal(wt$gradtol, 1e-4)
  expect_equal(wt$floor, 0)

  # The weight kept from the last backward pass meets the new settings on the
  # next one, not the ones it was recorded under.
  controls(x, NULL, "optionsReverse") <- list(gradtol = 1e-6, floor = 0.1)
  wt <- environment(k)$weightGet(NULL, .ctl_times)
  expect_equal(wt$gradtol, 1e-6)
  expect_equal(wt$floor, 0.1)
})

test_that("Xs keeps only user options in controls", {
  fx <- ctl_models()
  x <- Xs(fx$m)
  capture.output(nm <- controls(x))
  expect_setequal(nm[[1]], c("forcings", "names", "optionsOde", "optionsSens",
                             "optionsReverse"))
})


## ---- Pexpl --------------------------------------------------------------

test_that("Pexpl reads attach.input from controls in every entry", {
  fx <- ctl_models()
  pa <- fx$pa
  on.exit(controls(pa, NULL, "attach.input") <- FALSE, add = TRUE)
  pars <- c(A_log = 0, k = 0.5, s = 2)

  expect_false(controls(pa, "C1", "attach.input"))
  out <- pa(pars)
  expect_false(any(c("k", "s") %in% names(out$C1)))

  controls(pa, "C1", "attach.input") <- TRUE
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
  controls(pa, "C1", "attach.input") <- FALSE
  ct <- vjp(pars, NULL, c(A = 1, k = 3, s = 0))
  expect_equal(unname(ct["k", 1L]), 0)
})

test_that("attach.input set through controls reaches both sweeps", {
  fx <- ctl_models()
  pa <- fx$pa
  on.exit(controls(pa, NULL, "attach.input") <- FALSE, add = TRUE)
  controls(pa, NULL, "attach.input") <- TRUE

  tight <- list(atol = 1e-11, rtol = 1e-11)
  x <- Xs(fx$m, optionsOde = tight, optionsSens = tight)
  prd <- fx$g * x * pa
  pars <- c(A_log = 0, k = 0.5, s = 2)
  data <- as.datalist(data.frame(
    name = "y", time = c(1, 2, 4, 1, 2, 4),
    value = c(1.1, 0.8, 0.3, 2.3, 1.4, 0.6), sigma = 0.1,
    condition = rep(c("C1", "C2"), each = 3)), split.by = "condition")
  obj <- normL2(data, prd)

  # k and s reach the prediction only through attach.input, in the value pass
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


## ---- Pequil ---------------------------------------------------------------

test_that("Pequil does not answer a changed control from its memo", {
  fx <- ctl_models()
  pd <- fx$pd
  on.exit({
    controls(pd, NULL, "attach.input") <- TRUE
    controls(pd, NULL, "roottol") <- 1e-6
  }, add = TRUE)
  pars <- c(k_in = 1, k_out = 3)

  o1 <- pd(pars)[[1]]
  expect_true(all(c("k_in", "k_out") %in% names(o1)))

  # Same parameters, other controls: the memoised result no longer answers.
  controls(pd, NULL, "attach.input") <- FALSE
  o2 <- pd(pars)[[1]]
  expect_identical(names(o2), "A")
  expect_equal(as.numeric(o2["A"]), as.numeric(o1["A"]), tolerance = 1e-6)

  # roottol also sets the digits the root is rounded to.
  controls(pd, NULL, "roottol") <- 1e-2
  a3 <- as.numeric(pd(pars)[[1]]["A"])
  expect_equal(a3, round(a3, 3))
  expect_false(isTRUE(all.equal(a3, as.numeric(o1["A"]), tolerance = 0)))
})

test_that("Pequil reads its multistart and time window from controls", {
  fx <- ctl_models()
  pd <- fx$pd
  on.exit({
    controls(pd, NULL, "end.time") <- 1e10
    controls(pd, NULL, "controlsMS") <- list()
  }, add = TRUE)
  resetWarmStarts(pd, verbose = FALSE)
  pars <- c(k_in = 1, k_out = 3, A = 100)

  # Too short a window to reach the steady state from any start.
  controls(pd, NULL, "end.time") <- 1e-6
  controls(pd, NULL, "controlsMS") <- list(nStarts = 2L)
  expect_error(pd(pars), "after 2 integration attempt")
  controls(pd, NULL, "controlsMS") <- list()
  expect_error(pd(pars), "after 10 integration attempt")

  controls(pd, NULL, "end.time") <- 1e10
  expect_equal(as.numeric(pd(pars)[[1]]["A"]), 1 / 3, tolerance = 1e-5)
})

test_that("Pequil in totals mode reads keep.root and attach.input from controls", {
  fx <- ctl_models()
  erk <- fx$erk
  on.exit({
    controls(erk, NULL, "attach.input") <- TRUE
    controls(erk, NULL, "keep.root") <- TRUE
  }, add = TRUE)
  pars <- c(k1 = 1, k2 = 3, totalERK = 4, ERK = 1, pERK = 1)

  o1 <- erk(pars)[[1]]
  expect_true(all(c("k1", "k2", "totalERK") %in% names(o1)))
  controls(erk, NULL, "attach.input") <- FALSE
  o2 <- erk(pars)[[1]]
  expect_setequal(names(o2), c("ERK", "pERK"))
  expect_equal(as.numeric(o2["pERK"]), 1, tolerance = 1e-5)

  resetWarmStarts(erk, verbose = FALSE)
  cache <- environment(.ctl_kernel(erk))$reg$get(NULL)
  controls(erk, NULL, "keep.root") <- FALSE
  erk(pars)
  expect_null(cache$yini)
  expect_null(cache$last_result)
  controls(erk, NULL, "keep.root") <- TRUE
  erk(pars)
  expect_false(is.null(cache$yini))
  expect_false(is.null(cache$last_result))
})


## ---- Pimpl --------------------------------------------------------------

test_that("Pimpl reads keep.root and its solver options from controls", {
  fx <- ctl_models()
  pim <- fx$pim
  on.exit({
    controls(pim, NULL, "keep.root") <- TRUE
    controls(pim, NULL, "controlsNleqslv") <- list()
    controls(pim, NULL, "controlsMS") <- list(nStarts = 1L)
  }, add = TRUE)
  expect_identical(controls(pim, NULL, "controlsMS"), list(nStarts = 1L))

  resetWarmStarts(pim, verbose = FALSE)
  reg <- environment(.ctl_kernel(pim))$reg
  controls(pim, NULL, "keep.root") <- FALSE
  expect_equal(as.numeric(pim(c(a = 4, x = 1))[[1]]["x"]), 2, tolerance = 1e-3)
  expect_null(reg$get(NULL)$guess)
  controls(pim, NULL, "keep.root") <- TRUE
  pim(c(a = 4, x = 1))
  expect_false(is.null(reg$get(NULL)$guess))

  # One Newton step from far away does not reach the root, and with a single
  # start there is nothing else to try.
  resetWarmStarts(pim, verbose = FALSE)
  controls(pim, NULL, "controlsNleqslv") <- list(maxit = 1L)
  expect_error(pim(c(a = 4, x = 1e3)), "after 1 attempt")
  set.seed(1)
  controls(pim, NULL, "controlsMS") <- list(nStarts = 3L)
  expect_error(pim(c(a = 4, x = 1e3)), "after 4 attempt")

  # A partial replacement keeps the other defaults.
  controls(pim, NULL, "controlsNleqslv") <- list(ftol = 1e-10)
  expect_equal(as.numeric(pim(c(a = 4, x = 1e3))[[1]]["x"]), 2, tolerance = 1e-8)
})


## ---- The accessor ---------------------------------------------------------

test_that("condition = NULL sets every condition, and gets the first", {
  fx <- ctl_models()
  p <- fx$p
  on.exit(controls(p, NULL, "attach.input") <- FALSE, add = TRUE)

  controls(p, "C2", "attach.input") <- TRUE
  expect_false(controls(p, NULL, "attach.input"))
  expect_false(controls(p, "C1", "attach.input"))
  expect_true(controls(p, "C2", "attach.input"))
  expect_true(controls(p, 2, "attach.input"))

  controls(p, NULL, "attach.input") <- TRUE
  expect_true(controls(p, "C1", "attach.input"))
  expect_true(controls(p, "C2", "attach.input"))
})

test_that("the setter refuses unknown controls and conditions", {
  fx <- ctl_models()
  p <- fx$p
  expect_error(controls(p, NULL, "bogus") <- 1, "no function .* control 'bogus'")
  expect_error(controls(p, "C1", "bogus") <- 1, "Available: attach.input")
  expect_error(controls(p, "C3", "attach.input") <- TRUE, "unknown condition 'C3'")
  expect_error(controls(p, "C3", "attach.input"), "unknown condition 'C3'")
  expect_error(controls(p, 3, "attach.input") <- TRUE, "out of range")
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
    controls(g, NULL, "attach.input") <- FALSE
    controls(p, NULL, "attach.input") <- FALSE
  }, add = TRUE)
  gxp <- g * x * p
  pars <- c(A_log = 0, k_log = log(0.5))

  capture.output(nm <- controls(gxp))
  expect_identical(names(nm), c("obsfn", "prdfn", "parfn C1", "parfn C2"))
  expect_true("optionsOde" %in% nm$prdfn)
  capture.output(nm <- controls(gxp, "C2"))
  expect_identical(names(nm), c("obsfn", "prdfn", "parfn C2"))

  # Set on the composition, seen on the leaf and used by the composition.
  controls(gxp, NULL, "optionsOde") <- list(maxsteps = 2L, onFailure = "stop")
  expect_identical(controls(x, NULL, "optionsOde"),
                   list(maxsteps = 2L, onFailure = "stop"))
  expect_error(gxp(.ctl_times, pars, deriv = FALSE), "did not complete")
  controls(x, NULL, "optionsOde") <- list()
  expect_identical(controls(gxp, "C1", "optionsOde"), list())

  # A name several factors share is set on every one that answers the
  # condition: the observation function, which answers all of them, and the
  # C2 transformation, not the C1 one.
  controls(gxp, "C2", "attach.input") <- TRUE
  expect_true(controls(g, NULL, "attach.input"))
  expect_true(controls(p, "C2", "attach.input"))
  expect_false(controls(p, "C1", "attach.input"))
  out <- gxp(.ctl_times, pars, deriv = FALSE)
  expect_true("A" %in% colnames(out$C1))

  expect_error(controls(gxp, "C3", "attach.input") <- TRUE, "unknown condition")
  expect_error(controls(gxp, NULL, "bogus") <- TRUE, "no function")
})

test_that("controls reach through %.*% and objfn * parfn, which stay one term", {
  fx <- ctl_models()
  vali  <- datapointL2("A", 2, "newpoint", condition = "C1")
  vali2 <- datapointL2("A", 1, "otherpoint", condition = "C2")

  sobj <- 2 %.*% vali
  expect_setequal(controls(sobj), c("mu", "time", "sigma", "attr.name"))
  controls(sobj, "sigma") <- 2
  expect_identical(controls(vali, "sigma"), 2)
  expect_identical(controls(sobj, "sigma"), 2)
  expect_error(controls(sobj, "bogus") <- 1, "no control 'bogus'")

  op <- (vali + vali2) * fx$p
  controls(op, "sigma") <- 3
  expect_identical(controls(vali, "sigma"), 3)
  expect_identical(controls(vali2, "sigma"), 3)

  # A scaled or composed objective is one summand, not the objective inside.
  terms <- dMod2:::.objTerms(sobj)
  expect_length(terms, 1L)
  expect_identical(terms[[1]], sobj)
  expect_identical(dMod2:::.objTerms(op)[[1]], op)
  expect_length(dMod2:::.objTerms(sobj + vali), 2L)
  expect_length(dMod2:::.objTerms(sobj + op + vali), 3L)

  # A sum of wrappers reaches every objective inside.
  both <- sobj + 3 %.*% vali2
  controls(both, "time") <- 1.5
  expect_identical(controls(vali, "time"), 1.5)
  expect_identical(controls(vali2, "time"), 1.5)
})
