# Smoke tests for the cross-platform parallel-apply helper behind mstrust()
# and profile(), covering the PSOCK path a Windows runner would take.


test_that(".parallelLapply with cores=1 falls back to lapply", {
  out <- dMod2:::.parallelLapply(1:5, function(i) i * 2, cores = 1L)
  expect_identical(out, lapply(1:5, function(i) i * 2))
})


test_that(".parallelLapply with cores>1 returns correct results", {
  # Fork and PSOCK backends share one %dopar% loop and match the serial path.
  out_par <- dMod2:::.parallelLapply(1:8, function(i) i^2, cores = 2L)
  out_ser <- lapply(1:8, function(i) i^2)
  expect_identical(out_par, out_ser)
})


test_that(".parallelLapply preserves order across workers", {
  # R CMD check caps cores at 2 via _R_CHECK_LIMIT_CORES_; two workers over 20
  # items still exercise order preservation.
  n_cores <- if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) 2L else 4L
  set.seed(7L)
  X <- as.list(rnorm(20))
  out <- dMod2:::.parallelLapply(X, function(x) x + 1, cores = n_cores)
  expect_identical(out, lapply(X, function(x) x + 1))
})


test_that(".parallelLapply propagates worker errors", {
  # A worker error reaches the caller as a regular R error.
  err <- tryCatch(dMod2:::.parallelLapply(1:3, function(i) {
    if (i == 2L) stop("kaboom") else i
  }, cores = 2L), error = function(e) e)
  expect_s3_class(err, "error")
})


## ---- Two-axis core budget ------------------------------------------------

test_that(".splitCores reads a scalar as the outer axis only", {
  r <- dMod2:::.splitCores(4L, "fits")
  expect_identical(r$outer, 4L)
  expect_null(r$conditions)

  expect_identical(dMod2:::.splitCores(NULL)$outer, 1L)
})


test_that(".splitCores reads both axes, by name or by position", {
  r <- dMod2:::.splitCores(c(fits = 2L, conditions = 1L), "fits")
  expect_identical(r$outer, 2L)
  expect_identical(r$conditions, 1L)

  # positional: first is the outer axis
  r <- dMod2:::.splitCores(c(2L, 1L), "fits")
  expect_identical(r$outer, 2L)
  expect_identical(r$conditions, 1L)

  # the outer axis is named after the caller's axis
  r <- dMod2:::.splitCores(c(pars = 2L, conditions = 1L), "pars")
  expect_identical(r$outer, 2L)
})


test_that(".splitCores rejects nonsense and warns on over-subscription", {
  expect_error(dMod2:::.splitCores(0L), "positive")
  expect_error(dMod2:::.splitCores(c(1L, 2L, 3L)), "length 1 or 2")
  expect_warning(dMod2:::.splitCores(c(fits = 10000L, conditions = 10000L), "fits"),
                 "exceeds")
})


test_that(".sanitizeCores caps each axis, not the vector", {
  avail <- dMod2:::.availableCores()
  expect_identical(dMod2:::.sanitizeCores(c(fits = 1L, conditions = 2L)),
                   c(fits = 1L, conditions = min(2L, avail)))
  expect_identical(dMod2:::.sanitizeCores(c(fits = avail + 5L, conditions = 1L)),
                   c(fits = avail, conditions = 1L))
})


## ---- Layout of mstrust() on a compiled objective --------------------------

fx_register(multicond = TRUE)

# The four-condition fixture under normL2 and a weak prior, which makes the
# fit regular; starts from a fixed parframe. The objective names the joint
# shared object, so loadDLL() finds it in a fresh PSOCK worker.
.layout_fixture <- function() {
  fx <- fx_decay_multicond_compiled()
  data <- fx_decay_data_multi(setNames(rep(list(c(A = 1, k = 0.5)), 4), fx$conditions))
  obj <- normL2(data, fx$prd) + constraintL2(fx$outerpars, sigma = 3)
  modelname(obj) <- "fx_decay_all"
  starts <- msParframe(fx$outerpars, n = 4, sd = 0.3, seed = 3)
  list(obj = obj, starts = starts)
}

# trust() with the worker's PID and its dMod.cores in the fit.
.trustWhere <- function(objfun, parinit, ...) {
  fit <- trust(objfun, parinit, ...)
  fit$pid <- Sys.getpid()
  fit$dModCores <- getOption("dMod.cores", NA_integer_)
  fit
}

.alive <- function(pids) vapply(pids, tools::pskill, TRUE, signal = 0L)

.runLayout <- function(lf, cores) {
  mstrust(lf$obj, center = lf$starts, cores = cores, optmethod = .trustWhere,
          iterlim = 50)
}

.fitsOf <- function(fits)
  lapply(fits, `[`, c("value", "argument", "iterations", "nSwitch", "converged"))


test_that("a sum of objectives keeps the build info loadDLL() searches with", {
  fx <- fx_decay_multicond_compiled()
  data <- fx_decay_data_multi(setNames(rep(list(c(A = 1, k = 0.5)), 4), fx$conditions))
  l2 <- normL2(data, fx$prd)
  expect_identical(attr(l2 + constraintL2(fx$outerpars, sigma = 3), "compileInfo"),
                   attr(l2, "compileInfo"))
})


test_that("mstrust runs the layout asked for and leaves no workers or options behind", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if(dMod2:::.availableCores() < 4L)
  lf <- .layout_fixture()
  old <- options(dMod.cores = NULL); on.exit(options(old), add = TRUE)

  serial <- .runLayout(lf, 1L)
  expect_identical(attr(serial, "run")$cores, c(fits = 1L))
  expect_length(unique(vapply(serial, `[[`, 1L, "pid")), 1L)

  inner <- .runLayout(lf, c(fits = 1L, conditions = 2L))
  expect_identical(attr(inner, "run")$cores, c(fits = 1L, conditions = 2L))
  expect_true(all(vapply(inner, `[[`, 1L, "dModCores") == 2L))
  expect_null(getOption("dMod.cores"))

  fork <- .runLayout(lf, c(fits = 2L, conditions = 1L))
  expect_identical(attr(fork, "run")$cores, c(fits = 2L, conditions = 1L))
  expect_false(any(vapply(fork, `[[`, 1L, "pid") == Sys.getpid()))

  psock <- .runLayout(lf, c(2L, 2L))
  expect_identical(attr(psock, "run")$cores, c(fits = 2L, conditions = 2L))
  expect_length(setdiff(vapply(psock, `[[`, 1L, "pid"), Sys.getpid()), 2L)
  expect_true(all(vapply(psock, `[[`, 1L, "dModCores") == 2L))
  expect_null(getOption("dMod.cores"))

  # more workers than fits
  capped <- .runLayout(lf, c(fits = 10L, conditions = 1L))
  expect_identical(attr(capped, "run")$cores[["fits"]], 4L)

  for (fits in list(inner, fork, psock, capped))
    expect_equal(.fitsOf(fits), .fitsOf(serial))

  workers <- setdiff(vapply(psock, `[[`, 1L, "pid"), Sys.getpid())
  for (i in 1:50) {
    if (!any(.alive(workers))) break
    Sys.sleep(0.1)
  }
  expect_false(any(.alive(workers)))
})


test_that("profile gives the same profiles in every layout and closes its workers", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if(dMod2:::.availableCores() < 4L, "needs 4 cores")
  lf <- .layout_fixture()
  old <- options(dMod.cores = NULL); on.exit(options(old), add = TRUE)
  best <- trust(lf$obj, unlist(lf$starts[1, attr(lf$starts, "parameters")]))$argument
  conns <- nrow(showConnections())

  .prof <- function(cores)
    profile(lf$obj, best, names(best)[1:2], limits = c(-1, 1), cores = cores, verbose = FALSE)
  serial <- .prof(1L)
  expect_setequal(unique(serial$whichPar), names(best)[1:2])
  for (cores in list(c(pars = 1L, conditions = 2L), c(pars = 2L, conditions = 1L),
                     c(pars = 2L, conditions = 2L)))
    expect_equal(.prof(cores), serial)
  expect_null(getOption("dMod.cores"))
  expect_identical(nrow(showConnections()), conns)
})


test_that("profile leaves dMod.cores as it found it", {
  obj <- constraintL2(mu = c(a = 1, b = -0.5), sigma = 1)
  old <- options(dMod.cores = NULL); on.exit(options(old), add = TRUE)
  prof <- profile(obj, c(a = 1, b = -0.5), "a", limits = c(-2, 2),
                  cores = c(pars = 1L, conditions = 2L), verbose = FALSE)
  expect_true(nrow(prof) > 1L)
  expect_null(getOption("dMod.cores"))
})
