# Loading and unloading must not be able to break a model that is in use.
# `dyn.unload()` nulls every native symbol address already handed out, so
# nothing may keep one: entry points are dispatched by name and shared object.


# Two-condition decay model in its own workdir, so the reload tests cannot
# touch the shared fixture libraries. Two conditions are what puts the
# prediction on the batched cppDE entry point.
fx_loaddll <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    dir <- file.path(tempdir(), "dmod_loaddll")
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    oldwd <- setwd(dir); on.exit(setwd(oldwd), add = TRUE)

    el <- addReaction(eqnlist(), "A", "B", "k1*A")
    el <- addReaction(el, "B", "", "k2*B")
    m  <- odemodel(el, modelname = "ld_ode", compile = FALSE)
    x  <- Xs(m)
    g  <- Y(c(yB = "s*B"), f = x, compile = FALSE, modelname = "ld_obs")
    p1 <- P(eqnvec(A = "A0", B = "0", k1 = "k1", k2 = "k2", s = "s"),
            condition = "C1", compile = FALSE, modelname = "ld_p1")
    p2 <- P(eqnvec(A = "2*A0", B = "0", k1 = "k1", k2 = "k2", s = "s"),
            condition = "C2", compile = FALSE, modelname = "ld_p2")
    suppressMessages(compile(x, g, p1, p2, output = "ld_all", cores = 1))

    prd    <- (g * x) * (p1 + p2)
    pouter <- c(A0 = 1, k1 = 0.5, k2 = 0.3, s = 1)
    times  <- seq(0, 10, by = 1)
    pred   <- prd(times, pouter)
    data   <- as.datalist(lapply(pred, function(p)
      data.frame(name = "yB", time = p[, "time"], value = p[, "yB"] + 0.01,
                 sigma = 0.1)))

    cache <<- list(dir = dir, x = x, prd = prd, pouter = pouter,
                   obj = normL2(data, prd),
                   # dyn.unload() matches the path string dyn.load() was given.
                   # compile() loads a normalised path and names the library by
                   # the platform's extension, so both have to be reproduced.
                   so = normalizePath(
                     file.path(dir, paste0("ld_all", .Platform$dynlib.ext)),
                     winslash = "/"))
    cache
  }
})


test_that("loadDLL skips shared objects already loaded in this process", {
  fx <- fx_loaddll()
  oldwd <- setwd(fx$dir); on.exit(setwd(oldwd), add = TRUE)

  expect_length(suppressMessages(loadDLL(fx$x)), 0L)
  expect_silent(loadDLL(fx$obj))
})


test_that("an objective survives repeated loadDLL calls", {
  fx <- fx_loaddll()
  oldwd <- setwd(fx$dir); on.exit(setwd(oldwd), add = TRUE)

  before <- fx$obj(fx$pouter)$value
  suppressMessages(loadDLL(fx$obj))
  suppressMessages(loadDLL(fx$obj))
  expect_identical(fx$obj(fx$pouter)$value, before)
})


test_that("a reloaded shared object needs no cache flush", {
  fx <- fx_loaddll()
  oldwd <- setwd(fx$dir); on.exit(setwd(oldwd), add = TRUE)

  before <- fx$obj(fx$pouter)$value
  dyn.unload(fx$so)
  dyn.load(fx$so)
  # Both the scalar and the batched path, and the batch handle cached by Xs().
  expect_identical(fx$obj(fx$pouter)$value, before)
  expect_identical(fx$obj(fx$pouter)$value, before)
})


test_that("an unloaded model names the shared object it is missing", {
  fx <- fx_loaddll()
  oldwd <- setwd(fx$dir); on.exit(setwd(oldwd), add = TRUE)
  on.exit({ if (!fx$so %in% dMod2:::.loadedDLLPaths()) dyn.load(fx$so) }, add = TRUE)

  # Whichever entry point the chain reaches first, the message has to name it
  # and its library instead of failing on a null address.
  dyn.unload(fx$so)
  expect_error(fx$obj(fx$pouter), "ld_all")
})


test_that("loadDLL finds shared objects outside the working directory", {
  fx <- fx_loaddll()
  other <- file.path(tempdir(), "dmod_loaddll_elsewhere")
  dir.create(other, showWarnings = FALSE, recursive = TRUE)
  oldwd <- setwd(other); on.exit(setwd(oldwd), add = TRUE)

  dyn.unload(fx$so)
  loaded <- suppressMessages(loadDLL(fx$obj))
  expect_identical(loaded, fx$so)
  expect_true(is.finite(fx$obj(fx$pouter)$value))
})


test_that("mstrust converges when it reloads the objective per fit", {
  fx <- fx_loaddll()
  oldwd <- setwd(fx$dir); on.exit(setwd(oldwd), add = TRUE)

  set.seed(1L)
  fits <- mstrust(fx$obj, center = fx$pouter, fits = 3, sd = 0.1, cores = 1,
                  output = FALSE, name = "ld_mstrust", iterlim = 100)
  pf <- as.parframe(fits)
  expect_equal(nrow(pf), 3L)
  expect_true(all(pf$converged))
  # The parent session must still hold a usable objective afterwards.
  expect_true(is.finite(fx$obj(fx$pouter)$value))
})


test_that("recompiling into a loaded shared object serves the new model", {
  dir <- file.path(tempdir(), "dmod_recompile")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  oldwd <- setwd(dir); on.exit(setwd(oldwd), add = TRUE)

  build <- function(rate, nm) {
    m <- odemodel(addReaction(eqnlist(), "A", "", rate), modelname = nm,
                  compile = FALSE)
    x <- Xs(m)
    suppressMessages(compile(x, output = "rc_all", cores = 1))
    x
  }

  x1 <- build("k*A", "rc_one")
  v1 <- unname(x1(0:2, c(A = 1, k = 1))[[1]][3, "A"])
  expect_equal(v1, exp(-2), tolerance = 1e-4)

  # Same output library, different equations: R will not reload a shared object
  # it already holds, so the second model has to displace the first.
  x2 <- build("2*k*A", "rc_two")
  v2 <- unname(x2(0:2, c(A = 1, k = 1))[[1]][3, "A"])
  expect_equal(v2, exp(-4), tolerance = 1e-4)
})
