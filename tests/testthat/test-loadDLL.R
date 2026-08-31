# loadDLL() must not unload a shared object that is still in use: dyn.unload()
# nulls every native symbol pointer resolved so far, the one baked into a
# prepared batch handle included, and nothing resolves those again.


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
                   obj = normL2(data, prd), so = file.path(dir, "ld_all.so"))
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


test_that("a flushed symbol cache rebuilds the prepared batch handle", {
  fx <- fx_loaddll()
  oldwd <- setwd(fx$dir); on.exit(setwd(oldwd), add = TRUE)

  before <- fx$obj(fx$pouter)$value
  dyn.unload(fx$so)
  dyn.load(fx$so)
  dMod2:::.clearSymbols()
  expect_identical(fx$obj(fx$pouter)$value, before)
})
