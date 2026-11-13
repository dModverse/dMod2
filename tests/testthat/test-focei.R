# FOCEI: the C++ kernel of the inner problem against an R replica, the outer fit against the
# Theoph reference and nlmixr2, and the gradient against the slope of the objective.


# ---- Models ---------------------------------------------------------------

# Every model of this file, generated with compile = FALSE and linked into one shared object on
# first use. The random effects come from effect formulas on branched transformations.
.focei_models <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    withr::local_preserve_seed()
    dir <- file.path(tempdir(), "focei_models")
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    withr::local_dir(dir)
    trafo <- function(pars, tab) branch(define(NULL, "x <- x", x = pars), table = tab)
    # second derivatives forward over forward and forward over reverse
    dm <- c("forward", "reverse", "forward-reverse")

    # one-compartment model with first-order absorption, dose per subject
    th <- datasets::Theoph
    th$Subject <- as.character(th$Subject)
    th_subj <- sort(unique(th$Subject))
    doses <- vapply(th_subj, function(s) { r <- th[th$Subject == s, ][1, ]; r$Dose * r$Wt }, 0.0)
    reactions <- eqnlist()
    reactions <- addReaction(reactions, "Ag", "",  "Ka * Ag",     "absorption")
    reactions <- addReaction(reactions, "",   "Cc", "Ka * Ag / V", "appearance")
    reactions <- addReaction(reactions, "Cc", "",  "Cl/V * Cc",   "elimination")
    th_x <- Xs(odemodel(reactions, modelname = "theoph_cppreg", compile = FALSE,
                        backend = "cppDE", deriv2 = TRUE, derivMode = c(dm, "forward-forward")))
    th_g <- Y(c(y = "Cc"), th_x, modelname = "theoph_cppreg_obs", compile = FALSE, deriv2 = TRUE, derivMode = dm)
    th_err <- Y(eqnvec(y = "sigma_add"), th_g, attachInput = FALSE, compile = FALSE, deriv2 = TRUE,
                modelname = "theoph_cppreg_err", derivMode = dm)
    th_tr <- trafo(c("Ka", "V", "Cl", "Ag", "Cc", "sigma_add"),
                   data.frame(ID = th_subj, dose = doses, row.names = th_subj))
    th_tr <- insert(th_tr, "Ag <- dose", dose = dose)
    th_tr <- insert(th_tr, "Cc <- 0")
    th_tr <- insert(th_tr, "x <- exp(x)", x = c("Ka", "V", "Cl", "sigma_add"))
    th_tr <- insert(th_tr, Ka + V + Cl ~ 1 + (1 | ID))
    th_p <- P(th_tr, method = "explicit", compile = FALSE, modelname = "theoph_formula_p", deriv2 = TRUE,
              derivMode = dm)

    # one-compartment iv bolus with proportional error, three subjects
    se_subj <- c("A", "B", "C")
    reactions <- addReaction(eqnlist(), "Cc", "", "Cl/V * Cc", "elimination")
    se_ode <- odemodel(reactions, modelname = "sigeta_ode", compile = FALSE, backend = "cppDE",
                       deriv2 = TRUE, derivMode = c(dm, "forward-forward"))
    se_x <- Xs(se_ode)
    se_g <- Y(c(y = "Cc"), se_x, modelname = "sigeta_obs", compile = FALSE, deriv2 = TRUE, derivMode = dm)
    se_err <- Y(eqnvec(y = "sigma_prop * y"), se_g, attachInput = FALSE, compile = FALSE,
                deriv2 = TRUE, modelname = "sigeta_err", derivMode = dm)
    se_tr <- trafo(c("V", "Cl", "Cc", "sigma_prop"),
                   data.frame(ID = se_subj, dose = 100, row.names = se_subj))
    se_tr <- insert(se_tr, "Cc <- dose / V", dose = dose)
    se_tr <- insert(se_tr, "x <- exp(x)", x = c("V", "Cl", "sigma_prop"))
    se_tr <- insert(se_tr, V ~ 1 + (1 | ID))
    se_p <- P(se_tr, method = "explicit", compile = FALSE, modelname = "sigeta_p", deriv2 = TRUE, derivMode = dm)
    # the same with an observable and an error model that read the time
    se_gt <- Y(c(y = "Cc * exp(-b_t * time)"), se_x, modelname = "sigeta_obs_t", compile = FALSE, deriv2 = TRUE,
               derivMode = dm)
    se_errt <- Y(eqnvec(y = "sigma_prop * y + s_t * time"), se_gt, attachInput = FALSE, compile = FALSE,
                 deriv2 = TRUE, modelname = "sigeta_err_t", derivMode = dm)
    se_trt <- trafo(c("V", "Cl", "Cc", "sigma_prop", "b_t", "s_t"),
                    data.frame(ID = se_subj, dose = 100, row.names = se_subj))
    se_trt <- insert(se_trt, "Cc <- dose / V", dose = dose)
    se_trt <- insert(se_trt, "x <- exp(x)", x = c("V", "Cl", "sigma_prop", "b_t", "s_t"))
    se_trt <- insert(se_trt, V ~ 1 + (1 | ID))
    se_pt <- P(se_trt, method = "explicit", compile = FALSE, modelname = "sigeta_p_t", deriv2 = TRUE,
               derivMode = dm)

    # parent / metabolite chain A -> B -> 0, both observed, three subjects
    mo_subj <- c("S1", "S2", "S3")
    reactions <- eqnlist()
    reactions <- addReaction(reactions, "A", "B", "ka * A", "absorption")
    reactions <- addReaction(reactions, "B", "",  "ke * B", "elimination")
    mo_x <- Xs(odemodel(reactions, modelname = "mo_ode", compile = FALSE, backend = "cppDE",
                        deriv2 = FALSE))
    mo_g <- Y(c(yA = "A", yB = "B"), mo_x, modelname = "mo_obs", compile = FALSE, deriv2 = FALSE)
    mo_err <- Y(eqnvec(yA = "sigA", yB = "sigB"), mo_g, attachInput = FALSE, compile = FALSE,
                modelname = "mo_err")
    mo_tr <- trafo(c("ka", "ke", "A", "B", "sigA", "sigB"),
                   data.frame(ID = mo_subj, row.names = mo_subj))
    mo_tr <- insert(mo_tr, "A <- 10")
    mo_tr <- insert(mo_tr, "B <- 0")
    mo_tr <- insert(mo_tr, "x <- exp(x)", x = c("ka", "ke", "sigA", "sigB"))
    mo_tr <- insert(mo_tr, ka + ke ~ 1 + (1 | ID))
    mo_p <- P(mo_tr, method = "explicit", compile = FALSE, modelname = "mo_p", deriv2 = FALSE)

    compile(th_x, th_g, th_err, th_p, se_x, se_g, se_err, se_p, se_gt, se_errt, se_pt, mo_x, mo_g, mo_err, mo_p,
            output = "focei_models", cores = 2)
    cache <<- list(
      theoph = list(x = th_x, g = th_g, err = th_err, p = th_p, data = as.datalist(data.frame(
        name = "y", time = th$Time, value = th$conc, sigma = NA_real_, condition = th$Subject))),
      sigeta = list(x = se_x, g = se_g, err = se_err, p = se_p, ode = se_ode, gt = se_gt, errt = se_errt,
                    pt = se_pt),
      mo     = list(x = mo_x, g = mo_g, err = mo_err, p = mo_p))
    cache
  }
})

# The kernel description of the groups of `obj` and the full parameter vector it indexes.
.kernelMeta <- function(obj, p, outer, etas) {
  design <- dMod2:::.emDesign(obj, randomEffects(p))
  meta <- dMod2:::.foceiMeta(design, outer, NULL)
  full <- c(outer, etas)[meta$pars_full_names]
  list(design = design, meta = meta, full = full)
}


# ---- Fit -------------------------------------------------------------------

test_that("the outer SR1 correction recovers an indefinite Hessian from the steps of a quadratic", {
  A <- matrix(c(2, 1, 0, 1, -1, 0.5, 0, 0.5, 3), 3, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  gradient <- function(x) setNames(drop(A %*% x) + c(1, 0, -1), rownames(A))
  base <- diag(3)
  points <- list(c(0, 0, 0), c(0.3, 0, 0), c(0.3, -0.2, 0.1), c(0.1, 0.4, 0.5))
  learn <- dMod2:::.outerSecant(rownames(A), memory = 3L)
  for (x in points) learn$add(setNames(x, rownames(A)), gradient(x), function(v) base %*% v)
  expect_equal(unname(base + learn$dense()), unname(A), tolerance = 1e-10)
  # the low-rank factors are those of the dense correction
  f <- learn$factors()
  expect_equal(unname(f$U %*% (f$d * t(f$U))), unname(learn$dense()), tolerance = 1e-12)
  # only the last pairs enter: with a memory of 2 the first is forgotten
  short <- dMod2:::.outerSecant(rownames(A), memory = 2L)
  for (x in points) short$add(setNames(x, rownames(A)), gradient(x), function(v) base %*% v)
  expect_equal(ncol(short$factors()$U), 2L)
})

test_that("em() rejects an unknown method", {
  expect_error(em(obj = NULL, random = NULL, method = "doesNotExist"), "focei")
})

test_that("FOCEI reproduces the Theoph reference fit and agrees with nlmixr2", {
  skip_on_cran()
  ref <- readRDS("fixtures/focei_theoph_reference.rds")
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  mods <- .focei_models()$theoph
  obj <- normL2(mods$data, mods$g * mods$x * mods$p, errmodel = mods$err)
  map <- c(Ka = "tka", V = "tv", Cl = "tcl", sigma_add = "log_sigma_add",
           omega_Ka_Ka = "omega_Ka_Ka", omega_V_V = "omega_V_V", omega_Cl_Cl = "omega_Cl_Cl")
  fit <- em(obj, randomEffects(mods$p), init = setNames(ref$init[map], names(map)), verbose = FALSE,
            control = list(innerControl = list(iterlim = 30, ftol = 1e-7, mtol = 1e-7),
                           trustControl = list(iterlim = 100)))
  expect_true(fit$converged)
  expect_lt(abs(fit$value - ref$value), 1e-2)
  expect_lt(max(abs(fit$argument[names(map)] - ref$argument[map])), 5e-3)
  nl <- ref$nlmixr2
  inv <- setNames(names(map), map)
  expect_true(all(abs(fit$argument[inv[names(nl)]] - nl) < 0.25 * ref$se[names(nl)]))
  expect_lt(abs(fit$value - ref$nlmixr2Value), 0.5)
})

test_that("the FOCEI gradient with a random effect in the error model is the slope of the objective", {
  skip_on_cran()
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  mods <- .focei_models()$sigeta
  set.seed(5)
  # a tight solve, so that the slope is not that of the solver error
  model <- mods$g * Xs(mods$ode, options = list(atol = 1e-12, rtol = 1e-12)) * mods$p
  truth <- c(V = log(20), Cl = log(5), sigma_prop = log(0.1))
  eta <- c(eta_V_A = -0.2, eta_V_B = 0.1, eta_V_C = 0.05)
  sim <- model(c(0.5, 1, 2, 4, 8), c(truth, eta), deriv = FALSE)
  data <- as.datalist(do.call(rbind, lapply(names(sim), function(s) data.frame(
    name = "y", time = sim[[s]][, "time"], value = sim[[s]][, "y"] * (1 + 0.1 * rnorm(nrow(sim[[s]]))),
    sigma = NA_real_, condition = s))))
  obj <- normL2(data, model, errmodel = mods$err)
  init <- c(truth + 0.05, omega_V_V = log(0.3))
  ctl <- list(trustControl = list(iterlim = 0),
              innerControl = list(ftol = 1e-12, mtol = 1e-12, iterlim = 200))
  run <- function(p) em(obj, randomEffects(mods$p), setNames(p, names(init)), control = ctl, verbose = FALSE)
  at <- run(init)
  slope <- numDeriv::grad(function(p) run(p)$value, unname(init),
                          method.args = list(eps = 1e-4, d = 1e-3, r = 2))
  expect_equal(unname(at$gradient[names(init)]), slope, tolerance = 1e-3)
  # the sparse outer path: the same value, gradient and Hessian, and the same fit
  sp <- em(obj, randomEffects(mods$p), init, control = ctl, verbose = FALSE, sparse = TRUE)
  expect_equal(sp$value, at$value, tolerance = 1e-10)
  expect_equal(sp$gradient, at$gradient, tolerance = 1e-8)
  expect_equal(as.matrix(sp$hessian)[names(init), names(init)], at$hessian[names(init), names(init)],
               tolerance = 1e-8)
  ctl$trustControl$iterlim <- 100
  fd <- em(obj, randomEffects(mods$p), init, control = ctl, verbose = FALSE, sparse = FALSE)
  fs <- em(obj, randomEffects(mods$p), init, control = ctl, verbose = FALSE, sparse = TRUE)
  expect_equal(fs$value, fd$value, tolerance = 1e-6)
  expect_equal(fs$argument, fd$argument, tolerance = 1e-4)
  # forward over reverse, the default: the volume gradient and the mode shift of forward over
  # forward, and the same fit
  expect_identical(fd$derivMode, "forward-reverse")
  design <- .emDesign(obj, randomEffects(mods$p))
  E <- design$etaNames
  full <- c(init, setNames(c(-0.1, 0.05, 0.2), as.vector(E)))
  Hinv <- lapply(design$groups, function(g) matrix(0.4, 1, 1, dimnames = list(E[g, ], E[g, ])))
  ffV <- .foceiVolume(design, full, NULL, names(init), Hinv)
  frV <- .volumeReverse(design, full, NULL, names(init), Hinv)
  # the adjoint and the forward sensitivities differentiate discretizations that differ by the
  # tolerance of the solver
  expect_equal(frV$gradient, ffV$gradient, tolerance = 1e-3)
  expect_equal(frV$etaShift, ffV$etaShift, tolerance = 1e-3)
  # the sweeps through the native chain and through R
  V <- matrix(1, length(full), 1, dimnames = list(names(full), NULL))
  expect_false(is.null(.sweepChainForward(design, full, NULL, V, names(design$data), 1L, 1L)))
  frC <- .volumeReverse(design, full, NULL, names(init), Hinv, curvature = TRUE)
  frR <- withr::with_options(list(dMod2.evaluation = "R"),
                             .volumeReverse(design, full, NULL, names(init), Hinv, curvature = TRUE))
  expect_equal(frC, frR, tolerance = 1e-10)
  expect_equal(.schurReverse(design, full, NULL, names(init)),
               withr::with_options(list(dMod2.evaluation = "R"), .schurReverse(design, full, NULL, names(init))),
               tolerance = 1e-10)
  ff <- em(obj, randomEffects(mods$p), init, control = modifyList(ctl, list(derivMode = "forward-forward")),
            verbose = FALSE)
  expect_identical(ff$derivMode, "forward-forward")
  expect_equal(fd$value, ff$value, tolerance = 1e-6)
  expect_equal(fd$argument, ff$argument, tolerance = 1e-4)
  expect_equal(fd$hessian, ff$hessian, tolerance = 1e-4)
})


test_that("the native chain takes observables and error models that read the time", {
  skip_on_cran()
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  mods <- .focei_models()$sigeta
  model <- mods$gt * mods$x * mods$pt
  pars <- c(V = log(20), Cl = log(5), sigma_prop = log(0.1), b_t = log(0.05), s_t = log(0.01),
            eta_V_A = -0.2, eta_V_B = 0.1, eta_V_C = 0.05)
  sim <- model(c(0.5, 1, 2, 4, 8), pars, deriv = FALSE)
  set.seed(6)
  data <- as.datalist(do.call(rbind, lapply(names(sim), function(s) data.frame(
    name = "y", time = sim[[s]][, "time"], value = sim[[s]][, "y"] * (1 + 0.1 * rnorm(nrow(sim[[s]]))),
    sigma = NA_real_, condition = s))))
  obj <- normL2(data, model, errmodel = mods$errt)
  viaR <- function(expr) withr::with_options(list(dMod2.evaluation = "R"), expr)
  a <- obj(pars, deriv2 = TRUE); b <- viaR(obj(pars, deriv2 = TRUE))
  expect_equal(a$value, b$value, tolerance = 1e-12)
  expect_equal(a$gradient, b$gradient, tolerance = 1e-10)
  expect_equal(a$hessian, b$hessian, tolerance = 1e-10)
  entries <- environment(obj)$.chain$entries
  expect_false(any(vapply(ls(entries), function(k) isTRUE(get(k, envir = entries)[[1]]$none), TRUE)))
  # the adjoint gradient
  expect_equal(obj(pars, sweep = "reverse")$gradient, viaR(obj(pars, sweep = "reverse"))$gradient,
               tolerance = 1e-10)
  # the sweeps of forward over reverse
  design <- .emDesign(obj, randomEffects(mods$pt))
  outer <- c(pars[1:5], omega_V_V = log(0.3))
  full <- c(outer, pars[6:8])
  V <- matrix(1, length(full), 1, dimnames = list(names(full), NULL))
  expect_false(is.null(.sweepChainForward(design, full, NULL, V, names(design$data), 1L, 1L)))
  Hinv <- lapply(design$groups, function(g) matrix(0.4, 1, 1, dimnames = list(design$etaNames[g, ], design$etaNames[g, ])))
  expect_equal(.volumeReverse(design, full, NULL, names(outer), Hinv, curvature = TRUE),
               viaR(.volumeReverse(design, full, NULL, names(outer), Hinv, curvature = TRUE)), tolerance = 1e-10)
})


test_that("the volume rows give the exact joint blocks and d eta*/d theta is the slope of the modes", {
  skip_on_cran()
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  mods <- .focei_models()$sigeta
  set.seed(5)
  # a tight solve, so that the difference quotient of the modes is not that of the solver error
  model <- mods$g * Xs(mods$ode, options = list(atol = 1e-12, rtol = 1e-12)) * mods$p
  truth <- c(V = log(20), Cl = log(5), sigma_prop = log(0.1))
  sim <- model(c(0.5, 1, 2, 4, 8), c(truth, eta_V_A = -0.2, eta_V_B = 0.1, eta_V_C = 0.05), deriv = FALSE)
  data <- as.datalist(do.call(rbind, lapply(names(sim), function(s) data.frame(
    name = "y", time = sim[[s]][, "time"], value = sim[[s]][, "y"] * (1 + 0.1 * rnorm(nrow(sim[[s]]))),
    sigma = NA_real_, condition = s))))
  obj <- normL2(data, model, errmodel = mods$err)
  theta <- c(truth + 0.05, omega_V_V = log(0.3))
  ctl <- list(trustControl = list(iterlim = 0), innerControl = list(ftol = 1e-12, mtol = 1e-12, iterlim = 200))
  modes <- function(th) em(obj, randomEffects(mods$p), setNames(th, names(theta)), control = ctl,
                           verbose = FALSE)$etaModes[, 1]
  design <- dMod2:::.emDesign(obj, randomEffects(mods$p))
  full <- c(theta, setNames(modes(theta), design$etaNames[, 1]))
  unit <- dMod2:::.volumeUnits(design, full, NULL)[[1]]
  cn <- unit$conditions
  pred <- design$prdfn(times = dMod2:::.emGrid(design, cn), pars = full, deriv = TRUE, deriv2 = TRUE,
                       conditions = cn)
  frames <- dMod2:::.volumeFrameList(design, cn, setNames(as.list(pred), cn), NULL)
  A <- names(theta)
  rw <- dMod2:::.volumeRows(design, cn, frames[[cn]], unit$eta, A, FALSE)
  lr <- dMod2:::.rowLossDerivs(rw$p, rw$s, rw$y, rw$cens)
  jb <- dMod2:::.jointBlocks(rw, lr, unit, full, A, "omega_V_V")
  jh <- (design$objfn + dMod2:::.normalPrior(design))(pars = full, deriv = TRUE, deriv2 = TRUE)$hessian
  expect_equal(jb$Hee, jh[unit$eta, unit$eta, drop = FALSE], tolerance = 1e-10)
  expect_equal(jb$HeA, jh[unit$eta, A, drop = FALSE], tolerance = 1e-10)

  Hinv <- lapply(dMod2:::.volumeUnits(design, full, NULL), function(u) solve(jh[u$eta, u$eta, drop = FALSE]))
  shift <- dMod2:::.foceiVolume(design, full, NULL, names(theta), Hinv)$etaShift
  slope <- numDeriv::jacobian(modes, unname(theta), method.args = list(eps = 1e-4, d = 1e-3, r = 2))
  expect_equal(unname(do.call(rbind, shift)), slope, tolerance = 1e-4)
})


test_that("the inner search stops on the Newton step of the modes, not on a change of the value", {
  skip_on_cran()
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  mods <- .focei_models()$sigeta
  set.seed(5)
  model <- mods$g * mods$x * mods$p
  truth <- c(V = log(20), Cl = log(5), sigma_prop = log(0.1))
  sim <- model(c(0.5, 1, 2, 4, 8), c(truth, eta_V_A = -0.2, eta_V_B = 0.1, eta_V_C = 0.05), deriv = FALSE)
  data <- as.datalist(do.call(rbind, lapply(names(sim), function(s) data.frame(
    name = "y", time = sim[[s]][, "time"], value = sim[[s]][, "y"] * (1 + 0.1 * rnorm(nrow(sim[[s]]))),
    sigma = NA_real_, condition = s))))
  obj <- normL2(data, model, errmodel = mods$err)
  init <- c(truth + 0.05, omega_V_V = log(0.3))
  # a loose value tolerance ends only steps on the boundary of the region
  fit <- em(obj, randomEffects(mods$p), init, verbose = FALSE,
            control = list(innerControl = list(mtol = 1e-3), trustControl = list(iterlim = 0)))
  design <- .emDesign(obj, randomEffects(mods$p))
  En <- design$etaNames
  E <- as.vector(En)
  full <- c(init, setNames(as.vector(fit$etaModes[rownames(En), , drop = FALSE]), E))
  data <- obj(full, deriv = TRUE)
  prior <- .normalPrior(design)(pars = full, deriv = TRUE)
  step <- solve(data$hessian[E, E] + prior$hessian[E, E], data$gradient[E] + prior$gradient[E])
  expect_lt(max(abs(step)), 1e-8)
})


test_that("FOCEI does not end where failed trials shrank the region", {
  skip_on_cran()
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  mods <- .focei_models()$sigeta
  set.seed(5)
  model <- mods$g * mods$x * mods$p
  truth <- c(V = log(20), Cl = log(5), sigma_prop = log(0.1))
  sim <- model(c(0.5, 1, 2, 4, 8), c(truth, eta_V_A = -0.2, eta_V_B = 0.1, eta_V_C = 0.05), deriv = FALSE)
  data <- as.datalist(do.call(rbind, lapply(names(sim), function(s) data.frame(
    name = "y", time = sim[[s]][, "time"], value = sim[[s]][, "y"] * (1 + 0.1 * rnorm(nrow(sim[[s]]))),
    sigma = NA_real_, condition = s))))
  obj <- normL2(data, model, errmodel = mods$err)
  init <- c(truth + 0.3, omega_V_V = log(0.3))
  # the prediction fails away from the start, so only trials of a collapsed region succeed
  real <- .batchPredictor
  local_mocked_bindings(.batchPredictor = function(design) {
    f <- real(design)
    function(timesList, parsList, conditions, fixed) {
      if (any(vapply(parsList, function(p) max(abs(p[names(init)] - init)), 0) > 1e-7)) stop("fails away")
      f(timesList, parsList, conditions, fixed)
    }
  })
  fit <- em(obj, randomEffects(mods$p), init, verbose = FALSE, sparse = FALSE,
            control = list(derivMode = "forward-forward", trustControl = list(iterlim = 40)))
  expect_false(fit$converged)
})


# ---- C++ kernel parity -------------------------------------------------------

test_that("the inner kernel has the value, gradient and Gauss-Newton Hessian of the R replica for sigma(eta)", {
  skip_on_cran()
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  subjects <- c("A", "B", "C")
  times <- c(0.5, 1, 2, 4, 8)
  set.seed(123)
  pred_at <- function(t, V, Cl) (100 / V) * exp(-Cl / V * t)
  data_rows <- do.call(rbind, lapply(seq_along(subjects), function(i) {
    pred <- pred_at(times, 20 * exp(c(-0.2, 0.1, 0.05)[i]), 5)
    data.frame(name = "y", time = times, value = pred + rnorm(length(times), 0, 0.1 * pred),
               sigma = NA_real_, condition = subjects[i], stringsAsFactors = FALSE)
  }))
  mods <- .focei_models()$sigeta
  model <- mods$g * mods$x * mods$p
  err <- mods$err
  obj <- normL2(as.datalist(data_rows), model, errmodel = err)
  outer <- c(V = log(20), Cl = log(5), sigma_prop = log(0.1), omega_V_V = log(0.3))
  eta_vals <- c(eta_V_A = -0.15, eta_V_B = 0.05, eta_V_C = 0.1)
  km <- .kernelMeta(obj, mods$p, outer, eta_vals)
  Omega_inv <- matrix(1 / 0.3^2)
  Omega_log_det <- log(0.3^2)
  res <- dMod2:::focei_inner_trust(
    model_cb = model, err_cb = err, pars_full = km$full,
    eta_warmstart = matrix(eta_vals, ncol = 1), subject_meta = km$meta,
    Omega_inv_mat = Omega_inv, Omega_log_det = Omega_log_det, fixed = NULL,
    control = list(iterlim = 0))

  replica <- function(eta_block, s) {
    pars <- km$full
    pars[paste0("eta_V_", s)] <- eta_block
    prdf <- model(times = times, pars = pars, deriv = TRUE, conditions = s)[[1]]
    erm <- err(out = prdf, pars = attr(prdf, "parameters"), conditions = s)[[1]]
    d <- data_rows[data_rows$condition == s, ]
    tp <- match(d$time, prdf[, "time"]); te <- match(d$time, erm[, "time"])
    nm <- paste0("eta_V_", s)
    p <- prdf[tp, "y"]; sig <- erm[te, "y"]
    Jp <- attr(prdf, "deriv")[tp, "y", nm]; Js <- attr(erm, "deriv")[te, "y", nm]
    wr <- (p - d$value) / sig
    dwr <- (Jp - wr * Js) / sig
    value <- sum(wr^2 + log(2 * pi * sig^2)) + eta_block^2 * Omega_inv[1] + Omega_log_det
    grad <- 2 * sum(wr * dwr + Js / sig) + 2 * Omega_inv[1] * eta_block
    hess <- 2 * sum(dwr^2) + sum(-4 * wr * Jp * Js / sig^2) + sum(4 * wr^2 * Js^2 / sig^2) +
      sum(-2 * Js^2 / sig^2) + 2 * Omega_inv[1]
    list(value = value, grad = grad, hess = hess)
  }
  # the batched predictor of the lock-step trust
  batched <- dMod2:::focei_inner_trust(
    model_cb = model, err_cb = err, pars_full = km$full,
    eta_warmstart = matrix(eta_vals, ncol = 1), subject_meta = km$meta,
    Omega_inv_mat = Omega_inv, Omega_log_det = Omega_log_det, fixed = NULL,
    control = list(iterlim = 0), predict_cb = dMod2:::.batchPredictor(km$design))
  for (i in seq_along(subjects)) {
    s <- subjects[i]
    R <- replica(eta_vals[[i]], s)
    expect_equal(unname(res$value[[s]]), R$value, tolerance = 1e-12)
    expect_equal(as.numeric(res$gradient[s, ]), R$grad, tolerance = 1e-12)
    expect_equal(as.numeric(res$H_GN[[s]]), R$hess, tolerance = 1e-12)
    expect_equal(unname(batched$value[[s]]), R$value, tolerance = 1e-12)
    expect_equal(as.numeric(batched$gradient[s, ]), R$grad, tolerance = 1e-12)
    expect_equal(as.numeric(batched$H_GN[[s]]), R$hess, tolerance = 1e-12)
  }
})

test_that("a failed prediction at an inner trial point rejects the trial of that subject only", {
  skip_on_cran()
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  subjects <- c("A", "B", "C")
  times <- c(0.5, 1, 2, 4, 8)
  set.seed(123)
  pred_at <- function(t, V, Cl) (100 / V) * exp(-Cl / V * t)
  data_rows <- do.call(rbind, lapply(seq_along(subjects), function(i) {
    pred <- pred_at(times, 20 * exp(c(-0.2, 0.1, 0.05)[i]), 5)
    data.frame(name = "y", time = times, value = pred + rnorm(length(times), 0, 0.1 * pred),
               sigma = NA_real_, condition = subjects[i], stringsAsFactors = FALSE)
  }))
  mods <- .focei_models()$sigeta
  model <- mods$g * mods$x * mods$p
  obj <- normL2(as.datalist(data_rows), model, errmodel = mods$err)
  outer <- c(V = log(20), Cl = log(5), sigma_prop = log(0.1), omega_V_V = log(0.3))
  km <- .kernelMeta(obj, mods$p, outer, c(eta_V_A = 0, eta_V_B = 0, eta_V_C = 0))
  predict <- dMod2:::.batchPredictor(km$design)
  # the prediction of one subject fails at the first trial value of the subject eta, alone or in a batch
  bad <- NULL
  failing <- function(timesList, parsList, conditions, fixed) {
    for (q in which(conditions == "B")) {
      e <- parsList[[q]][["eta_V_B"]]
      if (e != 0 && is.null(bad)) bad <<- e
      if (identical(e, bad)) stop("1 of ", length(conditions), " conditions did not complete: B")
    }
    predict(timesList, parsList, conditions, fixed)
  }
  modes <- function(cb) dMod2:::focei_inner_trust(
    model_cb = model, err_cb = mods$err, pars_full = km$full, eta_warmstart = matrix(0, 3, 1),
    subject_meta = km$meta, Omega_inv_mat = matrix(1 / 0.3^2), Omega_log_det = log(0.3^2), fixed = NULL,
    control = list(iterlim = 50, ftol = 1e-12, mtol = 1e-12), predict_cb = cb)
  failure <- new.env(parent = emptyenv())
  ref <- modes(predict)
  got <- modes(dMod2:::.guardPredictor(failing, failure))
  expect_false(is.null(bad))
  expect_match(failure$message, "did not complete: B")
  expect_true(all(got$converged))
  # the other subjects take the same steps; the one with the failure ends within the inner tolerance
  expect_identical(got$eta_modes[c("A", "C"), ], ref$eta_modes[c("A", "C"), ])
  expect_lt(abs(got$eta_modes["B", 1] - ref$eta_modes["B", 1]), 1e-5)
})

test_that("an inner trial point whose derivatives are not finite rejects the trial of that subject only", {
  skip_on_cran()
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  subjects <- c("A", "B", "C")
  times <- c(0.5, 1, 2, 4, 8)
  set.seed(123)
  pred_at <- function(t, V, Cl) (100 / V) * exp(-Cl / V * t)
  data_rows <- do.call(rbind, lapply(seq_along(subjects), function(i) {
    pred <- pred_at(times, 20 * exp(c(-0.2, 0.1, 0.05)[i]), 5)
    data.frame(name = "y", time = times, value = pred + rnorm(length(times), 0, 0.1 * pred),
               sigma = NA_real_, condition = subjects[i], stringsAsFactors = FALSE)
  }))
  mods <- .focei_models()$sigeta
  model <- mods$g * mods$x * mods$p
  obj <- normL2(as.datalist(data_rows), model, errmodel = mods$err)
  outer <- c(V = log(20), Cl = log(5), sigma_prop = log(0.1), omega_V_V = log(0.3))
  km <- .kernelMeta(obj, mods$p, outer, c(eta_V_A = 0, eta_V_B = 0, eta_V_C = 0))
  predict <- dMod2:::.batchPredictor(km$design)
  # the sensitivities of one subject are NaN at the first trial value of the subject eta
  bad <- NULL
  poisoned <- function(timesList, parsList, conditions, fixed) {
    out <- predict(timesList, parsList, conditions, fixed)
    for (q in which(conditions == "B")) {
      e <- parsList[[q]][["eta_V_B"]]
      if (e != 0 && is.null(bad)) bad <<- e
      if (identical(e, bad)) attr(out$pred[[q]], "deriv")[] <- NaN
    }
    out
  }
  modes <- function(cb) dMod2:::focei_inner_trust(
    model_cb = model, err_cb = mods$err, pars_full = km$full, eta_warmstart = matrix(0, 3, 1),
    subject_meta = km$meta, Omega_inv_mat = matrix(1 / 0.3^2), Omega_log_det = log(0.3^2), fixed = NULL,
    control = list(iterlim = 50, ftol = 1e-12, mtol = 1e-12), predict_cb = cb)
  ref <- modes(predict)
  got <- modes(poisoned)
  expect_false(is.null(bad))
  expect_true(all(got$converged))
  expect_identical(got$eta_modes[c("A", "C"), ], ref$eta_modes[c("A", "C"), ])
  expect_lt(abs(got$eta_modes["B", 1] - ref$eta_modes["B", 1]), 1e-5)
})

test_that("the inner kernel matches the R replica for two observables per subject", {
  skip_on_cran()
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  subjects <- c("S1", "S2", "S3")
  times <- c(0.5, 1, 2, 4, 8)
  pred_at <- function(t, ka, ke) {
    A <- 10 * exp(-ka * t)
    B <- (ka * 10 / (ke - ka)) * (exp(-ka * t) - exp(-ke * t))
    list(A = A, B = B)
  }
  set.seed(321)
  data_rows <- do.call(rbind, lapply(seq_along(subjects), function(i) {
    pred <- pred_at(times, 0.7 * exp(c(-0.2, 0.15, 0.05)[i]), 0.3 * exp(c(0.1, -0.05, 0)[i]))
    rbind(data.frame(name = "yA", time = times, value = pred$A + rnorm(length(times), 0, 0.15),
                     sigma = NA_real_, condition = subjects[i]),
          data.frame(name = "yB", time = times, value = pred$B + rnorm(length(times), 0, 0.1),
                     sigma = NA_real_, condition = subjects[i]))
  }))
  mods <- .focei_models()$mo
  model <- mods$g * mods$x * mods$p
  err <- mods$err
  obj <- normL2(as.datalist(data_rows), model, errmodel = err)
  outer <- c(ka = log(0.7), ke = log(0.3), sigA = log(0.15), sigB = log(0.1),
             omega_ka_ka = log(0.3), omega_ke_ke = log(0.3))
  eta_vals <- c(eta_ka_S1 = -0.15, eta_ka_S2 = 0.10, eta_ka_S3 = 0,
                eta_ke_S1 = 0.05, eta_ke_S2 = -0.02, eta_ke_S3 = 0)
  km <- .kernelMeta(obj, mods$p, outer, eta_vals)
  K <- 2L
  Omega_inv <- diag(1 / 0.3^2, K)
  Omega_log_det <- 2 * log(0.3^2)
  res <- dMod2:::focei_inner_trust(
    model_cb = model, err_cb = err, pars_full = km$full,
    eta_warmstart = matrix(eta_vals, ncol = K), subject_meta = km$meta,
    Omega_inv_mat = Omega_inv, Omega_log_det = Omega_log_det, fixed = NULL,
    control = list(iterlim = 0))

  replica <- function(eta_block, s) {
    pars <- km$full
    nms <- paste0(c("eta_ka_", "eta_ke_"), s)
    pars[nms] <- eta_block
    prdf <- model(times = times, pars = pars, deriv = TRUE, conditions = s)[[1]]
    erm <- err(out = prdf, pars = attr(prdf, "parameters"), conditions = s)[[1]]
    d <- data_rows[data_rows$condition == s, ]
    tp <- match(d$time, prdf[, "time"]); te <- match(d$time, erm[, "time"])
    value <- 0; grad <- numeric(K); hess <- matrix(0, K, K)
    for (j in seq_len(nrow(d))) {
      o <- d$name[j]
      sig <- erm[te[j], o]
      Jp <- attr(prdf, "deriv")[tp[j], o, nms]
      Js <- if (!is.null(attr(erm, "deriv")) && o %in% dimnames(attr(erm, "deriv"))[[2]])
        attr(erm, "deriv")[te[j], o, nms] else numeric(K)
      wr <- (prdf[tp[j], o] - d$value[j]) / sig
      dwr <- (Jp - wr * Js) / sig
      value <- value + wr^2 + log(2 * pi * sig^2)
      grad <- grad + 2 * (wr * dwr + Js / sig)
      hess <- hess + 2 * tcrossprod(dwr) - 2 * wr / sig^2 * (tcrossprod(Jp, Js) + tcrossprod(Js, Jp)) +
        (4 * wr^2 - 2) / sig^2 * tcrossprod(Js)
    }
    list(value = value + drop(t(eta_block) %*% Omega_inv %*% eta_block) + Omega_log_det,
         grad = grad + 2 * drop(Omega_inv %*% eta_block), hess = hess + 2 * Omega_inv)
  }
  for (i in seq_along(subjects)) {
    s <- subjects[i]
    R <- replica(eta_vals[paste0(c("eta_ka_", "eta_ke_"), s)], s)
    expect_equal(unname(res$value[[s]]), unname(R$value), tolerance = 1e-10)
    expect_equal(as.numeric(res$gradient[s, ]), unname(R$grad), tolerance = 1e-10)
    expect_equal(as.numeric(res$H_GN[[s]]), as.numeric(R$hess), tolerance = 1e-10)
  }
})

test_that("an eta block with a vanishing variance is solved on the unit diagonal", {
  A <- matrix(c(5e18, 1, 0.5, 1, 4, 1, 0.5, 1, 3), 3, dimnames = rep(list(c("a", "b", "c")), 2))
  B <- matrix(1:6, 3)
  expect_error(solve(A, B), "singular")
  x <- dMod2:::.scaledSolve(A, B)
  expect_equal(A %*% x, B, ignore_attr = TRUE, tolerance = 1e-10)
  expect_equal(dMod2:::.scaledSolve(A)[2:3, 2:3], solve(A[2:3, 2:3]), tolerance = 1e-10)
})
