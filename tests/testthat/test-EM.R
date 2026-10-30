# em() on a linear Gaussian model with a closed-form marginal likelihood: y = a + b dose with a
# correlated random intercept and slope per group, two conditions per group and one condition
# without a group.

.lin <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    withr::local_preserve_seed()
    dir <- file.path(tempdir(), "nlmefit_lin")
    dir.create(dir, showWarnings = FALSE)
    withr::local_dir(dir)
    groups <- paste0("s", 1:8)
    tab <- data.frame(ID = c(rep(groups, each = 2), NA), dose = c(rep(c(1, 2), 8), 1),
                      row.names = c(paste0(rep(groups, each = 2), c("a", "b")), "ref"))
    trafo <- function(pars) {
      tr <- branch(define(NULL, "x <- x", x = pars), table = tab)
      tr <- insert(tr, a ~ 1 + (1 | p | ID), b ~ 1 + (1 | p | ID), conditionMatch = "^s")
      define(tr, "dose <- d", d = dose)
    }
    # second derivatives forward over forward and forward over reverse
    dm <- c("forward", "reverse", "forward-reverse")
    g <- Y(c(y = "a + b * dose"), f = NULL, parameters = c("a", "b", "dose"), compile = FALSE,
           modelname = "nlmefit_lin_g", deriv2 = TRUE, derivMode = dm)
    p <- P(trafo(c("a", "b")), method = "explicit", compile = FALSE, modelname = "nlmefit_lin_p", deriv2 = TRUE,
           derivMode = dm)
    # with the parameter s of an error model sigma = exp(s)
    ps <- P(trafo(c("a", "b", "s")), method = "explicit", compile = FALSE, modelname = "nlmefit_lin_ps",
            deriv2 = TRUE, derivMode = dm)
    err <- Y(c(y = "exp(s)"), g, attachInput = FALSE, compile = FALSE, modelname = "nlmefit_lin_err",
             deriv2 = TRUE, derivMode = dm)
    # second derivatives forward over forward only
    pf <- P(trafo(c("a", "b")), method = "explicit", compile = FALSE, modelname = "nlmefit_lin_pf", deriv2 = TRUE)
    compile(g, p, ps, err, pf, output = "nlmefit_lin", cores = 2)
    set.seed(3)
    eta <- cbind(rnorm(8, 0, 0.7), rnorm(8, 0, 0.4))
    rows <- lapply(rownames(tab), function(cn) {
      i <- match(tab[cn, "ID"], groups)
      e <- if (is.na(i)) c(0, 0) else eta[i, ]
      mu <- (1 + e[1]) + (0.5 + e[2]) * tab[cn, "dose"]
      data.frame(name = "y", time = c(0, 1, 2), value = mu + rnorm(3, 0, 0.3), sigma = 0.3,
                 condition = cn)
    })
    data <- as.datalist(do.call(rbind, rows))
    cache <<- list(g = g, p = p, ps = ps, pf = pf, err = err, data = data, tab = tab, groups = groups,
                   obj = normL2(data, g * Xt() * p))
    cache
  }
})

# -2 log of the closed-form marginal: per group a normal of its stacked rows with covariance
# sigma^2 I + Z Omega Z', the condition without a group with sigma^2 I.
.linMarginal <- function(pars, mod) {
  m <- mod
  Om <- local({
    L <- matrix(c(exp(pars[["omega_a_a"]]), pars[["omega_b_a"]], 0, exp(pars[["omega_b_b"]])), 2)
    tcrossprod(L)
  })
  out <- 0
  for (g in c(m$groups, "ref")) {
    cn <- if (g == "ref") "ref" else paste0(g, c("a", "b"))
    d <- do.call(rbind, lapply(cn, function(c) cbind(m$data[[c]], dose = m$tab[c, "dose"])))
    Z <- cbind(1, d$dose)
    C <- diag(d$sigma^2, nrow(d)) + if (g == "ref") 0 else Z %*% Om %*% t(Z)
    r <- d$value - (pars[["a"]] + pars[["b"]] * d$dose)
    out <- out + drop(t(r) %*% solve(C, r)) + as.numeric(determinant(C)$modulus) +
      nrow(d) * log(2 * pi)
  }
  out
}

test_that("the design collects the conditions of every group and the hyperparameters", {
  m <- .lin()
  d <- dMod2:::.emDesign(m$obj, randomEffects(m$p))
  expect_equal(d$groups, m$groups)
  expect_equal(d$conditions$s2, c("s2a", "s2b"))
  expect_equal(d$freeConditions, "ref")
  expect_equal(unname(d$etaNames["s3", ]), c("eta_a_s3", "eta_b_s3"))
  expect_equal(d$normal$cholPars, c("omega_a_a", "omega_b_a", "omega_b_b"))
  start <- dMod2:::.emStart(d, NULL, NULL, sd = 0.5)
  expect_equal(start[c("a", "b", "omega_a_a", "omega_b_a")], c(a = 0, b = 0, omega_a_a = log(0.5), omega_b_a = 0))
  expect_error(dMod2:::.emStart(d, c(a = 1, zz = 2)), "does not have: zz")
})

test_that("mu-referenced fixed effects: random intercepts on every condition, rows constant per group", {
  ids <- c("s1", "s2", "s3")
  tab <- data.frame(ID = rep(ids, each = 2), wt = rep(c(60, 70, 85), each = 2), dose = rep(1:2, 3),
                    row.names = paste0(rep(ids, each = 2), c("a", "b")))
  tr <- branch(define(NULL, "x <- x", x = c("ka", "cl", "v")), table = tab)
  tr <- insert(tr, ka ~ 1 + (1 | ID), cl ~ wt + (1 + dose | p | ID), v ~ dose + (1 | ID))
  re <- randomEffects(tr)
  design <- list(randomEffects = re, data = setNames(vector("list", 6), rownames(tab)))
  sub <- .groupingDesign(re, design$data, "ID")
  ref <- .muReference(design, sub)
  expect_named(ref, c("eta_ka", "eta_cl"))
  expect_equal(ref$eta_cl$X, cbind(cl = c(s1 = 1, s2 = 1, s3 = 1), cl_wt = c(60, 70, 85)))
  # a condition with data but without the formula leaves its symbol unreferenced
  design$data <- c(design$data, list(ref = NULL))
  design$randomEffects$etas["ref", ] <- NA
  expect_length(.muReference(design, sub), 0)

  set.seed(1)
  stems <- sub$normal$stems
  phi <- matrix(rnorm(12), 3, 4, dimnames = list(ids, stems))
  Omega <- crossprod(matrix(rnorm(16), 4)) + diag(4)
  loss <- function(b) sum(vapply(seq_along(ids), function(i) {
    r <- phi[i, ]
    r["eta_ka"] <- r["eta_ka"] - b[["ka"]]
    r["eta_cl"] <- r["eta_cl"] - b[["cl"]] - b[["cl_wt"]] * c(60, 70, 85)[i]
    drop(r %*% solve(Omega, r))
  }, 0))
  beta <- .glsFixedEffects(ref, phi, Omega, c(ka = 0, cl = 0, cl_wt = 0))
  expect_lt(max(abs(numDeriv::grad(loss, beta))), 1e-8)
  part <- .glsFixedEffects(ref, phi, Omega, c(ka = 0, cl = 0, cl_wt = 0.01), free = c("ka", "cl"))
  at <- function(b) loss(c(b, cl_wt = 0.01))
  expect_lt(max(abs(numDeriv::grad(at, part))), 1e-8)
})

test_that("the normal prior has the exact gradient and Hessian, and the M-step is its optimum", {
  m <- .lin()
  d <- dMod2:::.emDesign(m$obj, randomEffects(m$p))
  prior <- dMod2:::.normalPrior(d)
  set.seed(1)
  pars <- setNames(rnorm(length(attr(prior, "parameters")), 0, 0.5), attr(prior, "parameters"))
  out <- prior(pars = pars, deriv = TRUE)
  f <- function(v) prior(pars = setNames(v, names(pars)), deriv = FALSE)$value
  g <- function(v) prior(pars = setNames(v, names(pars)), deriv = TRUE)$gradient
  expect_equal(unname(out$gradient), numDeriv::grad(f, unname(pars)), tolerance = 1e-7)
  expect_equal(unname(out$hessian), numDeriv::jacobian(g, unname(pars)), tolerance = 1e-7)
  sp <- prior(pars = c(pars, other = 1), deriv = TRUE, hessian = "sparse")
  expect_s4_class(sp$hessian, "sparseMatrix")
  expect_equal(as.matrix(sp$hessian)[names(pars), names(pars)], out$hessian)
  expect_equal(unname(as.matrix(sp$hessian)["other", ]), numeric(length(pars) + 1))
  M <- lapply(1:8, function(i) crossprod(matrix(rnorm(6), 3, 2)))
  chol <- dMod2:::.updateNormal(d, M)
  # stationary: d/dchol [sum_g tr(Omega^-1 M_g) + G log|Omega|] = 0
  q <- function(v) {
    L <- dMod2:::.choleskyFactor(d$normal, setNames(v, names(chol)))
    O <- tcrossprod(L)
    sum(vapply(M, function(Mg) sum(diag(solve(O, Mg))), 0)) + 8 * as.numeric(determinant(O)$modulus)
  }
  expect_lt(max(abs(numDeriv::grad(q, unname(chol)))), 1e-6)
  # with a fixed entry the optimum over the free ones: a diagonal Omega of the mean moments
  held <- dMod2:::.updateNormal(d, M, fixed = c(omega_b_a = 0))
  expect_equal(unname(exp(held[c("omega_a_a", "omega_b_b")])), sqrt(diag(Reduce(`+`, M) / 8)), tolerance = 1e-8)
  expect_equal(held[["omega_b_a"]], 0)
})

test_that("FOCEI and quadrature reproduce the closed-form marginal with groups of several conditions", {
  m <- .lin()
  init <- c(a = 0.5, b = 0.2)
  fit <- em(m$obj, randomEffects(m$p), init, verbose = FALSE)
  expect_true(fit$converged)
  pars <- fit$argument
  expect_equal(fit$value, .linMarginal(pars, m), tolerance = 1e-6)
  best <- optim(pars, .linMarginal, mod = m, method = "BFGS", control = list(reltol = 1e-12))
  expect_equal(fit$value, best$value, tolerance = 1e-5)
  expect_equal(unname(fit$argument[c("a", "b")]), unname(best$par[c("a", "b")]), tolerance = 1e-3)
  # the Hessian of the marginal
  h <- numDeriv::jacobian(function(v) numDeriv::grad(function(u) .linMarginal(setNames(u, names(pars)), m), v),
                          unname(pars))
  expect_equal(unname(fit$hessian), h, tolerance = 1e-4)
  qu <- em(m$obj, randomEffects(m$p), init, method = "quadrature", verbose = FALSE)
  expect_equal(qu$value, best$value, tolerance = 1e-5)
  expect_equal(qu$etaSE, fit$etaSE, tolerance = 1e-5)
  expect_equal(nrow(fit$random), 16L)
  expect_equal(fit$omega$type, c("sd", "sd", "cor"))
})

test_that("rows with a recorded sigma keep it beside the error model, replicates included", {
  m <- .lin()
  rows <- do.call(rbind, Map(function(d, cn) data.frame(d[, c("name", "time", "value", "sigma")], condition = cn),
                             m$data, names(m$data)))
  rows$sigma[rows$time == 2] <- NA
  twin <- rows[rows$condition == "s1a" & rows$time == 1, ]
  twin$value <- twin$value + 0.4; twin$sigma <- 0.6
  data <- as.datalist(rbind(rows, twin))
  obj <- normL2(data, m$g * Xt() * m$ps, errmodel = m$err)
  fit <- em(obj, randomEffects(m$ps), c(a = 0.5, b = 0.2, s = log(0.3)), verbose = FALSE)
  closed <- m
  closed$data <- lapply(data, function(d) { d$sigma[is.na(d$sigma)] <- exp(fit$argument[["s"]]); d })
  expect_equal(fit$value, .linMarginal(fit$argument, closed), tolerance = 1e-6)
  ires <- residuals(fit, type = "IWRES")
  expect_equal(sort(unique(round(predict(fit)$sigma, 8))),
               sort(unique(round(c(0.3, 0.6, exp(fit$argument[["s"]])), 8))))
})

test_that("a fixed hyperparameter stays at its value", {
  m <- .lin()
  fit <- em(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), fixed = c(omega_b_a = 0), verbose = FALSE)
  expect_false("omega_b_a" %in% names(fit$argument))
  expect_equal(fit$value, .linMarginal(c(fit$argument, omega_b_a = 0), m), tolerance = 1e-6)
  expect_equal(fit$omega$type, c("sd", "sd"))
  expect_true(all(is.finite(fit$omega$se)))
})

test_that("a fit with every parameter fixed evaluates the marginal there; valueNonmem drops 2 pi", {
  m <- .lin()
  pars <- c(a = 0.9, b = 0.45, omega_a_a = log(0.6), omega_b_a = 0.1, omega_b_b = log(0.3))
  fit <- em(m$obj, randomEffects(m$p), fixed = pars, verbose = FALSE)
  expect_equal(fit$value, .linMarginal(pars, m), tolerance = 1e-6)
  expect_true(fit$converged)
  expect_length(fit$argument, 0L)
  expect_equal(fit$valueNonmem, fit$value - fit$nObs * log(2 * pi))
  expect_error(em(m$obj, randomEffects(m$p), fixed = pars, method = "quadrature", verbose = FALSE),
               "every parameter is fixed")
})

test_that("msem draws starts and seeds in the session, so the cores do not change the result", {
  skip_on_os("windows")
  m <- .lin()
  ctl <- list(saem = list(nBurnin = 10L, nEM = 10L))
  run <- function(cores) {
    set.seed(1)
    msem(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), method = "saem", control = ctl, fits = 2,
         cores = cores)$multistart$value
  }
  expect_equal(run(2), run(1))
})

test_that("an EM fit reports its predictions, residuals, log-likelihood and random effects", {
  m <- .lin()
  fit <- em(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), verbose = FALSE)
  pr <- predict(fit, times = c(0.5, 1.5))
  expect_true(all(c("condition", "group", "time", "name", "observed", "sigma", "sigmaPRED", "IPRED",
                    "PRED", "IRES", "IWRES", "PRES", "PWRES", "source") %in% names(pr)))
  obs <- pr[pr$source == "obs", ]
  expect_equal(unique(obs$group[obs$condition == "s2b"]), "s2")
  expect_true(is.na(unique(obs$group[obs$condition == "ref"])))
  expect_equal(obs$PWRES, obs$PRES / obs$sigmaPRED)
  expect_equal(residuals(fit, "IRES"), obs$IRES)
  # CWRES of a group whiten its stacked rows against the linearised population prediction
  cw <- residuals(fit, "CWRES")
  r <- obs$condition %in% c("s1a", "s1b")
  Z <- cbind(1, m$tab[obs$condition[r], "dose"])
  C <- Z %*% fit$Omega %*% t(Z) + diag(0.3^2, sum(r))
  expect_equal(cw[r], drop(backsolve(chol(C), obs$observed[r] - obs$PRED[r], transpose = TRUE)),
               tolerance = 1e-8)
  expect_equal(as.numeric(logLik(fit)), -fit$value / 2)
  expect_equal(BIC(fit), fit$value + log(51) * 5)
  expect_equal(fit$omega$estimate[1:2], unname(sqrt(diag(fit$Omega))), tolerance = 1e-12)
  expect_equal(fit$random$parameter[1:8], paste0("eta_a_", m$groups))
  expect_output(print(summary(fit)), "Random effects")
})

test_that("profile intervals bracket the estimate where -2 log L rises by the quantile", {
  m <- .lin()
  fit <- em(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), verbose = FALSE)
  ci <- confint(fit, parm = "a", method = "profile")
  at <- function(v) em(m$obj, randomEffects(m$p), fit$argument[-1], fixed = c(a = v), verbose = FALSE)$value
  expect_equal(at(ci$upper) - fit$value, qchisq(0.95, 1), tolerance = 0.02)
  expect_equal(at(ci$lower) - fit$value, qchisq(0.95, 1), tolerance = 0.02)
})

test_that("msem keeps the best start and lists every start", {
  m <- .lin()
  set.seed(1)
  fit <- msem(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), fits = 3, sd = 0.2)
  expect_s3_class(fit, "em")
  expect_equal(nrow(fit$multistart), 3L)
  expect_equal(fit$value, min(fit$multistart$value, na.rm = TRUE))
  expect_true(all(fit$multistart$iterations >= 1L))
  skip_if_not_installed("ggplot2")
  expect_equal(sort(attr(plotValues(fit), "data")$value), sort(fit$multistart$value))
  expect_error(plotValues(em(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), verbose = FALSE)), "msem")
})

test_that("FOCEI reports at most iterlim outer iterations", {
  m <- .lin()
  fit <- em(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), control = list(trustControl = list(iterlim = 2)),
            verbose = FALSE)
  expect_equal(fit$iterations, 2L)
  expect_false(fit$converged)
})

test_that("msem ends the other starts at the budget from the iterations of the first", {
  m <- .lin()
  run <- function(budget) {
    set.seed(1)
    msem(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), fits = 3, sd = 0.2, budget = budget)$multistart
  }
  free <- run(Inf)
  # a cap one below the iterations every other start needs
  cap <- min(free$iterations[-1]) - 1
  expect_gt(cap, 0)
  capped <- run((cap - 0.5) / free$iterations[1])
  expect_equal(capped[1, ], free[1, ])
  expect_equal(capped$iterations[-1], rep(cap, 2))
  expect_false(any(capped$converged[-1]))
  expect_error(run(0), "budget")
})

test_that("the diagnostic plots of an EM fit are ggplots", {
  skip_if_not_installed("ggplot2")
  m <- .lin()
  fit <- em(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), verbose = FALSE)
  expect_s3_class(plotCombined(fit), "ggplot")
  expect_s3_class(plot(fit), "ggplot")
  expect_s3_class(plotResiduals(fit), "ggplot")
  expect_true(inherits(plotRandomEffects(fit), "ggplot") || is.list(plotRandomEffects(fit)))
  expect_true(inherits(plotRandomEffects(fit, type = "modes"), "ggplot") ||
                is.list(plotRandomEffects(fit, type = "modes")))
  # the modes as they are, per group `draws` values from the conditional distribution
  modes <- dMod2:::.randomEffectValues(fit, "modes", 7L)
  draws <- dMod2:::.randomEffectValues(fit, "conditional", 7L)
  etas <- intersect(colnames(fit$etaModes), rownames(fit$Omega))
  expect_equal(modes$value, as.vector(fit$etaModes[, etas]))
  expect_equal(nrow(draws), 7L * nrow(modes))
  expect_true(all(is.finite(draws$value)))
  expect_error(plotTrace(fit), "stage trace")
})

test_that("forward over reverse gives the volume gradient and the mode shift of forward over forward", {
  m <- .lin()
  rows <- do.call(rbind, Map(function(d, cn) data.frame(d[, c("name", "time", "value", "sigma")], condition = cn),
                             m$data, names(m$data)))
  rows$sigma[rows$time == 2] <- NA
  rows$lloq <- ifelse(rows$condition == "s2a", 1.2, -Inf)
  rows$value[rows$condition == "s2a" & rows$time == 0] <- 1
  data <- as.datalist(rows)
  for (obj in list(normL2(data, m$g * Xt() * m$ps, errmodel = m$err, optBLOQ = "M3"),
                   normT(data, m$g * Xt() * m$ps, errmodel = m$err, df = 4))) {
    design <- .emDesign(obj, randomEffects(m$ps))
    start <- .emStart(design, c(a = 0.5, b = 0.2, s = log(0.3)))
    E <- design$etaNames
    full <- c(start, setNames(seq(-0.4, 0.4, length.out = length(E)), as.vector(E)))
    Hinv <- lapply(design$groups, function(g) matrix(c(0.5, 0.1, 0.1, 0.3), 2, dimnames = list(E[g, ], E[g, ])))
    fwd <- .foceiVolume(design, full, NULL, names(start), Hinv)
    rev <- .volumeReverse(design, full, NULL, names(start), Hinv)
    expect_equal(rev$gradient, fwd$gradient, tolerance = 1e-10)
    expect_equal(rev$etaShift, fwd$etaShift, tolerance = 1e-10)
    # a sweep with two seed blocks answers each as a sweep of its own
    units <- .volumeUnits(design, full, NULL)
    fw <- .sweepForward(design, full, NULL, units)
    seed <- function(f) lapply(fw$rows, function(r) list(pred = f(r$pd, r$p), sigma = f(r$sd, r$s)))
    one <- seed(function(d, v) cbind(v, d))
    two <- seed(function(d, v) cbind(v^2, d * v))
    both <- .sweepBackward(fw, Map(function(a, b) list(pred = cbind(a$pred, b$pred), sigma = cbind(a$sigma, b$sigma)),
                                   one, two), blocks = 2L)
    w <- fw$K + 1L
    for (b in 1:2) {
      single <- .sweepBackward(fw, list(one, two)[[b]])
      cols <- (b - 1L) * w + seq_len(w)
      expect_equal(Map(function(ct, s) ct[rownames(s), cols, drop = FALSE], both, single), single, tolerance = 1e-12)
      expect_true(all(unlist(Map(function(ct, s) ct[setdiff(rownames(ct), rownames(s)), cols], both, single)) == 0))
    }
  }
  # the model is linear in the etas: the volume curvature from dH/dtheta is the whole one, and the
  # Hessian forward over reverse reports is the exact one
  fit <- em(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), verbose = FALSE)
  expect_identical(fit$derivMode, "forward-reverse")
  design <- .emDesign(m$obj, randomEffects(m$p))
  H <- .reportedReverse(design, design$objfn + .normalPrior(design), c(fit$argument, .etaVector(fit)), NULL,
                        names(fit$argument), lapply(fit$etaCov, `/`, 2))
  pars <- fit$argument
  exact <- numDeriv::jacobian(function(v) numDeriv::grad(function(u) .linMarginal(setNames(u, names(pars)), m), v),
                              unname(pars))
  expect_equal(unname(H), exact, tolerance = 1e-4)
  ff <- em(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), verbose = FALSE,
           control = list(derivMode = "forward-forward"))
  expect_equal(fit$value, ff$value, tolerance = 1e-8)
  expect_equal(fit$argument, ff$argument, tolerance = 1e-4)
  expect_error(em(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), control = list(derivMode = "reverse")),
               "derivMode")
})

test_that("em() and msem() say once that the models give second derivatives forward over forward only", {
  m <- .lin()
  messages <- function(expr) {
    out <- character()
    value <- withCallingHandlers(expr, message = function(w) {
      out <<- c(out, conditionMessage(w))
      invokeRestart("muffleMessage")
    })
    list(value = value, messages = out)
  }
  objF <- normL2(m$data, m$g * Xt() * m$pf)
  ff <- messages(em(objF, randomEffects(m$pf), c(a = 0.5, b = 0.2), verbose = FALSE))
  expect_identical(ff$value$derivMode, "forward-forward")
  expect_length(grep("forward-reverse", ff$messages), 1L)
  ms <- messages(msem(objF, randomEffects(m$pf), c(a = 0.5, b = 0.2), fits = 3, cores = 1))
  expect_length(grep("forward-reverse", ms$messages), 1L)
  expect_null(ms$value$derivModeNote)
  fr <- messages(em(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), verbose = FALSE))
  expect_identical(fr$value$derivMode, "forward-reverse")
  expect_length(grep("derivMode", fr$messages), 0L)
  expect_equal(ff$value$value, fr$value$value, tolerance = 1e-8)
})

test_that("the batched density of the groups equals the kernel value, censored and Student-t rows included", {
  m <- .lin()
  rows <- do.call(rbind, Map(function(d, cn) data.frame(d[, c("name", "time", "value", "sigma")], condition = cn),
                             m$data, names(m$data)))
  rows$sigma[rows$time == 2] <- NA
  rows$lloq <- ifelse(rows$condition == "s2a", 1.2, -Inf)
  rows$value[rows$condition == "s2a" & rows$time == 0] <- 1
  data <- as.datalist(rows)
  for (obj in list(normL2(data, m$g * Xt() * m$ps, errmodel = m$err, optBLOQ = "M3"),
                   normT(data, m$g * Xt() * m$ps, errmodel = m$err, df = 4))) {
    design <- .emDesign(obj, randomEffects(m$ps))
    start <- .emStart(design, c(a = 0.5, b = 0.2, s = log(0.3)))
    meta <- .foceiMeta(design, start, NULL)
    E <- design$etaNames
    full <- setNames(numeric(length(meta$pars_full_names)), meta$pars_full_names)
    full[names(start)] <- start
    eta <- matrix(seq(-0.4, 0.4, length.out = length(E)), nrow(E), dimnames = dimnames(E))
    L <- .choleskyFactor(design$normal, start)
    Oinv <- chol2inv(t(L)); logDet <- 2 * sum(log(diag(L)))
    kernel <- vapply(seq_len(nrow(E)), function(i) {
      p <- full; p[E[i, ]] <- eta[i, ]
      focei_eval_one_subject(model_cb = design$prdfn, err_cb = .errorModel(design), pars_full = p, fixed = NULL,
                             meta_i = meta$fast_meta[[i]], eta_block = eta[i, ], Omega_inv = Oinv,
                             Omega_log_det = logDet)$value
    }, 0)
    expect_equal(unname(.groupDensity(design, meta, NULL)(full, eta, Oinv, logDet)), kernel, tolerance = 1e-12)
  }
})

test_that("laplaceEM and SAEM reach the closed-form marginal", {
  m <- .lin()
  best <- optim(c(a = 0.5, b = 0.2, omega_a_a = log(0.5), omega_b_a = 0, omega_b_b = log(0.5)),
                .linMarginal, mod = m, method = "BFGS", control = list(reltol = 1e-12))
  le <- em(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), method = "laplaceEM", verbose = FALSE)
  expect_equal(le$value, best$value, tolerance = 1e-4)
  # a and b enter the condition without a group: the centred trust step
  set.seed(2)
  sa <- em(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), method = "saem", verbose = FALSE)
  expect_true(sa$converged)
  expect_lt(sa$value - best$value, 0.15)
  expect_setequal(unique(sa$stageTrace$phase), c("burnin", "averaging"))
  expect_s3_class(plotTrace(le), "ggplot")
  expect_s3_class(plotTrace(sa), "ggplot")
  # without it a and b come by GLS; FOCEI is exact on this model
  grouped <- normL2(m$data[setdiff(names(m$data), "ref")], m$g * Xt() * m$p)
  fo <- em(grouped, randomEffects(m$p), c(a = 0.5, b = 0.2), verbose = FALSE)
  set.seed(2)
  sg <- em(grouped, randomEffects(m$p), c(a = 0.5, b = 0.2), method = "saem", verbose = FALSE)
  expect_lt(sg$value - fo$value, 0.15)
  # the random numbers come from the session stream
  short <- function() {
    set.seed(5)
    em(grouped, randomEffects(m$p), c(a = 0.5, b = 0.2), method = "saem", verbose = FALSE,
       control = list(saem = list(nBurnin = 5, nEM = 5)))$argument
  }
  expect_identical(short(), short())
})

test_that("the SAEM step on the averaged quadratic models zeroes the mean gradient", {
  # -2 log L of a residual scale s over n residuals with sum of squares r; the draws alternate
  # between two r of different curvature, the score of their mean vanishes at exp(2 s) = mean(r) / n
  n <- 10; r <- c(10, 40)
  draw <- function(s, r) list(gradient = c(s = 2 * n - 2 * r * exp(-2 * s)),
                              hessian = matrix(4 * r * exp(-2 * s), 1, 1, dimnames = list("s", "s")))
  theta <- c(s = 0); model <- NULL
  for (k in 1:4000) {
    model <- .averageQuadratic(model, draw(theta[["s"]], r[k %% 2 + 1]), theta, k^-0.7)
    theta <- .quadraticStep(model, theta)
  }
  expect_equal(theta[["s"]], log(mean(r) / n) / 2, tolerance = 1e-2)
})

test_that("anova tests nested fits and simulate draws data that refit to the estimates", {
  m <- .lin()
  full <- em(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), verbose = FALSE)
  diag <- em(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), fixed = c(omega_b_a = 0), verbose = FALSE)
  tab <- anova(full, diag)
  expect_equal(tab$npar, c(4, 5))
  expect_equal(tab$chisq[2], diag$value - full$value)
  expect_equal(tab$p[2], pchisq(diag$value - full$value, 1, lower.tail = FALSE))
  sims <- simulate(full, nsim = 2, seed = 1)
  expect_length(sims, 2L)
  expect_s3_class(sims[[1]], "datalist")
  expect_equal(names(sims[[1]]), names(m$data))
  expect_false(isTRUE(all.equal(sims[[1]][["s1a"]]$value, sims[[2]][["s1a"]]$value)))
  # many groups simulated from the fit refit near it
  refit <- em(normL2(sims[[1]], m$g * Xt() * m$p), randomEffects(m$p), full$argument[c("a", "b")],
              verbose = FALSE)
  se <- sqrt(diag(vcov(full)))[c("a", "b")]
  expect_true(all(abs(refit$argument[c("a", "b")] - full$argument[c("a", "b")]) < 4 * se))
})


# crossed groupings, several replicates per condition
.crossed <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    withr::local_preserve_seed()
    dir <- file.path(tempdir(), "nlmefit_crossed")
    dir.create(dir, showWarnings = FALSE)
    withr::local_dir(dir)
    set.seed(7)
    pat <- paste0("p", 1:8); cen <- paste0("c", 1:3)
    tab <- do.call(rbind, lapply(seq_along(pat), function(i)
      data.frame(ID = pat[i], centre = cen[c(i %% 3 + 1, (i + 1) %% 3 + 1)])))
    rownames(tab) <- paste(tab$ID, tab$centre, sep = "_")
    tr <- branch(define(NULL, "x <- x", x = "a"), table = tab)
    tr <- insert(tr, a ~ 1 + (1 | ID) + (1 | centre))
    dm <- c("forward", "reverse", "forward-reverse")
    g <- Y(c(y = "a"), f = NULL, parameters = "a", compile = FALSE, modelname = "crossed_g", deriv2 = TRUE,
           derivMode = dm)
    p <- P(tr, method = "explicit", compile = FALSE, modelname = "crossed_p", deriv2 = TRUE, derivMode = dm)
    compile(g, p, output = "nlmefit_crossed", cores = 2)
    ep <- setNames(rnorm(8, 0, 0.6), pat); ec <- setNames(rnorm(3, 0, 0.4), cen)
    rows <- lapply(rownames(tab), function(cn) data.frame(name = "y", time = 1:3,
      value = 1 + ep[[tab[cn, "ID"]]] + ec[[tab[cn, "centre"]]] + rnorm(3, 0, 0.3), sigma = 0.3, condition = cn))
    data <- as.datalist(do.call(rbind, rows))
    cache <<- list(p = p, data = data, tab = tab, obj = normL2(data, g * Xt() * p))
    cache
  }
})

# -2 log of the closed-form marginal of all data: sigma^2 I + omega_p^2 Zp Zp' + omega_c^2 Zc Zc'.
.crossedMarginal <- function(pars, mod) {
  d <- do.call(rbind, lapply(names(mod$data), function(cn) cbind(mod$data[[cn]], cond = cn)))
  Zp <- outer(mod$tab[d$cond, "ID"], unique(mod$tab$ID), `==`) * 1
  Zc <- outer(mod$tab[d$cond, "centre"], unique(mod$tab$centre), `==`) * 1
  C <- diag(0.3^2, nrow(d)) + exp(2 * pars[["omega_a_a"]]) * tcrossprod(Zp) +
    exp(2 * pars[["omega_a_centre_a_centre"]]) * tcrossprod(Zc)
  r <- d$value - pars[["a"]]
  drop(t(r) %*% solve(C, r)) + as.numeric(determinant(C)$modulus) + nrow(d) * log(2 * pi)
}

test_that("crossed groupings: FOCEI over the connected component hits the closed-form marginal", {
  m <- .crossed()
  re <- randomEffects(m$p)
  expect_equal(re$terms$stem, c("eta_a", "eta_a_centre"))
  d <- dMod2:::.emDesign(m$obj, re)
  expect_equal(names(d$groupings), c("ID", "centre"))
  expect_length(d$components, 1L)
  fit <- em(m$obj, re, c(a = 0.5), verbose = FALSE)
  expect_equal(fit$value, .crossedMarginal(fit$argument, m), tolerance = 1e-6)
  best <- optim(fit$argument, .crossedMarginal, mod = m, method = "BFGS", control = list(reltol = 1e-12))
  expect_equal(fit$value, best$value, tolerance = 1e-6)
  h <- numDeriv::jacobian(function(v) numDeriv::grad(function(u) .crossedMarginal(setNames(u, names(fit$argument)), m), v),
                          unname(fit$argument))
  expect_equal(unname(fit$hessian), h, tolerance = 1e-4)
  expect_equal(unique(fit$omega$grouping), c("ID", "centre"))
  expect_equal(nrow(fit$random), 11L)
  expect_true(all(is.finite(residuals(fit, "CWRES"))))
  # a fixed structural parameter enters the joint modes
  at <- em(m$obj, re, NULL, fixed = c(a = 1.1), verbose = FALSE)
  expect_equal(at$value, .crossedMarginal(c(at$argument, a = 1.1), m), tolerance = 1e-6)
})

test_that("one grouping through the component path equals the kernel path", {
  m <- .lin()
  fit <- em(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), verbose = FALSE)
  d <- dMod2:::.emDesign(m$obj, randomEffects(m$p))
  d$groupings <- list(ID = d[c("grouping", "groups", "conditions", "etaNames", "blocks", "normal")])
  d$components <- dMod2:::.components(d$groupings)
  joint <- d$objfn + dMod2:::.normalPrior(d$groupings[[1]])
  at <- dMod2:::.fitFoceiJoint(d, joint, fit$argument, control = list(trustControl = list(iterlim = 0)))
  expect_equal(at$value, fit$value, tolerance = 1e-8)
  expect_equal(unname(at$hessian), unname(fit$hessian), tolerance = 1e-6)
})

test_that("sparse FOCEI equals the dense fit, for one grouping and for crossed groupings", {
  m <- .lin()
  re <- randomEffects(m$p)
  init <- c(a = 0.5, b = 0.2)
  dense <- em(m$obj, re, init, verbose = FALSE, sparse = FALSE)
  sp <- em(m$obj, re, init, verbose = FALSE, sparse = TRUE)
  expect_s4_class(sp$hessian, "sparseMatrix")
  expect_true(sp$converged)
  # the dense outer model has the BFGS correction, the sparse one not: equal to the tolerance of the
  # trust region
  expect_equal(sp$value, dense$value, tolerance = 1e-8)
  expect_equal(sp$argument, dense$argument, tolerance = 1e-4)
  expect_equal(as.matrix(sp$hessian)[names(dense$argument), names(dense$argument)], dense$hessian,
               tolerance = 1e-4)
  expect_equal(sp$etaModes, dense$etaModes, tolerance = 1e-4)
  expect_equal(summary(sp)$population, summary(dense)$population, tolerance = 1e-4)
  expect_equal(sp$omega$se, dense$omega$se, tolerance = 1e-4)
  # auto takes the sparse path above the thresholds
  auto <- withr::with_options(list(dMod2.sparseAuto = list(parameters = 0, fill = 1)),
                              em(m$obj, re, init, verbose = FALSE))
  expect_s4_class(auto$hessian, "sparseMatrix")
  expect_true(is.matrix(em(m$obj, re, init, verbose = FALSE, control = list(trustControl = list(iterlim = 0)))$hessian))

  m <- .crossed()
  re <- randomEffects(m$p)
  dense <- em(m$obj, re, c(a = 0.5), verbose = FALSE, sparse = FALSE)
  sp <- em(m$obj, re, c(a = 0.5), verbose = FALSE, sparse = TRUE)
  expect_s4_class(sp$hessian, "sparseMatrix")
  expect_identical(dense$derivMode, "forward-reverse")
  expect_equal(sp$value, dense$value, tolerance = 1e-8)
  expect_equal(sp$argument, dense$argument, tolerance = 1e-4)
  expect_equal(as.matrix(sp$hessian)[names(dense$argument), names(dense$argument)], dense$hessian,
               tolerance = 1e-4)
  ff <- em(m$obj, re, c(a = 0.5), verbose = FALSE, control = list(derivMode = "forward-forward"))
  expect_equal(dense$value, ff$value, tolerance = 1e-8)
  expect_equal(dense$hessian, ff$hessian, tolerance = 1e-4)
})
