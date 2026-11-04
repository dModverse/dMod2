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
    # uncorrelated random effects
    pu <- P(define(insert(branch(define(NULL, "x <- x", x = c("a", "b")), table = tab), a + b ~ 1 + (1 | ID),
                          conditionMatch = "^s"), "dose <- d", d = dose),
            method = "explicit", compile = FALSE, modelname = "nlmefit_lin_pu", deriv2 = TRUE, derivMode = dm)
    compile(g, p, ps, err, pf, pu, output = "nlmefit_lin", cores = 2)
    # random effects at quantiles and noise of mean 0 in every condition: the means of the
    # conditions, and with them the maximum likelihood of Omega, are the same for any draw and in
    # the interior; drawn group means leave it on the boundary in about one draw of four
    set.seed(3)
    q <- qnorm(ppoints(8))
    eta <- cbind(0.7 * q, 0.6 * q[c(5, 2, 7, 4, 1, 8, 3, 6)])
    rows <- lapply(rownames(tab), function(cn) {
      i <- match(tab[cn, "ID"], groups)
      e <- if (is.na(i)) c(0, 0) else eta[i, ]
      mu <- (1 + e[1]) + (0.5 + e[2]) * tab[cn, "dose"]
      noise <- rnorm(3, 0, 0.3)
      data.frame(name = "y", time = c(0, 1, 2), value = mu + noise - mean(noise), sigma = 0.3,
                 condition = cn)
    })
    data <- as.datalist(do.call(rbind, rows))
    cache <<- list(g = g, p = p, ps = ps, pf = pf, pu = pu, err = err, data = data, tab = tab, groups = groups,
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
  # with a fixed entry the optimum over the free ones, by trust() to its tolerance: a diagonal
  # Omega of the mean moments
  held <- dMod2:::.updateNormal(d, M, fixed = c(omega_b_a = 0))
  expect_equal(unname(exp(held[c("omega_a_a", "omega_b_b")])), sqrt(diag(Reduce(`+`, M) / 8)), tolerance = 1e-6)
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

test_that("a variance the data leave at 0 is set to 0, and the fit without its random effect is returned", {
  m <- .lin()
  withr::local_seed(5)
  # the intercepts vary between the groups, the slopes do not for any draw: the noise of a group is
  # the same at both doses, and the condition without a group is the mean of the first doses
  ea <- setNames(rnorm(8, 0, 0.7), m$groups)
  noise <- lapply(ea, function(e) rnorm(3, 0, 0.3))
  grouped <- setdiff(rownames(m$tab), "ref")
  values <- lapply(setNames(nm = grouped), function(cn) {
    g <- m$tab[cn, "ID"]
    1 + ea[[g]] + 0.5 * m$tab[cn, "dose"] + noise[[g]]
  })
  first <- grouped[m$tab[grouped, "dose"] == 1]
  values$ref <- Reduce(`+`, values[first]) / length(first)
  d0 <- as.datalist(do.call(rbind, lapply(names(values), function(cn)
    data.frame(name = "y", time = c(0, 1, 2), value = values[[cn]], sigma = 0.3, condition = cn))))
  obj0 <- normL2(d0, m$g * Xt() * m$pu)
  random <- randomEffects(m$pu)
  off <- em(obj0, random, c(a = 0.5, b = 0.2), control = list(boundary = FALSE), verbose = FALSE)
  expect_lt(off$argument[["omega_b_b"]], log(0.1))
  fit <- em(obj0, random, c(a = 0.5, b = 0.2), verbose = FALSE)
  expect_identical(fit$boundary, "omega_b_b")
  expect_true(fit$converged)
  expect_false("omega_b_b" %in% names(fit$argument))
  expect_identical(unname(fit$fixed[paste0("eta_b_", m$groups)]), rep(0, 8))
  expect_lte(fit$value, off$value + 0.1)
  expect_equal(fit$argument[c("a", "b", "omega_a_a")], off$argument[c("a", "b", "omega_a_a")], tolerance = 1e-2)
  # a fit as the start goes on with the random effects it kept
  again <- em(obj0, random, fit, control = list(trustControl = list(iterlim = 0)), verbose = FALSE)
  expect_false("omega_b_b" %in% names(again$argument))
  expect_equal(again$value, fit$value, tolerance = 1e-8)
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

test_that("msem draws a new start where a fit fails, and keeps the first attempts", {
  m <- .lin()
  fitEm <- dMod2:::.emFit
  local_mocked_bindings(.emFit = function(design, start, ...) {
    if (start[["b"]] > 0.2) stop("no stable root at the start")
    fitEm(design, start, ...)
  })
  run <- function(retry) {
    set.seed(2)
    msem(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), fits = 8, sd = 0.2, retry = retry)$multistart
  }
  again <- run(TRUE)
  once <- run(FALSE)
  expect_true(all(is.na(again$error)))
  expect_true(any(again$tries > 1L))
  expect_equal(!is.na(once$error), again$tries > 1L)
  expect_equal(once$value[again$tries == 1L], again$value[again$tries == 1L])
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
  # quadrature caps its FOCEI stage by the iterations of that stage; SAEM is not capped
  first <- list(iterations = 3L, trustIterations = 12L)
  expect_equal(dMod2:::.msemCap(first, "quadrature", 2), 24L)
  expect_equal(dMod2:::.msemCap(first["iterations"], "focei", 2), 6L)
  expect_null(dMod2:::.msemCap(first, "saem", 2))
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

test_that("objectives of single conditions evaluated together equal their evaluations one by one", {
  m <- .lin()
  dataE <- as.datalist(do.call(rbind, lapply(names(m$data), function(cn)
    transform(m$data[[cn]], sigma = NA_real_, condition = cn))))
  obj <- normL2(dataE, m$g * Xt() * m$ps, errmodel = m$err)
  objs <- dMod2:::.conditionObjectives(dMod2:::.emDesign(obj, randomEffects(m$ps)))
  pars <- attr(obj, "parameters")
  set.seed(1)
  parsList <- lapply(seq_along(objs), function(j) setNames(rnorm(length(pars), 0, 0.3), pars))
  for (deriv in c(TRUE, FALSE)) {
    together <- dMod2:::.objEvalMany(objs, parsList, deriv = deriv)
    single <- lapply(seq_along(objs), function(j) objs[[j]](parsList[[j]], deriv = deriv))
    for (j in seq_along(objs)) {
      expect_equal(together[[j]]$value, single[[j]]$value, tolerance = 1e-12)
      expect_equal(together[[j]]$gradient, single[[j]]$gradient, tolerance = 1e-12)
      expect_equal(together[[j]]$hessian, single[[j]]$hessian, tolerance = 1e-12)
    }
  }
})

test_that("laplaceEM and SAEM reach the closed-form marginal", {
  m <- .lin()
  best <- optim(c(a = 0.5, b = 0.2, omega_a_a = log(0.5), omega_b_a = 0, omega_b_b = log(0.5)),
                .linMarginal, mod = m, method = "BFGS", control = list(reltol = 1e-12))
  le <- em(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), method = "laplaceEM", verbose = FALSE)
  expect_equal(le$value, best$value, tolerance = 1e-4)
  # SAEM stops where its Monte Carlo gap falls below 0.05 or at the end of its 500 iterations,
  # -2 log L at its estimate above the optimum by a few times that gap
  gap <- function(fit) max(utils::tail(fit$stageTrace$mcGap, 1), 0.05)
  # a and b enter the condition without a group: the centred trust step
  set.seed(2)
  sa <- em(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), method = "saem", verbose = FALSE)
  expect_true(sa$converged || nrow(sa$stageTrace) == 500L)
  expect_lt(sa$value - best$value, 5 * gap(sa))
  expect_setequal(unique(sa$stageTrace$phase), c("burnin", "averaging"))
  expect_s3_class(plotTrace(le), "ggplot")
  expect_s3_class(plotTrace(sa), "ggplot")
  # without it a and b come by GLS; FOCEI is exact on this model
  grouped <- normL2(m$data[setdiff(names(m$data), "ref")], m$g * Xt() * m$p)
  fo <- em(grouped, randomEffects(m$p), c(a = 0.5, b = 0.2), verbose = FALSE)
  set.seed(2)
  sg <- em(grouped, randomEffects(m$p), c(a = 0.5, b = 0.2), method = "saem", verbose = FALSE)
  expect_lt(sg$value - fo$value, 5 * gap(sg))
  # the random numbers come from the session stream
  short <- function() {
    set.seed(5)
    em(grouped, randomEffects(m$p), c(a = 0.5, b = 0.2), method = "saem", verbose = FALSE,
       control = list(saem = list(nBurnin = 5, nEM = 5)))$argument
  }
  expect_identical(short(), short())
})

test_that("SAEM goes on where the step of the trust parameters fails at the draws", {
  m <- .lin()
  start <- c(a = 0.5, b = 0.2)
  inner <- dMod2:::.innerTrust
  failedAt <- NULL
  after <- 0L
  # the first trust step that starts away from the start of the fit fails
  local_mocked_bindings(.innerTrust = function(args) {
    p <- args$parinit
    if (!is.null(failedAt)) after <<- after + 1L
    if (is.null(failedAt) && setequal(names(p), names(start)) && !isTRUE(all.equal(p[names(start)], start))) {
      failedAt <<- p
      stop("parinit not feasible: objfun failed: 1 of 2 conditions did not complete")
    }
    inner(args)
  }, .package = "dMod2")
  set.seed(2)
  sa <- em(m$obj, randomEffects(m$p), start, method = "saem", verbose = FALSE,
           control = list(saem = list(nBurnin = 5, nEM = 5)))
  expect_false(is.null(failedAt))
  expect_gt(after, 0L)
  expect_true(is.finite(sa$value))
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
    # the same groupings, nonlinear in the random effects with a proportional error
    trn <- branch(define(NULL, "x <- x", x = c("a", "s")), table = tab)
    trn <- insert(trn, "s <- exp(s)")
    trn <- insert(trn, a ~ 1 + (1 | ID) + (1 | centre))
    gn <- Y(c(y = "exp(a)"), f = NULL, parameters = "a", compile = FALSE, modelname = "crossed_gn", deriv2 = TRUE,
            derivMode = dm)
    en <- Y(c(y = "s * y"), gn, attachInput = FALSE, compile = FALSE, modelname = "crossed_en", deriv2 = TRUE,
            derivMode = dm)
    pn <- P(trn, method = "explicit", compile = FALSE, modelname = "crossed_pn", deriv2 = TRUE, derivMode = dm)
    compile(g, p, gn, en, pn, output = "nlmefit_crossed", cores = 2)
    # the effects of the three centres spread by construction: three drawn ones leave the maximum
    # likelihood of their variance at 0 in about one draw of ten
    ep <- setNames(rnorm(8, 0, 0.6), pat); ec <- setNames(c(-0.4, 0, 0.4), cen)
    rows <- lapply(rownames(tab), function(cn) data.frame(name = "y", time = 1:3,
      value = 1 + ep[[tab[cn, "ID"]]] + ec[[tab[cn, "centre"]]] + rnorm(3, 0, 0.3), sigma = 0.3, condition = cn))
    data <- as.datalist(do.call(rbind, rows))
    rows <- lapply(rownames(tab), function(cn) data.frame(name = "y", time = 1:3,
      value = exp(1 + ep[[tab[cn, "ID"]]] + ec[[tab[cn, "centre"]]]) * (1 + rnorm(3, 0, 0.2)), sigma = NA_real_,
      condition = cn))
    datan <- as.datalist(do.call(rbind, rows))
    cache <<- list(p = p, data = data, tab = tab, obj = normL2(data, g * Xt() * p),
                   pn = pn, objn = normL2(datan, gn * Xt() * pn, errmodel = en))
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
  h <- numDeriv::hessian(function(u) .crossedMarginal(setNames(u, names(fit$argument)), m), unname(fit$argument))
  expect_equal(unname(fit$hessian), h, tolerance = 1e-4)
  expect_equal(unique(fit$omega$grouping), c("ID", "centre"))
  expect_equal(nrow(fit$random), 11L)
  expect_true(all(is.finite(residuals(fit, "CWRES"))))
  # a fixed structural parameter enters the joint modes
  at <- em(m$obj, re, NULL, fixed = c(a = 1.1), verbose = FALSE)
  expect_equal(at$value, .crossedMarginal(c(at$argument, a = 1.1), m), tolerance = 1e-6)
})

test_that("the lock-step search of the modes stops on its Newton step, not on a change of the value", {
  m <- .crossed()
  re <- randomEffects(m$pn)
  d <- dMod2:::.emDesign(m$objn, re)
  # a loose value tolerance ends only steps on the boundary of the region
  fit <- em(m$objn, re, c(a = 0.5, s = log(0.3)), verbose = FALSE,
            control = list(innerControl = list(mtol = 1e-4), trustControl = list(iterlim = 0)))
  subs <- dMod2:::.subDesigns(d)
  etas <- unlist(lapply(seq_along(subs), function(k) {
    En <- subs[[k]]$etaNames
    setNames(as.vector(fit$etaModes[[k]][rownames(En), , drop = FALSE]), as.vector(En))
  }))
  full <- c(fit$argument, etas)
  E <- names(etas)
  data <- m$objn(full, deriv = TRUE)
  prior <- Reduce(`+`, lapply(subs, function(sd) dMod2:::.normalPrior(sd)(pars = full, deriv = TRUE)))
  step <- solve(data$hessian[E, E] + prior$hessian[E, E], data$gradient[E] + prior$gradient[E])
  expect_lt(max(abs(step)), 1e-8)
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
  # both to the optimum: the outer models learn the volume curvature by SR1 along different steps,
  # so at the default tolerances they stop at different points within them
  tight <- list(trustControl = list(mtol = 1e-10, ftol = 1e-12))
  dense <- em(m$obj, re, init, verbose = FALSE, sparse = FALSE, control = tight)
  sp <- em(m$obj, re, init, verbose = FALSE, sparse = TRUE, control = tight)
  expect_s4_class(sp$hessian, "sparseMatrix")
  expect_true(sp$converged)
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
  dense <- em(m$obj, re, c(a = 0.5), verbose = FALSE, sparse = FALSE, control = tight)
  sp <- em(m$obj, re, c(a = 0.5), verbose = FALSE, sparse = TRUE, control = tight)
  expect_s4_class(sp$hessian, "sparseMatrix")
  expect_identical(dense$derivMode, "forward-reverse")
  expect_equal(sp$value, dense$value, tolerance = 1e-8)
  expect_equal(sp$argument, dense$argument, tolerance = 1e-4)
  expect_equal(as.matrix(sp$hessian)[names(dense$argument), names(dense$argument)], dense$hessian,
               tolerance = 1e-4)
  ff <- em(m$obj, re, c(a = 0.5), verbose = FALSE, control = c(tight, list(derivMode = "forward-forward")))
  expect_equal(dense$value, ff$value, tolerance = 1e-8)
  expect_equal(dense$hessian, ff$hessian, tolerance = 1e-4)
})

test_that("FOCEI names second derivatives that are not finite at the start", {
  m <- .lin()
  many <- dMod2:::.predictMany
  local_mocked_bindings(.predictMany = function(...) {
    out <- many(...)
    if (isTRUE(list(...)$deriv2)) for (i in seq_along(out)) attr(out[[i]], "deriv2")[] <- NaN
    out
  }, .package = "dMod2")
  expect_error(em(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), verbose = FALSE,
                  control = list(derivMode = "forward-forward")),
               "second derivatives of the model are not finite")
})

test_that("NPDE of the linear model are its residuals decorrelated by the exact marginal, a censored reading below its fraction", {
  m <- .lin()
  fit <- em(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), verbose = FALSE)
  npde <- residuals(fit, "NPDE", nsim = 4000, seed = 1)
  pr <- predict(fit)
  pr <- pr[pr$source == "obs", ]
  expect_length(npde, nrow(pr))
  expect_equal(attr(npde, "tests")$test, c("mean", "variance", "normality", "global"))
  expect_true(all(attr(npde, "tests")$p >= 0 & attr(npde, "tests")$p <= 1))
  mu <- function(cn) fit$argument[["a"]] + fit$argument[["b"]] * m$tab[cn, "dose"]
  for (g in c("s3", "s6")) {
    r <- which(pr$condition %in% paste0(g, c("a", "b")))
    Z <- cbind(1, m$tab[pr$condition[r], "dose"])
    C <- Z %*% fit$Omega %*% t(Z) + diag(0.3^2, length(r))
    expect_lt(max(abs(npde[r] - forwardsolve(t(chol(C)), pr$observed[r] - mu(pr$condition[r])))), 0.15)
  }
  # rows outside every group are independent; a censored reading lies below the quantile of its
  # censored fraction
  d <- m$data
  d[["ref"]]$lloq <- c(d[["ref"]]$value[1] + 0.1, -Inf, -Inf)
  cf <- em(normL2(d, m$g * Xt() * m$p), randomEffects(m$p), fit$argument[c("a", "b")], verbose = FALSE)
  nc <- residuals(cf, "NPDE", nsim = 4000, seed = 2)
  pc <- predict(cf)
  pc <- pc[pc$source == "obs", ]
  r <- which(pc$condition == "ref")[1]
  limit <- pnorm((d[["ref"]]$lloq[1] - (cf$argument[["a"]] + cf$argument[["b"]] * m$tab["ref", "dose"])) / 0.3)
  expect_lt(nc[r], qnorm(limit) + 0.1)
  expect_true(all(is.finite(nc)))
})

test_that("NPDE of data simulated from the fit itself, with censored readings, are standard normal", {
  m <- .lin()
  d <- m$data
  for (cn in names(d)) d[[cn]]$lloq <- 1.4
  cf <- em(normL2(d, m$g * Xt() * m$p), randomEffects(m$p), c(a = 0.5, b = 0.2), verbose = FALSE)
  sims <- simulate(cf, nsim = 20, seed = 3)
  expect_gt(mean(unlist(lapply(sims, function(s) lapply(s, function(x) x$value <= x$lloq)))), 0.2)
  z <- unlist(lapply(sims, function(s) {
    f <- cf
    for (cn in names(f$design$data)) f$design$data[[cn]]$value <- s[[cn]]$value
    residuals(f, "NPDE", nsim = 500, seed = 4)
  }))
  # between draws the mean of the 1020 values scatters by some 0.03, their variance from 0.95 to 1.27
  expect_lt(abs(mean(z)), 0.15)
  expect_lt(abs(var(z) - 1), 0.4)
})

test_that("the visual predictive check, plain, prediction corrected, stratified and with censored readings", {
  skip_if_not_installed("ggplot2")
  m <- .lin()
  fit <- em(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), verbose = FALSE)
  p <- plotVPC(fit, nsim = 200, bins = 3, seed = 1)
  expect_s3_class(p, "ggplot")
  v <- attr(p, "vpc")
  expect_equal(nrow(v), 3L * 3L)
  expect_setequal(v$percentile, c("5%", "50%", "95%"))
  expect_true(all(v$lo <= v$mid & v$mid <= v$hi))
  expect_gte(mean(v$observed >= v$lo & v$observed <= v$hi), 2 / 3)
  # population predictions that differ between conditions: the correction moves the readings
  pcv <- attr(plotVPC(fit, nsim = 200, bins = 3, pc = TRUE, seed = 1), "vpc")
  expect_false(isTRUE(all.equal(pcv$observed, v$observed)))
  # readings whose population prediction is 0 leave the prediction-corrected check
  zf <- fit
  zf$argument[["a"]] <- -zf$argument[["b"]]
  zv <- attr(plotVPC(zf, nsim = 50, bins = 3, pc = TRUE, seed = 1), "vpc")
  expect_true(all(is.finite(as.matrix(zv[c("observed", "lo", "mid", "hi")]))))
  fit$design$data <- structure(fit$design$data, condition.grid = data.frame(
    condition = names(m$data), dose = m$tab[names(m$data), "dose"], row.names = names(m$data)))
  sv <- attr(plotVPC(fit, nsim = 100, bins = 3, strata = "dose", seed = 1), "vpc")
  expect_setequal(unique(sv$stratum), c("1", "2"))
  expect_error(plotVPC(fit, nsim = 10, strata = "weight"), "condition grid")
  d <- m$data
  for (cn in names(d)) d[[cn]]$lloq <- 1
  cf <- em(normL2(d, m$g * Xt() * m$p), randomEffects(m$p), fit$argument[c("a", "b")], verbose = FALSE)
  cv <- attr(plotVPC(cf, nsim = 100, bins = 3, seed = 1), "vpc")
  obs <- predict(cf)
  obs <- obs[obs$source == "obs", ]
  expect_equal(sum(cv$blq[cv$percentile == "50%"] * table(findInterval(obs$time, sort(unique(
    quantile(obs$time, seq(0, 1, length.out = 4), names = FALSE))), rightmost.closed = TRUE, all.inside = TRUE))),
    sum(obs$observed <= 1))
  expect_true(all(cv$blqLo <= cv$blqHi))
  # the prediction correction leaves the censored fractions and their panel
  cpv <- plotVPC(cf, nsim = 100, bins = 3, pc = TRUE, seed = 1)
  expect_equal(attr(cpv, "vpc")[c("blq", "blqLo", "blqMid", "blqHi")], cv[c("blq", "blqLo", "blqMid", "blqHi")])
  # on a log scale the correction shifts by the difference of the population predictions
  lv <- attr(plotVPC(fit, nsim = 200, bins = 3, pc = TRUE, logScale = TRUE, seed = 1), "vpc")
  pr <- predict(fit)
  pr <- pr[pr$source == "obs", ]
  bin <- findInterval(pr$time, sort(unique(quantile(pr$time, seq(0, 1, length.out = 4), names = FALSE))),
                      rightmost.closed = TRUE, all.inside = TRUE)
  shifted <- pr$observed + ave(pr$PRED, bin, FUN = median) - pr$PRED
  mid <- lv[lv$percentile == "50%", ]
  expect_equal(mid$observed, as.vector(tapply(shifted, bin, median)[as.character(mid$bin)]))
})

test_that("the residual plots of a fit draw its residuals, covariates come from the condition grid", {
  skip_if_not_installed("ggplot2")
  m <- .lin()
  fit <- em(m$obj, randomEffects(m$p), c(a = 0.5, b = 0.2), verbose = FALSE)
  p <- plotResiduals(fit, type = c("IWRES", "CWRES", "NPDE"), nsim = 200, seed = 1)
  expect_equal(p$data$y[p$data$type == "CWRES" & p$data$panel == "vs PRED"], residuals(fit, "CWRES"))
  cov <- data.frame(dose = m$tab$dose, arm = ifelse(seq_len(nrow(m$tab)) %% 4 < 2, "x", "y"),
                    row.names = rownames(m$tab))
  expect_s3_class(plotRandomEffects(fit, covariates = cov["arm"]), "ggplot")
  expect_error(plotRandomEffects(fit, covariates = "weight"), "condition grid")
})

test_that("FOCEI stops with an error where -2 log L is not finite at the start", {
  m <- .lin()
  d <- m$data
  d[["s1a"]]$sigma[1] <- 0
  obj <- normL2(d, m$g * Xt() * m$p)
  expect_error(em(obj, randomEffects(m$p), c(a = 0.5, b = 0.2), verbose = FALSE), "not finite")
})

test_that("a failed prediction gives the FOCEI kernel frames it refuses, and keeps the message", {
  failure <- new.env(parent = emptyenv())
  guarded <- dMod2:::.guardPredictor(function(timesList, parsList, conditions, fixed) stop("no prediction"), failure)
  out <- guarded(list(1), list(c(a = 1)), "s1a", NULL)
  expect_identical(failure$message, "no prediction")
  expect_identical(out$pred, list("failed"))
})

test_that("a rejected outer step of dense FOCEI leaves the modes of the point it stands on", {
  m <- .lin()
  start <- c(a = 3, b = -2)
  at <- function(iterlim, rinit)
    em(m$obj, randomEffects(m$p), start, verbose = FALSE, sparse = FALSE,
       control = list(trustControl = list(iterlim = iterlim, rinit = rinit)))
  # an oversized first step is rejected, so the fit stays at the start
  one <- at(1, 100)
  expect_equal(one$argument, at(0, 1)$argument)
  expect_equal(one$etaModes, at(0, 1)$etaModes, tolerance = 1e-10)
})
