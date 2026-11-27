# Spike-and-slab random effects with a spike at -Inf: selection of reaction rates by em().

.spkDir <- file.path(tempdir(), "spikeSubsets")
dir.create(.spkDir, showWarnings = FALSE)
.spk <- withr::with_dir(.spkDir, {
  times <- seq(0, 10, by = 0.5)
  decay <- function(x0, k0, k) k0 / k + (x0 - k0 / k) * exp(-k * times)
  set.seed(1)
  noisy <- function(v) data.frame(name = "x", time = times, value = v + rnorm(length(times), 0, 0.05), sigma = 0.05)
  # data with and without the production term, one condition or several
  one <- datalist(C1 = noisy(decay(1, 0, 0.5)))
  two <- datalist(A = noisy(decay(1, 1, 0.5)), B = noisy(decay(1, 0, 0.5)))
  counts <- datalist(C1 = data.frame(name = "x", time = times, sigma = NA,
                                     value = rpois(length(times), 100 * exp(-0.5 * times))))
  shift <- datalist(A = noisy(decay(1, 1, 0.5)), B = noisy(decay(1, 1.6, 0.5)))
  ids <- paste0("i", 1:6)
  k1i <- 0.5 * exp(c(0.2, -0.15, 0.1, -0.2, 0.15, -0.1))
  individuals <- function(k0) do.call(datalist, setNames(lapply(1:6, function(i) noisy(decay(1, k0[i], k1i[i]))), ids))
  mixed <- individuals(rep(1:0, each = 3))
  pure <- individuals(rep(0, 6))

  ff <- "forward-forward"
  x <- Xs(odemodel(eqnvec(x = "k0 - k1*x"), modelname = "spk_x", compile = FALSE))
  g <- Y(eqnvec(x = "x"), x, modelname = "spk_g", compile = FALSE)
  err <- Y(c(x = "sqrt(x)"), g, attachInput = FALSE, modelname = "spk_err", compile = FALSE)
  xFF <- Xs(odemodel(eqnvec(x = "k0 - k1*x"), modelname = "spk_xff", derivMode = ff, compile = FALSE))
  gFF <- Y(eqnvec(x = "x"), xFF, modelname = "spk_gff", derivMode = ff, compile = FALSE)
  base <- function(rows) insert(branch(define(NULL, "p <- p", p = c("x", "k0", "k1")),
                                       table = data.frame(row.names = rows, ID = rows)),
                                "k <- exp(k)", k = c("k0", "k1"))
  slab <- function(...) spikeSlab(tau = 1, spike = -Inf, ...)
  sel <- insert(base("C1"), k ~ 0 + (1 | gr(dist = slab(prob = 0.5, mean = log(0.5)))), k = c("k0", "k1"))
  selProb <- insert(base("C1"), k ~ 0 + (1 | gr(dist = slab(mean = log(0.5)))), k = c("k0", "k1"))
  grp <- insert(base(c("A", "B")), k0 ~ 1 + (1 | gr(ID, dist = slab(prob = 0.5))))
  def <- insert(base(c("A", "B")), k0 ~ 1 + (1 | gr(ID, dist = spikeSlab(prob = 0.5, tau = 1))))
  mix <- insert(base(ids), k1 ~ 1 + (1 | ID), k0 ~ 0 + (1 | gr(ID, dist = slab(prob = 0.5))))
  mixOff <- insert(base(ids), k1 ~ 1 + (1 | ID), k0 ~ 0 + (1 | gr(ID, dist = slab(prob = 1e-6))))
  norm <- insert(branch(define(NULL, "p <- p", p = c("x", "k0", "k1")),
                        table = data.frame(row.names = ids, ID = ids)), "k0 <- 0")
  norm <- insert(insert(norm, "k <- exp(k)", k = "k1"), k1 ~ 1 + (1 | ID))
  linear <- insert(branch(define(NULL, "p <- p", p = c("x", "k0", "k1")),
                          table = data.frame(row.names = "C1", one = 1)),
                   k0 ~ 0 + (1 | gr(dist = spikeSlab(spike = -Inf))))
  pSel <- P(sel, modelname = "spk_psel", compile = FALSE)
  pProb <- P(selProb, modelname = "spk_pprob", compile = FALSE)
  pGrp <- P(grp, modelname = "spk_pgrp", compile = FALSE)
  pDef <- P(def, modelname = "spk_pdef", compile = FALSE)
  pDefFF <- P(def, modelname = "spk_pdefff", derivMode = ff, compile = FALSE)
  pMix <- P(mix, modelname = "spk_pmix", derivMode = ff, compile = FALSE)
  pMixOff <- P(mixOff, modelname = "spk_pmixoff", derivMode = ff, compile = FALSE)
  pNorm <- P(norm, modelname = "spk_pnorm", derivMode = ff, compile = FALSE)
  pLin <- P(linear, modelname = "spk_plin", compile = FALSE)
  ten <- insert(branch(define(NULL, "p <- p", p = c("x", "k0", "k1")), table = data.frame(row.names = "C1", one = 1)),
                "k <- exp10(k)", k = c("k0", "k1"))
  ten <- insert(ten, k ~ 0 + (1 | gr(dist = spikeSlab(prob = 0.5, mean = log10(0.5), tau = 1 / log(10), spike = -Inf))),
                k = c("k0", "k1"))
  pTen <- P(ten, modelname = "spk_pten", compile = FALSE)
  # two independent species with four candidate rates in one group
  ks <- c("a0", "a1", "b0", "b1")
  x2 <- Xs(odemodel(eqnvec(x = "a0 - a1*x", y = "b0 - b1*y"), modelname = "spk_x2", compile = FALSE))
  g2 <- Y(eqnvec(x = "x", y = "y"), x2, modelname = "spk_g2", compile = FALSE)
  ec <- insert(branch(define(NULL, "p <- p", p = c("x", "y", ks)), table = data.frame(row.names = "C1", one = 1)),
               "k <- exp(k)", k = ks)
  ec <- insert(ec, k ~ 0 + (1 | gr(dist = slab(prob = 0.5, mean = log(0.5)))), k = ks)
  pEC <- P(ec, modelname = "spk_pec", compile = FALSE)
  compile(g, x, err, xFF, gFF, pSel, pProb, pGrp, pDef, pDefFF, pMix, pMixOff, pNorm, pLin, pTen, g2, x2, pEC,
          output = "spk_all", cores = test_cores())
  list(one = one, two = two, counts = counts, shift = shift, mixed = mixed, pure = pure, x = x, g = g, err = err,
       xFF = xFF, gFF = gFF, pSel = pSel, pProb = pProb, pGrp = pGrp, pDef = pDef, pDefFF = pDefFF, pMix = pMix,
       pMixOff = pMixOff, pNorm = pNorm, pLin = pLin, pTen = pTen, x2 = x2, g2 = g2, pEC = pEC, dir = .spkDir)
})

.spkFit <- function(data, p, ..., errmodel = NULL, method = "focei") {
  withr::local_dir(.spk$dir)
  obj <- normL2(data, .spk$g * .spk$x * p, errmodel = errmodel)
  list(fit = em(obj, randomEffects(p), method = method, verbose = FALSE, ...), obj = obj)
}

# fits that several tests share, each computed once
.spkFits <- new.env()
.spkShared <- function(key, expr) {
  if (is.null(.spkFits[[key]])) .spkFits[[key]] <- expr
  .spkFits[[key]]
}
.spkCounts <- function(method) .spkShared(paste0("counts_", method),
  .spkFit(.spk$counts, .spk$pSel, init = c(x = 100), errmodel = .spk$err, method = method))
.spkMixed <- function(method) .spkShared(paste0("mixed_", method), withr::with_dir(.spk$dir, {
  set.seed(1)
  em(normL2(.spk$mixed, .spk$gFF * .spk$xFF * .spk$pMix), randomEffects(.spk$pMix),
     init = c(x = 1, k1 = log(0.5), omega_k1_k1 = log(0.2)), method = method, verbose = FALSE)
}))

test_that("a production term the data do not need is switched off, the decay rate stays", {
  r <- .spkShared("one_sel", .spkFit(.spk$one, .spk$pSel, init = c(x = 1)))
  fit <- r$fit
  expect_true(fit$converged)
  p0 <- setNames(fit$random$p0, fit$random$eta)
  expect_gt(p0[["k0"]], 0.99)
  expect_lt(p0[["k1"]], 1e-6)
  expect_equal(fit$random$mode[fit$random$eta == "k0"], -Inf)
  expect_equal(exp(fit$random$mode[fit$random$eta == "k1"]), 0.5, tolerance = 0.1)
  expect_equal(unname(fit$argument["x"]), 1, tolerance = 0.05)

  # p0 from the Laplace integrals of the subsets of non-negligible weight, by trust() on objective plus prior
  withr::local_dir(.spk$dir)
  modes <- setNames(fit$random$mode, fit$random$parameter)
  tau <- 1
  m2logZ <- function(on) {
    prior <- constraintL2(setNames(rep(log(0.5), length(on)), on), sigma = tau)
    off <- setdiff(names(modes), on)
    start <- modes[on]
    start[!is.finite(start)] <- log(0.5)
    tr <- trust(r$obj + prior, start, fixed = c(fit$argument["x"], setNames(rep(-Inf, length(off)), off)),
                hessianMethod = "gn")
    tr$value + length(on) * log(tau^2) + as.numeric(determinant(tr$hessian / 2)$modulus)
  }
  z1 <- m2logZ("eta_k1_all")
  z2 <- m2logZ(c("eta_k0_all", "eta_k1_all"))
  expect_equal(p0[["k0"]], 1 / (1 + exp(-(z2 - z1) / 2)), tolerance = 1e-4)
})

test_that("exp10() with the slab scaled by 1 / log(10) gives the fit of exp()", {
  ten <- .spkFit(.spk$one, .spk$pTen, init = c(x = 1))$fit
  e <- .spkShared("one_sel", .spkFit(.spk$one, .spk$pSel, init = c(x = 1)))$fit
  expect_equal(ten$random$p0, e$random$p0, tolerance = 1e-4)
  expect_equal(ten$value, e$value, tolerance = 1e-6)
  expect_equal(10^ten$random$mode, exp(e$random$mode), tolerance = 1e-4)
})

test_that("with counts and an error model of the prediction the production term is switched off as well", {
  fit <- .spkCounts("focei")$fit
  expect_true(fit$converged)
  p0 <- setNames(fit$random$p0, fit$random$eta)
  expect_gt(p0[["k0"]], 0.7)
  expect_lt(p0[["k1"]], 1e-6)
  expect_equal(exp(fit$random$mode[fit$random$eta == "k1"]), 0.5, tolerance = 0.15)
  expect_equal(unname(fit$argument["x"]), 100, tolerance = 0.1)
})

test_that("EC with an error model of the prediction selects as the exact sum", {
  exact <- .spkFit(.spk$counts, .spk$pSel, init = c(x = 100), errmodel = .spk$err,
                   control = list(sparseFit = list(phases = "exact")))$fit
  ec <- .spkFit(.spk$counts, .spk$pSel, init = c(x = 100), errmodel = .spk$err,
                control = list(sparseFit = list(phases = "ec")))$fit
  expect_true(ec$converged)
  expect_equal(ec$random$p0 > 0.5, exact$random$p0 > 0.5)
  expect_equal(ec$value, exact$value, tolerance = 0.01)
})

# -2 log Z of the subset `on` of the etas of `fit` (slab N(log(0.5), tau^2), the others at -Inf) at the estimate x of
# `fit`, by a product rule of 20 Gauss-Hermite nodes per coordinate around the mode
.spkTensorM2logZ <- function(obj, fit, on, tau = 1) {
  withr::local_dir(.spk$dir)
  etas <- fit$random$parameter
  gh <- dMod2:::.gaussHermite(20)
  prior <- constraintL2(setNames(rep(log(0.5), length(on)), on), sigma = tau)
  fx <- c(fit$argument["x"], setNames(rep(-Inf, length(setdiff(etas, on))), setdiff(etas, on)))
  tr <- trust(obj + prior, setNames(rep(log(0.5), length(on)), on), fixed = fx, hessianMethod = "gn")
  R <- chol(tr$hessian / 2)
  z <- as.matrix(expand.grid(rep(list(gh$x), length(on))))
  w <- apply(as.matrix(expand.grid(rep(list(gh$w), length(on)))), 1, prod)
  lt <- vapply(seq_len(nrow(z)), function(i) {
    eta <- setNames(tr$argument + sqrt(2) * backsolve(R, z[i, ]), on)
    -(obj(eta, fixed = fx, deriv = FALSE)$value + prior(eta, deriv = FALSE)$value) / 2 + sum(z[i, ]^2)
  }, 0) + log(w)
  top <- max(lt)
  -2 * (top + log(sum(exp(lt - top))) + length(on) / 2 * log(2) - sum(log(diag(R)))) +
    length(on) * log(2 * pi * tau^2)
}

test_that("quadrature per subset agrees with a tensor Gauss-Hermite rule", {
  r <- .spkCounts("quadrature")
  fit <- r$fit
  expect_true(fit$converged)
  expect_equal(fit$method, "quadrature")
  p0 <- setNames(fit$random$p0, fit$random$eta)
  z1 <- .spkTensorM2logZ(r$obj, fit, "eta_k1_all")
  z2 <- .spkTensorM2logZ(r$obj, fit, c("eta_k0_all", "eta_k1_all"))
  expect_lt(abs(p0[["k0"]] - 1 / (1 + exp(-(z2 - z1) / 2))), 5e-3)
  laplace <- .spkCounts("focei")$fit
  expect_false(isTRUE(all.equal(laplace$random$p0, fit$random$p0, tolerance = 1e-3)))
})

test_that("quadrature takes the default spike at 0 as well, near the Laplace sum", {
  init <- c(x = 1, k0 = 0, k1 = log(0.5))
  quad <- .spkFit(.spk$shift, .spk$pDef, init = init, method = "quadrature")$fit
  withr::local_dir(.spk$dir)
  laplace <- suppressMessages(em(normL2(.spk$shift, .spk$gFF * .spk$xFF * .spk$pDefFF),
                                 randomEffects(.spk$pDefFF), init = init, verbose = FALSE))
  expect_true(quad$converged)
  expect_equal(quad$random$p0, laplace$random$p0, tolerance = 0.05)
  expect_equal(quad$value, laplace$value, tolerance = 0.01)
})

test_that("normal random effects beside a spike at -Inf: focei and quadrature select the same", {
  fl <- .spkMixed("focei")
  fq <- .spkMixed("quadrature")
  for (f in list(fl, fq)) {
    expect_true(f$converged)
    p0 <- setNames(f$random$p0[f$random$eta == "k0"], f$random$group[f$random$eta == "k0"])
    expect_true(all(p0[c("i1", "i2", "i3")] < 0.05))
    expect_true(all(p0[c("i4", "i5", "i6")] > 0.5))
    expect_equal(exp(unname(f$argument["k1"])), 0.5, tolerance = 0.15)
  }
  expect_equal(fq$value, fl$value, tolerance = 0.01)
})

test_that("SAEM with a spike at -Inf beside normal random effects agrees with quadrature", {
  fq <- .spkMixed("quadrature")
  fs <- .spkMixed("saem")
  expect_true(fs$converged)
  k0 <- fs$random$eta == "k0"
  expect_equal(fs$random$p0[k0], fq$random$p0[k0], tolerance = 0.05)
  expect_equal(fs$random$mode[k0][4:6], rep(-Inf, 3))
  expect_equal(fs$argument[c("x", "k1")], fq$argument[c("x", "k1")], tolerance = 0.03)
  expect_equal(exp(fs$argument[["omega_k1_k1"]]), exp(fq$argument[["omega_k1_k1"]]), tolerance = 0.2)
})

test_that("SAEM switches a production term off without groups", {
  set.seed(2)
  fit <- .spkFit(.spk$one, .spk$pSel, init = c(x = 1), method = "saem")$fit
  p0 <- setNames(fit$random$p0, fit$random$eta)
  expect_gt(p0[["k0"]], 0.9)
  expect_lt(p0[["k1"]], 0.01)
  expect_equal(exp(fit$random$mode[fit$random$eta == "k1"]), 0.5, tolerance = 0.1)
})

test_that("with prob near 0 and no production the subset quadrature is the quadrature of the normal model", {
  withr::local_dir(.spk$dir)
  init <- c(x = 1, k1 = log(0.5), omega_k1_k1 = log(0.2))
  ref <- suppressMessages(em(normL2(.spk$pure, .spk$gFF * .spk$xFF * .spk$pNorm), randomEffects(.spk$pNorm),
                             init = init, method = "quadrature", verbose = FALSE))
  sub <- em(normL2(.spk$pure, .spk$gFF * .spk$xFF * .spk$pMixOff), randomEffects(.spk$pMixOff), init = init,
            method = "quadrature", verbose = FALSE)
  # the weight of the spikes enters -2 log L
  expect_equal(sub$value + 2 * 6 * log1p(-1e-6), ref$value, tolerance = 1e-3)
  expect_equal(sub$argument[names(init)], ref$argument[names(init)], tolerance = 1e-3)
})

test_that("an estimated prob is the mean inclusion probability", {
  fit <- .spkFit(.spk$one, .spk$pProb, init = c(x = 1))$fit
  expect_true(fit$converged)
  expect_equal(unname(fit$argument["prob_all"]), mean(1 - fit$random$p0), tolerance = 1e-4)
})

test_that("a spike per group switches the production term off in the group without production", {
  fit <- .spkFit(.spk$two, .spk$pGrp, init = c(x = 1, k0 = 0, k1 = log(0.5)))$fit
  p0 <- setNames(fit$random$p0, fit$random$group)
  expect_lt(p0[["A"]], 0.1)
  expect_gt(p0[["B"]], 0.9)
  expect_equal(exp(unname(fit$argument["k0"])), 1, tolerance = 0.2)
})

test_that("a spike at -Inf needs an exponential, and laplaceEM refuses it", {
  withr::local_dir(.spk$dir)
  objLin <- normL2(.spk$one, .spk$g * .spk$x * .spk$pLin)
  expect_error(em(objLin, randomEffects(.spk$pLin), init = c(x = 1, k1 = 0.5)),
               "summand in the exponent")
  obj <- normL2(.spk$one, .spk$g * .spk$x * .spk$pSel)
  expect_error(em(obj, randomEffects(.spk$pSel), init = c(x = 0.5), method = "laplaceEM"), "method = \"focei\"")
})

test_that("EC beyond the budget agrees with the sum over all subsets", {
  withr::local_dir(.spk$dir)
  times <- seq(0, 10, by = 0.5)
  set.seed(2)
  noisy <- function(n, v) data.frame(name = n, time = times, value = v + rnorm(length(times), 0, 0.05),
                                     sigma = 0.05)
  data <- datalist(C1 = rbind(noisy("x", 2 - 1.5 * exp(-0.5 * times)), noisy("y", 2 * exp(-0.3 * times))))
  fit <- function(...) em(normL2(data, .spk$g2 * .spk$x2 * .spk$pEC), randomEffects(.spk$pEC),
                          init = c(x = 0.5, y = 2), verbose = FALSE, ...)
  exact <- fit()
  ec <- fit(control = list(sparseFit = list(termBudget = 8)))
  expect_true(ec$converged)
  expect_equal(ec$random$p0, exact$random$p0, tolerance = 0.1)
  expect_equal(ec$random$p0 > 0.5, exact$random$p0 > 0.5)
  expect_equal(ec$value, exact$value, tolerance = 0.01)
  expect_equal(ec$argument, exact$argument, tolerance = 1e-3)
})

# data of the four candidate rates; em() of the fixture pEC with sparseFit settings
.spkECData <- local({
  times <- seq(0, 10, by = 0.5)
  set.seed(2)
  noisy <- function(n, v) data.frame(name = n, time = times, value = v + rnorm(length(times), 0, 0.05),
                                     sigma = 0.05)
  datalist(C1 = rbind(noisy("x", 2 - 1.5 * exp(-0.5 * times)), noisy("y", 2 * exp(-0.3 * times))))
})
.spkECObj <- function() withr::with_dir(.spk$dir, normL2(.spkECData, .spk$g2 * .spk$x2 * .spk$pEC))
.spkECFit <- function(init, ...) withr::with_dir(.spk$dir, em(.spkECObj(), randomEffects(.spk$pEC), init = init,
                                                          verbose = FALSE, control = list(sparseFit = list(...))))

test_that("mean field selects as the exact sum", {
  near <- c(x = 0.5, y = 2)
  exact <- .spkShared("ec_exact", .spkECFit(near, phases = "exact"))
  mf <- .spkECFit(near, phases = "meanfield")
  expect_true(mf$converged)
  expect_equal(mf$random$p0 > 0.5, exact$random$p0 > 0.5)
  expect_equal(mf$random$p0, exact$random$p0, tolerance = 0.01)
  expect_equal(mf$value, exact$value, tolerance = 1e-3)
  expect_equal(mf$phases$approx, "meanfield")
  # the default phases end with the exact sum
  def <- .spkECFit(near)
  expect_equal(def$phases$approx, c("meanfield", "ec", "exact"))
  expect_equal(def$value, exact$value, tolerance = 1e-6)
})

test_that("from a far start EC alone fails, after a burn-in and mean field it converges", {
  far <- c(x = 20, y = 20)
  expect_error(.spkECFit(far, phases = "ec"),
               "^em: the expectation consistent \\(EC\\) approximation of group all did not converge")
  ref <- .spkECFit(c(x = 0.5, y = 2), phases = "ec")
  for (start in c("spike", "slab")) {
    fit <- .spkECFit(far, burnin = list(start = start), phases = list(meanfield = list(), ec = list(maxECM = 50)))
    expect_true(fit$converged)
    expect_equal(fit$phases$approx, c("burnin", "meanfield", "ec"))
    expect_true(all(fit$phases$converged))
    expect_equal(fit$phases$value[3], fit$value)
    expect_equal(fit$value, ref$value, tolerance = 1e-6)
    expect_equal(fit$argument, ref$argument, tolerance = 1e-4)
    expect_equal(fit$random$p0, ref$random$p0, tolerance = 1e-4)
    expect_equal(unique(fit$stageTrace$phase), 1:3)
  }
  expect_equal(fit$burnin, ref$random$parameter)
})

test_that("a burn-in from a given subset, and msem() with random subsets", {
  far <- c(x = 20, y = 20)
  ref <- .spkECFit(c(x = 0.5, y = 2), phases = "ec")
  given <- c("eta_a0_all", "eta_a1_all", "eta_b1_all")
  fit <- .spkECFit(far, burnin = list(start = given), phases = c("meanfield", "ec"))
  expect_equal(fit$burnin, given)
  expect_equal(fit$value, ref$value, tolerance = 1e-6)
  expect_error(.spkECFit(far, burnin = list(start = c("eta_a0_all", "eta_z_all"))),
               "names no spike-and-slab random effect: eta_z_all")

  withr::local_dir(.spk$dir)
  ctl <- list(sparseFit = list(burnin = TRUE, phases = c("meanfield", "ec")))
  run <- function() {
    set.seed(3)
    msem(.spkECObj(), randomEffects(.spk$pEC), far, control = ctl, fits = 4, sd = 0.3)
  }
  ms <- run()
  tab <- ms$multistart
  expect_true(all(tab$converged))
  expect_true(all(is.na(tab$phase)))
  expect_true(all(unlist(strsplit(tab$subset, ",")) %in% ref$random$parameter))
  expect_gt(length(unique(tab$subset)), 1)
  expect_equal(run()$multistart$subset, tab$subset)
  expect_equal(ms$value, ref$value, tolerance = 1e-6)
  set.seed(1)
  expect_error(msem(.spkECObj(), randomEffects(.spk$pEC), far, control = list(sparseFit = list(phases = "ec")),
                    fits = 2, sd = 0.3, retry = FALSE),
               "every fit failed.*did not converge")
})

test_that("burnin, phases and termBudget are checked, a failed phase is named", {
  near <- c(x = 0.5, y = 2)
  e <- tryCatch(.spkECFit(c(x = 20, y = 20), phases = c("ec", "exact")), error = function(e) e)
  expect_s3_class(e, "dMod2PhaseError")
  expect_equal(e$phase, 1L)
  expect_match(conditionMessage(e), "^em: phase 1 \\(ec\\): the expectation consistent")
  expect_error(.spkECFit(near, phases = c("ec", "laplace")), "`control\\$sparseFit\\$phases` names each of")
  expect_error(.spkECFit(near, phases = c("ec", "ec")), "at most once")
  expect_error(.spkECFit(near, phases = list(ec = list(maxECM = 0))), "`maxECM` \\(a whole number of at least 1\\)")
  expect_error(.spkECFit(near, burnin = "random"), "`control\\$sparseFit\\$burnin` is FALSE, TRUE or a list")
  expect_error(.spkECFit(near, burnin = list(start = "spike", iterlim = 0)), "`iterlim` \\(a whole number")
  expect_error(.spkECFit(near, budget = 8), "unknown entries of `control\\$sparseFit`: budget")
  expect_error(.spkECFit(near, approx = "ec"), "unknown entries of `control\\$sparseFit`: approx")
  # beyond termBudget the phase "exact" is skipped after EC, and stops without it
  expect_error(.spkECFit(near, phases = c("meanfield", "exact"), termBudget = 8),
               "16 subsets, more than `termBudget` \\(8\\); add \"ec\"")
  ec <- .spkECFit(near, phases = "ec")
  skipped <- .spkECFit(near, phases = c("ec", "exact"), termBudget = 8)
  expect_equal(skipped$value, ec$value, tolerance = 1e-8)
  expect_equal(skipped$phases$iterations[2], 0L)
  expect_true(is.na(skipped$phases$converged[2]))
  expect_true(skipped$converged)
})

# the random effects `re` with `prob` and `tau` of every spike-and-slab term replaced: a number or NA, or a hyperprior
.spkReprior <- function(re, prob = NA, tau = NULL) {
  slab <- re$terms$dist == "spikeSlab"
  if (inherits(prob, "hyperprior")) { re$terms$prob[slab] <- NA; re$terms$probPrior[slab] <- list(prob) }
  else re$terms$prob[slab] <- prob
  if (inherits(tau, "hyperprior")) { re$terms$tau[slab] <- NA; re$terms$tauPrior[slab] <- list(tau) }
  else if (!is.null(tau)) re$terms$tau[slab] <- tau
  re
}
# the four candidate rates with the data of y read at sigma `sy`
.spkBetaObj <- function(sy) {
  d <- .spkECData
  y <- d$C1$name == "y"
  d$C1$sigma[y] <- sy
  set.seed(5)
  if (sy > 0.05) d$C1$value[y] <- 2 * exp(-0.3 * seq(0, 10, by = 0.5)) + rnorm(sum(y), 0, sy)
  withr::with_dir(.spk$dir, normL2(d, .spk$g2 * .spk$x2 * .spk$pEC))
}
.spkBetaFit <- function(obj, random, init, ...) withr::with_dir(.spk$dir,
  em(obj, random, init = init, verbose = FALSE, control = list(sparseFit = list(...))))
# -2 log of the integral over prob in (0, 1) of the beta density times exp(-V(prob) / 2), by Gauss-Legendre
.spkOverProb <- function(V, a = 1, b = 1, n = 12) {
  off <- seq_len(n - 1L) / sqrt(4 * seq_len(n - 1L)^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(seq_len(n - 1L), 2:n)] <- J[cbind(2:n, seq_len(n - 1L))] <- off
  e <- eigen(J, symmetric = TRUE)
  rho <- (e$values + 1) / 2
  w <- e$vectors[1, ]^2
  v <- vapply(rho, V, 0)
  min(v) - 2 * log(sum(w * stats::dbeta(rho, a, b) * exp(-(v - min(v)) / 2)))
}

test_that("a beta prior on prob integrates the marginal of a given prob over the beta density", {
  obj <- .spkBetaObj(0.05)
  re <- randomEffects(.spk$pEC)
  init <- c(x = 0.5, y = 2)
  # -2 log L at the start, before the first M-step
  value <- function(r) .spkBetaFit(obj, r, init, phases = list(exact = list(maxECM = 1)))$phases$value[1]
  beta <- value(.spkReprior(re, betaPrior(2, 3)))
  expect_equal(beta, .spkOverProb(function(p) value(.spkReprior(re, p)), 2, 3), tolerance = 1e-8)
})

test_that("the groups of a block under a beta prior are summed jointly over their count", {
  withr::local_dir(.spk$dir)
  obj <- normL2(.spk$two, .spk$g * .spk$x * .spk$pGrp)
  re <- randomEffects(.spk$pGrp)
  init <- c(x = 1, k0 = 0, k1 = log(0.5))
  value <- function(r) em(obj, r, init = init, verbose = FALSE,
                          control = list(sparseFit = list(phases = list(exact = list(maxECM = 1)))))$phases$value[1]
  beta <- value(.spkReprior(re, betaPrior(1, 1)))
  expect_equal(beta, .spkOverProb(function(p) value(.spkReprior(re, p))), tolerance = 1e-8)
})

test_that("where empirical Bayes takes prob to 1, a beta prior keeps it inside", {
  obj <- .spkBetaObj(2)
  re <- randomEffects(.spk$pEC)
  init <- c(x = 0.5, y = 2)
  # an M-step may end at the resolution of the parameters, with a warning of trust()
  eb <- suppressWarnings(.spkBetaFit(obj, .spkReprior(re, NA), c(init, prob_all = 0.5), phases = "exact"))
  beta <- .spkBetaFit(obj, .spkReprior(re, betaPrior(1, 1)), init, phases = "exact")
  p0 <- function(f) setNames(f$random$p0, f$random$eta)
  expect_gt(eb$slab$estimate[1], 0.999)
  expect_lt(p0(eb)[["b0"]], 1e-3)
  expect_lt(beta$slab$estimate[1], 0.9)
  expect_gt(p0(beta)[["b0"]], 0.2)
  expect_equal(beta$slab$prior[1], "Beta(1, 1), integrated")
  expect_true(beta$slab$lower[1] < beta$slab$estimate[1] && beta$slab$estimate[1] < beta$slab$upper[1])
  expect_false("prob_all" %in% names(beta$argument))
})

test_that("mean field and EC with a beta prior sum as the exact phase", {
  re <- .spkReprior(randomEffects(.spk$pEC), betaPrior(1, 1))
  init <- c(x = 0.5, y = 2)
  for (sy in c(0.05, 1)) {
    obj <- .spkBetaObj(sy)
    exact <- .spkBetaFit(obj, re, init, phases = "exact")
    mf <- .spkBetaFit(obj, re, init, phases = "meanfield")
    expect_equal(mf$random$p0, exact$random$p0, tolerance = 0.01)
    expect_equal(mf$value, exact$value, tolerance = 1e-3)
  }
  obj <- .spkBetaObj(0.05)
  ec <- .spkBetaFit(obj, re, init, phases = "ec")
  exact <- .spkBetaFit(obj, re, init, phases = "exact")
  expect_equal(ec$random$p0 > 0.5, exact$random$p0 > 0.5)
  expect_equal(ec$value, exact$value, tolerance = 0.01)
})

test_that("tau under a scale prior is the posterior mode of the M-step", {
  m2Q <- function(tau, num, den, prior) den * log(tau^2) + num / tau^2 + dMod2:::.spikeTauPenalty(tau, prior)
  priors <- list(halfCauchyPrior(1), halfNormalPrior(0.5), invGammaPrior(2, 1))
  for (prior in priors) for (nd in list(c(3, 2), c(0.1, 3))) {
    mode <- dMod2:::.spikeTauMode(nd[1], nd[2], prior, 1)
    ref <- if (prior$dist == "invGamma")
      sqrt(stats::optimize(function(v) m2Q(sqrt(v), nd[1], nd[2], prior), c(1e-6, 1e3), tol = 1e-12)$minimum)
    else stats::optimize(function(t) m2Q(t, nd[1], nd[2], prior), c(1e-6, 1e3), tol = 1e-12)$minimum
    expect_equal(mode, ref, tolerance = 1e-6)
  }

  obj <- .spkBetaObj(0.05)
  re <- .spkReprior(randomEffects(.spk$pEC), betaPrior(1, 1), halfCauchyPrior(1))
  fit <- .spkBetaFit(obj, re, c(x = 0.5, y = 2, tau_all = 1), phases = "exact")
  expect_true(fit$converged)
  tau <- fit$argument[["tau_all"]]
  expect_equal(fit$hyperprior, dMod2:::.spikeTauPenalty(tau, halfCauchyPrior(1)))
  expect_equal(fit$slab$prior, c("Beta(1, 1), integrated", "half-Cauchy(1), mode"))
  free <- .spkBetaFit(obj, .spkReprior(randomEffects(.spk$pEC), betaPrior(1, 1), NA_real_),
                      c(x = 0.5, y = 2, tau_all = 1), phases = "exact")
  expect_lt(tau, free$argument[["tau_all"]])
})

test_that("quadrature takes a beta prior on prob and a scale prior on tau", {
  withr::local_dir(.spk$dir)
  obj <- normL2(.spk$counts, .spk$g * .spk$x * .spk$pSel, errmodel = .spk$err)
  re <- .spkReprior(randomEffects(.spk$pSel), betaPrior(1, 1))
  fit <- function(re, init, method) em(obj, re, init = init, method = method, verbose = FALSE)
  quad <- fit(re, c(x = 100), "quadrature")
  expect_true(quad$converged)
  # with k1 in the slab, k0 off or on: one or two of the two candidates, prior weights B(2, 2) and B(3, 1)
  z1 <- .spkTensorM2logZ(obj, quad, "eta_k1_all")
  z2 <- .spkTensorM2logZ(obj, quad, c("eta_k0_all", "eta_k1_all"))
  p0 <- setNames(quad$random$p0, quad$random$eta)
  expect_lt(abs(p0[["k0"]] - 1 / (1 + beta(3, 1) / beta(2, 2) * exp(-(z2 - z1) / 2))), 5e-3)

  reTau <- .spkReprior(randomEffects(.spk$pSel), betaPrior(1, 1), halfNormalPrior(0.5))
  tau <- fit(reTau, c(x = 100, tau_all = 0.5), "quadrature")
  expect_true(tau$converged)
  expect_equal(tau$hyperprior, dMod2:::.spikeTauPenalty(tau$argument[["tau_all"]], halfNormalPrior(0.5)))
  free <- fit(.spkReprior(randomEffects(.spk$pSel), betaPrior(1, 1), NA_real_), c(x = 100, tau_all = 0.5), "quadrature")
  expect_lt(tau$argument[["tau_all"]], free$argument[["tau_all"]])
})

test_that("the phase exact integrates tau over log tau under a scale prior with integrate = TRUE", {
  obj <- .spkBetaObj(0.05)
  re <- randomEffects(.spk$pEC)
  fx <- c(x = 0.5, y = 2)
  prior <- halfCauchyPrior(1, integrate = TRUE)
  fit <- withr::with_dir(.spk$dir, em(obj, .spkReprior(re, betaPrior(1, 1), prior), init = c(tau_all = 1), fixed = fx,
                                      verbose = FALSE, control = list(sparseFit = list(phases = "exact"))))
  expect_true(fit$converged)
  expect_equal(fit$slab$prior[2], "half-Cauchy(1), integrated")
  # the exact marginal at a given tau integrated over log tau by Simpson's rule
  V <- function(tau) suppressWarnings(withr::with_dir(.spk$dir, em(obj, .spkReprior(re, betaPrior(1, 1), tau),
    fixed = fx, verbose = FALSE, control = list(sparseFit = list(phases = "exact")))))$value
  u <- seq(-3, 2, length.out = 31)
  lf <- -(vapply(exp(u), V, 0) + vapply(exp(u), dMod2:::.spikeTauPenalty, 0, prior = prior)) / 2 + u
  w <- c(1, rep(c(4, 2), length.out = length(u) - 2), 1) * diff(u)[1] / 3
  expect_equal(fit$value, -2 * (max(lf) + log(sum(w * exp(lf - max(lf))))), tolerance = 1e-4)
  counts <- withr::with_dir(.spk$dir, normL2(.spk$counts, .spk$g * .spk$x * .spk$pSel, errmodel = .spk$err))
  expect_error(withr::with_dir(.spk$dir, em(counts, .spkReprior(randomEffects(.spk$pSel), 0.5, prior),
                                            init = c(x = 100, tau_all = 1), method = "quadrature", verbose = FALSE)),
               "takes the posterior mode of tau; drop `integrate = TRUE`")
})

test_that("SAEM, msem() and print() take the hyperpriors", {
  re <- .spkReprior(randomEffects(.spk$pEC), betaPrior(1, 1))
  expect_output(print(re), "prob Beta\\(1, 1\\), integrated")
  expect_output(print(spikeSlab(prob = betaPrior(2, 3), tau = halfCauchyPrior(1))),
                "spikeSlab\\(prob = betaPrior\\(2, 3\\), tau = halfCauchyPrior\\(1\\)")
  obj <- .spkBetaObj(2)
  init <- c(x = 0.5, y = 2)
  exact <- .spkBetaFit(obj, re, init, phases = "exact")
  set.seed(2)
  saem <- withr::with_dir(.spk$dir, em(obj, re, init = init, method = "saem", verbose = FALSE))
  expect_lt(max(abs(saem$random$p0 - exact$random$p0)), 0.1)
  expect_equal(saem$value, exact$value, tolerance = 1e-3)
  expect_equal(saem$slab$prior[1], "Beta(1, 1), integrated")

  set.seed(1)
  ms <- suppressWarnings(withr::with_dir(.spk$dir, msem(obj, re, init, control = list(sparseFit = list(phases = "exact")),
                                                        fits = 3, sd = 0.2, keepFits = TRUE)))
  expect_true(all(ms$multistart$converged))
  expect_true(all(vapply(ms$fits, function(f) !"prob_all" %in% names(f$argument), TRUE)))
  expect_equal(ms$slab$prior[1], "Beta(1, 1), integrated")

  withr::local_dir(.spk$dir)
  grp <- .spkReprior(randomEffects(.spk$pGrp), betaPrior(1, 1))
  expect_error(em(normL2(.spk$two, .spk$g * .spk$x * .spk$pGrp), grp, init = c(x = 1, k0 = 0, k1 = log(0.5)),
                  method = "saem", verbose = FALSE),
               "takes betaPrior\\(\\) for spike-and-slab random effects that reach")
  expect_error(spikeSlab(prob = halfCauchyPrior()), "the prior of `prob` is betaPrior")
  expect_error(betaPrior(0, 1), "positive numbers")
})
