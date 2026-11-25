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
  # two independent species with four candidate rates in one group
  ks <- c("a0", "a1", "b0", "b1")
  x2 <- Xs(odemodel(eqnvec(x = "a0 - a1*x", y = "b0 - b1*y"), modelname = "spk_x2", compile = FALSE))
  g2 <- Y(eqnvec(x = "x", y = "y"), x2, modelname = "spk_g2", compile = FALSE)
  ec <- insert(branch(define(NULL, "p <- p", p = c("x", "y", ks)), table = data.frame(row.names = "C1", one = 1)),
               "k <- exp(k)", k = ks)
  ec <- insert(ec, k ~ 0 + (1 | gr(dist = slab(prob = 0.5, mean = log(0.5)))), k = ks)
  pEC <- P(ec, modelname = "spk_pec", compile = FALSE)
  compile(g, x, err, xFF, gFF, pSel, pProb, pGrp, pDef, pDefFF, pMix, pMixOff, pNorm, pLin, g2, x2, pEC,
          output = "spk_all", cores = test_cores())
  list(one = one, two = two, counts = counts, shift = shift, mixed = mixed, pure = pure, x = x, g = g, err = err,
       xFF = xFF, gFF = gFF, pSel = pSel, pProb = pProb, pGrp = pGrp, pDef = pDef, pDefFF = pDefFF, pMix = pMix,
       pMixOff = pMixOff, pNorm = pNorm, pLin = pLin, x2 = x2, g2 = g2, pEC = pEC, dir = .spkDir)
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
  r <- .spkFit(.spk$one, .spk$pSel, init = c(x = 1))
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

test_that("with counts and an error model of the prediction the production term is switched off as well", {
  fit <- .spkCounts("focei")$fit
  expect_true(fit$converged)
  p0 <- setNames(fit$random$p0, fit$random$eta)
  expect_gt(p0[["k0"]], 0.7)
  expect_lt(p0[["k1"]], 1e-6)
  expect_equal(exp(fit$random$mode[fit$random$eta == "k1"]), 0.5, tolerance = 0.15)
  expect_equal(unname(fit$argument["x"]), 100, tolerance = 0.1)
})

test_that("quadrature per subset agrees with a tensor Gauss-Hermite rule", {
  r <- .spkCounts("quadrature")
  fit <- r$fit
  expect_true(fit$converged)
  expect_equal(fit$method, "quadrature")
  p0 <- setNames(fit$random$p0, fit$random$eta)

  # -2 log Z of a subset by a product rule of 20 nodes per coordinate around the mode
  withr::local_dir(.spk$dir)
  tau <- 1
  etas <- fit$random$parameter
  gh <- dMod2:::.gaussHermite(20)
  m2logZ <- function(on) {
    prior <- constraintL2(setNames(rep(log(0.5), length(on)), on), sigma = tau)
    fx <- c(fit$argument["x"], setNames(rep(-Inf, length(setdiff(etas, on))), setdiff(etas, on)))
    tr <- trust(r$obj + prior, setNames(rep(log(0.5), length(on)), on), fixed = fx, hessianMethod = "gn")
    R <- chol(tr$hessian / 2)
    z <- as.matrix(expand.grid(rep(list(gh$x), length(on))))
    w <- apply(as.matrix(expand.grid(rep(list(gh$w), length(on)))), 1, prod)
    lt <- vapply(seq_len(nrow(z)), function(i) {
      eta <- setNames(tr$argument + sqrt(2) * backsolve(R, z[i, ]), on)
      -(r$obj(eta, fixed = fx, deriv = FALSE)$value + prior(eta, deriv = FALSE)$value) / 2 + sum(z[i, ]^2)
    }, 0) + log(w)
    top <- max(lt)
    -2 * (top + log(sum(exp(lt - top))) + length(on) / 2 * log(2) - sum(log(diag(R)))) +
      length(on) * log(2 * pi * tau^2)
  }
  z1 <- m2logZ("eta_k1_all")
  z2 <- m2logZ(c("eta_k0_all", "eta_k1_all"))
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
  ec <- fit(control = list(sparseFit = list(budget = 8)))
  expect_true(ec$converged)
  expect_equal(ec$random$p0, exact$random$p0, tolerance = 0.1)
  expect_equal(ec$random$p0 > 0.5, exact$random$p0 > 0.5)
  expect_equal(ec$value, exact$value, tolerance = 0.01)
  expect_equal(ec$argument, exact$argument, tolerance = 1e-3)
})
