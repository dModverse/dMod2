# SAEM with spike-and-slab random effects: the indicators against the exact posterior of a linear
# model with a normal intercept, the selection and the estimates against FOCEI, beside normal and Dirichlet-process effects.

# The decay model y = A exp(-k dose) of test-EMsparse.R, five groups with two conditions each, the
# linear model y = a + b dose with a known sigma, and a decay model with 14 spike-and-slab random
# effects per group, three of 56 not 0.
.saemSparse <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    withr::local_preserve_seed()
    dir <- file.path(tempdir(), "saemsparse_models")
    dir.create(dir, showWarnings = FALSE)
    withr::local_dir(dir)
    groups <- paste0("L", 1:5)
    tab <- data.frame(ID = rep(groups, each = 2), dose = rep(c(0.5, 2), 5),
                      row.names = paste0(rep(groups, each = 2), c("a", "b")))
    base <- function() {
      tr <- branch(define(NULL, "x <- x", x = c("k", "A", "s")), table = tab)
      tr <- define(tr, "dose <- d", d = dose)
      insert(tr, "x <- exp(x)", x = c("k", "A", "s"))
    }
    # second derivatives forward over reverse
    .dm <- c("forward", "reverse", "forward-reverse")
    pOf <- function(tr, name) P(tr, method = "explicit", compile = FALSE, modelname = name, derivMode = .dm)
    g <- Y(c(y = "A * exp(-k * dose)"), f = NULL, parameters = c("A", "k", "dose"), compile = FALSE,
           modelname = "ss_g", derivMode = .dm)
    err <- Y(c(y = "s"), g, attachInput = FALSE, compile = FALSE, modelname = "ss_err", derivMode = .dm)
    ps <- pOf(insert(base(), k + A ~ 1 + (1 | gr(ID, dist = spikeSlab()))), "ss_slab")
    pdp <- pOf(insert(base(), k ~ 1 + (1 | gr(ID, dist = dirichletProcess())), A ~ 1 + (1 | gr(ID, dist = spikeSlab(prob = 0.3, tau = 0.5)))),
               "ss_dpslab")
    tr <- branch(define(NULL, "x <- x", x = c("a", "b")), table = tab)
    tr <- define(tr, "dose <- d", d = dose)
    gl <- Y(c(y = "a + b * dose"), f = NULL, parameters = c("a", "b", "dose"), compile = FALSE, modelname = "ss_glin",
            derivMode = .dm)
    pls <- pOf(insert(tr, a ~ 1 + (1 | ID), b ~ 1 + (1 | gr(ID, dist = spikeSlab()))), "ss_linslab")
    # seven observables y_j = A_j exp(-(k_j + k_(j+1)) dose) on a cycle, every parameter with a
    # spike-and-slab random effect: 14 per group
    J <- 7
    ks <- paste0("k", 1:J); As <- paste0("A", 1:J)
    trm <- branch(define(NULL, "x <- x", x = c(ks, As, "s")), table = tab[tab$ID != "L5", ])
    trm <- define(trm, "dose <- d", d = dose)
    trm <- insert(trm, "x <- exp(x)", x = c(ks, As, "s"))
    pm <- pOf(insert(trm, as.formula(paste(paste(c(ks, As), collapse = " + "), "~ 1 + (1 | gr(ID, dist = spikeSlab()))"))),
              "ss_many")
    gm <- Y(setNames(paste0(As, " * exp(-(", ks, " + ", c(ks[-1], ks[1]), ") * dose)"), paste0("y", 1:J)), f = NULL,
            parameters = c(ks, As, "dose"), compile = FALSE, modelname = "ss_gmany", derivMode = .dm)
    errm <- Y(setNames(rep("s", J), paste0("y", 1:J)), gm, attachInput = FALSE, compile = FALSE, modelname = "ss_errmany",
              derivMode = .dm)
    # ten observables y_j = A_j exp(-(k_j + k_(j+1)) dose) on a chain, k_11 = 0, 20 spike-and-slab random
    # effects per group of five, three conditions without a group; in the groups seven observables
    # with a noise far above the signal, so most effects are hardly informed
    J2 <- 10
    ks2 <- paste0("k", 1:J2); As2 <- paste0("A", 1:J2)
    groups2 <- paste0("C", 1:5)
    tab2 <- data.frame(ID = c(rep(groups2, each = 2), NA, NA, NA), dose = c(rep(c(0.5, 2), 5), 0.5, 1, 2),
                       row.names = c(paste0(rep(groups2, each = 2), c("a", "b")), "R1", "R2", "R3"))
    trc <- branch(define(NULL, "x <- x", x = c(ks2, As2)), table = tab2)
    trc <- define(trc, "dose <- d", d = dose)
    trc <- insert(trc, "x <- exp(x)", x = c(ks2, As2))
    pc <- pOf(insert(trc, as.formula(paste(paste(c(ks2, As2), collapse = " + "), "~ 1 + (1 | gr(ID, dist = spikeSlab()))"))),
              "ss_cmp")
    gc <- Y(setNames(paste0(As2, " * exp(-(", ks2, c(paste(" +", ks2[-1]), ""), ") * dose)"), paste0("y", 1:J2)), f = NULL,
            parameters = c(ks2, As2, "dose"), compile = FALSE, modelname = "ss_gcmp", derivMode = .dm)
    compile(g, err, ps, pdp, gl, pls, gm, errm, pm, gc, pc, output = "saemsparse_models", cores = test_cores())
    set.seed(7)
    dev <- matrix(0, 4, 2 * J, dimnames = list(groups[1:4], c(ks, As)))
    dev["L3", "k1"] <- 1; dev["L2", "k2"] <- -1; dev["L4", "A3"] <- 0.7
    rows <- lapply(rownames(tab)[tab$ID != "L5"], function(cn) {
      id <- tab[cn, "ID"]
      do.call(rbind, lapply(1:J, function(j) {
        k <- 0.6 * (exp(dev[id, ks[j]]) + exp(dev[id, ks[j %% J + 1]]))
        data.frame(name = paste0("y", j), time = 0:3, value = 10 * exp(dev[id, As[j]] - k * tab[cn, "dose"]) + rnorm(4, 0, 1),
                   sigma = NA_real_, condition = cn)
      }))
    })
    many <- list(g = gm, err = errm, p = pm, data = as.datalist(do.call(rbind, rows)), deviates = dev != 0,
                 init = c(setNames(rep(log(0.6), J), ks), setNames(rep(log(10), J), As), s = 0))
    set.seed(11)
    dev2 <- matrix(0, 5, 2 * J2, dimnames = list(groups2, c(ks2, As2)))
    dev2["C2", "k3"] <- 1; dev2["C4", "A2"] <- -0.8; dev2["C5", "k1"] <- -0.9; dev2["C1", "A1"] <- 0.7
    dev2["C3", "A3"] <- 0.6; dev2["C1", "k4"] <- 0.8; dev2["C3", "k2"] <- -0.7; dev2["C5", "A2"] <- 0.5
    A0 <- c(rep(10, 3), rep(1, 7)); sig <- c(rep(0.5, 3), rep(50, 7))
    rows <- lapply(rownames(tab2), function(cn) {
      id <- tab2[cn, "ID"]
      d <- if (is.na(id)) setNames(numeric(2 * J2), c(ks2, As2)) else dev2[id, ]
      do.call(rbind, lapply(1:J2, function(j) {
        k <- 0.6 * (exp(d[[ks2[j]]]) + if (j < J2) exp(d[[ks2[j + 1]]]) else 0)
        sj <- if (is.na(id)) 0.05 * A0[j] else sig[j]
        data.frame(name = paste0("y", j), time = 0:3, value = A0[j] * exp(d[[As2[j]]] - k * tab2[cn, "dose"]) +
                     rnorm(4, 0, sj), sigma = sj, condition = cn)
      }))
    })
    compounds <- list(g = gc, p = pc, data = as.datalist(do.call(rbind, rows)),
                      init = c(setNames(rep(log(0.6), J2), ks2), setNames(log(A0), As2), prob_ID = 0.2, tau_ID = 0.5))
    set.seed(4)
    dk <- c(L1 = 0, L2 = 0, L3 = 0.8, L4 = 0.8, L5 = -0.6)
    rows <- lapply(rownames(tab), function(cn) {
      k <- 0.6 * exp(dk[[tab[cn, "ID"]]])
      data.frame(name = "y", time = c(0, 1, 2, 3), value = 10 * exp(-k * tab[cn, "dose"]) + rnorm(4, 0, 0.1),
                 sigma = NA_real_, condition = cn)
    })
    data <- as.datalist(do.call(rbind, rows))
    ea <- setNames(rnorm(5, 0, 0.5), groups)
    eb <- c(L1 = 0, L2 = 0, L3 = 0.4, L4 = -0.3, L5 = 0)
    rows <- lapply(rownames(tab), function(cn) {
      id <- tab[cn, "ID"]
      data.frame(name = "y", time = c(0, 1, 2, 3),
                 value = 1 + ea[[id]] + (0.5 + eb[[id]]) * tab[cn, "dose"] + rnorm(4, 0, 0.3),
                 sigma = 0.3, condition = cn)
    })
    cache <<- list(g = g, err = err, ps = ps, pdp = pdp, gl = gl, pls = pls, tab = tab,
                   data = data, linear = as.datalist(do.call(rbind, rows)),
                   init = c(k = log(0.6), A = log(10), s = log(0.1)), many = many, compounds = compounds)
    cache
  }
})

# Per group of the linear model the exact posterior probability of the spike, and -2 log L: y of
# a group is normal with covariance sigma^2 I + omega^2 11', plus tau^2 dose dose' in the slab.
.linearSpike <- function(theta, prob, tau) {
  m <- .saemSparse()
  ids <- unique(m$tab$ID)
  terms <- vapply(ids, function(id) {
    cn <- rownames(m$tab)[m$tab$ID == id]
    y <- unlist(lapply(m$linear[cn], `[[`, "value"))
    dose <- rep(m$tab[cn, "dose"], each = 4)
    C <- diag(0.3^2, length(y)) + exp(theta[["omega_a_a"]])^2
    logN <- function(V) {
      r <- y - theta[["a"]] - theta[["b"]] * dose
      -0.5 * (sum(r * solve(V, r)) + as.numeric(determinant(2 * pi * V)$modulus))
    }
    c(log1p(-prob) + logN(C), log(prob) + logN(C + tau^2 * tcrossprod(dose)))
  }, numeric(2))
  top <- apply(terms, 2, max)
  logI <- top + log(colSums(exp(sweep(terms, 2, top))))
  list(p0 = exp(terms[1, ] - logI), value = -2 * sum(logI))
}

test_that("SAEM with spike-and-slab random effects: the indicators follow the exact posterior", {
  m <- .saemSparse()
  set.seed(5)
  theta <- c(a = 1.1, b = 0.4, omega_a_a = log(0.4))
  obj <- normL2(m$linear, m$gl * Xt() * m$pls)
  fit <- em(obj, randomEffects(m$pls), fixed = c(theta, prob_ID = 0.4, tau_ID = 0.5), method = "saem",
            control = list(saem = list(nEM = 600L)), verbose = FALSE)
  exact <- .linearSpike(theta, 0.4, 0.5)
  b <- fit$random[fit$random$eta == "b", ]
  expect_lt(max(abs(setNames(b$p0, b$group)[names(exact$p0)] - exact$p0)), 0.08)
  # with one sparse coordinate per group of a linear model the reported marginal is exact
  expect_equal(fit$value, exact$value, tolerance = 1e-6)
  expect_null(fit$valueSE)
})

test_that("the collapsed indicator move samples the exact posterior of the subsets of a Gaussian model", {
  set.seed(6)
  J <- 3; G <- 4
  H <- matrix(c(6, 2, 1, 2, 4, 1.5, 1, 1.5, 3), J)
  grad <- c(-3, 1, -0.5); ref <- c(0.2, 0, 0)
  prob <- c(0.3, 0.3, 0.3); tau <- c(0.8, 0.8, 0.8)
  quad <- function(v) drop(sum(grad * (v - ref)) + (v - ref) %*% H %*% (v - ref) / 2)
  dataAt <- function(V) if (is.list(V)) vapply(V, function(v) apply(v, 1, quad), numeric(G)) else apply(V, 1, quad)
  model <- list(ref = matrix(ref, G, J, byrow = TRUE), grad = matrix(grad, G, J, byrow = TRUE), H = rep(list(H), G))
  # the exact posterior of every subset with the slab values integrated out
  subsets <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), J)))
  logPost <- apply(subsets, 1, function(z) {
    s <- which(z)
    out <- sum(ifelse(z, log(prob), log1p(-prob)))
    if (!length(s)) return(out)
    Lam <- H[s, s, drop = FALSE] / 2 + diag(1 / tau[s]^2, length(s))
    h <- ((H %*% ref)[s] - grad[s]) / 2
    out - sum(log(tau[s])) - as.numeric(determinant(Lam)$modulus) / 2 + sum(h * solve(Lam, h)) / 2
  })
  exact <- exp(logPost - max(logPost)); exact <- exact / sum(exact)
  inclusion <- colSums(subsets * exact)
  V <- matrix(0, G, J); Z <- matrix(FALSE, G, J); value <- dataAt(V)
  gh <- dMod2:::.hermiteNodes(3L)
  n <- 1500
  visits <- matrix(0, G, nrow(subsets)); rb <- matrix(0, G, J)
  for (it in seq_len(n)) {
    for (r in sample.int(J)) {
      fl <- dMod2:::.slabCollapsed(V, Z, value, dataAt, model, prob, tau, rep(r, G), 1.5, gh)
      V <- fl$V; Z <- fl$Z; value <- fl$value
      rb[, r] <- rb[, r] + fl$pinc
    }
    key <- drop(Z %*% 2^(seq_len(J) - 1)) + 1
    visits[cbind(seq_len(G), key)] <- visits[cbind(seq_len(G), key)] + 1
  }
  expect_equal(value, dataAt(ifelse(Z, V, 0)))
  # subsets ordered as expand.grid, the first coordinate fastest: the binary key above
  expect_lt(max(abs(colSums(visits) / (G * n) - exact)), 0.03)
  # on a Gaussian model the Gauss-Hermite sum is exact: the Rao-Blackwellized inclusion has less noise
  expect_lt(max(abs(colMeans(rb) / n - inclusion)), 0.015)
})

# The estimates and the -2 log L of a SAEM fit and of the FOCEI fit of the same model.
# Without the refit of a normal random effect whose variance tends to 0.
.againstFocei <- function(p, nEM = 150L) {
  m <- .saemSparse()
  obj <- normL2(m$data, m$g * Xt() * p, errmodel = m$err)
  focei <- suppressMessages(em(obj, randomEffects(p), m$init, control = list(boundary = FALSE), verbose = FALSE))
  saem <- suppressMessages(em(obj, randomEffects(p), m$init, method = "saem",
                              control = list(boundary = FALSE, saem = list(nEM = nEM)), verbose = FALSE))
  list(focei = focei, saem = saem, diff = saem$argument[names(focei$argument)] - focei$argument,
       value = saem$value - focei$value, k = saem$random[saem$random$eta == "k", ],
       A = saem$random[saem$random$eta == "A", ])
}

test_that("SAEM selects the groups whose decay rate deviates, with the estimates of FOCEI", {
  m <- .saemSparse()
  set.seed(1)
  # spike and slab on k and A in one block
  r <- .againstFocei(m$ps)
  expect_identical(r$k$group[r$k$cluster > 0], c("L3", "L4", "L5"))
  expect_true(all(r$k$p0[r$k$group %in% c("L1", "L2")] > 0.9))
  expect_true(all(r$k$p0[r$k$group %in% c("L3", "L4", "L5")] < 0.05))
  expect_equal(r$k$mode, r$focei$random$mode[r$focei$random$eta == "k"], tolerance = 0.05)
  expect_lt(max(abs(r$diff[c("k", "A", "s")])), 0.05)
  expect_lt(max(abs(r$diff[c("prob_ID", "tau_ID")])), 0.1)
  expect_lt(abs(r$value), 0.3)
  expect_null(r$saem$valueSE)
  expect_true(all(r$A$cluster == 0 & r$A$mode == 0))
  # importance sampling at the SAEM estimate starts from the inner modes of the fit, so the
  # approximation it reports is the value of the fit
  expect_identical(dim(r$saem$innerModes), dim(r$saem$etaModes))
  l <- logLik(r$saem, "importance", draws = 500L)
  expect_equal(attr(l, "approximation"), -r$saem$value / 2, tolerance = 1e-6)
})

test_that("SAEM with Dirichlet-process and spike-and-slab random effects in one fit", {
  m <- .saemSparse()
  set.seed(2)
  r <- .againstFocei(m$pdp, nEM = 300L)
  cl <- setNames(r$k$cluster, r$k$group)
  expect_identical(cl[["L1"]], 0)
  expect_identical(cl[["L2"]], 0)
  expect_identical(cl[["L3"]], cl[["L4"]])
  expect_false(cl[["L5"]] %in% cl[c("L1", "L3")])
  # the slab of A has fixed hyperparameters: no group deviates in A, they are not identified
  expect_setequal(names(r$saem$argument), c("k", "A", "s", "concentration_ID_1", "tau_ID_1"))
  expect_true(all(r$A$cluster == 0))
  expect_lt(max(abs(r$diff[c("k", "A", "s", "tau_ID_1")])), 0.05)
  expect_lt(abs(r$value), 0.3)
  # at the FOCEI estimate -2 log L by importance sampling over the partitions against the exact sum
  set.seed(2)
  is <- suppressMessages(em(normL2(m$data, m$g * Xt() * m$pdp, errmodel = m$err), randomEffects(m$pdp),
                            fixed = r$focei$argument, method = "saem", verbose = FALSE,
                            control = list(saem = list(nBurnin = 5L, nEM = 5L, integral = "importance"))))
  expect_true(is.finite(is$valueSE) && is$valueSE < 0.2)
  expect_lt(abs(is$value - r$focei$value), 3 * is$valueSE + 0.02)
  expect_error(em(normL2(m$data, m$g * Xt() * m$ps, errmodel = m$err), randomEffects(m$ps), m$init,
                  method = "laplaceEM"), "\"focei\" or \"saem\"")
})

test_that("SAEM with many spike-and-slab random effects per group, most of them 0", {
  m <- .saemSparse()$many
  obj <- normL2(m$data, m$g * Xt() * m$p, errmodel = m$err)
  focei <- suppressMessages(em(obj, randomEffects(m$p), m$init, verbose = FALSE))
  set.seed(1)
  saem <- suppressMessages(em(obj, randomEffects(m$p), m$init, method = "saem", verbose = FALSE,
                              control = list(saem = list(nEM = 100L))))
  selected <- function(fit) setNames(fit$random$p0 < 0.5, fit$random$parameter)
  truth <- setNames(as.vector(m$deviates), paste0("eta_", rep(colnames(m$deviates), each = 4), "_",
                                                  rownames(m$deviates)))
  expect_identical(selected(saem)[names(truth)], truth)
  expect_identical(selected(focei)[names(truth)], truth)
  expect_lt(abs(saem$argument[["prob_ID"]] - focei$argument[["prob_ID"]]), 0.03)
  expect_lt(abs(saem$value - focei$value), 2)
  # every move type is tried and accepted at some rate
  tr <- saem$stageTrace[saem$stageTrace$phase == "averaging", ]
  rates <- colMeans(tr[c("acceptBirth", "acceptDeath", "acceptWalk")], na.rm = TRUE)
  expect_true(all(rates > 0 & rates < 1))
  # a FOCEI fit as init: the chain starts at the selected effects of the fit and keeps them
  set.seed(2)
  again <- suppressMessages(em(obj, randomEffects(m$p), focei, method = "saem", verbose = FALSE,
                               control = list(saem = list(nBurnin = 30L, nEM = 60L))))
  expect_identical(selected(again)[names(truth)], truth)
  expect_lt(abs(again$value - focei$value), 2)
})

test_that("SAEM with 20 spike-and-slab random effects per group, most of them hardly informed", {
  m <- .saemSparse()$compounds
  obj <- normL2(m$data, m$g * Xt() * m$p)
  set.seed(1)
  saem <- suppressMessages(em(obj, randomEffects(m$p), m$init, method = "saem", verbose = FALSE,
                              control = list(saem = list(nBurnin = 100L, nEM = 100L))))
  # the FOCEI fit of this fixture, 1 min: -2 log L, prob, the 9 effects with p0 below 1/2 and those
  # of them with p0 below 0.05
  focei <- list(value = 2980.078, prob = 0.4289,
                clear = c("eta_k1_C5", "eta_A1_C1", "eta_A2_C4", "eta_A2_C5", "eta_A3_C1", "eta_A3_C3"))
  selected <- saem$random$parameter[saem$random$p0 < 0.5]
  expect_true(all(focei$clear %in% selected))
  expect_lte(abs(length(selected) - 9), 3)
  expect_lt(abs(saem$value - focei$value), 2)
  # prob is weakly determined: -2 log L of FOCEI changes by about 0.2 over 0.06
  expect_lt(abs(saem$argument[["prob_ID"]] - focei$prob), 0.1)
  # the hardly informed effects keep prob where the informed ones put it, no run to every effect in the slab
  tr <- saem$stageTrace
  expect_lt(max(tr$included), 0.8 * nrow(saem$random))
})

test_that("msem runs SAEM fits of spike-and-slab random effects, the parts of a run join into one fit", {
  m <- .saemSparse()
  obj <- normL2(m$data, m$g * Xt() * m$ps, errmodel = m$err)
  run <- function(which = NULL) {
    set.seed(3)
    suppressMessages(msem(obj, randomEffects(m$ps), m$init, method = "saem", fits = 2L, sd = 0.2, which = which,
                          control = list(saem = list(nBurnin = 20L, nEM = 20L))))
  }
  whole <- run()
  expect_equal(whole$method, "saem")
  expect_true(all(is.finite(whole$multistart$value)))
  expect_true(all(c("cluster", "p0") %in% names(whole$random)))
  expect_setequal(whole$slab$name, c("prob_ID", "tau_ID"))
  joined <- dMod2:::.joinMsem(lapply(1:2, run))
  expect_equal(joined$multistart$value, whole$multistart$value)
  expect_equal(joined$argument, whole$argument)
})
