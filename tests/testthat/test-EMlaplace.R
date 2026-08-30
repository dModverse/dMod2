# ============================================================================
# EM / sparsify: marginal-likelihood L1/Laplace model selection.
#
# Sections:
#   * normalLaplace          - 1-D normal-Laplace marginal + moments vs a
#                              numerical integral (no compilation).
#   * .laplaceSubjectMarginal    - GN-exact per-subject FOCE marginal on a
#                              linear-Gaussian (quadratic) objfun (no compilation).
#   * EM + sparsify     - end-to-end on a tiny decay ODE, reference
#                              encoding, 2 cell lines: recover the individual /
#                              shared split and the parsimonious support.
#
# The math anchors live at the normal-Laplace section of R/nlmeLaplace.R.
# ============================================================================

## Context: "EM + sparsify (Laplace-NLME model selection)"  (context() is deprecated in testthat 3e; kept as a note)


# ---- normalLaplace primitive ---------------------------------------------

test_that("normalLaplace marginal + moments match a numerical integral", {
  # I(a,m,lambda) = int (lambda/2) e^{-lambda|eta|} e^{-(a/2)(eta-m)^2} deta,
  # its posterior moments E[eta], E|eta|, E[(eta-m)^2].
  cases <- list(c(a = 2.0, m = 0.7, lambda = 1.5),
                c(a = 8.0, m = -1.2, lambda = 0.4),
                c(a = 0.5, m = 3.0, lambda = 3.0))
  for (cs in cases) {
    a <- cs[["a"]]; m <- cs[["m"]]; lambda <- cs[["lambda"]]
    nl <- normalLaplace(a, m, lambda)

    kern  <- function(eta) (lambda / 2) * exp(-lambda * abs(eta)) *
      exp(-(a / 2) * (eta - m)^2)
    I     <- integrate(kern, -Inf, Inf, rel.tol = 1e-10)$value
    Eeta  <- integrate(function(e) e * kern(e), -Inf, Inf,
                       rel.tol = 1e-10)$value / I
    Eabs  <- integrate(function(e) abs(e) * kern(e), -Inf, Inf,
                       rel.tol = 1e-10)$value / I
    Ecen2 <- integrate(function(e) (e - m)^2 * kern(e), -Inf, Inf,
                       rel.tol = 1e-10)$value / I

    expect_equal(nl$logI,  log(I), tolerance = 1e-8)
    expect_equal(nl$Eeta,  Eeta,   tolerance = 1e-7)
    expect_equal(nl$Eabs,  Eabs,   tolerance = 1e-7)
    expect_equal(nl$Ecen2, Ecen2,  tolerance = 1e-7)
  }
})

test_that("normalLaplace recycles its arguments to a common length", {
  nl <- normalLaplace(a = c(1, 2, 4), m = 0.5, lambda = 1)
  expect_length(nl$logI, 3L)
  expect_length(nl$Eabs, 3L)
})


# ---- .laplaceSubjectMarginal on a linear-Gaussian objfun ---------------------

test_that(".laplaceSubjectMarginal is GN-exact on a quadratic objfun", {
  # A quadratic data objfun in -2logL convention:
  #   value(eta) = 0.5 Hd (eta - m0)^2 + v0,  grad = Hd (eta - m0),  hess = Hd.
  # Then the FOCE marginal -2logL_i = v0 - 2 log I with
  #   I = int (lambda/2) e^{-lambda|z|} e^{-(Hd/4)(z - (m0 - target))^2} dz,
  #   z = eta - target  (a = Hd/2, m_shift = m0 - target).
  Hd <- 3.0; m0 <- 1.5; v0 <- 4.0; target <- 0.2; lambda <- 2.0
  nm <- "eta_x"
  objfun <- function(p, ...) {
    eta <- unname(p[nm])
    list(value    = 0.5 * Hd * (eta - m0)^2 + v0,
         gradient = stats::setNames(Hd * (eta - m0), nm),
         hessian  = matrix(Hd, 1, 1, dimnames = list(nm, nm)))
  }
  sm <- .laplaceSubjectMarginal(objfun,
                            eta0    = stats::setNames(0, nm),
                            targets = stats::setNames(target, nm),
                            lambda  = lambda)

  kern <- function(z) (lambda / 2) * exp(-lambda * abs(z)) *
    exp(-(Hd / 4) * (z - (m0 - target))^2)
  I <- integrate(kern, -Inf, Inf, rel.tol = 1e-10)$value

  expect_equal(sm$value, v0 - 2 * log(I), tolerance = 1e-6)
  expect_true(sm$converged)
  # closed-form data-alone mode is m0 regardless of the penalised mode
  expect_equal(unname(sm$m), m0, tolerance = 1e-6)
})


# ---- end-to-end EM + sparsify on a tiny decay ODE -------------------

# One-state decay dA/dt = -k A, observed as log(A+1); 2 cell lines A/B with the
# reference encoding (line A = baseline, line B carries the deviations). Truth:
# k is individual (line B differs), A0 is shared.
.build_reg_smoke <- function(tag = "reg", seed = 1L) {
  set.seed(seed)
  reactions <- addReaction(NULL, from = "A", to = "", rate = "k * A",
                           description = "decay")
  m <- odemodel(reactions, modelname = paste0("regsmoke_", tag),
                backend = "cppDE", compile = TRUE, deriv2 = TRUE)
  x <- Xs(m, compile = TRUE)
  g <- Y(eqnvec(y = "log(A + 1)"), x, modelname = paste0("regsmoke_obs_", tag),
         compile = TRUE, deriv2 = TRUE, attach.input = FALSE)
  e <- Y(eqnvec(y = "sigma"), g, modelname = paste0("regsmoke_err_", tag),
         compile = TRUE, deriv2 = TRUE, attach.input = FALSE)
  trafo <- eqnvec(k = "exp(log_k + eta_k)", A = "exp(log_A0 + eta_A0)",
                  sigma = "exp(log_sigma)")
  lines <- c("A", "B")
  cl_table <- data.frame(eta_k  = c("0", "eta_k_B"),
                         eta_A0 = c("0", "eta_A0_B"),
                         row.names = lines, stringsAsFactors = FALSE)
  p <- P(branch(trafo, table = cl_table, apply = "insert"), method = "explicit",
         modelname = paste0("regsmoke_p_", tag), compile = TRUE, deriv2 = TRUE)
  prd <- g * x * p

  st     <- c(log_k = log(0.3), log_A0 = log(10), log_sigma = log(0.05))
  d_true <- c(eta_k_B = 0.8, eta_A0_B = 0)
  times  <- c(0.5, 1, 2, 3, 4, 6, 8, 10)
  sim <- prd(times, c(st, d_true), deriv = FALSE)
  dl_rows <- do.call(rbind, lapply(lines, function(s) {
    pr <- sim[[s]]
    data.frame(name = "y", time = pr[, "time"],
               value = pr[, "y"] + rnorm(nrow(pr), 0, 0.05),
               sigma = NA_real_, condition = s, stringsAsFactors = FALSE)
  }))
  dlist <- as.datalist(dl_rows)
  list(prd = prd, e = e, dlist = dlist, st = st)
}

test_that("EM recovers the individual/shared split on a decay ODE", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  s <- .build_reg_smoke(tag = "fit")

  pen    <- penaltyL1(c("k", "A0"), subjects = "B")
  obj    <- normL2(s$dlist, s$prd, errmodel = s$e) +
    constraintL1(pen)
  center <- emInit(c(log_k = s$st[["log_k"]], log_A0 = s$st[["log_A0"]],
                      log_sigma = log(0.1)), pen, lambda = 1)

  fit <- suppressWarnings(suppressMessages(EM(obj, center, method = "focei",
                                 control = list(cm1 = list(iterlim = 100L)),
                                 verbose = FALSE)))

  expect_s3_class(fit, "em")
  expect_true(fit$converged)
  expect_true(is.finite(fit$value))
  expect_true(fit$lambda > 0)
  # k deviates (individual), A0 does not (shared)
  expect_true(abs(fit$etaModes["B", "eta_k"])  > 0.3)
  expect_true(abs(fit$etaModes["B", "eta_A0"]) < 0.1)
})

test_that("the FOCEI correction closes the error-parameter gradient gap", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  s <- .build_reg_smoke(tag = "grad")

  pen   <- penaltyL1(c("k", "A0"), subjects = "B")
  obj   <- normL2(s$dlist, s$prd, errmodel = s$e) +
    constraintL1(pen)
  rec   <- .laplaceReconstruct(obj)
  conds <- names(rec$data)
  resObjList <- stats::setNames(lapply(conds, function(cc) {
    d <- rec$data[cc]; class(d) <- class(rec$data)
    normL2(d, rec$prdfn, errmodel = rec$errfn)
  }), conds)
  tgt   <- .laplaceFactorisingTargets(pen)
  theta <- c(log_k = s$st[["log_k"]] + 0.15, log_A0 = s$st[["log_A0"]] - 0.1,
             log_sigma = log(0.08), lambda = 1.5)

  ctrl <- list(iterlim = 300L)
  fval <- function(th) .laplaceOuterEval(th, rec, resObjList, tgt, ctrl)
  gF <- .laplaceOuterEval(theta, rec, resObjList, tgt, ctrl,
                      grad = TRUE, correction = FALSE)$gradient
  gI <- .laplaceOuterEval(theta, rec, resObjList, tgt, ctrl,
                      grad = TRUE, correction = TRUE)$gradient
  # moving-mode central FD of the marginal value = the ground-truth gradient
  gFD <- vapply(names(theta), function(q) {
    h  <- 1e-6 * max(abs(theta[[q]]), 1)
    tp <- theta; tp[[q]] <- tp[[q]] + h
    tm <- theta; tm[[q]] <- tm[[q]] - h
    (fval(tp) - fval(tm)) / (2 * h)
  }, numeric(1))
  relF <- abs(gF[names(theta)] - gFD) / pmax(abs(gFD), 1e-8)
  relI <- abs(gI[names(theta)] - gFD) / pmax(abs(gFD), 1e-8)

  # FOCE misses the interaction (sigma, several %) AND the structural + lambda
  # mode shift (a few tenths of a %); the full FOCEI total derivative closes all
  # of them to FD accuracy.
  expect_gt(relF[["log_sigma"]], 0.02)
  expect_gt(max(relF[c("log_k", "log_A0", "lambda")]), 0.002)
  expect_lt(max(relI[c("log_k", "log_A0", "log_sigma")]), 1e-3)
  expect_lt(relI[["lambda"]], 5e-3)
})

test_that("EM(method='saem') agrees with the FOCEI backend", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  s <- .build_reg_smoke(tag = "saem")

  pen    <- penaltyL1(c("k", "A0"), subjects = "B")
  obj    <- normL2(s$dlist, s$prd, errmodel = s$e) +
    constraintL1(pen)
  center <- emInit(c(log_k = s$st[["log_k"]], log_A0 = s$st[["log_A0"]],
                      log_sigma = log(0.1)), pen, lambda = 1)

  f_focei <- suppressMessages(EM(obj, center, method = "focei",
                                     control = list(cm1 = list(iterlim = 50L)),
                                     verbose = FALSE))
  set.seed(4)
  f_saem <- suppressMessages(EM(
    obj, center, method = "saem",
    control = list(saem = list(nBurnin = 120L, nEM = 120L, nMcmc = 12L)),
    verbose = FALSE))

  expect_s3_class(f_saem, "em")
  expect_equal(f_saem$method, "saem")
  expect_true(f_saem$converged)
  # the two independent backends land at the same structural estimate + lambda
  sn <- c("log_k", "log_A0", "log_sigma")
  expect_lt(max(abs(f_focei$argument[sn] - f_saem$argument[sn])), 0.2)
  expect_gt(f_saem$lambda / f_focei$lambda, 0.6)
  expect_lt(f_saem$lambda / f_focei$lambda, 1.6)
  # both recover the individual k deviation
  expect_gt(f_saem$etaModes["B", "eta_k"], 0.4)
  expect_lt(abs(f_saem$etaModes["B", "eta_A0"]), 0.2)
})

test_that("sparsify picks the parsimonious support {k}", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  s <- .build_reg_smoke(tag = "sel")

  pen    <- penaltyL1(c("k", "A0"), subjects = "B")
  obj    <- normL2(s$dlist, s$prd, errmodel = s$e) +
    constraintL1(pen)
  center <- emInit(c(log_k = s$st[["log_k"]], log_A0 = s$st[["log_A0"]],
                      log_sigma = log(0.1)), pen, lambda = 1)

  set.seed(2)
  sel <- suppressMessages(sparsify(obj, center, method = "focei",
                                    control = list(cm1 = list(iterlim = 50L)),
                                    fits = 6, cores = 1, sd = 0.5,
                                    verbose = FALSE))

  expect_s3_class(sel, "sparsify")
  expect_identical(sel$support, "k")
  # k ranks above A0
  expect_gt(sel$ranking[["k"]], sel$ranking[["A0"]])
  # the selected support is the marginal-minimiser and beats the full support
  expect_true(sel$chain$selected[sel$chain$support == "k"])
  m_k   <- sel$chain$m2ll_marg[sel$chain$support == "k"]
  m_all <- sel$chain$m2ll_marg[sel$chain$support == "k,A0"]
  expect_lt(m_k, m_all)
})


# ============================================================================
# Clustered inner solver (complete-graph fusion): solveFusedComplete (sort +
# weighted PAVA), .clusterSolve (joint prox-linear driver), and the runnable
# MAP-EM fit that recovers the grouping. Math anchors: R/L1Clustering.R.
# ============================================================================

## Context: "EM clustered inner solver (complete-graph fusion)"  (context() is deprecated in testthat 3e; kept as a note)

# ---- solveFusedComplete vs brute force (no compilation) ------------------

test_that("solveFusedComplete matches exhaustive brute force + limits", {
  fobj  <- function(u, m, w, lambda) 0.5 * sum(w * (u - m)^2) + lambda * sum(dist(u))
  ## exact reference: min objective over ordered contiguous partitions
  brute <- function(m, w, lambda) {
    n <- length(m); ord <- order(m); ms <- m[ord]; ws <- w[ord]
    v <- ms - lambda * (2 * seq_len(n) - n - 1) / ws
    best <- NULL
    for (mask in 0:(2^(n - 1) - 1)) {
      br <- if (n > 1) which(as.logical(intToBits(mask)[seq_len(n - 1)])) else integer(0)
      starts <- c(1, br + 1); ends <- c(br, n)
      us <- numeric(n); prev <- -Inf; ok <- TRUE
      for (gg in seq_along(starts)) {
        idx <- starts[gg]:ends[gg]; wg <- sum(ws[idx] * v[idx]) / sum(ws[idx])
        us[idx] <- wg; if (wg < prev - 1e-9) ok <- FALSE; prev <- wg
      }
      if (!ok) next
      u <- numeric(n); u[ord] <- us; o <- fobj(u, m, w, lambda)
      if (is.null(best) || o < best - 1e-12) { best <- o; bu <- u }
    }
    list(o = best, u = bu)
  }
  set.seed(11); maxdiff <- 0
  for (t in 1:120) {
    n <- sample(3:6, 1)
    m <- round(rnorm(n, 0, 2), 3); w <- round(runif(n, 0.3, 4), 3)
    lambda <- sample(c(0, 0.1, 0.5, 2, 10), 1)
    sol <- solveFusedComplete(setNames(m, paste0("s", seq_len(n))), w, lambda)
    bf  <- brute(m, w, lambda)
    maxdiff <- max(maxdiff, max(abs(sol$u - bf$u)))
  }
  expect_lt(maxdiff, 1e-7)
  ## limits
  m <- c(a = 0, b = 1, c = 5, d = 6); w <- c(1, 2, 1, 3)
  s0 <- solveFusedComplete(m, w, 1e-9)
  expect_equal(unname(s0$u), unname(m), tolerance = 1e-6)  # lambda ~ 0 -> singletons
  expect_equal(s0$G, 4L)
  sI <- solveFusedComplete(m, w, 1e7)
  expect_equal(sI$G, 1L)                                    # lambda large -> one cluster
  expect_equal(unname(sI$u[1]), sum(w * m) / sum(w), tolerance = 1e-5)
})

# ---- .clusterSolve + MAP-EM on a compiled ODE -------------------------

# 4 cell lines A/B/C/D, one clustered parameter k; truth two clusters
# {A,B} (eta 0.35) and {C,D} (eta -0.35). Reuses the decay ODE.
.build_reg_cluster <- function(tag = "cl", seed = 3L) {
  set.seed(seed)
  reactions <- addReaction(NULL, from = "A", to = "", rate = "k * A",
                           description = "decay")
  m <- odemodel(reactions, modelname = paste0("regcl_", tag),
                backend = "cppDE", compile = TRUE, deriv2 = TRUE)
  x <- Xs(m, compile = TRUE)
  g <- Y(eqnvec(y = "log(A + 1)"), x, modelname = paste0("regcl_obs_", tag),
         compile = TRUE, deriv2 = TRUE, attach.input = FALSE)
  e <- Y(eqnvec(y = "sigma"), g, modelname = paste0("regcl_err_", tag),
         compile = TRUE, deriv2 = TRUE, attach.input = FALSE)
  lines <- c("A", "B", "C", "D")
  trafo <- eqnvec(k = "exp(log_k + eta_k)", A = "exp(log_A0)",
                  sigma = "exp(log_sigma)")
  cl_table <- data.frame(eta_k = paste0("eta_k_", lines), row.names = lines,
                         stringsAsFactors = FALSE)
  p <- P(branch(trafo, table = cl_table, apply = "insert"), method = "explicit",
         modelname = paste0("regcl_p_", tag), compile = TRUE, deriv2 = TRUE)
  prd <- g * x * p
  st     <- c(log_k = log(0.3), log_A0 = log(10), log_sigma = log(0.05))
  d_true <- c(eta_k_A = 0.35, eta_k_B = 0.35, eta_k_C = -0.35, eta_k_D = -0.35)
  times  <- c(0.5, 1, 2, 3, 4, 6, 8, 10)
  sim <- prd(times, c(st, d_true), deriv = FALSE)
  dl <- do.call(rbind, lapply(lines, function(s) {
    pr <- sim[[s]]
    data.frame(name = "y", time = pr[, "time"],
               value = pr[, "y"] + rnorm(nrow(pr), 0, 0.05),
               sigma = NA_real_, condition = s, stringsAsFactors = FALSE)
  }))
  list(prd = prd, e = e, dlist = as.datalist(dl), st = st, lines = lines)
}

test_that(".clusterSolve finds the joint clustered mode (vs optim)", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  s <- .build_reg_cluster(tag = "solve")
  pen <- penaltyL1("k", method = "clustered", subjects = s$lines)
  obj <- normL2(s$dlist, s$prd, errmodel = s$e) +
    constraintL1(pen)
  rec <- .laplaceReconstruct(obj); conds <- names(rec$data)
  resObjList <- stats::setNames(lapply(conds, function(cc) {
    d <- rec$data[cc]; class(d) <- class(rec$data)
    normL2(d, rec$prdfn, errmodel = rec$errfn)
  }), conds)
  all_eta <- as.vector(pen$subjectEtas)
  Fobj <- function(eta4, lam) {
    ev <- stats::setNames(rep(0, length(all_eta)), all_eta)
    ev[paste0("eta_k_", s$lines)] <- eta4
    dv <- sum(vapply(s$lines, function(ss)
      resObjList[[ss]](c(s$st, ev), deriv = FALSE)$value, 0.0))
    dv + lam * sum(dist(eta4))
  }
  # away from the fusion boundary the exact solver matches Nelder-Mead closely
  # (near/at fusion the general optimiser is the looser reference, not the solver)
  for (lam in c(0, 50, 200)) {
    sol <- .clusterSolve(resObjList, s$st, pen, lam, maxit = 150L, tol = 1e-11)
    eh  <- sol$etahat[s$lines, "eta_k"]
    opt <- optim(rep(0, 4), Fobj, lam = lam, method = "Nelder-Mead",
                 control = list(reltol = 1e-13, maxit = 15000))
    expect_lt(max(abs(eh - opt$par)), 5e-3)
  }
  ## intermediate lambda -> the correct two clusters {A,B},{C,D} emerge
  sol <- .clusterSolve(resObjList, s$st, pen, 150, maxit = 150L, tol = 1e-11)
  expect_equal(sol$G, 2L)
  grp <- lapply(sol$clusters[["eta_k"]], sort)
  has <- function(gg) any(vapply(grp, function(z) identical(z, gg), NA))
  expect_true(has(c("A", "B")) && has(c("C", "D")))
})

test_that("EM clustered (MAP-EM) recovers the grouping at fixed sigma", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  s <- .build_reg_cluster(tag = "map")
  pen <- penaltyL1("k", method = "clustered", subjects = s$lines)
  obj <- normL2(s$dlist, s$prd, errmodel = s$e) +
    constraintL1(pen)
  fixed <- c(log_sigma = log(0.05))                  # fixed sigma breaks the feedback
  init  <- emInit(c(log_k = log(0.3), log_A0 = log(10)), pen, lambda = 300)
  fit <- suppressMessages(EM(obj, init, fixed = fixed, method = "focei",
                                 control = list(cm1 = list(iterlim = 50L),
                                                maxOuter = 40L), verbose = FALSE))
  expect_s3_class(fit, "em")
  expect_true(isTRUE(fit$clustered))
  expect_true(fit$converged)
  expect_equal(fit$G, 2L)                            # {A,B} and {C,D}
  # anchoring: deviations mean-zero per parameter
  expect_lt(abs(mean(fit$etaModes[, "eta_k"])), 1e-8)
  # the recovered clusters are {A,B} and {C,D}
  grp <- lapply(fit$clusters[["eta_k"]], sort)
  has <- function(g) any(vapply(grp, function(x) identical(x, g), NA))
  expect_true(has(c("A", "B")))
  expect_true(has(c("C", "D")))
})


# ---- exact clustered marginal (order-cone / AGHQ) ------------------------

test_that(".clusterParamMarginal is exact vs a direct integral + Fisher moments", {
  gh <- .gaussHermite(24)
  directMarg <- function(H, m, lambda) {
    n <- length(H); B <- .anchorBasis(n)
    s <- seq(-9, 9, length.out = if (n == 3) 251 else 61); dz <- s[2] - s[1]
    grid <- as.matrix(expand.grid(rep(list(s), n - 1L)))
    Dp <- function(e) sum(dist(e)); etas <- B %*% t(grid)
    ph  <- 0.25 * colSums(H * (etas - m)^2) + 0.5 * lambda * apply(etas, 2, Dp)
    phZ <- 0.5 * lambda * apply(etas, 2, Dp)
    mI <- max(-ph);  I <- log(sum(exp(-ph  - mI))) + mI + (n - 1) * log(dz)
    mZ <- max(-phZ); Z <- log(sum(exp(-phZ - mZ))) + mZ + (n - 1) * log(dz)
    -2 * I + 2 * Z
  }
  set.seed(5)
  for (n in c(3, 4)) {
    H <- runif(n, 1, 5); m <- rnorm(n, 0, 1); lambda <- 1.5
    pm <- .clusterParamMarginal(H, m, lambda, gh)
    expect_lt(abs(pm$value - directMarg(H, m, lambda)), 0.1)   # exact corner (quad accuracy)
  }
  # Fisher: d(-2logL)/dm_i = E_post[-H_i(eta_i - m_i)] = -H_i (Eeta_i - m_i). Hold the
  # quadrature centre fixed so the FD sees only the integrand's m-dependence.
  H <- c(3, 2, 4); m <- c(0.6, -0.2, 0.9); lambda <- 2
  B  <- MASS::Null(matrix(1, 3, 1)); Az <- crossprod(B, H * B)
  ctr <- as.numeric(solve(Az, crossprod(B, H * m)))
  pm <- .clusterParamMarginal(H, m, lambda, gh, center = ctr)
  grad_m <- -H * (pm$Eeta - m)
  fd <- vapply(seq_len(3), function(i) { h <- 1e-4
    mp <- m; mp[i] <- mp[i] + h; mm <- m; mm[i] <- mm[i] - h
    (.clusterParamMarginal(H, mp, lambda, gh, center = ctr)$value -
       .clusterParamMarginal(H, mm, lambda, gh, center = ctr)$value) / (2 * h) }, 0.0)
  expect_lt(max(abs(grad_m - fd)), 1e-2)
})

test_that("the exact-marginal comparison recovers the true grouping", {
  # linear-Gaussian n=4, precise data, true clusters {1,2},{3,4}. The exact
  # marginal must reject full fusion (data cost) AND full separation (Occam).
  gh <- .gaussHermite(24)
  Hd <- c(50, 45, 52, 48); m <- c(0.55, 0.52, -0.48, -0.51)
  lamgrid <- 10^seq(-3, 3, length.out = 40)
  paths <- unique(lapply(lamgrid, function(lam)
    lapply(solveFusedComplete(stats::setNames(m, as.character(1:4)), Hd, lam)$clusters,
           function(g) sort(as.integer(g)))))
  scores <- vapply(paths, function(P) .clusterGroupScore(P, Hd, m, gh)$value, 0.0)
  best <- lapply(paths[[which.min(scores)]], sort)
  has <- function(g) any(vapply(best, function(z) identical(z, g), NA))
  expect_true(has(c(1L, 2L)) && has(c(3L, 4L)))   # picks {1,2},{3,4}
  expect_equal(length(best), 2L)
})


test_that("sparsify recovers the grouping on a 4-subject ODE", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  s <- .build_reg_cluster(tag = "sel")
  pen  <- penaltyL1("k", method = "clustered", subjects = s$lines)
  obj  <- normL2(s$dlist, s$prd, errmodel = s$e) +
    constraintL1(pen)
  init <- emInit(c(log_k = log(0.3), log_A0 = log(10), log_sigma = log(0.1)),
                  pen, lambda = 1)
  res <- suppressMessages(sparsify(
    obj, init, control = list(cm1 = list(iterlim = 50L), maxOuter = 30L),
    verbose = FALSE))
  expect_s3_class(res, "sparsify")
  expect_equal(res$perParam[["eta_k"]]$G, 2L)          # two clusters
  grp <- lapply(res$perParam[["eta_k"]]$clusters, sort)
  has <- function(gg) any(vapply(grp, function(z) identical(z, gg), NA))
  expect_true(has(c("A", "B")) && has(c("C", "D")))    # {A,B},{C,D}
})


test_that("the joint MH sampler matches the AGHQ clustered marginal moments", {
  # SAEM/MCMC exact cross-check: the coupled sampler and the quadrature target the
  # same posterior, so E[eta] and E||D eta||_1 must agree (to Monte-Carlo error).
  set.seed(7)
  gh <- .gaussHermite(24)
  H <- c(5, 4, 6, 5); m <- c(0.5, 0.4, -0.45, -0.55); lambda <- 2
  aq <- .clusterParamMarginal(H, m, lambda, gh)
  mh <- .clusterMHmoments(H, m, lambda, nsamp = 40000L, burn = 5000L)
  expect_lt(max(abs(aq$Eeta - mh$Eeta)), 0.05)
  expect_lt(abs(aq$Efus - mh$Efus), 0.20)   # sum over pairs -> higher MC variance
  # weighted reduced (grouped) model too
  aq2 <- .clusterParamMarginal(c(16, 17), c(0.55, -0.52), lambda, gh, sizes = c(2, 2))
  mh2 <- .clusterMHmoments(c(16, 17), c(0.55, -0.52), lambda, sizes = c(2, 2),
                           nsamp = 40000L, burn = 5000L)
  expect_lt(max(abs(aq2$Eeta - mh2$Eeta)), 0.05)
})


test_that("EM(method='saem') fits the clustered model on the true ODE posterior", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  s <- .build_reg_cluster(tag = "csaem")
  pen  <- penaltyL1("k", method = "clustered", subjects = s$lines)
  obj  <- normL2(s$dlist, s$prd, errmodel = s$e) +
    constraintL1(pen)
  init <- emInit(c(log_k = log(0.3), log_A0 = log(10), log_sigma = log(0.1)),
                  pen, lambda = 5)
  set.seed(9)
  fit <- suppressMessages(EM(obj, init, method = "saem",
    control = list(saem = list(nBurnin = 100L, nEM = 100L, nMcmc = 8L)),
    verbose = FALSE))
  expect_s3_class(fit, "em")
  expect_equal(fit$method, "saem")
  # SAEM samples the TRUE (nonlinear) posterior -> cross-checks the FOCE path
  expect_lt(abs(fit$argument[["log_k"]] - log(0.3)), 0.1)
  # posterior means reveal the two groups: A,B > 0, C,D < 0, well separated
  E <- fit$etaModes[, "eta_k"]
  expect_gt(min(E[c("A", "B")]), 0.15)
  expect_lt(max(E[c("C", "D")]), -0.15)
})


# ---- WS2: per-parameter lambda_k + many-lines ----------------------------

test_that("penalty machinery supports per-parameter and merged strengths", {
  # global (backward compatible): one shared strength name
  pg <- penaltyL1(c("k", "A0"), subjects = paste0("s", 1:3))
  expect_equal(pg$lambdaName, "lambda")
  expect_equal(unname(pg$lambdaByParam[c("k", "A0")]), c("lambda", "lambda"))

  # per-parameter: one strength name each
  pp <- penaltyL1(c("k", "A0"), subjects = paste0("s", 1:3), perParam = TRUE)
  expect_setequal(pp$lambdaName, c("lambda_k", "lambda_A0"))
  expect_equal(pp$lambdaByParam[["k"]], "lambda_k")

  # merged specs keep their distinct (per-parameter) strength names
  m <- penaltyL1("k", subjects = "s1", perParam = TRUE) +
       penaltyL1("A0", subjects = "s1", perParam = TRUE)
  expect_setequal(m$lambdaName, c("lambda_k", "lambda_A0"))

  # clustered rejects per-parameter strengths
  expect_error(penaltyL1("k", method = "clustered", perParam = TRUE,
                            subjects = c("a", "b", "c")), "single global")

  # constraintL1 accumulates each term's magnitude into its own strength slot
  o    <- constraintL1(pp)
  etas <- setNames(c(0.5, -0.3, 0.1, 0, 0.2, 0.05), as.vector(pp$subjectEtas))
  allp <- c(etas, lambda_k = 2, lambda_A0 = 5)
  r    <- o(allp)
  expect_equal(r$value, 2 * (0.5 + 0.3 + 0.1) + 5 * (0 + 0.2 + 0.05))
  expect_equal(unname(r$gradient["lambda_k"]),  0.9)
  expect_equal(unname(r$gradient["lambda_A0"]), 0.25)
})

# 4-subject fixture: parameter k is individual (two groups +-0.6), A0 is shared.
.build_reg_manyline <- function(tag = "ml", seed = 3L) {
  set.seed(seed)
  reactions <- addReaction(NULL, from = "A", to = "", rate = "k * A",
                           description = "decay")
  m <- odemodel(reactions, modelname = paste0("regml_", tag),
                backend = "cppDE", compile = TRUE, deriv2 = TRUE)
  x <- Xs(m, compile = TRUE)
  g <- Y(eqnvec(y = "log(A + 1)"), x, modelname = paste0("regml_obs_", tag),
         compile = TRUE, deriv2 = TRUE, attach.input = FALSE)
  e <- Y(eqnvec(y = "sigma"), g, modelname = paste0("regml_err_", tag),
         compile = TRUE, deriv2 = TRUE, attach.input = FALSE)
  subj  <- paste0("s", 1:4)
  trafo <- eqnvec(k = "exp(log_k + eta_k)", A = "exp(log_A0 + eta_A0)",
                  sigma = "exp(log_sigma)")
  tab <- data.frame(eta_k  = paste0("eta_k_",  subj),
                    eta_A0 = paste0("eta_A0_", subj),
                    row.names = subj, stringsAsFactors = FALSE)
  p <- P(branch(trafo, table = tab, apply = "insert"), method = "explicit",
         modelname = paste0("regml_p_", tag), compile = TRUE, deriv2 = TRUE)
  prd <- g * x * p
  st  <- c(log_k = log(0.3), log_A0 = log(10), log_sigma = log(0.05))
  dev <- c(s1 = 0.6, s2 = -0.6, s3 = 0.55, s4 = -0.55)
  times <- c(0.5, 1, 2, 3, 4, 6, 8, 10)
  pars  <- c(st, setNames(dev[subj], paste0("eta_k_", subj)),
             setNames(rep(0, 4), paste0("eta_A0_", subj)))
  sim <- prd(times, pars, deriv = FALSE)
  dl_rows <- do.call(rbind, lapply(subj, function(s) {
    pr <- sim[[s]]
    data.frame(name = "y", time = pr[, "time"],
               value = pr[, "y"] + rnorm(nrow(pr), 0, 0.05),
               sigma = NA_real_, condition = s, stringsAsFactors = FALSE)
  }))
  list(prd = prd, e = e, dlist = as.datalist(dl_rows), st = st, subj = subj)
}

test_that("per-parameter FOCEI gradient matches the moving-mode FD", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  s <- .build_reg_manyline(tag = "grad")
  pen <- penaltyL1(c("k", "A0"), subjects = s$subj, perParam = TRUE)
  obj <- normL2(s$dlist, s$prd, errmodel = s$e) +
    constraintL1(pen)
  rec <- .laplaceReconstruct(obj); conds <- names(rec$data)
  resObjList <- setNames(lapply(conds, function(cc) {
    d <- rec$data[cc]; class(d) <- class(rec$data)
    normL2(d, rec$prdfn, errmodel = rec$errfn)
  }), conds)
  tgt   <- .laplaceFactorisingTargets(pen)
  theta <- c(log_k = s$st[["log_k"]] + 0.1, log_A0 = s$st[["log_A0"]] - 0.05,
             log_sigma = log(0.06), lambda_k = 2, lambda_A0 = 6)
  ctrl  <- list(iterlim = 300L)
  fv  <- function(t) .laplaceOuterEval(t, rec, resObjList, tgt, ctrl)
  gFD <- vapply(names(theta), function(q) {
    h  <- 1e-6 * max(abs(theta[[q]]), 1)
    tp <- theta; tp[[q]] <- tp[[q]] + h
    tm <- theta; tm[[q]] <- tm[[q]] - h
    (fv(tp) - fv(tm)) / (2 * h)
  }, numeric(1))
  gI <- .laplaceOuterEval(theta, rec, resObjList, tgt, ctrl,
                      grad = TRUE, correction = TRUE)$gradient
  # compare on an absolute scale (some coordinates have ~0 gradient here)
  expect_lt(max(abs(gI[names(theta)] - gFD)) / max(abs(gFD)), 1e-3)
})

test_that("per-parameter EM gives soft sparsity (shared par -> capped lambda)", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  s <- .build_reg_manyline(tag = "fit")
  pen <- penaltyL1(c("k", "A0"), subjects = s$subj, perParam = TRUE)
  obj <- normL2(s$dlist, s$prd, errmodel = s$e) +
    constraintL1(pen)
  center <- emInit(c(log_k = s$st[["log_k"]], log_A0 = s$st[["log_A0"]],
                      log_sigma = log(0.05)), pen, lambda = 1)
  fit <- suppressWarnings(suppressMessages(EM(obj, center, method = "focei",
                                 control = list(cm1 = list(iterlim = 100L)),
                                 verbose = FALSE)))
  expect_s3_class(fit, "em")
  expect_true(fit$converged)
  expect_length(fit$lambda, 2L)
  # shared parameter A0: its deviations collapse, its strength is frozen high
  expect_lt(mean(abs(fit$etaModes[, "eta_A0"])), 0.05)
  expect_gt(fit$lambda[["lambda_A0"]], 50 * fit$lambda[["lambda_k"]])
  # individual parameter k: deviations recovered near the truth (+-0.6)
  expect_gt(mean(abs(fit$etaModes[, "eta_k"])), 0.4)
})

test_that("single global lambda stays scalar (backward compatible)", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  s <- .build_reg_manyline(tag = "compat")
  pen <- penaltyL1(c("k", "A0"), subjects = s$subj)
  obj <- normL2(s$dlist, s$prd, errmodel = s$e) +
    constraintL1(pen)
  center <- emInit(c(log_k = s$st[["log_k"]], log_A0 = s$st[["log_A0"]],
                      log_sigma = log(0.05)), pen, lambda = 1)
  fit <- suppressWarnings(suppressMessages(EM(obj, center, method = "focei",
                                 control = list(cm1 = list(iterlim = 100L)),
                                 verbose = FALSE)))
  expect_length(fit$lambda, 1L)
  expect_equal(names(fit$lambda), "lambda")
  expect_true(fit$lambda > 0)
})


# ---- WS1: adaptive / kink-robust quadrature for the clustered marginal ----

test_that("adaptive-tensor marginal equals the dense tensor at small d", {
  gh <- .gaussHermite(24L)
  set.seed(1)
  for (i in 1:3) {
    H <- runif(4, 5, 40); m <- rnorm(4, 0, 0.7); lam <- 10^runif(1, -0.5, 1.5)
    va <- .clusterParamMarginal(H, m, lam, gh, rule = "auto")
    vt <- .clusterParamMarginal(H, m, lam, gh, rule = "tensor")
    # d = G-1 = 3, adaptiveNq(3) = 24, so auto and tensor(24) coincide exactly
    expect_equal(va$value, vt$value)
    expect_equal(va$Efus,  vt$Efus)
  }
})

test_that("adaptive-tensor marginal stays feasible and finite at larger G", {
  # d = G-1 = 5: the dense 24^5 (~8M nodes) tensor is infeasible; the adaptive
  # rule shrinks nq (here 7 -> 7^5 ~ 1.7e4) and never produces a signed-sum NaN.
  expect_lte(.adaptiveNq(5L)^5, 3e4)
  set.seed(2); H <- runif(6, 8, 30); m <- c(0.8, 0.75, 0.7, -0.7, -0.75, -0.8)
  for (lam in c(1, 5, 30)) {
    v <- .clusterParamMarginal(H, m, lam, rule = "auto")
    expect_true(is.finite(v$value))
    expect_true(all(is.finite(v$Eeta)))
  }
})

test_that("sparse rule never escapes a NaN (falls back to the dense tensor)", {
  gh <- .gaussHermite(24L)
  set.seed(1); H <- runif(4, 5, 40); m <- rnorm(4, 0, 0.7)
  for (lam in c(1, 7.5, 20)) {
    v <- suppressWarnings(.clusterParamMarginal(H, m, lam, gh, rule = "sparse"))
    expect_true(is.finite(v$value))
  }
})

test_that("grouping selection recovers the truth at n=6 with the adaptive rule", {
  Hd <- rep(50, 6); m <- c(0.5, 0.55, 0.52, -0.5, -0.55, -0.52)  # two groups
  paths <- list(as.list(1:6), list(1:3, 4:6), list(1:6))
  sc <- vapply(paths, function(P) .clusterGroupScore(P, Hd, m, rule = "auto")$value, 0.0)
  expect_equal(which.min(sc), 2L)                 # the G=2 grouping wins
  expect_lt(sc[2], sc[1]); expect_lt(sc[2], sc[3])
})


# ---- WS5: C++ marginal kernels are bit-comparable to the R reference ------

test_that("normalLaplaceCpp reproduces the R normalLaplace", {
  set.seed(1)
  a <- runif(8, 0.5, 50); m <- rnorm(8, 0, 1.5); lam <- c(0.3, 1, 5)
  rR <- normalLaplace(a, m, lam); rC <- normalLaplaceCpp(a, m, lam)
  for (nm in c("logI", "dlogI_dm", "dlogI_da", "dlogI_dlambda",
               "Eeta", "Eabs", "Ecen2", "Ppos"))
    expect_lt(max(abs(rR[[nm]] - rC[[nm]])), 1e-10)
})

test_that("clusterMarginalCpp reproduces the R .clusterParamMarginal", {
  gh <- .gaussHermite(24L)
  withr::local_options(laplaceUseCpp = TRUE)
  for (G in c(3, 4, 5)) {
    set.seed(G); H <- runif(G, 5, 40); m <- rnorm(G, 0, 0.7)
    vC <- .clusterParamMarginal(H, m, 3, gh, rule = "auto")            # C++
    vR <- local({ op <- options(laplaceUseCpp = FALSE); on.exit(options(op))
      .clusterParamMarginal(H, m, 3, gh, rule = "auto") })            # R reference
    expect_lt(abs(vR$value - vC$value), 1e-9)
    expect_lt(max(abs(vR$Eeta  - vC$Eeta)),  1e-9)
    expect_lt(max(abs(vR$Ecen2 - vC$Ecen2)), 1e-9)
    expect_lt(abs(vR$Efus - vC$Efus), 1e-9)
  }
})

test_that("the C++ toggle does not change the clustered grouping selection", {
  Hd <- rep(40, 5); m <- c(0.5, 0.55, 0.52, -0.5, -0.55)
  paths <- list(as.list(1:5), list(1:3, 4:5), list(1:5))
  score <- function(useCpp) {
    op <- options(laplaceUseCpp = useCpp); on.exit(options(op))
    vapply(paths, function(P) .clusterGroupScore(P, Hd, m, rule = "auto")$value, 0.0)
  }
  expect_lt(max(abs(score(TRUE) - score(FALSE))), 1e-7)
})

test_that("the inner trustL1 lock-step reproduces the per-subject solves", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  s <- .build_reg_manyline(tag = "lock")
  pen <- penaltyL1(c("k", "A0"), subjects = s$subj)
  obj <- normL2(s$dlist, s$prd, errmodel = s$e) +
    constraintL1(pen)
  rec <- .laplaceReconstruct(obj); conds <- names(rec$data)
  resObjList <- setNames(lapply(conds, function(cc) {
    d <- rec$data[cc]; class(d) <- class(rec$data)
    normL2(d, rec$prdfn, errmodel = rec$errfn)
  }), conds)
  tgt  <- .laplaceFactorisingTargets(pen)
  ctrl <- list(iterlim = 300L)
  run  <- function(on, theta, correction)
    withr::with_options(list(dMod.laplace.lockstep = on, dMod.cores = 1L),
      .laplaceOuterEval(theta, rec, resObjList, tgt, ctrl, grad = TRUE,
                        hess = TRUE, correction = correction, collect = TRUE))

  for (lam in c(0.5, 3, 20)) for (correction in c(FALSE, TRUE)) {
    theta <- c(log_k = s$st[["log_k"]] + 0.1, log_A0 = s$st[["log_A0"]] - 0.05,
               log_sigma = log(0.06), lambda = lam)
    a <- run(FALSE, theta, correction)
    b <- run(TRUE,  theta, correction)
    # Bit-identical, and per subject: a differing mode or convergence flag means
    # the propose/accept split is not behaviour-preserving, which the summed
    # value alone would not show.
    expect_identical(a$value, b$value)
    expect_equal(a$gradient, b$gradient, tolerance = 0)
    expect_equal(a$hessian,  b$hessian,  tolerance = 0)
    expect_equal(lapply(a$perSubject, `[[`, "etahat"),
                 lapply(b$perSubject, `[[`, "etahat"), tolerance = 0)
    expect_equal(lapply(a$perSubject, `[[`, "Hd"),
                 lapply(b$perSubject, `[[`, "Hd"), tolerance = 0)
    expect_identical(vapply(a$perSubject, `[[`, TRUE, "converged"),
                     vapply(b$perSubject, `[[`, TRUE, "converged"))
  }
})

test_that("the lock-step issues one batched request per round, not one per subject", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  s <- .build_reg_manyline(tag = "lockn")
  pen <- penaltyL1(c("k", "A0"), subjects = s$subj)
  obj <- normL2(s$dlist, s$prd, errmodel = s$e) +
    constraintL1(pen)
  rec <- .laplaceReconstruct(obj); conds <- names(rec$data)
  n <- 0L
  resObjList <- setNames(lapply(conds, function(cc) {
    d <- rec$data[cc]; class(d) <- class(rec$data)
    f <- normL2(d, rec$prdfn, errmodel = rec$errfn)
    w <- function(...) { n <<- n + 1L; f(...) }
    attributes(w) <- attributes(f); w
  }), conds)
  theta <- c(log_k = s$st[["log_k"]] + 0.1, log_A0 = s$st[["log_A0"]] - 0.05,
             log_sigma = log(0.06), lambda = 3)
  batches <- 0L
  orig <- get(".predictMany", envir = asNamespace("dMod2"))
  local_mocked_bindings(
    .predictMany = function(...) { batches <<- batches + 1L; orig(...) },
    .package = "dMod2")
  withr::local_options(dMod.laplace.lockstep = TRUE, dMod.cores = 1L)
  .laplaceOuterEval(theta, rec, resObjList, .laplaceFactorisingTargets(pen),
                    list(iterlim = 300L), grad = TRUE, hess = TRUE)
  expect_gt(n, 0L)
  expect_gt(batches, 0L)
  # every round covers several subjects, so requests must be well below evaluations
  expect_lt(batches, n / 2)
})
