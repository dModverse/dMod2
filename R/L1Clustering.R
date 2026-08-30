## Clustered (complete-graph fusion) inner solver for .fitLaplace.
##
## The clustered penalty penaltyL1(method = "clustered") penalises the
## complete-graph pairwise differences sum_{i<j} |eta_{i,k} - eta_{j,k}| across
## the n subjects, per candidate parameter k, so individuals crystallise into
## clusters that share a value. Unlike penalty()/groupPenalty("fused") (which
## factorise over subjects and are solved coordinate-wise by trustL1), the
## clustered penalty COUPLES the subjects, so the inner conditional-mode solve is
## joint over all n subjects.
##
## The GN-linearised per-parameter subproblem is the weighted complete-graph
## fused lasso, which reduces EXACTLY to sort-by-mode + weighted isotonic
## regression (PAVA); see solveFusedComplete(). .clusterSolve() wraps it in a
## joint prox-linear Gauss-Newton driver over all subjects. .clusterMap()
## is the runnable MAP-EM fit (grouping + a MAP lambda); the exact marginal-ML
## FOCEI value/gradient is a later step (see notes/laplace_nlme_theory.Rmd, the
## clustered marginal 2 f + log|H_red| + 2 log Z(lambda) - G log 2pi).


## Weighted complete-graph fused lasso in one dimension.
##
## Solves, for n subjects with data precisions w_i > 0 and data-alone modes m_i,
##   argmin_u  (1/2) sum_i w_i (u_i - m_i)^2  +  lambda sum_{i<j} |u_i - u_j| .
## Because the complete-graph total-variation term is linear in the order
## statistics (sum_{i<j}|u_i-u_j| = sum_i (2 rank(u_i) - n - 1) u_(i)), the
## solution is EXACTLY the weighted isotonic regression (pool-adjacent-violators)
## of the rank-shifted values v_i = m_i - lambda (2 i - n - 1)/w_i, taken in the
## order of m. The PAVA pool blocks are the clusters; each block value is the
## w-weighted mean of its v_i. O(n log n), exact, dependency-free. Limits:
## lambda = 0 -> u = m (singletons); lambda -> Inf -> one block at the w-weighted
## grand mean of m (all fused), independent of lambda (since sum(2i-n-1) = 0).
##
## Returns list(u = <named per-subject values, original order>,
##   clusters = <list of subject-id character groups>, values = <per-block value>,
##   G = <number of clusters>).
solveFusedComplete <- function(m, w, lambda) {
  n  <- length(m)
  nm <- names(m); if (is.null(nm)) nm <- as.character(seq_len(n))
  if (length(w) != n) stop("`w` must have the same length as `m`.")
  if (any(w <= 0))    stop("`w` (data precisions) must be positive.")
  if (lambda < 0)     stop("`lambda` must be non-negative.")
  if (n == 1L)
    return(list(u = stats::setNames(as.numeric(m), nm), clusters = list(nm),
                values = as.numeric(m), G = 1L))

  ord <- order(m)                                   # sort ascending by data mode
  ms  <- as.numeric(m)[ord]; ws <- as.numeric(w)[ord]; nms <- nm[ord]
  v   <- ms - lambda * (2 * seq_len(n) - n - 1) / ws  # rank-shifted targets

  ## weighted PAVA: stack of blocks (value, weight, member positions in sorted order)
  blkVal <- numeric(0); blkW <- numeric(0); blkMem <- list()
  for (i in seq_len(n)) {
    cv <- v[i]; cw <- ws[i]; cmem <- i
    while (length(blkVal) > 0L && blkVal[length(blkVal)] >= cv) {
      j    <- length(blkVal)
      newW <- blkW[j] + cw
      cv   <- (blkVal[j] * blkW[j] + cv * cw) / newW
      cw   <- newW
      cmem <- c(blkMem[[j]], cmem)
      blkVal <- blkVal[-j]; blkW <- blkW[-j]; blkMem <- blkMem[-j]
    }
    blkVal <- c(blkVal, cv); blkW <- c(blkW, cw); blkMem <- c(blkMem, list(cmem))
  }

  u_sorted <- numeric(n); clusters <- vector("list", length(blkVal))
  for (g in seq_along(blkVal)) {
    idx <- blkMem[[g]]
    u_sorted[idx] <- blkVal[g]
    clusters[[g]] <- nms[idx]
  }
  u <- numeric(n); u[ord] <- u_sorted; names(u) <- nm
  list(u = u, clusters = clusters, values = blkVal, G = length(blkVal))
}


## Per-subject PROFILE (conditional) linearisation at a given deviation matrix
## `Eta` (subjects x K). For each subject it evaluates the data objfn (deriv), takes
## the K x K joint Hessian over the candidate parameters, and returns, per parameter,
## the profile precision 1/diag(H^{-1}) and the joint Newton data mode
## eta - H^{-1} g. The profile precision (unlike the raw diagonal H_kk) accounts for
## the correlation between candidate parameters, so a per-subject mode that only
## differs from its cluster because of cross-attribution is not over-confident and
## the marginal can fuse it. For K = 1 it reduces to the diagonal. A tiny ridge +
## singular fallback keep it stable for degenerate directions.
.clusterProfileLinearise <- function(resObjList, outer_struct, se, Eta) {
  subjects <- rownames(se); eta_cols <- colnames(se)
  n <- length(subjects); K <- length(eta_cols)
  all_eta  <- as.vector(se)
  base_eta <- stats::setNames(rep(0, length(all_eta)), all_eta)
  Hd <- m <- grad <- matrix(0, n, K, dimnames = list(subjects, eta_cols))
  fulls <- lapply(seq_len(n), function(si) {
    f <- c(outer_struct, base_eta); f[as.character(se[subjects[si], ])] <- Eta[si, ]; f
  })
  ## One batched pass over the subjects; on a batch failure fall back to the
  ## per-subject loop, which is what carries the tolerated-failure branch below.
  os <- tryCatch(.objEvalMany(resObjList[subjects], fulls, deriv = TRUE),
                 error = function(e) NULL)
  if (is.null(os)) os <- lapply(seq_len(n), function(si)
    tryCatch(resObjList[[subjects[si]]](fulls[[si]], deriv = TRUE),
             error = function(e) NULL))
  for (si in seq_len(n)) {
    en <- as.character(se[subjects[si], ])
    o  <- os[[si]]
    if (is.null(o)) { Hd[si, ] <- 1e-8; m[si, ] <- Eta[si, ]; next }
    gj <- o$gradient[en]; Hj <- o$hessian[en, en, drop = FALSE]; grad[si, ] <- gj
    ## Invert the joint Hessian exactly as dMod's trust() treats its Newton system:
    ## when the Hessian is positive definite (all eigenvalues > 0, the normal case
    ## for a Gauss-Newton information matrix) take the plain Newton inverse, no
    ## regularisation. A degenerate direction (eigenvalue <= 0, e.g. the kon<->koff
    ## binding degeneracy) carries no data information; trust() would handle it via
    ## its trust radius, and here, with no radius, we simply floor those
    ## non-positive eigenvalues so the direction gets negligible precision and is
    ## correctly left shared (never over-split). Identifiable directions are
    ## untouched, so the profile precision / joint Newton mode are exact.
    eg  <- eigen(Hj, symmetric = TRUE)
    lam <- eg$values
    if (any(lam <= 0)) lam <- pmax(lam, 1e-8 * max(abs(lam)))
    Cj   <- eg$vectors %*% (t(eg$vectors) / lam)
    prof <- 1 / diag(Cj)
    low  <- prof < max(1e-8, 1e-8 * max(prof))
    mr   <- Eta[si, ] - as.numeric(Cj %*% gj); mr[low] <- Eta[si, low]
    Hd[si, ] <- pmax(prof, 1e-8); m[si, ] <- mr
  }
  list(Hd = Hd, m = m, grad = grad)
}


## Joint clustered inner solve (the "trustGL" role).
##
## Finds the joint penalised conditional mode of all `n` subjects' deviations for
## the clustered penalty, by proximal Gauss-Newton: at the current eta, linearise
## each subject's data (Gauss-Newton, separable per parameter -> data precision
## Hd_{i,k} and data mode m_{i,k}), solve each parameter's complete-graph fused
## lasso exactly with solveFusedComplete(), then backtrack-line-search the joint
## step on the true objective sum_i value_data_i(eta) + lambda * ||D eta||_1.
## Pure R (a C++ trustGL kernel is a later performance step).
##
## `resObjList` is .fitLaplace's per-condition normL2 list; only the penalised
## `subjects` (rownames(penalty$subjectEtas)) are solved here (reference
## conditions carry no eta and are handled by the outer/structural step).
##
## Returns the joint bundle the clustered marginal / MAP-EM consume:
##   etahat    subjects x K matrix of the fused deviations (colnames = eta_<par>)
##   clusters  list (length K) of per-parameter subject-id partitions
##   Gk, G     per-parameter and total cluster counts (G = sum(Gk))
##   fusionL1  ||D etahat||_1 = sum_k sum_{i<j} |etahat_{i,k}-etahat_{j,k}|
##   Hd, m     subjects x K data precisions / data modes at the mode
##   value     the penalised inner objective at the mode
##   converged logical
.clusterSolve <- function(resObjList, outer_struct, penalty, lambda,
                             eta0 = NULL, maxit = 50L, tol = 1e-8,
                             maxls = 20L, verbose = FALSE) {
  se       <- penalty$subjectEtas
  subjects <- rownames(se)
  eta_cols <- colnames(se)                       # eta_<par>
  n <- length(subjects); K <- length(eta_cols)
  all_eta  <- as.vector(se)
  base_eta <- stats::setNames(rep(0, length(all_eta)), all_eta)

  ## current deviations as a subjects x K matrix
  Eta <- matrix(0, n, K, dimnames = list(subjects, eta_cols))
  if (!is.null(eta0)) {
    have <- intersect(names(eta0), all_eta)
    for (nm in have) { rc <- which(se == nm, arr.ind = TRUE); Eta[rc[1], rc[2]] <- eta0[[nm]] }
  }

  fullFor <- function(s, eta_row) {
    full <- c(outer_struct, base_eta)
    full[as.character(se[s, ])] <- eta_row
    full
  }
  ## True objective: summed subject data -2logL + lambda * complete-graph fusion.
  ## ROBUST: an ODE-solver failure or non-finite value is treated as +Inf so the
  ## step is rejected (mirrors how a trust-region method rejects a bad model
  ## evaluation), instead of throwing and aborting the whole fit.
  Fval <- function(Eta) {
    ## Any failing subject makes the whole value +Inf, so a failing batch and a
    ## failing single solve have the same outcome.
    os <- tryCatch(.objEvalMany(resObjList[subjects],
                                lapply(subjects, function(s) fullFor(s, Eta[s, ])),
                                deriv = FALSE),
                   error = function(e) NULL)
    if (is.null(os)) return(Inf)
    dv <- 0
    for (o in os) {
      if (!is.finite(o$value)) return(Inf)
      dv <- dv + o$value
    }
    fus <- 0
    for (k in seq_len(K)) fus <- fus + sum(stats::dist(Eta[, k]))
    dv + lambda * fus
  }
  ## per-subject GN linearisation (gradient + Hessian diagonal), robust to a
  ## solver failure at the current point.
  linearise <- function(Eta) {
    G_ <- H_ <- matrix(0, n, K, dimnames = list(subjects, eta_cols))
    os <- tryCatch(.objEvalMany(resObjList[subjects],
                                lapply(seq_len(n), function(si) fullFor(subjects[si], Eta[si, ])),
                                deriv = TRUE),
                   error = function(e) NULL)
    if (is.null(os)) return(list(grad = G_, Hd = H_, ok = FALSE))
    for (si in seq_len(n)) {
      o <- os[[si]]
      en <- as.character(se[subjects[si], ])
      G_[si, ] <- o$gradient[en]
      H_[si, ] <- pmax(diag(o$hessian[en, en, drop = FALSE]), 0)  # Gauss-Newton Hd >= 0
    }
    list(grad = G_, Hd = H_, ok = TRUE)
  }

  Fcur <- Fval(Eta); conv <- FALSE
  clusters <- vector("list", K); Gk <- integer(K)
  Hd <- m <- grad <- matrix(0, n, K, dimnames = list(subjects, eta_cols))
  mu <- 0; mu_init <- 1e-3                         # Levenberg-Marquardt damping
  for (it in seq_len(maxit)) {
    lin <- linearise(Eta); if (!lin$ok) { conv <- FALSE; break }
    grad <- lin$grad; Hd <- lin$Hd
    ## LM-damped prox-linear step: the damped precision w = Hd + mu bounds the data
    ## mode m = eta - g/(Hd + mu), so an ill-conditioned or unidentified direction
    ## (Hd -> 0) cannot propose an absurd, ODE-breaking mode; for a well-informed
    ## direction (Hd >> mu) the damping is negligible. mu is adapted like a trust
    ## radius: raised when a step fails (worse objective OR a solver failure, both
    ## +Inf via Fval), lowered when it succeeds. No hard identifiability threshold.
    accepted <- FALSE; dF <- 0; step <- 0
    for (tries in seq_len(maxls)) {
      w  <- Hd + mu
      md <- Eta - grad / w
      Eta_try <- Eta
      for (k in seq_len(K)) {
        sol <- solveFusedComplete(stats::setNames(md[, k], subjects), w[, k], lambda)
        Eta_try[, k] <- sol$u[subjects]; clusters[[k]] <- sol$clusters; Gk[k] <- sol$G
      }
      Ftry <- Fval(Eta_try)
      if (is.finite(Ftry) && Ftry <= Fcur + 1e-12) {
        step <- max(abs(Eta_try - Eta)); dF <- Fcur - Ftry
        Eta <- Eta_try; Fcur <- Ftry; accepted <- TRUE; mu <- mu / 2; break
      }
      mu <- if (mu <= 0) mu_init else mu * 4
    }
    if (verbose)
      cat(sprintf("  regCluster it %d: F=%.6f step=%.2e mu=%.2e\n", it, Fcur, step, mu))
    if (!accepted) { conv <- TRUE; break }         # damping exhausted -> stationary
    if (step < tol || abs(dF) < tol * (abs(Fcur) + 1)) { conv <- TRUE; break }
  }
  ## Final linearisation ingredients for the marginal. Use the PROFILE (conditional)
  ## precision: the diagonal of the inverse per-subject joint Hessian over the K
  ## candidate parameters, inverted back, not the raw diagonal. With several
  ## candidates the raw diagonal Hd_kk is over-confident: a per-subject mode that
  ## only appears to differ from its group because of cross-attribution with the
  ## OTHER (jointly estimated) candidates would then look statistically significant
  ## and the marginal would over-split. The profile precision deflates exactly that
  ## correlation, so genuinely-shared subjects fuse. For K = 1 it equals the diagonal
  ## (no correlation), so the single-candidate path is unchanged. The mode is the
  ## joint Newton step m = eta - H^{-1} g (data-alone conditional mode).
  ## An unidentified direction (precision -> 0) still overflows m, so floor it
  ## relative to the row scale and hold that coordinate's mode at eta.
  lp <- .clusterProfileLinearise(resObjList, outer_struct, se, Eta)
  Hd <- lp$Hd; m <- lp$m; grad <- lp$grad
  names(clusters) <- eta_cols
  list(etahat = Eta, clusters = clusters, Gk = Gk, G = sum(Gk),
       fusionL1 = sum(vapply(seq_len(K), function(k) sum(stats::dist(Eta[, k])), 0.0)),
       Hd = Hd, m = m, grad = grad, value = Fcur, converged = conv)
}


## Runnable clustered fit via MAP-EM (penalised point estimate).
##
## The clustered analogue of .fitLaplace's ECM, but a MAP (Laplace/point) EM rather
## than the marginal ML: the inner clustered mode is the E-step surrogate, and the
## structural + lambda M-steps are closed form. It recovers the grouping and a MAP
## lambda; the exact marginal-ML FOCEI value/gradient (with the log|H_red| volume
## term) is a later step, so `value` here is the MAP penalised objective, NOT a
## marginal -2logL, and clustered fits are not (yet) comparable across models.
##
## Loop per outer iteration:
##   inner  .clusterSolve at (mu, lambda) -> joint clustered mode eta_hat.
##   anchor recentre eta_hat_{.,k} to mean 0, fold the mean into the structural
##          mu_k = log_<par> (the shift eta->eta+c, mu->mu-c leaves prediction AND
##          ||D eta|| unchanged, so mu / eta-level is unidentified without it).
##   CM-2   lambda = K(n-1) / ||D eta_hat||_1  (complete-graph normalizer rho=K(n-1)).
##   CM-1   trust() on the complete-data objective at the frozen anchored eta_hat.
.clusterMap <- function(rec, resObjList, free, fixed, lam_name,
                              control, verbose, eta_init = NULL) {
  penalty  <- rec$penalty
  kinds    <- vapply(penalty$blocks, `[[`, "", "kind")
  if (!all(kinds == "clustered"))
    stop(".fitLaplace: the clustered MAP-EM path needs an all-clustered penalty ",
         "(mixing 'clustered' with 'single'/'fused' is not supported).",
         call. = FALSE)
  se       <- penalty$subjectEtas
  subjects <- rownames(se); eta_cols <- colnames(se); params <- penalty$params
  n <- length(subjects); K <- length(eta_cols)
  all_eta  <- as.vector(se)
  base_eta <- stats::setNames(rep(0, length(all_eta)), all_eta)
  struct_names   <- setdiff(names(free), lam_name)
  ref_conditions <- setdiff(names(rec$data), subjects)
  rho    <- K * (n - 1L)
  mu_of  <- stats::setNames(paste0("log_", params), eta_cols)   # eta_<par> -> log_<par>

  cm1 <- .trustControl(list(rinit = 1, rmax = 10, iterlim = 30L,
                            ftol = 1e-6, mtol = 1e-6),
                       control$cm1, label = "control$cm1")
  clc <- modifyList(list(maxit = 60L, tol = 1e-9),
                    if (is.null(control$cluster)) list() else control$cluster)
  ## lambda is HELD FIXED by default. The closed-form MAP M-step
  ## lambda = rho/||D eta_hat||_1 degenerates here (estimation is not selection):
  ## with precise data ||D eta_hat||_1 is dominated by the large between-cluster
  ## gaps, so the point estimate of lambda is far too small to fuse anything, and
  ## the fixed point sits at no fusion. Consistent lambda selection needs the
  ## marginal E||D eta||_1 (a later step) or a lambda-scan + BIC. Opt into the
  ## (degenerate) MAP update with control$estimateLambda = TRUE only to inspect it.
  estimateLambda <- isTRUE(control$estimateLambda)
  maxOuter  <- if (is.null(control$maxOuter))  50L  else as.integer(control$maxOuter)
  epsPar    <- if (is.null(control$epsPar))    1e-4 else control$epsPar
  epsOfvRel <- if (is.null(control$epsOfvRel)) 1e-5 else control$epsOfvRel

  struct_cur <- free[struct_names]
  lambda_cur <- as.numeric(free[[lam_name]])
  eta_warm   <- if (is.null(eta_init)) NULL else eta_init[all_eta]

  ## CM-1 objective: structural at frozen etas (data only; fusion const).
  cm1_obj <- function(theta_struct, Eta) {
    os <- c(stats::setNames(as.numeric(theta_struct), struct_names), fixed)
    total <- 0
    gg <- stats::setNames(numeric(length(struct_names)), struct_names)
    HH <- matrix(0, length(struct_names), length(struct_names),
                 dimnames = list(struct_names, struct_names))
    ## robust to a solver failure at a trust trial point: return +Inf so trust()
    ## rejects the step and shrinks its radius (the trust-region way).
    conds <- c(subjects, ref_conditions)
    fulls <- lapply(conds, function(cc) {
      full <- c(os, base_eta)
      if (cc %in% subjects) full[as.character(se[cc, ])] <- Eta[match(cc, subjects), ]
      full
    })
    ## Any failure already collapses the whole objective to +Inf, so catching the
    ## batch as a whole matches the per-condition catch it replaces.
    got <- tryCatch(.objEvalMany(resObjList[conds], fulls, deriv = TRUE),
                    error = function(e) NULL)
    if (is.null(got)) return(objlist(value = Inf, gradient = gg, hessian = HH))
    for (ci in seq_along(conds)) {
      o <- got[[ci]]
      if (!is.finite(o$value))
        return(objlist(value = Inf, gradient = gg, hessian = HH))
      total <- total + o$value
      sn <- intersect(struct_names, names(o$gradient))
      gg[sn] <- gg[sn] + o$gradient[sn]
      HH[sn, sn] <- HH[sn, sn] + o$hessian[sn, sn, drop = FALSE]
    }
    objlist(value = total, gradient = gg, hessian = HH)
  }

  prev_par <- c(struct_cur, stats::setNames(lambda_cur, lam_name)); prev_ofv <- NA_real_
  rows <- list(); conv <- FALSE; nOuter <- 0L; sol <- NULL
  for (it in seq_len(maxOuter)) {
    nOuter <- it
    outer_struct <- c(struct_cur, fixed)
    ## inner clustered mode
    sol <- .clusterSolve(resObjList, outer_struct, penalty, lambda_cur,
                            eta0 = eta_warm, maxit = clc$maxit, tol = clc$tol)
    Eta <- sol$etahat
    ## anchor: recentre each parameter column, fold the mean into log_<par>
    for (kk in seq_len(K)) {
      c_k <- mean(Eta[, kk])
      Eta[, kk] <- Eta[, kk] - c_k
      mu_nm <- mu_of[[eta_cols[kk]]]
      if (mu_nm %in% names(struct_cur)) struct_cur[[mu_nm]] <- struct_cur[[mu_nm]] + c_k
    }
    eta_warm <- stats::setNames(as.vector(Eta), as.vector(se))
    ## CM-2: closed-form MAP lambda (opt-in; degenerates, see the note above).
    fusion <- sum(vapply(seq_len(K), function(kk) sum(stats::dist(Eta[, kk])), 0.0))
    if (estimateLambda && fusion > 0) lambda_cur <- rho / fusion
    ## CM-1: structural M-step at the frozen anchored etas
    cm1_fit <- suppressMessages(do.call(trust, modifyList(cm1, list(
      objfun = function(th, ...) cm1_obj(th, Eta), parinit = struct_cur))))
    struct_cur <- cm1_fit$argument
    ofv <- cm1_fit$value + lambda_cur * fusion         # MAP penalised objective

    cur_par     <- c(struct_cur, stats::setNames(lambda_cur, lam_name))
    deltaPar    <- max(abs(cur_par - prev_par))
    deltaOfvRel <- if (is.na(prev_ofv) || abs(prev_ofv) < .Machine$double.eps)
                     Inf else abs(ofv - prev_ofv) / abs(prev_ofv)
    rows[[length(rows) + 1L]] <- data.frame(
      outerIter = it, OFV = ofv, lambda = lambda_cur, fusion = fusion,
      G = sol$G, deltaPar = deltaPar)
    if (verbose)
      cat(sprintf("  EM(clustered) outer %d : MAP=%.4f lambda=%.5g G=%d |dpar|=%.2e\n",
                  it, ofv, lambda_cur, sol$G, deltaPar))
    prev_par <- cur_par; prev_ofv <- ofv
    if (deltaPar < epsPar || deltaOfvRel < epsOfvRel) { conv <- TRUE; break }
  }

  Eta <- sol$etahat
  for (kk in seq_len(K)) Eta[, kk] <- Eta[, kk] - mean(Eta[, kk])   # report anchored
  colnames(Eta) <- eta_cols
  out <- list(
    argument   = c(struct_cur, stats::setNames(lambda_cur, lam_name)),
    value      = prev_ofv,
    valueType  = "MAP penalised objective (not a marginal -2logL)",
    lambda     = lambda_cur,
    etaModes   = Eta,
    clusters   = sol$clusters,
    Gk         = sol$Gk,
    G          = sol$G,
    clustered  = TRUE,
    converged  = conv,
    iterations = nOuter,
    stageTrace = do.call(rbind, rows),
    method     = "focei-map",
    prior      = "penaltyL1",
    penalty    = penalty,
    prdfn      = rec$prdfn,
    data       = rec$data,
    errfn      = rec$errfn)
  class(out) <- c("em", "list")
  out
}


## ---- exact clustered marginal via adaptive Gauss-Hermite quadrature ----------
##
## The complete-graph fusion penalty is piecewise-linear over the n! order cones
## (in each cone sum|eta_i-eta_j| = c_sigma' eta is linear), so the FOCE-linearised
## integrand exp(-1/2[V_lin + lambda ||D eta||_1]) is a Gaussian truncated to each
## cone, the multi-dimensional analogue of the 1-D normal-Laplace (normalLaplace).
## Rather than the intractable closed-form order-cone probabilities, the exact
## marginal is computed by mode-centred adaptive Gauss-Hermite quadrature over the
## anchored (sum eta = 0) space, and its gradient by Fisher's identity (posterior
## expectations over the quadrature nodes). This is EXACT (to quadrature accuracy),
## unlike the conditional-on-partition Gaussian-volume approximation whose corner
## error is O(nats) at strong fusion.

## physicists' Gauss-Hermite nodes/weights (weight e^{-x^2}) via Golub-Welsch.
.gaussHermite <- function(nq) {
  i <- seq_len(nq - 1L); b <- sqrt(i / 2)
  J <- matrix(0, nq, nq); J[cbind(i, i + 1L)] <- b; J[cbind(i + 1L, i)] <- b
  ev <- eigen(J, symmetric = TRUE); ord <- order(ev$values)
  list(x = ev$values[ord], w = sqrt(pi) * (ev$vectors[1, ord])^2)
}

## Dimension-adaptive dense-tensor node count: keep the total node count nq^d
## bounded (~3e4) by shrinking the per-dimension nodes as the anchored dimension
## d = G-1 grows, never exceeding `cap` (the user's node budget). All-positive
## weights (unlike the signed Smolyak sparse grid, it never produces a negative
## partial sum on the sharp fusion kink), so it is the robust default for the
## clustered marginal at larger G.
.adaptiveNq <- function(d, cap = 24L) {
  if (d <= 2L) return(cap)
  max(4L, min(cap, as.integer(floor(exp(log(3e4) / d)))))
}

## Helmert orthonormal basis of {sum eta = 0}: n x (n-1), columns orthonormal, sum 0.
.anchorBasis <- function(n) {
  B <- matrix(0, n, n - 1L)
  for (k in seq_len(n - 1L)) {
    v <- c(rep(1, k), -k, rep(0, n - k - 1L)); B[, k] <- v / sqrt(sum(v^2))
  }
  B
}

## Mode-centred Gauss-Hermite log-integral of exp(-phi) with local precision
## `Prec` at the centre `zhat`, over the change of variables z = zhat + sqrt(2) L u
## (L = chol of the local covariance). Two node sets:
##   sparse = FALSE : full tensor product of the 1-D rule `gh` (nq^d nodes,
##                    all-positive weights). Exact for small d, blows up for d>=4.
##   sparse = TRUE  : Smolyak sparse grid (`sparseGridGH`, the C++ rule reused
##                    from the Gaussian nlme path) at Smolyak depth `level`
##                    (default d + 3). Node count grows polynomially in d, so it
##                    stays feasible for larger clusters (more groups). Weights are
##                    signed, so the reduction is a signed log-sum-exp.
## Returns the log-integral, the (possibly signed) posterior weights over nodes,
## and the node coordinates (d x M) for moment accumulation. The tensor path is
## bit-identical to the previous implementation.
.aghqLogI <- function(phi, zhat, Prec, gh = NULL, sparse = FALSE, level = NULL) {
  d   <- length(zhat)
  L   <- t(chol(solve(Prec)))                      # z = zhat + sqrt(2) L u
  cst <- 0.5 * d * log(2) + sum(log(diag(L)))
  if (!sparse) {
    grid <- as.matrix(expand.grid(rep(list(seq_along(gh$x)), d)))
    M    <- nrow(grid)
    logw <- rowSums(matrix(log(gh$w[grid]), M, d)) # log prod GH weights
    U    <- matrix(gh$x[grid], M, d)               # M x d node coords
    z2   <- rowSums(U^2)
    sgn  <- rep(1, M)
  } else {
    if (is.null(level)) level <- d + 3L
    g2   <- .getSparseGrid(d, as.integer(level))   # signed Smolyak weights
    U    <- g2$zNodes; M <- g2$B                   # M x d
    logw <- g2$logAbsW; z2 <- g2$z2Sum; sgn <- g2$signs
  }
  Z <- matrix(0, d, M); logterms <- numeric(M)
  for (q in seq_len(M)) {
    z <- as.numeric(zhat + sqrt(2) * (L %*% U[q, ])); Z[, q] <- z
    logterms[q] <- logw[q] + z2[q] - phi(z)        # z2 undoes the exp(-u^2) weight
  }
  mx    <- max(logterms)
  terms <- sgn * exp(logterms - mx)                # signed; sgn == 1 for tensor
  S     <- sum(terms)
  ## a well-resolved positive integrand gives S > 0 even with signed weights; if
  ## the sparse rule under-resolves it (rare), refine the Smolyak depth once.
  if (sparse && (!is.finite(S) || S <= 0)) {
    g3   <- .getSparseGrid(d, as.integer(level + 2L))
    U    <- g3$zNodes; M <- g3$B; logw <- g3$logAbsW; z2 <- g3$z2Sum; sgn <- g3$signs
    Z <- matrix(0, d, M); logterms <- numeric(M)
    for (q in seq_len(M)) {
      z <- as.numeric(zhat + sqrt(2) * (L %*% U[q, ])); Z[, q] <- z
      logterms[q] <- logw[q] + z2[q] - phi(z)
    }
    mx <- max(logterms); terms <- sgn * exp(logterms - mx); S <- sum(terms)
  }
  post <- terms / S
  ## S <= 0 only when the signed sparse rule fails on a kink; return NaN (without
  ## a log-of-negative warning) so the caller falls back to the dense tensor.
  logI <- if (is.finite(S) && S > 0) cst + mx + log(S) else NaN
  list(logI = logI, post = post, Z = Z)
}

## Exact clustered marginal for ONE parameter over the n subjects.
##   H, m   : n-vectors, data precisions Hd and data modes (FOCE linearisation)
##   lambda : penalty strength; etahat: penalised mode (n-vector) for centring
## Returns: value (-2logL contribution, prior-normalised via Z), the posterior
## moments Eeta (n), Ecen2 (n, = E[(eta-m)^2]) for the Fisher outer gradient, and
## Efus (= E_post sum_{i<j}|eta_i-eta_j|) for the closed-form lambda M-step.
## `sizes` (default all 1) gives the WEIGHTED complete-graph fusion for a REDUCED
## (grouped) model: if `H`/`m` are the G group precisions / data modes and `sizes`
## the group sizes n_g, the fusion is sum_{g<h} n_g n_h |w_g - w_h| and the
## anchoring is sum_g n_g w_g = 0 (so the full model is `sizes = 1`). This is what
## the grouping comparison (structure selection) consumes.
## `center` (anchored coords, length G-1) overrides the quadrature centre; NULL
## uses the data-mode projection. Fisher's identity gives the EXACT marginal
## gradient regardless of the centre, but a centre that moves with the argument
## adds a quadrature-grid-motion term to a finite difference of the value; hold it
## fixed to check the analytic gradient against FD.
## `rule` selects the quadrature over the anchored dimension d = G-1:
##   "tensor" : dense nq^d product using the supplied `gh` (fixed nq).
##   "auto"   : dense tensor with a DIMENSION-ADAPTIVE nq (`.adaptiveNq`, shrinks
##              with d to bound the node count): the robust default; all-positive
##              weights, no signed-cancellation failure on the fusion kink.
##   "sparse" : Smolyak sparse grid at depth `level` (default d + 3). Cheapest at
##              large d for SMOOTH integrands, but its signed weights can fail at
##              strong fusion (kink); use only when the fusion is weak.
## The value is EXACT to quadrature accuracy in every case.
.clusterParamMarginal <- function(H, m, lambda, gh = NULL, sizes = rep(1, length(H)),
                                  center = NULL, rule = c("auto", "tensor", "sparse"),
                                  level = NULL) {
  rule <- match.arg(rule)
  G <- length(H)
  if (G == 1L) return(list(value = 0, Eeta = m, Ecen2 = 0, Efus = 0))
  B  <- MASS::Null(matrix(sizes, G, 1))             # G x (G-1) orthonormal basis of {sizes' w = 0}
  d  <- G - 1L
  sparse <- identical(rule, "sparse")
  ## "auto": dimension-adaptive dense tensor, capped at the caller's node budget
  ## (the length of the supplied `gh`, i.e. quadNodes; 24 when none is given).
  if (identical(rule, "auto"))
    gh <- .gaussHermite(.adaptiveNq(d, cap = if (!is.null(gh)) length(gh$x) else 24L))
  if (!sparse && is.null(gh))
    stop(".clusterParamMarginal: the tensor rule needs a Gauss-Hermite rule `gh`.")
  Wp <- outer(sizes, sizes); ut <- upper.tri(Wp)
  Dp <- function(w) { dd <- as.matrix(stats::dist(w)); sum(dd[ut] * Wp[ut]) }
  Az <- crossprod(B, H * B)                         # B' diag(H) B  (d x d)
  zc <- if (is.null(center)) as.numeric(solve(Az, crossprod(B, H * m))) else center
  phi  <- function(z) { e <- as.numeric(B %*% z)
    0.25 * sum(H * (e - m)^2) + 0.5 * lambda * Dp(e) }
  phiZ <- function(z) 0.5 * lambda * Dp(as.numeric(B %*% z))
  ## C++ node loop (dense/adaptive tensor only; all-positive weights). Builds the
  ## same tensor node set the R .aghqLogI does and runs the hot sweep in C++.
  if (!sparse && .laplaceUseCpp()) {
    grid <- as.matrix(expand.grid(rep(list(seq_along(gh$x)), d)))
    Mc   <- nrow(grid)
    U    <- matrix(gh$x[grid], Mc, d)
    logw <- rowSums(matrix(log(gh$w[grid]), Mc, d))
    z2   <- rowSums(U^2)
    Lc   <- t(chol(solve(0.5 * Az)))
    cm   <- clusterMarginalCpp(H, m, lambda, sizes, B, Lc, zc, U, logw, z2,
                               rep(1, Mc))
    if (is.finite(cm$value))
      return(list(value = cm$value, Eeta = cm$Eeta, Ecen2 = cm$Ecen2,
                  Efus = cm$Efus))
  }
  aI <- .aghqLogI(phi,  zc,        0.5 * Az, gh, sparse = sparse, level = level)
  aZ <- .aghqLogI(phiZ, rep(0, d), 0.5 * Az, gh, sparse = sparse, level = level)
  value <- -2 * aI$logI + 2 * aZ$logI
  ## sparse safety net: the signed Smolyak rule can fail (negative partial sum ->
  ## NaN) on a sharp fusion kink; fall back to the robust adaptive dense tensor.
  if (sparse && !is.finite(value)) {
    gh2 <- .gaussHermite(.adaptiveNq(d))
    aI  <- .aghqLogI(phi,  zc,        0.5 * Az, gh2)
    aZ  <- .aghqLogI(phiZ, rep(0, d), 0.5 * Az, gh2)
    value <- -2 * aI$logI + 2 * aZ$logI
  }
  ws    <- B %*% aI$Z                               # G x M reconstructed group values
  Eeta  <- as.numeric(ws %*% aI$post)
  Ecen2 <- as.numeric(((ws - m)^2) %*% aI$post)
  Efus  <- sum(aI$post * apply(ws, 2, Dp))
  list(value = value, Eeta = Eeta, Ecen2 = Ecen2, Efus = Efus)
}


## Exact-marginal score of a grouping P (list of member-index vectors) given the
## FOCE linearisation (per-subject data precisions Hd, data modes m). The grouping
## constrains members to share a value: -2logL(P) = within-group data cost sum_g C_g
## (C_g = sum_{i in g} 1/2 Hd_i (m_i - m_g)^2, the cost of forcing the group equal)
## plus the reduced weighted clustered marginal on the G group values, evaluated at
## the grouping's own marginal-ML lambda (M-step fixed point lambda = (G-1)/E||Dw||).
## The null grouping (G = 1) has no random effect: -2logL = C only.
.clusterGroupScore <- function(P, Hd, m, gh = NULL, lamIter = 30L,
                                  rule = c("auto", "tensor", "sparse"),
                                  level = NULL) {
  rule <- match.arg(rule)
  G  <- length(P)
  Hg <- vapply(P, function(g) sum(Hd[g]), 0.0)
  mg <- vapply(P, function(g) sum(Hd[g] * m[g]) / sum(Hd[g]), 0.0)
  sz <- vapply(P, length, 0L)
  Cg <- sum(vapply(P, function(g) sum(0.5 * Hd[g] * (m[g] - sum(Hd[g] * m[g]) / sum(Hd[g]))^2), 0.0))
  if (G == 1L) return(list(value = Cg, lambda = NA_real_))
  rho <- G - 1L; lam <- 1
  for (it in seq_len(lamIter)) {                    # marginal-ML lambda fixed point
    ef <- .clusterParamMarginal(Hg, mg, lam, gh, sizes = sz,
                                rule = rule, level = level)$Efus
    lam <- if (ef > 0) rho / ef else lam
  }
  list(value = Cg + .clusterParamMarginal(Hg, mg, lam, gh, sizes = sz,
                                          rule = rule, level = level)$value,
       lambda = lam)
}


## Joint agglomerative clustering across ALL candidate parameters at once.
##
## Per subject the linearisation supplies a joint data mode `M[i, ]` (K-vector) and
## a joint precision `Omega[[i]]` (K x K, the data Hessian block over the K
## candidates). A clustering-combination C (a partition of subjects per parameter)
## is scored by the JOINT cost
##   score(C) = D(C) + log(N) * sum_k G_k,
## where D(C) = 0.5 sum_i (A_i(C) w_hat - M_i)^T Omega_i (A_i w_hat - M_i) is the
## cross-parameter data cost at the constrained generalised-least-squares group
## values w_hat = (sum_i A_i^T Omega_i A_i)^{-1} sum_i A_i^T Omega_i M_i, and A_i(C)
## maps the group values to subject i's K deviations. Because D uses the FULL
## Omega_i, a subject whose deviation only appears to differ from its group because
## of cross-attribution with the OTHER candidates is explained by those instead and
## is fused; the per-parameter selector cannot see this. Minimised by greedy
## agglomeration from the all-separate partition (merge the (parameter, group-pair)
## that most lowers the score, until none does). O(n^2 K) merges, each a small
## dense solve on the LINEARISED data (no ODE), cheap.
.clusterJoint <- function(M, Omega, subjects, N) {
  n <- nrow(M); K <- ncol(M)
  score <- function(cl) {
    gid  <- lapply(cl, function(z) { u <- sort(unique(z)); match(z, u) })
    offs <- c(0L, cumsum(vapply(gid, max, 0L)))
    P    <- offs[K + 1L]
    Amats <- vector("list", n)
    Hred  <- matrix(0, P, P); b <- numeric(P)
    for (i in seq_len(n)) {
      Ai <- matrix(0, K, P)
      for (k in seq_len(K)) Ai[k, offs[k] + gid[[k]][i]] <- 1
      Amats[[i]] <- Ai
      Hred <- Hred + crossprod(Ai, Omega[[i]] %*% Ai)
      b    <- b + crossprod(Ai, Omega[[i]] %*% M[i, ])
    }
    w <- solve(Hred + diag(1e-8, P), b)
    D <- 0
    for (i in seq_len(n)) {
      r <- as.numeric(Amats[[i]] %*% w) - M[i, ]
      D <- D + 0.5 * as.numeric(crossprod(r, Omega[[i]] %*% r))
    }
    D + log(N) * P
  }
  cl <- lapply(seq_len(K), function(k) seq_len(n))     # all-separate
  cur <- score(cl); trace <- list()
  repeat {
    best <- NULL
    for (k in seq_len(K)) {
      grps <- sort(unique(cl[[k]])); if (length(grps) < 2L) next
      for (a in seq_along(grps)) for (bb in seq_len(a - 1L)) {
        cl2 <- cl; cl2[[k]][cl2[[k]] == grps[a]] <- grps[bb]
        sc <- score(cl2)
        if (is.null(best) || sc < best$sc) best <- list(sc = sc, cl = cl2)
      }
    }
    if (is.null(best) || best$sc >= cur - 1e-9) break
    cl <- best$cl; cur <- best$sc
    trace[[length(trace) + 1L]] <- list(G = sum(vapply(cl, function(z) length(unique(z)), 0L)),
                                        score = cur)
  }
  clusters <- lapply(seq_len(K), function(k) {
    u <- sort(unique(cl[[k]])); lapply(u, function(gg) sort(subjects[cl[[k]] == gg]))
  })
  list(clusters = clusters, score = cur, trace = trace)
}


#' Recover the grouping of individuals for a clustered regularised NLME fit
#'
#' @description
#' The clustered counterpart of [sparsify]. For a `penaltyL1(method =
#' "clustered")` model it recovers, per candidate parameter, how the individuals
#' group into clusters that share a value, by **exact-marginal model selection**
#' over the cluster path: a single all-separate structural fit provides the FOCE
#' linearisation, the complete-graph fused lasso over a `lambda` grid enumerates the O(n)
#' nested groupings, and each grouping's exact clustered marginal (the multi-
#' dimensional normal-Laplace, integrated by adaptive Gauss-Hermite quadrature) is
#' compared. The grouping with the smallest marginal `-2 log L` wins: the marginal
#' rejects over-fusion (the data cost of forcing unequal individuals together) and
#' over-splitting (the Occam cost of extra random-effect dimensions). This is the
#' consistent selector; reading a hard grouping off a point estimate is not (see
#' `notes/laplace_nlme_theory.Rmd`, "estimation is not selection").
#'
#' @param obj Composed objective `normL2(...) + constraintL1(groupPenalty(...,
#'   kind = "clustered"))`.
#' @param init Named numeric start (structural + error + `lambda`); see [emInit].
#'   The structural start is held fixed (only the deviation starts are perturbed
#'   across `fits`), so it is an `init`, not a perturbation `center` as in
#'   [sparsify].
#' @param fixed Optional named numeric held fixed.
#' @param control List forwarded to the structural fit (`cm1`, `maxOuter`, ...).
#' @param quadNodes Node budget per anchored dimension for the Gauss-Hermite
#'   quadrature: the fixed node count for `quadRule = "tensor"`, and the upper
#'   cap for the dimension-adaptive `"auto"` rule (which shrinks the count as the
#'   anchored dimension grows). Default 24.
#' @param quadRule Quadrature for the clustered marginal: `"auto"` (the default:
#'   dense tensor when the anchored dimension is small, sparse Smolyak grid once
#'   it grows, so larger clusters stay feasible), `"tensor"`, or `"sparse"`.
#' @param level Smolyak depth for the sparse rule (default: anchored dimension
#'   plus 3). Ignored by the tensor rule.
#' @param lamGrid `lambda` grid that generates the candidate groupings (the
#'   cluster path); wide enough to span all-separate to all-fused.
#' @param fits,cores,sd Multistart of the all-separate structural fit: `fits`
#'   perturbed deviation starts (Normal `sd` on the log scale) run over `cores`,
#'   keeping the lowest-objective fit. The joint deviation mode is non-convex, so
#'   with several correlated candidate parameters a single start can settle in a
#'   local optimum that mis-attributes variation and mars the clustering; a few
#'   starts give a cleaner attribution. Default `fits = 1` (single start).
#' @param refine Integer; rounds of iterative refinement (block coordinate
#'   descent). After the first clustering, each candidate is re-linearised with
#'   every OTHER candidate fused to its current grouping, so a per-subject mode no
#'   longer absorbs cross-attribution from the freely-varying others; then it is
#'   re-clustered. Repeats until the groupings stabilise. Default `2`; set `0` for
#'   the single-pass profile linearisation only. Ignored when `joint = TRUE`.
#' @param joint Logical (default `TRUE`); with more than one candidate parameter,
#'   select the groupings JOINTLY by greedy agglomeration on a score that carries
#'   the full cross-parameter covariance (the per-subject joint Hessian), so a
#'   subject that only appears individual-specific in one parameter because of
#'   cross-attribution with the others is fused rather than split. This is the
#'   accurate multi-parameter selector; the per-parameter marginal (with `profile`
#'   precision + `refine`) is used when `FALSE` or for a single candidate.
#' @param verbose Logical; print the recovered grouping.
#' @return An object of class `c("sparsify", "list")` with `type = "grouping"`,
#'   `perParam` (a list, per candidate parameter, of the selected `clusters`, `G`,
#'   and the full scored `chain`), the fitted `structural` parameters, the
#'   per-subject deviation modes `etaModes` (subjects x candidates, as in [EM]),
#'   and `subjects`. A [plot][plot.sparsify] method visualises the recovered
#'   grouping and the score chain.
#' @seealso [sparsify], [EM], [penaltyL1]
#' @noRd
.selectCluster <- function(obj, init, fixed = NULL, control = list(),
                             quadNodes = 24L,
                             quadRule = c("auto", "tensor", "sparse"),
                             level = NULL,
                             lamGrid = 10^seq(-3, 6, length.out = 80L),
                             fits = 1L, cores = 1L, sd = 0.5, refine = 2L,
                             joint = TRUE, verbose = TRUE) {
  quadRule <- match.arg(quadRule)
  rec     <- .laplaceReconstruct(obj)
  penalty <- rec$penalty
  if (!all(vapply(penalty$blocks, `[[`, "", "kind") == "clustered"))
    stop("sparsify: needs an all-clustered penalty ",
         "(penaltyL1(method = 'clustered')).", call. = FALSE)
  se <- penalty$subjectEtas; subjects <- rownames(se); eta_cols <- colnames(se)
  params <- penalty$params; lam_name <- penalty$lambdaName
  if (!lam_name %in% names(init))
    stop("sparsify: `init` is missing the penalty parameter '", lam_name,
         "' (assemble it with emInit()).", call. = FALSE)
  all_conditions <- names(rec$data)
  resObjList <- stats::setNames(lapply(all_conditions, function(s) {
    d <- rec$data[s]; class(d) <- class(rec$data)
    normL2(d, rec$prdfn, errmodel = rec$errfn)
  }), all_conditions)
  free <- init[setdiff(names(init), names(fixed))]

  ## 1. all-separate structural fit (lambda -> 0) for the FOCE linearisation.
  ## No fusion is forced here, so sigma is well estimated (no fusion feedback).
  ## MULTISTART (fits > 1): the joint deviation mode is non-convex (the nonlinear
  ## ODE + the L1 fusion), so with many candidate parameters the single fit can
  ## settle in a local optimum that mis-attributes variation across correlated
  ## parameters and mars the clustering. Run from `fits` perturbed eta starts and
  ## keep the lowest-objective fit: a cleaner attribution, cleaner modes.
  fits <- as.integer(fits)
  all_eta <- as.vector(se)
  free0   <- free; free0[[lam_name]] <- 1e-6
  ctrl0   <- modifyList(control, list(estimateLambda = FALSE))
  fitOne  <- function(i) {
    e0 <- if (i == 1L) NULL
          else stats::setNames(stats::rnorm(length(all_eta), 0, sd), all_eta)
    try(.clusterMap(rec, resObjList, free0, fixed, lam_name, ctrl0,
                          verbose = FALSE, eta_init = e0), silent = TRUE)
  }
  if (verbose)
    cat(sprintf("sparsify: all-separate structural fit (%d start%s) ...\n",
                fits, if (fits > 1L) "s" else ""))
  fitsList <- if (cores > 1L && fits > 1L)
                .parallelLapply(seq_len(fits), fitOne, cores = .splitCores(cores, "fits")$outer,
                  coresConditions = .splitCores(cores, "fits")$conditions)
              else lapply(seq_len(fits), fitOne)
  okFits <- Filter(function(f) !inherits(f, "try-error") && is.finite(f$value), fitsList)
  if (!length(okFits)) stop("sparsify: all structural fits failed.",
                            call. = FALSE)
  bf <- okFits[[which.min(vapply(okFits, function(f) f$value, 0.0))]]
  struct_hat <- bf$argument[setdiff(names(bf$argument), lam_name)]
  eta_warm   <- stats::setNames(as.vector(bf$etaModes), as.vector(se))

  ## 2. linearise at the fitted structural: per-(subject, parameter) Hd, m.
  sol <- .clusterSolve(resObjList, c(struct_hat, fixed), penalty, 1e-6,
                          eta0 = eta_warm, maxit = 60L, tol = 1e-9)
  Eta_cur <- sol$etahat

  ## per-parameter grouping selection from a linearisation's (m, Hd).
  gh <- .gaussHermite(as.integer(quadNodes))
  clusterFrom <- function(mMat, HdMat) {
    pp <- lapply(seq_along(eta_cols), function(k) {
      Hk <- HdMat[subjects, k]; mk <- mMat[subjects, k]
      paths <- unique(lapply(lamGrid, function(lam)
        lapply(solveFusedComplete(stats::setNames(mk, subjects), Hk, lam)$clusters, sort)))
      sc <- vapply(paths, function(P) .clusterGroupScore(
        lapply(P, function(gg) match(gg, subjects)), Hk, mk, gh,
        rule = quadRule, level = level)$value, 0.0)
      jstar <- which.min(sc)
      chain <- data.frame(
        grouping = vapply(paths, function(P)
          paste(vapply(P, paste, "", collapse = "="), collapse = " | "), ""),
        size = vapply(paths, length, 0L), score = sc,
        selected = seq_along(paths) == jstar, stringsAsFactors = FALSE)
      list(param = params[k], clusters = paths[[jstar]], G = length(paths[[jstar]]),
           chain = chain[order(sc), , drop = FALSE])
    })
    names(pp) <- eta_cols; pp
  }
  perParam <- clusterFrom(sol$m, sol$Hd)

  if (joint && length(eta_cols) > 1L) {
    ## 3a. JOINT agglomerative clustering: score groupings with the full
    ## cross-parameter covariance (per-subject K x K joint Hessian) so a subject
    ## that only differs because of cross-attribution is fused, not split. This is
    ## the multi-parameter selector; the per-parameter chains above are kept for
    ## inspection. (For K = 1 there is no cross-parameter effect, so the single
    ## per-parameter marginal is used unchanged.)
    Ndata <- sum(vapply(names(rec$data),
      function(cc) nrow(as.data.frame(rec$data[[cc]])), 0L))
    base_eta <- stats::setNames(rep(0, length(as.vector(se))), as.vector(se))
    fulls <- lapply(subjects, function(s) {
      full <- c(struct_hat, fixed, base_eta)
      full[as.character(se[s, ])] <- sol$etahat[s, ]
      full
    })
    hs <- .objEvalMany(resObjList[subjects], fulls, deriv = TRUE)
    Omega <- lapply(seq_along(subjects), function(si) {
      en <- as.character(se[subjects[si], ])
      hs[[si]]$hessian[en, en, drop = FALSE]
    })
    jc <- .clusterJoint(sol$m[subjects, , drop = FALSE], Omega, subjects, Ndata)
    perParam <- stats::setNames(lapply(seq_along(eta_cols), function(k)
      list(param = params[k], clusters = jc$clusters[[k]],
           G = length(jc$clusters[[k]]), chain = perParam[[k]]$chain)), eta_cols)
    if (verbose) cat("sparsify: joint agglomerative clustering\n")
  } else if (as.integer(refine) > 0L) {
    ## 3b. Iterative refinement (per-parameter path): re-linearise with every OTHER
    ## candidate fused to its current clustering, so its per-subject freedom no
    ## longer leaks cross-attribution into the parameter being clustered; re-cluster.
    ## Repeat until the groupings stop changing.
    sig <- function(pp) paste(vapply(pp, function(z)
      paste(vapply(lapply(z$clusters, sort), paste, "", collapse = ","), collapse = "|"), ""),
      collapse = ";")
    prevSig <- sig(perParam)
    for (r in seq_len(as.integer(refine))) {
      Eta_cl <- Eta_cur
      for (k in seq_along(eta_cols))
        for (cl in perParam[[k]]$clusters) {
          idx <- match(cl, subjects); Eta_cl[idx, k] <- mean(Eta_cur[idx, k])
        }
      lp <- .clusterProfileLinearise(resObjList, c(struct_hat, fixed), se, Eta_cl)
      ppNew <- clusterFrom(lp$m, lp$Hd)
      Eta_cur <- Eta_cl; perParam <- ppNew
      newSig <- sig(ppNew)
      if (verbose) cat(sprintf("sparsify: refinement round %d\n", r))
      if (identical(newSig, prevSig)) break
      prevSig <- newSig
    }
  }

  out <- list(type = "grouping", perParam = perParam, structural = struct_hat,
              etaModes = sol$m, subjects = subjects, penalty = penalty)
  class(out) <- c("sparsify", "list")
  if (verbose) print(out)
  out
}


#' Print a sparsify (structural model selection) result
#'
#' @param x A `sparsify` object from [sparsify].
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.sparsify <- function(x, ...) {
  if (identical(x$type, "grouping")) .print_sparsify_grouping(x, ...)
  else .print_sparsify_support(x, ...)
  invisible(x)
}

## Grouping branch of print.sparsify.
.print_sparsify_grouping <- function(x, ...) {
  cat("sparsify (grouping; exact-marginal)\n")
  cat("  structural :\n"); print(round(x$structural, 4))
  for (pp in x$perParam)
    cat(sprintf("  %-6s -> {%s}  (G = %d)\n", pp$param,
                paste(vapply(pp$clusters, paste, "", collapse = ","),
                      collapse = "} {"), pp$G))
  invisible(x)
}


#' Diagnostic plots for a clustered sparsify result
#'
#' @description
#' Visualises a clustered (grouping) [sparsify] result. `type = "grouping"`
#' (default) draws, per candidate parameter, each subject's deviation mode
#' \eqn{\hat\eta} with the recovered clusters coloured and their shared centroids
#' drawn as a line: the recovered group structure at a glance. `type = "chain"`
#' draws the per-parameter cluster-path score chain (marginal \eqn{-2\log L}
#' against the number of clusters \eqn{G}), with its marginal-optimal grouping
#' highlighted, so the Occam trade-off is explicit; under the default joint
#' selection the final grouping (the "grouping" panel) may refine this
#' per-parameter one.
#'
#' @param x A clustered (grouping) [sparsify] object.
#' @param type One of `"grouping"` (recovered clusters, default) or `"chain"`
#'   (the per-parameter scored cluster path).
#' @param ... Ignored.
#' @return A ggplot object.
#' @seealso [sparsify]
#' @export
plot.sparsify <- function(x, type = c("grouping", "chain"), ...) {
  if (!identical(x$type, "grouping"))
    stop("plot() is available for clustered (grouping) sparsify results; ",
         "a support selection exposes its `chain` data frame directly.",
         call. = FALSE)
  type     <- match.arg(type)
  subjects <- x$subjects
  eta_cols <- names(x$perParam)

  if (identical(type, "grouping")) {
    if (is.null(x$etaModes))
      stop("plot(type = 'grouping') needs `etaModes`; refit with sparsify().",
           call. = FALSE)
    rows <- lapply(eta_cols, function(ec) {
      pp    <- x$perParam[[ec]]
      cl_of <- stats::setNames(rep(NA_integer_, length(subjects)), subjects)
      for (g in seq_along(pp$clusters)) cl_of[pp$clusters[[g]]] <- g
      mode  <- x$etaModes[subjects, ec]
      cent  <- vapply(subjects, function(s) mean(mode[subjects[cl_of == cl_of[[s]]]]), 0.0)
      data.frame(parameter = pp$param, subject = subjects, mode = as.numeric(mode),
                 cluster = factor(cl_of), centroid = as.numeric(cent),
                 stringsAsFactors = FALSE)
    })
    df <- do.call(rbind, rows)
    df$subject <- factor(df$subject, levels = subjects)
    return(
      ggplot2::ggplot(df, ggplot2::aes(x = subject, y = mode, color = cluster)) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
        ggplot2::geom_line(ggplot2::aes(y = centroid, group = cluster)) +
        ggplot2::geom_point(size = 2.4) +
        ggplot2::facet_wrap(~ parameter, scales = "free_y") +
        scale_color_dMod() +
        ggplot2::labs(x = "subject", y = "eta_hat (deviation)", color = "cluster",
                      title = "sparsify: recovered grouping of individuals") +
        theme_dMod(base_size = 11))
  }

  ## type == "chain": the scored cluster path per parameter
  rows <- lapply(eta_cols, function(ec) {
    ch <- x$perParam[[ec]]$chain
    data.frame(parameter = x$perParam[[ec]]$param, G = ch$size, score = ch$score,
               selected = ch$selected, stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  ggplot2::ggplot(df, ggplot2::aes(x = G, y = score)) +
    ggplot2::geom_line(color = "grey60") +
    ggplot2::geom_point(ggplot2::aes(color = selected), size = 2.6) +
    ggplot2::facet_wrap(~ parameter, scales = "free") +
    ggplot2::scale_color_manual(values = c(`FALSE` = "grey40", `TRUE` = "#c0392b"),
                                labels = c(`FALSE` = "candidate", `TRUE` = "selected"),
                                name = NULL) +
    ggplot2::labs(x = "number of clusters G", y = "marginal -2 log L",
                  title = "sparsify: cluster-path score chain") +
    theme_dMod(base_size = 11)
}


## ---- SAEM/MCMC exact cross-check of the clustered marginal (P5) --------------
##
## A joint random-walk Metropolis sampler of the (anchored, FOCE-linearised)
## clustered posterior for one parameter, the coupled analogue of the SAEM
## E-step (prior-agnostic, value-only, no subgradient for the L1 kink). It targets
## the SAME posterior exp(-phi(z)) that .clusterParamMarginal integrates by
## quadrature, so their posterior moments (E[eta], E||D eta||_1) must agree: an
## independent MCMC check that the deterministic exact marginal is correct
## (`sizes` gives the weighted reduced/grouped model, as in .clusterParamMarginal).
## Returns list(Eeta, Efus, accept).
.clusterMHmoments <- function(H, m, lambda, sizes = rep(1, length(H)),
                              nsamp = 40000L, burn = 5000L, step = NULL) {
  G <- length(H); B <- MASS::Null(matrix(sizes, G, 1)); d <- G - 1L
  Wp <- outer(sizes, sizes); ut <- upper.tri(Wp)
  Dp  <- function(w) { dd <- as.matrix(stats::dist(w)); sum(dd[ut] * Wp[ut]) }
  phi <- function(z) { e <- as.numeric(B %*% z)
    0.25 * sum(H * (e - m)^2) + 0.5 * lambda * Dp(e) }
  Az  <- crossprod(B, H * B)
  if (is.null(step)) step <- 0.6 * sqrt(diag(2 * solve(Az)))
  z   <- as.numeric(solve(Az, crossprod(B, H * m)))  # start at the data-mode projection
  cur <- phi(z); Esum <- numeric(G); Efus <- 0; nkeep <- 0L; acc <- 0L
  for (it in seq_len(nsamp)) {
    zp <- z + stats::rnorm(d, 0, step); pp <- phi(zp)
    if (log(stats::runif(1)) < cur - pp) { z <- zp; cur <- pp; acc <- acc + 1L }
    if (it > burn) { e <- as.numeric(B %*% z)
      Esum <- Esum + e; Efus <- Efus + Dp(e); nkeep <- nkeep + 1L }
  }
  list(Eeta = Esum / nkeep, Efus = Efus / nkeep, accept = acc / nsamp)
}


## Full SAEM fit for the clustered model on the TRUE (nonlinear) ODE posterior.
##
## The stochastic-approximation EM cross-check of the FOCE/exact-marginal path:
## unlike .selectCluster (which linearises the data at the mode), SAEM samples
## the true coupled posterior, so agreement of its structural / posterior-mean
## estimates validates the FOCE linearisation itself. E-step: per-subject-row
## random-walk Metropolis over all subjects' etas, coupled by the complete-graph
## penalty (value-only, no subgradient), with an incremental log-density (only the
## proposed subject's data + its fusion terms change). CM-2: lambda = 2 rho /
## E||D eta|| (rho = K(n-1)); CM-1: trust() on the complete-data objective at the
## drawn etas; anchoring by re-centring each iteration. Reports the posterior mean
## etaModes (which reveals the group structure); the hard grouping comes from the
## marginal selector .selectCluster, not from a point estimate.
.clusterSaem <- function(rec, resObjList, free, fixed, lam_name,
                               control, verbose) {
  penalty <- rec$penalty
  if (!all(vapply(penalty$blocks, `[[`, "", "kind") == "clustered"))
    stop(".fitLaplace: clustered SAEM needs an all-clustered penalty.", call. = FALSE)
  se <- penalty$subjectEtas; subjects <- rownames(se); eta_cols <- colnames(se)
  params <- penalty$params; n <- length(subjects); K <- length(eta_cols)
  all_eta <- as.vector(se); base_eta <- stats::setNames(rep(0, length(all_eta)), all_eta)
  struct_names   <- setdiff(names(free), lam_name)
  ref_conditions <- setdiff(names(rec$data), subjects)
  rho    <- K * (n - 1L)
  mu_of  <- stats::setNames(paste0("log_", params), eta_cols)

  sc <- modifyList(list(nBurnin = 150L, nEM = 150L, nMcmc = 8L, stepsize = 0.15,
                        cm1 = list()),
                   if (is.null(control$saem)) list() else control$saem)
  cm1 <- .trustControl(list(rinit = 1, rmax = 10, iterlim = 30L,
                            ftol = 1e-6, mtol = 1e-6),
                       sc$cm1, label = "control$saem$cm1")
  nBurnin <- as.integer(sc$nBurnin); nEM <- as.integer(sc$nEM); nMcmc <- as.integer(sc$nMcmc)

  struct_cur <- free[struct_names]; lambda_cur <- as.numeric(free[[lam_name]])
  fullFor <- function(os, s, er) { f <- c(os, base_eta); f[as.character(se[s, ])] <- er; f }
  fus_of  <- function(E) sum(vapply(seq_len(K), function(k) sum(stats::dist(E[, k])), 0.0))
  cm1_obj <- function(th, E) {
    os <- c(stats::setNames(th, struct_names), fixed); tot <- 0
    gg <- stats::setNames(numeric(length(struct_names)), struct_names)
    HH <- matrix(0, length(struct_names), length(struct_names),
                 dimnames = list(struct_names, struct_names))
    conds <- c(subjects, ref_conditions)
    got <- .objEvalMany(resObjList[conds],
                        c(lapply(seq_len(n), function(si) fullFor(os, subjects[si], E[si, ])),
                          rep(list(c(os, base_eta)), length(ref_conditions))),
                        deriv = TRUE)
    for (si in seq_len(n)) {
      o <- got[[si]]
      tot <- tot + o$value; gg <- gg + o$gradient[struct_names]
      HH <- HH + o$hessian[struct_names, struct_names, drop = FALSE]
    }
    for (ri in seq_along(ref_conditions)) {
      o <- got[[n + ri]]; sn <- intersect(struct_names, names(o$gradient))
      tot <- tot + o$value; gg[sn] <- gg[sn] + o$gradient[sn]; HH[sn, sn] <- HH[sn, sn] + o$hessian[sn, sn, drop = FALSE]
    }
    objlist(value = tot, gradient = gg, hessian = HH)
  }

  Eta  <- matrix(0, n, K, dimnames = list(subjects, eta_cols))
  os0  <- c(struct_cur, fixed)
  dvc  <- vapply(.objEvalMany(resObjList[subjects],
                              lapply(seq_len(n), function(si) fullFor(os0, subjects[si], Eta[si, ])),
                              deriv = FALSE), `[[`, 0.0, "value")
  names(dvc) <- subjects
  step <- rep(sc$stepsize, n)
  S_fus <- NULL; EtaMean <- matrix(0, n, K, dimnames = list(subjects, eta_cols)); nAcc <- 0L
  rows <- list(); prev_par <- c(struct_cur, stats::setNames(lambda_cur, lam_name))
  for (it in seq_len(nBurnin + nEM)) {
    gamma <- if (it <= nBurnin) 1.0 else 1.0 / (it - nBurnin)
    os <- c(struct_cur, fixed); acc <- integer(n)
    for (mm in seq_len(nMcmc)) for (si in seq_len(n)) {
      er  <- Eta[si, ] + stats::rnorm(K, 0, step[si])
      dvp <- resObjList[[subjects[si]]](fullFor(os, subjects[si], er), deriv = FALSE)$value
      Ep  <- Eta; Ep[si, ] <- er
      dlp <- -0.5 * ((dvp - dvc[si]) + lambda_cur * (fus_of(Ep) - fus_of(Eta)))
      if (log(stats::runif(1)) < dlp) { Eta <- Ep; dvc[si] <- dvp; acc[si] <- acc[si] + 1L }
    }
    ## smooth stepsize adaptation toward a 0.3 acceptance target
    a <- acc / nMcmc; step <- pmin(pmax(step * exp(0.4 * (a - 0.3)), 0.02), 5)
    ## anchor: recentre each parameter column, fold the mean into log_<par>
    for (k in seq_len(K)) { c_k <- mean(Eta[, k]); Eta[, k] <- Eta[, k] - c_k
      mn <- mu_of[[eta_cols[k]]]; if (mn %in% names(struct_cur)) struct_cur[[mn]] <- struct_cur[[mn]] + c_k }
    os <- c(struct_cur, fixed)
    dvc <- vapply(subjects, function(s) resObjList[[s]](fullFor(os, s, Eta[match(s, subjects), ]), deriv = FALSE)$value, 0.0)
    ## CM-2: closed-form lambda (2 rho / E||D eta||)
    fc <- fus_of(Eta); S_fus <- if (is.null(S_fus)) fc else S_fus + gamma * (fc - S_fus)
    lambda_cur <- if (S_fus > 0) 2 * rho / S_fus else lambda_cur
    ## CM-1: structural M-step at the drawn etas, RM-damped in convergence
    cf <- suppressMessages(do.call(trust, modifyList(cm1, list(
      objfun = function(th, ...) cm1_obj(th, Eta), parinit = struct_cur))))
    struct_cur <- if (it <= nBurnin) cf$argument else struct_cur + gamma * (cf$argument - struct_cur)
    if (it > nBurnin) { EtaMean <- EtaMean + Eta; nAcc <- nAcc + 1L }
    cur_par <- c(struct_cur, stats::setNames(lambda_cur, lam_name))
    rows[[length(rows) + 1L]] <- data.frame(iter = it, lambda = lambda_cur,
      deltaPar = max(abs(cur_par - prev_par)))
    if (verbose && (it %% 50L == 0L || it == nBurnin + nEM))
      cat(sprintf("  EM(saem/clustered) iter %d lambda=%.4g\n", it, lambda_cur))
    prev_par <- cur_par
  }
  EtaMean <- EtaMean / nAcc
  for (k in seq_len(K)) EtaMean[, k] <- EtaMean[, k] - mean(EtaMean[, k])
  theta_hat <- c(struct_cur, stats::setNames(lambda_cur, lam_name))
  val <- sum(vapply(subjects, function(s)
    resObjList[[s]](fullFor(c(struct_cur, fixed), s, EtaMean[match(s, subjects), ]), deriv = FALSE)$value, 0.0))
  out <- list(argument = theta_hat, value = val,
              valueType = "data -2logL at the posterior mean (not a marginal)",
              lambda = lambda_cur, etaModes = EtaMean, clusters = NULL, clustered = TRUE,
              converged = nEM > 0L && mean(utils::tail(vapply(rows, `[[`, 0.0, "deltaPar"), 20L)) < 5e-2,
              iterations = nBurnin + nEM, stageTrace = do.call(rbind, rows),
              method = "saem", prior = "penaltyL1", penalty = penalty,
              prdfn = rec$prdfn, data = rec$data, errfn = rec$errfn)
  class(out) <- c("em", "list")
  out
}
