# Expectation consistent (EC) approximation of the integral of one group over spike-and-slab random effects with a
# given spike and slab mean, in the coordinate u in which the spike is 0 (methods article Sparsity priors, section
# Spike at a given value).


## Gauss-Hermite nodes and weights for the weight exp(-x^2), `n` nodes (Golub-Welsch).
.gaussHermite <- function(n) {
  if (n == 1L) return(list(x = 0, w = sqrt(pi)))
  J <- matrix(0, n, n)
  off <- sqrt(seq_len(n - 1L) / 2)
  J[cbind(seq_len(n - 1L), 2:n)] <- off
  J[cbind(2:n, seq_len(n - 1L))] <- off
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = sqrt(pi) * e$vectors[1, ]^2)
}

## The tilted distribution of one site, the Gaussian cavity N(u; mc, vc) times the slab part of the
## prior with weight `prob`: log of the slab mass, E[u], E[u^2] and E[(eta - mean)^2] under the slab;
## `lb` = log b for the log-normal slab of ln u = lb eta, NA for a normal slab in u = eta - spike.
.ecSlabMoments <- function(mc, vc, mean, tau, lb, spike, gh) {
  if (is.na(lb)) {
    # u = eta - spike, slab N(mean - spike, tau^2): the product of two Gaussians in closed form
    ms <- mean - spike
    v <- 1 / (1 / vc + 1 / tau^2)
    m <- v * (mc / vc + ms / tau^2)
    logZ <- stats::dnorm(mc, ms, sqrt(vc + tau^2), log = TRUE)
    return(list(logZ = logZ, m1 = m, m2 = v + m^2, dev2 = v + (m + spike - mean)^2))
  }
  # in w = ln u: log N(e^w; mc, vc) + log N(w; lb mean, (lb tau)^2), adaptive Gauss-Hermite at the mode
  mw <- lb * mean; sw2 <- (lb * tau)^2
  h <- function(w) -(exp(w) - mc)^2 / (2 * vc) - (w - mw)^2 / (2 * sw2)
  dh <- function(w) -(exp(w) - mc) * exp(w) / vc - (w - mw) / sw2
  d2h <- function(w) -(2 * exp(2 * w) - mc * exp(w)) / vc - 1 / sw2
  w <- if (mc > 0) min(log(mc), mw + 6 * sqrt(sw2)) else mw
  for (i in seq_len(50L)) {
    c2 <- d2h(w)
    step <- if (c2 < 0) -dh(w) / c2 else sign(dh(w)) * sqrt(sw2)
    step <- max(min(step, 2), -2)
    w <- w + step
    if (abs(step) < 1e-10) break
  }
  s <- 1 / sqrt(max(-d2h(w), 1e-12))
  nodes <- w + sqrt(2) * s * gh$x
  lf <- h(nodes) + log(gh$w) + gh$x^2 + log(sqrt(2) * s) - 0.5 * log(2 * pi * vc) - 0.5 * log(2 * pi * sw2)
  top <- max(lf)
  q <- exp(lf - top)
  Z <- sum(q)
  u <- exp(nodes)
  list(logZ = top + log(Z), m1 = sum(q * u) / Z, m2 = sum(q * u^2) / Z,
       dev2 = sum(q * (nodes / lb - mean)^2) / Z)
}

## EC over the slab candidates `cand` of group `g`, the normal etas `nrm` exact in the Gaussian, expanded first at
## `start`: -2 log I, inclusion probabilities, E[(eta - mean)^2 | slab], E[eta eta'] of the normal etas and the
## subset of inclusion above 1/2; NULL where EC or the relinearization does not converge.
.spikeECGroup <- function(design, objs, g, theta, fixed, etas, cand, nrm, pr, alpha, prob, lb, start,
                          control = list()) {
  conds <- design$conditions[[g]]
  J <- length(cand)
  co <- c(cand, nrm)
  N <- length(co)
  mu <- pr$mu; tau <- pr$tau
  lbc <- c(lb[cand], setNames(rep(NA_real_, length(nrm)), nrm))
  gh <- .gaussHermite(control$nodes %||% 40L)
  # u = 0 is the spike of a slab coordinate, a normal coordinate is the eta itself
  rate <- !is.na(lbc)
  etaOf <- function(u) {
    out <- u + alpha[co]
    out[rate] <- log(u[rate]) / lbc[rate]
    out
  }
  floorU <- c(ifelse(is.na(lb[cand]), -Inf, exp(lb[cand] * (mu[cand] - 4 * tau[cand]))), rep(-Inf, length(nrm)))
  um <- ifelse(is.na(lbc), start[co] - alpha[co], exp(lbc * start[co]))
  lam <- gam <- numeric(N)
  # the normal prior: precision Omega^-1 in the Gaussian, the normalizer of the prior in the constant
  P0 <- matrix(0, N, N)
  if (length(nrm)) P0[J + seq_along(nrm), J + seq_along(nrm)] <- pr$Oinv
  c0prior <- if (length(nrm)) -0.5 * pr$logdetO - length(nrm) / 2 * log(2 * pi) else 0
  relinearized <- FALSE
  for (outer in seq_len(control$maxRelinearize %||% 30L)) {
    um <- pmax(um, floorU)
    eta <- setNames(etaOf(um), co)
    out <- Reduce(`+`, .objEvalMany(objs[conds], rep(list(eta), length(conds)), deriv = TRUE,
                                    fixed = c(theta, fixed, alpha[.setdiff(etas, co)])))
    if (!is.finite(out$value)) return(NULL)
    D <- ifelse(rate, 1 / (um * lbc), 1)
    gu <- D * out$gradient[co]
    # Gauss-Newton in u: the residual Jacobian times D
    Hu <- (D * as.matrix(out$hessian)[co, co]) * rep(D, each = N)
    A <- Hu / 2 + P0
    a <- drop((Hu / 2) %*% um) - gu / 2
    c0 <- -out$value / 2 + sum(gu * um) / 2 - drop(um %*% Hu %*% um) / 4 + c0prior
    if (outer == 1L) {
      # sites from the prior moments of every slab coordinate
      for (j in seq_len(J)) {
        s <- .ecSlabMoments(0, 1e12, mu[cand[j]], tau[cand[j]], lb[cand[j]], alpha[cand[j]], gh)
        m1 <- prob[cand[j]] * s$m1; m2 <- prob[cand[j]] * s$m2
        lam[j] <- 1 / max(m2 - m1^2, 1e-12); gam[j] <- m1 * lam[j]
      }
    }
    ok <- FALSE
    prevM <- prevV <- NULL
    for (sweep in seq_len(control$maxSweeps %||% 500L)) {
      Sigma <- tryCatch(chol2inv(chol(A + diag(lam, N))), error = function(e) NULL)
      if (is.null(Sigma)) return(NULL)
      m <- drop(Sigma %*% (a + gam))
      sd <- sqrt(diag(Sigma))
      # converged when the marginal means move by less than tolSweep standard deviations and the
      # standard deviations by less than tolSweep relative
      if (!is.null(prevM) && max(abs(m - prevM) / sd, abs(sd - prevV) / sd) < (control$tolSweep %||% 1e-5)) {
        ok <- TRUE; break
      }
      prevM <- m; prevV <- sd
      # one site at a time, each against the posterior after the update before
      for (j in seq_len(J)) {
        if (j > 1L) {
          Sigma <- tryCatch(chol2inv(chol(A + diag(lam, N))), error = function(e) NULL)
          if (is.null(Sigma)) return(NULL)
          m <- drop(Sigma %*% (a + gam))
        }
        v <- Sigma[j, j]
        pc <- 1 / v - lam[j]
        if (pc <= 0) next
        vc <- 1 / pc
        mc <- vc * (m[j] / v - gam[j])
        e <- cand[j]
        s <- .ecSlabMoments(mc, vc, mu[[e]], tau[[e]], lb[[e]], alpha[[e]], gh)
        lSpike <- log1p(-prob[[e]]) + stats::dnorm(0, mc, sqrt(vc), log = TRUE)
        lSlab <- log(prob[[e]]) + s$logZ
        top <- max(lSpike, lSlab)
        logZj <- top + log(exp(lSpike - top) + exp(lSlab - top))
        pin <- exp(lSlab - logZj)
        mt <- pin * s$m1
        vt <- max(pin * s$m2 - mt^2, 1e-14)
        # a site of negative precision is replaced by the one that matches the tilted mean at the cavity variance
        newLam <- max(1 / vt - 1 / vc, 0)
        newGam <- mt * (1 / vc + newLam) - mc / vc
        lam[j] <- 0.5 * lam[j] + 0.5 * newLam
        gam[j] <- 0.5 * gam[j] + 0.5 * newGam
      }
    }
    if (!ok) return(NULL)
    Sigma <- chol2inv(chol(A + diag(lam, N)))
    m <- drop(Sigma %*% (a + gam))
    # the tilted moments and site constants of every site at the converged posterior
    logC <- incl <- dev2 <- numeric(J)
    for (j in seq_len(J)) {
      v <- Sigma[j, j]
      pc <- 1 / v - lam[j]
      if (pc <= 0) return(NULL)
      vc <- 1 / pc
      mc <- vc * (m[j] / v - gam[j])
      e <- cand[j]
      s <- .ecSlabMoments(mc, vc, mu[[e]], tau[[e]], lb[[e]], alpha[[e]], gh)
      lSpike <- log1p(-prob[[e]]) + stats::dnorm(0, mc, sqrt(vc), log = TRUE)
      lSlab <- log(prob[[e]]) + s$logZ
      top <- max(lSpike, lSlab)
      logZj <- top + log(exp(lSpike - top) + exp(lSlab - top))
      incl[j] <- exp(lSlab - logZj)
      dev2[j] <- s$dev2
      # the site constant: the tilted mass over the mass of cavity times site
      p <- 1 / vc + lam[j]; r <- mc / vc + gam[j]
      logC[j] <- logZj - (-0.5 * log(vc * p) + 0.5 * r^2 / p - 0.5 * mc^2 / vc)
    }
    logI <- c0 + 0.5 * as.numeric(determinant(2 * pi * Sigma)$modulus) + 0.5 * sum((a + gam) * m) + sum(logC)
    step <- max(abs(pmax(m, floorU) - um) / pmax(abs(um), 1e-8))
    um <- m
    if (step < (control$tolRelinearize %||% 1e-4)) { relinearized <- TRUE; break }
  }
  if (!relinearized) return(NULL)
  k <- J + seq_along(nrm)
  normalMoment <- if (length(nrm)) Sigma[k, k, drop = FALSE] + tcrossprod(m[k])
  list(m2logI = -2 * logI, inclusion = setNames(incl, cand), dev2 = setNames(dev2, cand),
       normalMoment = normalMoment, on = cand[incl > 0.5])
}
