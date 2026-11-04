# Spike-and-slab random effects, each 0 or from N(0, tau^2) with probability prob, and
# Dirichlet-process random effects, a partition of the groups under the Ewens prior with one
# N(0, tau^2) value per block but that of the population value; the guide chapter Group-specific effects derives both.


## The integral of a cluster symbol on the expansion (a, m) of its groups: the exact sum over all
## set partitions of the groups and the anchor under the Ewens prior (dirichletProcess_kernel.cpp), its
## Fisher derivatives, posterior moments, M-step statistics and the most probable partition.
.slabPartition <- function(a, m, concentration, tau) {
  r <- dirichlet_process_partitions(as.numeric(a), as.numeric(m), concentration, tau)
  c(r[c("logI", "dm", "da", "mean", "cen2", "spike", "dconcentration", "dtau", "labels", "values")],
    list(stats = c(r$free, length(a), r$slab2)))
}

## The largest number of spike-and-slab coordinates of a group integrated jointly, 2^K subsets.
.slabJoint <- 12L

## The integral of the spike-and-slab coordinates of a group on their full curvature, a sum over
## the subsets F of slab coordinates of closed-form Gaussian integrals at eta = 0 off F; Fisher
## derivatives, posterior moments and the most probable subset.
.slabSubsets <- function(P, m, prob, tau) {
  K <- length(m)
  if (K > .slabJoint)
    stop(sprintf("em: at most %d spike-and-slab coordinates integrated jointly.", .slabJoint), call. = FALSE)
  b <- drop(P %*% m); c0 <- sum(m * b)
  subsets <- unname(as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), K))))
  n <- nrow(subsets)
  logw <- numeric(n)
  means <- matrix(0, n, K)
  second <- array(0, c(n, K, K))
  for (s in seq_len(n)) {
    F <- subsets[s, ]
    logw[s] <- sum(log(prob[F])) + sum(log1p(-prob[!F])) - c0 / 2
    if (!any(F)) next
    R <- chol(P[F, F, drop = FALSE] + diag(1 / tau[F]^2, sum(F)))
    mu <- backsolve(R, forwardsolve(t(R), b[F]))
    logw[s] <- logw[s] + sum(b[F] * mu) / 2 - sum(log(diag(R))) - sum(log(tau[F]))
    means[s, F] <- mu
    second[s, F, F] <- chol2inv(R) + tcrossprod(mu)
  }
  mx <- max(logw)
  logI <- mx + log(sum(exp(logw - mx)))
  w <- exp(logw - logI)
  mean <- drop(w %*% means)
  M2 <- matrix(colSums(w * matrix(second, n)), K, K)
  cen <- M2 - tcrossprod(mean, m) - tcrossprod(m, mean) + tcrossprod(m)
  incl <- drop(w %*% subsets)
  slab2 <- diag(M2)
  pick <- which.max(logw)
  list(logI = logI, dm = drop(P %*% (mean - m)), dP = -cen / 2, mean = mean, cov = M2 - tcrossprod(mean),
       incl = incl, slab2 = slab2, dprob = incl / prob - (1 - incl) / (1 - prob),
       dtau = -incl / tau + slab2 / tau^3, spike = 1 - incl, selected = as.numeric(subsets[pick, ]),
       values = means[pick, ])
}

## log I of the Dirichlet-process stems on the expansion (precision P[[g]], centre m[g, ]) by
## sequential importance sampling over the seating of the groups, each by prior times predictive;
## the product of the normalisers is unbiased for I, the standard error of log I beside it.
.dirichletImportance <- function(P, m, concentration, tau, particles = 500L) {
  G <- nrow(m); K <- ncol(m)
  # log of the integral over the free values with precision Lam (data and prior), linear term h
  logIntegral <- function(Lam, h, logTau) {
    if (!length(h)) return(0)
    R <- chol(Lam)
    -sum(logTau) - sum(log(diag(R))) + sum(backsolve(R, h, transpose = TRUE)^2) / 2
  }
  logZ <- numeric(particles)
  for (pt in seq_len(particles)) {
    sizes <- rep(list(1L), K)
    stemOf <- integer(0); clusterOf <- integer(0)
    Lam <- matrix(0, 0, 0); h <- numeric(0); logTau <- numeric(0); c0 <- 0; current <- 0
    for (g in seq_len(G)) {
      choices <- as.matrix(expand.grid(lapply(sizes, function(n) seq_len(length(n) + 1L))))
      Pg <- P[[g]]; mg <- m[g, ]
      states <- vector("list", nrow(choices))
      logw <- numeric(nrow(choices))
      for (r in seq_len(nrow(choices))) {
        ch <- choices[r, ]
        new <- which(ch > lengths(sizes))
        d <- length(h) + length(new)
        A <- matrix(0, K, d)
        sOf <- c(stemOf, new); cOf <- c(clusterOf, ch[new])
        for (s in seq_len(K)) if (ch[s] > 1L) A[s, which(sOf == s & cOf == ch[s])] <- 1
        L2 <- matrix(0, d, d); L2[seq_along(h), seq_along(h)] <- Lam
        lt <- c(logTau, log(tau[new]))
        diag(L2)[length(h) + seq_along(new)] <- 1 / tau[new]^2
        L2 <- L2 + crossprod(A, Pg %*% A)
        h2 <- c(h, numeric(length(new))) + drop(crossprod(A, Pg %*% mg))
        cc <- c0 - sum(mg * (Pg %*% mg)) / 2
        prior <- sum(vapply(seq_len(K), function(s)
          log(if (ch[s] > length(sizes[[s]])) concentration[s] else sizes[[s]][ch[s]]) - log(g + concentration[s]), 0))
        value <- cc + logIntegral(L2, h2, lt)
        logw[r] <- prior + value - current
        states[[r]] <- list(Lam = L2, h = h2, logTau = lt, c0 = cc, value = value, stemOf = sOf, clusterOf = cOf)
      }
      mx <- max(logw)
      logZ[pt] <- logZ[pt] + mx + log(sum(exp(logw - mx)))
      pick <- sample.int(length(logw), 1L, prob = exp(logw - mx))
      st <- states[[pick]]; ch <- choices[pick, ]
      for (s in seq_len(K)) {
        if (ch[s] > length(sizes[[s]])) sizes[[s]] <- c(sizes[[s]], 1L) else sizes[[s]][ch[s]] <- sizes[[s]][ch[s]] + 1L
      }
      Lam <- st$Lam; h <- st$h; logTau <- st$logTau; c0 <- st$c0; current <- st$value
      stemOf <- st$stemOf; clusterOf <- st$clusterOf
    }
  }
  mx <- max(logZ)
  w <- exp(logZ - mx)
  list(logI = mx + log(mean(w)), se = stats::sd(w) / sqrt(particles) / mean(w))
}
