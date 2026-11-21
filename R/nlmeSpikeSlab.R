# Spike-and-slab random effects, 0 or from N(0, tau^2) with probability prob, and Dirichlet-process ones, a
# partition of the groups under the Ewens prior with one N(0, tau^2) value per block but that of the population
# value; derived in the methods article Sparsity priors.


## The most nodes the joint sum over the subsets of the spike-and-slab coordinates of a group
## visits (src/sparse_kernel.cpp, SlabSum) before the group goes to the expectation consistent
## (EC) approximation.
.slabBudget <- 2^24

## log I of the Dirichlet-process stems on the expansion (precision P[[g]], center m[g, ]) by
## sequential importance sampling over the seatings, with the standard error of log I; the stems of
## `slab` are spike-and-slab ones, summed over subsets, at most `budget` seatings and subsets.
.dirichletImportance <- function(P, m, concentration, tau, particles = 500L, slab = NULL, prob = NULL,
                                 budget = 4096L) {
  G <- nrow(m); K <- ncol(m)
  slab <- slab %||% rep(FALSE, K)
  dpIdx <- which(!slab)
  # log of the integral over the free values with precision Lam (data and prior), linear term h
  logIntegral <- function(Lam, h, logTau) {
    if (!length(h)) return(0)
    R <- chol(Lam)
    -sum(logTau) - sum(log(diag(R))) + sum(backsolve(R, h, transpose = TRUE)^2) / 2
  }
  logZ <- numeric(particles)
  for (pt in seq_len(particles)) {
    sizes <- rep(list(1L), length(dpIdx))
    stemOf <- integer(0); clusterOf <- integer(0)
    Lam <- matrix(0, 0, 0); h <- numeric(0); logTau <- numeric(0); c0 <- 0; current <- 0
    for (g in seq_len(G)) {
      Pg <- P[[g]]; mg <- m[g, ]
      open <- which(slab & rowSums(abs(Pg)) > 0)
      seats <- if (length(dpIdx)) as.matrix(expand.grid(lapply(sizes, function(n) seq_len(length(n) + 1L))))
               else matrix(0L, 1, 0)
      subsets <- if (length(open)) as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), length(open))))
                 else matrix(FALSE, 1, 0)
      if (nrow(seats) * nrow(subsets) > budget)
        stop("em: importance sampling sums over at most ", budget, " seatings and spike-and-slab subsets of ",
             "a group.", call. = FALSE)
      combos <- expand.grid(seat = seq_len(nrow(seats)), subset = seq_len(nrow(subsets)))
      states <- vector("list", nrow(combos))
      logw <- numeric(nrow(combos))
      for (r in seq_len(nrow(combos))) {
        ch <- seats[combos$seat[r], ]
        on <- open[as.logical(subsets[combos$subset[r], ])]
        new <- which(ch > lengths(sizes))
        nNew <- length(new) + length(on)
        d <- length(h) + nNew
        A <- matrix(0, K, d)
        sOf <- c(stemOf, dpIdx[new]); cOf <- c(clusterOf, ch[new])
        for (s in seq_along(dpIdx)) if (ch[s] > 1L) A[dpIdx[s], which(sOf == dpIdx[s] & cOf == ch[s])] <- 1
        for (k in seq_along(on)) A[on[k], length(h) + length(new) + k] <- 1
        L2 <- matrix(0, d, d); L2[seq_along(h), seq_along(h)] <- Lam
        lt <- c(logTau, log(tau[dpIdx[new]]), log(tau[on]))
        diag(L2)[length(h) + seq_len(nNew)] <- 1 / exp(2 * lt[length(h) + seq_len(nNew)])
        L2 <- L2 + crossprod(A, Pg %*% A)
        h2 <- c(h, numeric(nNew)) + drop(crossprod(A, Pg %*% mg))
        cc <- c0 - sum(mg * (Pg %*% mg)) / 2
        prior <- sum(vapply(seq_along(dpIdx), function(s) {
          a <- concentration[dpIdx[s]]
          log(if (ch[s] > length(sizes[[s]])) a else sizes[[s]][ch[s]]) - log(g + a)
        }, 0)) + sum(ifelse(open %in% on, log(prob[open]), log1p(-prob[open])))
        value <- cc + logIntegral(L2, h2, lt)
        logw[r] <- prior + value - current
        # the slab values of the group are integrated out: the Schur complement keeps the clusters
        keep <- seq_len(length(h) + length(new))
        if (length(on)) {
          drop <- setdiff(seq_len(d), keep)
          Sinv <- solve(L2[drop, drop, drop = FALSE])
          h2 <- h2[keep] - drop(L2[keep, drop, drop = FALSE] %*% Sinv %*% h2[drop])
          L2 <- L2[keep, keep, drop = FALSE] - L2[keep, drop, drop = FALSE] %*% Sinv %*% L2[drop, keep, drop = FALSE]
          lt <- lt[keep]
          cc <- value - logIntegral(L2, h2, lt)
        }
        states[[r]] <- list(Lam = L2, h = h2, logTau = lt, c0 = cc, value = value, stemOf = sOf, clusterOf = cOf)
      }
      mx <- max(logw)
      logZ[pt] <- logZ[pt] + mx + log(sum(exp(logw - mx)))
      pick <- sample.int(length(logw), 1L, prob = exp(logw - mx))
      st <- states[[pick]]; ch <- seats[combos$seat[pick], ]
      for (s in seq_along(dpIdx)) {
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
