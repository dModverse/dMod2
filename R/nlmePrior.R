# Priors of the random effects and their closed-form M-steps. The normal prior is
# sum_g eta_g' Omega^-1 eta_g + G log|Omega|, Omega = L L' block diagonal by the correlation
# blocks with log on the diagonal of L, the -2 log density without the 2 pi constant.


## Omega^-1, the first and second derivatives of Omega^-1 by the Cholesky parameters, and the
## Hessian of tr(Omega^-1 RR) in them. log|Omega| is linear in the Cholesky parameters, so
## `hcc` is the Hessian of tr(Omega^-1 RR) + G log|Omega| as well.
.choleskyDerivatives <- function(layout, L, RR = NULL) {
  K <- layout$K
  loc <- unname(layout$loc)
  dg <- unname(layout$isDiag)
  nc <- length(layout$cholPars)
  Oinv <- crossprod(forwardsolve(L, diag(K)))
  dL <- lapply(seq_len(nc), function(m) {
    d <- matrix(0, K, K)
    d[loc[m, 1L], loc[m, 2L]] <- if (dg[m]) L[loc[m, , drop = FALSE]] else 1
    d
  })
  dO <- lapply(dL, function(d) d %*% t(L) + L %*% t(d))
  dOinv <- lapply(dO, function(d) -Oinv %*% d %*% Oinv)
  d2Oinv <- array(0, c(nc, nc, K, K))
  hcc <- matrix(0, nc, nc, dimnames = list(layout$cholPars, layout$cholPars))
  for (a in seq_len(nc)) for (b in seq_len(a)) {
    d2O <- dL[[a]] %*% t(dL[[b]]) + dL[[b]] %*% t(dL[[a]])
    if (a == b && dg[a]) d2O <- d2O + dO[[a]]
    d2 <- Oinv %*% (dO[[a]] %*% Oinv %*% dO[[b]] + dO[[b]] %*% Oinv %*% dO[[a]] - d2O) %*% Oinv
    d2Oinv[a, b, , ] <- d2Oinv[b, a, , ] <- d2
    if (!is.null(RR)) hcc[a, b] <- hcc[b, a] <- sum(d2 * RR)
  }
  list(Oinv = Oinv, dOinv = dOinv, d2Oinv = d2Oinv, hcc = hcc)
}


## The normal prior of the design as an objective over the etas of every group and the
## Cholesky parameters, with gradient and exact Hessian. It is 0 while a parameter is missing.
.normalPrior <- function(design) {
  lay <- design$normal
  K <- lay$K
  E <- design$etaNames[, lay$stems, drop = FALSE]
  G <- nrow(E)
  etaNames <- as.vector(E)
  chol <- lay$cholPars
  parnames <- c(etaNames, chol)
  loc <- unname(lay$loc)
  dg <- unname(lay$isDiag)

  prior <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, hessian = NULL,
                    conditions = NULL, env = NULL, cores = getOption("dMod2.cores", 1L)) {
    p <- list(...)[[match.fnargs(list(...), "pars")]]
    cv <- .resolveCurvature(deriv, deriv2, hessian, "forward")
    allp <- c(p, fixed)
    result <- function(value, gradient = NULL, hess = NULL) {
      gr <- setNames(numeric(length(p)), names(p))
      free <- .intersect(names(p), parnames)
      if (!is.null(gradient)) gr[free] <- gradient[free]
      hs <- if (cv$sparse) .sparseOver(hess, names(p))
            else if (cv$hessian) matrix(0, length(p), length(p), dimnames = list(names(p), names(p)))
      if (!is.null(hess) && !cv$sparse) hs[free, free] <- hess[free, free]
      out <- objlist(value = value, gradient = gr, hessian = hs)
      attr(out, "prior") <- value
      attr(out, "env") <- env
      out
    }
    if (!all(parnames %in% names(allp))) return(result(0))

    # whitened residuals Z = L^-1 eta_g, one column per group
    L <- .choleskyFactor(lay, allp)
    R <- t(matrix(allp[etaNames], G, K))
    Z <- forwardsolve(L, R)
    W <- backsolve(t(L), Z)
    value <- sum(Z * Z) + 2 * G * sum(log(diag(L)))
    if (!deriv) return(result(value))

    gc <- -2 * (W %*% t(Z))[loc]
    gc[dg] <- gc[dg] * L[loc][dg] + 2 * G
    gradient <- setNames(c(2 * as.vector(t(W)), gc), parnames)
    if (!cv$hessian) return(result(value, gradient))

    cd <- .choleskyDerivatives(lay, L, tcrossprod(R))
    hec <- matrix(vapply(cd$dOinv, function(d) 2 * as.vector(t(d %*% R)), numeric(G * K)),
                  G * K, length(chol))
    hess <- if (cv$sparse)
      rbind(cbind(2 * Matrix::kronecker(cd$Oinv, Matrix::Diagonal(G)), Matrix::Matrix(hec, sparse = TRUE)),
            cbind(Matrix::Matrix(t(hec), sparse = TRUE), Matrix::Matrix(cd$hcc, sparse = TRUE)))
    else rbind(cbind(2 * kronecker(cd$Oinv, diag(G)), hec), cbind(t(hec), cd$hcc))
    dimnames(hess) <- list(parnames, parnames)
    result(value, gradient, hess)
  }
  class(prior) <- c("objfn", "fn")
  attr(prior, "conditions") <- NULL
  attr(prior, "parameters") <- parnames
  prior
}


## M-step of the normal blocks: Omega restricted to the block structure that maximizes the
## expected prior given the second moments E[eta_g eta_g'] per group, as Cholesky parameters; in
## closed form per block, by trust over the free entries when some are fixed.
.updateNormal <- function(design, moments, fixed = NULL) {
  lay <- design$normal
  S <- Reduce(`+`, moments) / length(moments)
  dimnames(S) <- list(lay$stems, lay$stems)
  L <- matrix(0, lay$K, lay$K)
  start <- 0L
  for (b in design$blocks[vapply(design$blocks, function(b) b$dist == "normal", TRUE)]) {
    k <- start + seq_along(b$stems)
    L[k, k] <- t(chol(S[k, k, drop = FALSE]))
    start <- start + length(k)
  }
  v <- L[unname(lay$loc)]
  v[lay$isDiag] <- log(v[lay$isDiag])
  chol <- setNames(v, lay$cholPars)
  held <- .intersect(names(fixed), lay$cholPars)
  if (!length(held)) return(chol)
  chol[held] <- fixed[held]
  free <- .setdiff(lay$cholPars, held)
  # tr(Omega^-1 S) + log|Omega|, log|Omega| = 2 sum of the diagonal parameters
  expected <- function(p, ...) {
    all <- chol; all[free] <- p
    cd <- .choleskyDerivatives(lay, .choleskyFactor(lay, all), S)
    m <- match(free, lay$cholPars)
    objlist(value = sum(cd$Oinv * S) + 2 * sum(all[lay$isDiag]),
            gradient = setNames(vapply(m, function(j) sum(cd$dOinv[[j]] * S) + 2 * lay$isDiag[[j]], 0), free),
            hessian = cd$hcc[free, free, drop = FALSE])
  }
  if (!length(free)) return(chol)
  fit <- trust(expected, chol[free], rinit = 1, rmax = 10, iterlim = 100)
  chol[free] <- fit$argument
  chol
}
