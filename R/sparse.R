## Sparse Hessians: assembly from local blocks, realignment, sums, and the
## trust-region step on a sparse Cholesky factor.

# Sparse matrix over `pnames` from the nonzero entries of the local Hessian blocks of the
# conditions; entries of the same pair are summed. The names of all blocks are matched at once.
.sparseFromBlocks <- function(blocks, pnames) {
  blocks <- Filter(function(h) !is.null(h) && length(h), blocks)
  n <- length(pnames)
  nm <- lapply(blocks, rownames)
  g <- split(match(unlist(nm, use.names = FALSE), pnames),
             factor(rep(seq_along(blocks), lengths(nm)), levels = seq_along(blocks)))
  parts <- Map(function(h, gi) {
    k <- which(!is.na(gi))
    if (!length(k)) return(NULL)
    x <- as.vector(h[k, k, drop = FALSE])
    nz <- which(x != 0)
    list(i = rep(gi[k], times = length(k))[nz], j = rep(gi[k], each = length(k))[nz], x = x[nz])
  }, blocks, g)
  parts <- Filter(Negate(is.null), parts)
  i <- unlist(lapply(parts, `[[`, "i")); j <- unlist(lapply(parts, `[[`, "j"))
  x <- unlist(lapply(parts, `[[`, "x"))
  Matrix::drop0(Matrix::sparseMatrix(i = as.integer(i), j = as.integer(j), x = as.numeric(x),
                                     dims = c(n, n), dimnames = list(pnames, pnames)))
}

# A Hessian, sparse or dense, as a sparse matrix over `pnames`; rows and
# columns outside `pnames` are dropped, missing ones are zero.
.sparseOver <- function(H, pnames) {
  n <- length(pnames)
  if (is.null(H) || !length(H))
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                dims = c(n, n), dimnames = list(pnames, pnames)))
  # general storage: a symmetric one holds one triangle only
  T <- methods::as(methods::as(methods::as(H, "CsparseMatrix"), "generalMatrix"), "TsparseMatrix")
  i <- match(rownames(H), pnames)[T@i + 1L]
  j <- match(colnames(H), pnames)[T@j + 1L]
  k <- !is.na(i) & !is.na(j)
  Matrix::sparseMatrix(i = i[k], j = j[k], x = T@x[k], dims = c(n, n),
                       dimnames = list(pnames, pnames))
}

.isSparse <- function(H) inherits(H, "sparseMatrix")

# A parameter-local term restricted to the inner parameters `inner` and the outer
# columns of the chain `dP` (and `dP2`) they depend on; without `dP` the chain is the
# identity on them. `outer` names the columns of the unrestricted result.
.sparseRestrict <- function(p, dP, dP2, inner) {
  inner <- intersect(names(p), inner)
  if (is.null(dP)) {
    dPs <- diag(1, length(inner))
    dimnames(dPs) <- list(inner, inner)
    return(list(pars = p[inner], dP = dPs, dP2 = NULL, outer = names(p)))
  }
  rows <- intersect(inner, rownames(dP))
  hit <- colSums(abs(dP[rows, , drop = FALSE])) > 0
  rows2 <- if (!is.null(dP2)) intersect(inner, dimnames(dP2)[[1]]) else character(0)
  if (length(rows2))
    hit <- hit | apply(abs(dP2[rows2, , , drop = FALSE]), 2, sum) > 0
  cols <- colnames(dP)[hit]
  list(pars = p[inner], dP = dP[rows, cols, drop = FALSE],
       dP2 = if (length(rows2)) dP2[rows2, cols, cols, drop = FALSE],
       outer = colnames(dP))
}

# The objlist of a restricted term over all of `outer`, its Hessian sparse.
.sparseLift <- function(value, gradient, hessian, outer) {
  g <- setNames(numeric(length(outer)), outer)
  if (length(gradient)) g[names(gradient)] <- gradient
  objlist(value = value, gradient = g, hessian = .sparseOver(hessian, outer))
}

# Whether a sparse Cholesky factor is numerically singular: a squared pivot below `tol` times
# the largest, which CHOLMOD leaves for a semidefinite matrix where solve() would stop.
.cholSingular <- function(ch, tol = 1e-14) {
  d <- Matrix::diag(methods::as(ch, "Matrix"))^2
  !length(d) || !all(is.finite(d)) || min(d) <= tol * max(d)
}

# The rule of .cholSingular() for a sparse LU factor of an indefinite matrix, by its pivots;
# a factorisation that failed is singular.
.luSingular <- function(lu, tol = 1e-14) {
  if (!methods::is(lu, "sparseLU")) return(TRUE)
  d <- abs(Matrix::diag(lu@U))
  !length(d) || !all(is.finite(d)) || min(d) <= tol * max(d)
}

# The columns of a singular sparse M at the pivots .luSingular() finds, one per null
# direction. An exactly singular M is factorised with a diagonal shift below `tol`.
.luSingularColumns <- function(M, tol = 1e-14) {
  lu <- Matrix::lu(M, errSing = FALSE)
  if (!methods::is(lu, "sparseLU"))
    lu <- Matrix::lu(M + Matrix::Diagonal(nrow(M), tol / 10 * max(abs(M@x))), errSing = FALSE)
  if (!methods::is(lu, "sparseLU")) return(integer(0))
  d <- abs(Matrix::diag(lu@U))
  q <- if (length(lu@q)) lu@q + 1L else seq_along(d)
  sort(q[d <= tol * max(d)])
}

# Trust-region step: minimises g'p + p'Bp/2 subject to |p| <= radius with
# sparse Cholesky factors of B + lambda I (More and Sorensen 1983), B = H plus
# the low-rank part U diag(d) U' of `lowrank`, solved by Woodbury.
.sparseTrustStep <- function(H, g, radius, lambda0 = 0, lowrank = NULL) {
  n <- length(g)
  I <- Matrix::Diagonal(n)
  U <- lowrank$U
  d <- lowrank$d
  # NULL where B + lambda I is not positive definite: with H + lambda I positive definite, B + lambda I
  # has as many negative eigenvalues as diag(1/d) + U'(H + lambda I)^-1 U has positive ones beyond d;
  # a solve that is not finite counts as not positive definite
  fac <- function(lambda) {
    ch <- tryCatch(suppressWarnings(Matrix::Cholesky(Matrix::forceSymmetric(H + lambda * I),
                                    perm = TRUE, LDL = FALSE)), error = function(e) NULL)
    if (is.null(ch) || !length(d)) return(list(ch = ch)[!is.null(ch)])
    W <- as.matrix(Matrix::solve(ch, U))
    K <- diag(1 / d, length(d)) + crossprod(U, W)
    if (!all(is.finite(K))) return(NULL)
    ev <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
    if (any(ev == 0) || sum(ev > 0) != sum(d > 0)) return(NULL)
    list(ch = ch, W = W, K = K)
  }
  # the limit of the step for lambda to infinity
  steepest <- function() list(p = -radius * g / sqrt(sum(g^2)), lambda = Inf, boundary = TRUE)
  solveB <- function(f, b) {
    x <- as.numeric(Matrix::solve(f$ch, b))
    if (is.null(f$W)) x else x - drop(f$W %*% solve(f$K, crossprod(f$W, b)))
  }
  scale <- max(1, max(abs(Matrix::diag(H))))
  lambda <- lambda0
  f <- fac(lambda)
  while (!length(f)) {
    lambda <- if (lambda <= 0) 1e-10 * scale else 4 * lambda
    if (!is.finite(lambda)) return(steepest())
    f <- fac(lambda)
  }
  p <- -solveB(f, g)
  np <- sqrt(sum(p^2))
  if (lambda == 0 && np <= radius) return(list(p = p, lambda = 0, boundary = FALSE))
  for (it in 1:50) {
    if (abs(np - radius) <= 1e-3 * radius) break
    q2 <- sum(p * solveB(f, p))
    lnew <- lambda + (np / radius - 1) * np^2 / q2
    if (lnew <= 0) lnew <- lambda / 10
    f2 <- fac(lnew)
    while (!length(f2) && is.finite(lnew)) { lnew <- 2 * lnew + 1e-10 * scale; f2 <- fac(lnew) }
    if (!length(f2)) return(steepest())
    lambda <- lnew; f <- f2
    p <- -solveB(f, g); np <- sqrt(sum(p^2))
    if (np < radius && lambda < 1e-12 * scale) break
  }
  list(p = p, lambda = lambda, boundary = TRUE)
}

# The diagonal of the inverse of half a sparse Hessian over `which`, by one Cholesky factor
# solved in blocks of columns; a singular one through vcov().
.sparseVarianceDiagonal <- function(H, which, block = 256L) {
  which <- intersect(which, colnames(H))
  Hw <- Matrix::forceSymmetric(0.5 * H[which, which, drop = FALSE])
  ch <- tryCatch(suppressWarnings(Matrix::Cholesky(Hw, perm = TRUE, LDL = FALSE)), error = function(e) NULL)
  if (is.null(ch) || .cholSingular(ch)) return(diag(.vcovSparse(H, NULL, which)))
  n <- length(which)
  out <- numeric(n)
  for (start in seq(1L, n, by = block)) {
    j <- start:min(n, start + block - 1L)
    E <- Matrix::sparseMatrix(i = j, j = seq_along(j), x = 1, dims = c(n, length(j)))
    out[j] <- as.matrix(Matrix::solve(ch, E))[cbind(j, seq_along(j))]
  }
  setNames(out, which)
}

# `sparse` of a fitting function: TRUE, FALSE or "auto".
.checkSparse <- function(sparse, fn) {
  if (!(isTRUE(sparse) || isFALSE(sparse) || identical(sparse, "auto")))
    stop(fn, ": `sparse` is TRUE, FALSE or \"auto\".", call. = FALSE)
  sparse
}

# The thresholds of trust(sparse = "auto"), from getOption("dMod2.sparseAuto").
.sparseThresholds <- function()
  utils::modifyList(list(parameters = 500, fill = 0.05), as.list(getOption("dMod2.sparseAuto")))

# Whether the objective takes the sparse path: at least `parameters` free parameters and a
# sparse Hessian at `parinit` that fills at most `fill`; `first` is that evaluation.
.sparseChoice <- function(objfun, parinit, dots) {
  th <- .sparseThresholds()
  P <- length(parinit)
  if (P < th$parameters) return(list(sparse = FALSE))
  out <- tryCatch(do.call(objfun, c(list(parinit, deriv = TRUE, hessian = "sparse"), dots)),
                  error = function(e) NULL)
  if (is.null(out) || !is.finite(out$value) || !.isSparse(out$hessian)) return(list(sparse = FALSE))
  fill <- Matrix::nnzero(.sparseOver(out$hessian, names(parinit))) / P^2
  if (fill > th$fill) return(list(sparse = FALSE))
  list(sparse = TRUE, first = out)
}

# The terms of a sum that return a dense Hessian at `pars`, for the error of trust(sparse = TRUE).
.denseTerms <- function(objfun, pars, dots) {
  terms <- if (inherits(objfun, "objfn")) .objTerms(objfun) else list(objfun)
  dense <- vapply(terms, function(f) {
    out <- tryCatch(do.call(f, c(list(pars, deriv = TRUE, hessian = "sparse"), dots)), error = function(e) NULL)
    !is.null(out) && !is.null(out$hessian) && !.isSparse(out$hessian)
  }, TRUE)
  which(dense)
}

# Trust-region minimisation with sparse Hessians. objfun(pars, deriv, hessian = "sparse", ...)
# returns value, gradient and a sparse Hessian; `first` is its evaluation at `parinit`.
.trustSparse <- function(objfun, parinit, rinit, rmax, iterlim, fterm, mterm, dots, printIter = FALSE,
                         first = NULL) {
  pnames <- names(parinit)
  failure <- "its value is not finite"
  prepare <- function(out) {
    if (is.null(out) || !is.finite(out$value)) return(NULL)
    if (!.isSparse(out$hessian)) {
      dense <- .denseTerms(objfun, setNames(x, pnames), dots)
      stop("trust(sparse = TRUE): the objective returns a dense Hessian",
           if (length(dense) && inherits(objfun, "objfn"))
             paste0(", from term ", paste(dense, collapse = ", "), " of its sum"),
           "; use sparse = FALSE.", call. = FALSE)
    }
    out$gradient <- out$gradient[pnames]
    out$hessian <- .sparseOver(out$hessian, pnames)
    out
  }
  call <- function(x) {
    out <- tryCatch(do.call(objfun, c(list(setNames(x, pnames), deriv = TRUE, hessian = "sparse"), dots)),
                    error = function(e) { failure <<- conditionMessage(e); NULL })
    prepare(out)
  }
  x <- as.numeric(parinit)
  f <- if (!is.null(first)) prepare(first) else call(x)
  if (is.null(f)) stop("trust(sparse = TRUE): the objective fails at `parinit`: ", failure, call. = FALSE)
  r <- rinit
  converged <- FALSE
  it <- 0L
  while (it < iterlim) {
    it <- it + 1L
    g <- as.numeric(f$gradient)
    st <- .sparseTrustStep(f$hessian, g, r)
    pred <- -(sum(g * st$p) + 0.5 * sum(st$p * as.numeric(f$hessian %*% st$p)))
    if (pred <= mterm) { converged <- TRUE; break }
    fn <- call(x + st$p)
    rho <- if (is.null(fn)) -Inf else (f$value - fn$value) / pred
    if (printIter) cat(sprintf("%4d  value %.8g  radius %.3g  rho %.3g\n", it, f$value, r, rho))
    if (rho < 0.25) r <- r / 4 else if (rho > 0.75 && st$boundary) r <- min(2 * r, rmax)
    if (rho > 0) {
      done <- abs(f$value - fn$value) <= fterm
      x <- x + st$p; f <- fn
      if (done) { converged <- TRUE; break }
    }
  }
  list(value = f$value, gradient = setNames(f$gradient, pnames), hessian = f$hessian,
       argument = setNames(x, pnames), converged = converged, iterations = it)
}
