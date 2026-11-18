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

# The sparse matrix over `pnames` from the 0-based slots list(i, p, x) of a column-compressed one.
.sparseFromSlots <- function(s, pnames) {
  n <- length(pnames)
  methods::new("dgCMatrix", i = s$i, p = s$p, x = s$x, Dim = c(n, n), Dimnames = list(pnames, pnames))
}

# A Hessian, sparse or dense, as a sparse matrix over `pnames`; rows and
# columns outside `pnames` are dropped, missing ones are zero.
.sparseOver <- function(H, pnames) {
  n <- length(pnames)
  if (is.null(H) || !length(H))
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                dims = c(n, n), dimnames = list(pnames, pnames)))
  if (methods::is(H, "dgCMatrix") && identical(rownames(H), pnames) && identical(colnames(H), pnames))
    return(H)
  # general storage: a symmetric one holds one triangle only
  T <- methods::as(methods::as(methods::as(H, "CsparseMatrix"), "generalMatrix"), "TsparseMatrix")
  i <- match(rownames(H), pnames)[T@i + 1L]
  j <- match(colnames(H), pnames)[T@j + 1L]
  k <- !is.na(i) & !is.na(j)
  Matrix::sparseMatrix(i = i[k], j = j[k], x = T@x[k], dims = c(n, n),
                       dimnames = list(pnames, pnames))
}

.isSparse <- function(H) inherits(H, "sparseMatrix")

# The low-rank part `lr` of a Hessian, list(U, d), with the rows of U over `pnames`: missing ones
# are zero, others are dropped. NULL for none.
.lowrankOver <- function(lr, pnames) {
  if (is.null(lr) || !length(lr$d)) return(NULL)
  U <- matrix(0, length(pnames), length(lr$d), dimnames = list(pnames, NULL))
  hit <- intersect(rownames(lr$U), pnames)
  U[hit, ] <- lr$U[hit, , drop = FALSE]
  list(U = U, d = lr$d)
}

# The low-rank parts of two Hessians as one, over the union of their rows.
.lowrankBind <- function(l1, l2) {
  if (is.null(l1) || !length(l1$d)) return(.lowrankOver(l2, rownames(l2$U)))
  if (is.null(l2) || !length(l2$d)) return(.lowrankOver(l1, rownames(l1$U)))
  rn <- .union(rownames(l1$U), rownames(l2$U))
  a <- .lowrankOver(l1, rn)
  b <- .lowrankOver(l2, rn)
  list(U = cbind(a$U, b$U), d = c(a$d, b$d))
}

# Two low-rank parts over the same rows as one.
.lowrankJoin <- function(l1, l2) {
  if (is.null(l1) || !length(l1$d)) return(if (length(l2$d)) l2)
  if (is.null(l2) || !length(l2$d)) return(l1)
  list(U = cbind(l1$U, l2$U), d = c(l1$d, l2$d))
}

# H x plus the low-rank part U diag(d) U' x, U over the rows of H.
.lowrankApply <- function(H, lr, x) {
  y <- as.numeric(H %*% x)
  if (is.null(lr) || !length(lr$d)) return(y)
  y + drop(lr$U %*% (lr$d * crossprod(lr$U, x)))
}

# A sparse H plus the low-rank part U diag(d) U' as one sparse matrix, dense where U is.
.lowrankDense <- function(H, lr) {
  if (is.null(lr) || !length(lr$d)) return(H)
  H + Matrix::Matrix(lr$U %*% (lr$d * t(lr$U)), sparse = TRUE)
}

# The Hessian of an objlist over `pnames` as `hessian`, sparse, and `lowrank` over `pnames`, where the
# low-rank part has no rows outside `keep`; else that part is added to `hessian`. An operation
# on the rows outside `keep`, such as a Schur complement over them, then leaves `lowrank` as it is.
.hessianParts <- function(out, pnames, keep = pnames) {
  H <- .sparseOver(out$hessian, pnames)
  lr <- .lowrankOver(out$lowrank, pnames)
  if (is.null(lr)) return(list(hessian = H, lowrank = NULL))
  out <- setdiff(pnames, keep)
  if (length(out) && any(lr$U[out, ] != 0)) return(list(hessian = .lowrankDense(H, lr), lowrank = NULL))
  list(hessian = H, lowrank = .lowrankOver(lr, keep))
}

# A factor P'LDL'P of the sparse symmetric H (D = I by Cholesky, LDL' where H is indefinite) for
# B = H + U diag(d) U' by Woodbury, with Y = L^-1 P U, Kinv and log det B; NULL where B is not
# positive definite: it has neg(H) + pos(K) - pos(d) negative eigenvalues, K = diag(1/d) + Y'D^-1 Y.
.lowrankCholesky <- function(H, lr = NULL) {
  Hs <- Matrix::forceSymmetric(H)
  ch <- tryCatch(suppressWarnings(Matrix::Cholesky(Hs, perm = TRUE, LDL = FALSE)), error = function(e) NULL)
  Dinv <- NULL
  if (!is.null(ch) && !.cholSingular(ch)) logdet <- 2 * sum(log(Matrix::diag(methods::as(ch, "Matrix"))))
  else {
    if (is.null(lr) || !any(lr$d > 0)) return(NULL)
    ch <- tryCatch(suppressWarnings(Matrix::Cholesky(Hs, perm = TRUE, LDL = TRUE, super = FALSE)),
                   error = function(e) NULL)
    if (is.null(ch)) return(NULL)
    D <- Matrix::diag(Matrix::expand1(ch, "D"))
    if (!all(is.finite(D)) || !all(abs(D) > 1e-14 * max(abs(D)))) return(NULL)
    Dinv <- 1 / D
    logdet <- sum(log(abs(D)))
  }
  negH <- sum(Dinv < 0)
  f <- list(ch = ch, Dinv = Dinv, logdet = logdet)
  if (is.null(lr) || !length(lr$d)) return(if (!negH) f)
  Y <- as.matrix(.lowrankZ(f, lr$U))
  K <- diag(1 / lr$d, length(lr$d)) + crossprod(Y, .lowrankD(f, Y))
  ev <- eigen((K + t(K)) / 2, symmetric = TRUE, only.values = TRUE)$values
  if (!all(is.finite(ev)) || min(abs(ev)) <= 1e-14 * max(abs(ev)) || negH + sum(ev > 0) != sum(lr$d > 0))
    return(NULL)
  c(f, list(Y = Y, Kinv = solve(K), logdet = logdet + sum(log(abs(lr$d))) + sum(log(abs(ev)))))
}

# z = L^-1 P X of a factor from .lowrankCholesky(), sparse for a sparse X; D^-1 z; and the
# quadratic forms of B^-1 between the columns of X from z: z' D^-1 z - W' Kinv W, W = Y' D^-1 z,
# the diagonal or the matrix.
.lowrankZ <- function(f, X) Matrix::solve(f$ch, Matrix::solve(f$ch, X, system = "P"), system = "L")
.lowrankD <- function(f, Z) if (is.null(f$Dinv)) Z else f$Dinv * Z
.lowrankQuadratic <- function(f, Z, diagonal = TRUE) {
  DZ <- .lowrankD(f, Z)
  W <- if (!is.null(f$Y)) as.matrix(crossprod(f$Y, DZ))
  if (diagonal) {
    q <- Matrix::colSums(Z * DZ)
    if (!is.null(W)) q <- q - colSums(W * (f$Kinv %*% W))
    return(as.numeric(q))
  }
  Q <- as.matrix(Matrix::crossprod(Z, DZ))
  if (!is.null(W)) Q <- Q - crossprod(W, f$Kinv %*% W)
  (Q + t(Q)) / 2
}

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

# The objlist of a restricted term over all of `outer`, with a sparse Hessian.
.sparseLift <- function(value, gradient, hessian, outer) {
  g <- setNames(numeric(length(outer)), outer)
  if (length(gradient)) g[names(gradient)] <- gradient
  objlist(value = value, gradient = g, hessian = .sparseOver(hessian, outer))
}

# Whether a Cholesky factor, sparse or dense, is numerically singular: a squared pivot below `tol`
# times the largest, which CHOLMOD leaves for a semidefinite matrix where solve() would stop.
.cholSingular <- function(ch, tol = 1e-14) {
  d <- Matrix::diag(ch)
  !length(d) || !all(is.finite(d)) || min(d) <= tol * max(d)
}

# The rule of .cholSingular() for a sparse LU factor of an indefinite matrix, by the pivots;
# a factorization that failed is singular.
.luSingular <- function(lu, tol = 1e-14) {
  if (!methods::is(lu, "sparseLU")) return(TRUE)
  d <- abs(Matrix::diag(lu@U))
  !length(d) || !all(is.finite(d)) || min(d) <= tol * max(d)
}

# The columns of a singular sparse M at the pivots .luSingular() finds, one per null
# direction. An exactly singular M is factorized with a diagonal shift below `tol`.
.luSingularColumns <- function(M, tol = 1e-14) {
  lu <- Matrix::lu(M, errSing = FALSE)
  if (!methods::is(lu, "sparseLU"))
    lu <- Matrix::lu(M + Matrix::Diagonal(nrow(M), tol / 10 * max(abs(M@x))), errSing = FALSE)
  if (!methods::is(lu, "sparseLU")) return(integer(0))
  d <- abs(Matrix::diag(lu@U))
  q <- if (length(lu@q)) lu@q + 1L else seq_along(d)
  sort(q[d <= tol * max(d)])
}

# The Euclidean norm of `x`, free of under- and overflow.
.norm2 <- function(x) {
  m <- max(abs(x))
  if (m == 0 || !is.finite(m)) m else m * sqrt(sum((x / m)^2))
}

# The diagonal of the inverse of half a sparse Hessian plus the low-rank part over `which`, by one
# Cholesky factor solved in blocks of columns; a singular one through vcov().
.sparseVarianceDiagonal <- function(H, which, block = 256L, lowrank = NULL) {
  which <- intersect(which, colnames(H))
  lr <- .lowrankOver(lowrank, which)
  if (!is.null(lr)) lr$d <- lr$d / 2
  f <- .lowrankCholesky(0.5 * H[which, which, drop = FALSE], lr)
  if (is.null(f)) return(diag(.vcovSparse(H, NULL, which, lowrank)))
  n <- length(which)
  out <- numeric(n)
  for (start in seq(1L, n, by = block)) {
    j <- start:min(n, start + block - 1L)
    E <- Matrix::sparseMatrix(i = j, j = seq_along(j), x = 1, dims = c(n, length(j)))
    out[j] <- .lowrankQuadratic(f, .lowrankZ(f, E))
  }
  setNames(out, which)
}

# `sparse` of a fitting function: TRUE, FALSE or "auto".
.checkSparse <- function(sparse, fn) {
  if (!(isTRUE(sparse) || isFALSE(sparse) || identical(sparse, "auto")))
    stop(fn, ": `sparse` is TRUE, FALSE or \"auto\".", call. = FALSE)
  sparse
}

# The settings of sparse = "auto", from getOption("dMod2.sparseAuto"): below `parameters` free
# parameters the dense path, above it the cheaper one by .sparseCost(); a `fill` given there
# replaces the cost by the rule fill of the sparse Hessian <= `fill`.
.sparseThresholds <- function()
  utils::modifyList(list(parameters = 100, fill = NULL), as.list(getOption("dMod2.sparseAuto")))

# The flops of a trust-region step on a sparse Hessian H plus a low-rank part of m columns against
# the dense one: about eight factorizations of H + lambda I, each with m + 2 solves, and the
# Woodbury blocks, against the eigendecomposition of the dense model, 4 P^3.
.sparseCost <- function(H, m) {
  P <- nrow(H)
  T <- methods::as(H, "TsparseMatrix")
  up <- T@i <= T@j
  U <- Matrix::sparseMatrix(i = T@i[up] + 1L, j = T@j[up] + 1L, x = 1, dims = c(P, P))
  f <- sparse_factor_cost(U@p, U@i, P)
  c(sparse = 8 * (f[["flops"]] + 4 * f[["lnz"]] * (m + 2) + 2 * m^2 * P), dense = 4 * P^3)
}

# Whether the objective takes the sparse path: at least `parameters` free parameters and a
# sparse Hessian at `parinit` whose steps cost less than dense ones; `first` is that evaluation.
.sparseChoice <- function(objfun, parinit, dots) {
  th <- .sparseThresholds()
  P <- length(parinit)
  if (P < th$parameters) return(list(sparse = FALSE))
  out <- tryCatch(do.call(objfun, c(list(parinit, deriv = TRUE, hessian = "sparse"), dots)),
                  error = function(e) NULL)
  if (is.null(out) || !is.finite(out$value)) return(list(sparse = FALSE))
  if (!.isSparse(out$hessian)) return(list(sparse = FALSE, first = out))
  H <- .sparseOver(out$hessian, names(parinit))
  cheaper <- if (!is.null(th$fill)) Matrix::nnzero(H) / P^2 <= th$fill
             else { cost <- .sparseCost(H, length(out$lowrank$d)); cost[["sparse"]] < cost[["dense"]] }
  if (cheaper) return(list(sparse = TRUE, first = out))
  # the dense path takes the same evaluation, with the Hessian over the names of the gradient
  pn <- names(out$gradient) %||% names(parinit)
  out$hessian <- as.matrix(.lowrankDense(.sparseOver(out$hessian, pn), .lowrankOver(out$lowrank, pn)))
  out$lowrank <- NULL
  list(sparse = FALSE, first = out)
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
