# Batched matrix product through BLAS, the batch index first: [B,M,K] x [K,N], [M,K] x [B,K,N] and
# [B,M,K] x [B,K,N], each to [B,M,N]
`%bmm%` <- function(A, B) {
  da <- dim(A)
  db <- dim(B)

  if (length(da) == 3 && length(db) == 2) {
    # [B,M,K] x [K,N] -> [B,M,N]
    stopifnot(da[3] == db[1])
    bmm_lb(A, B, da[1], da[2], da[3], db[2])

  } else if (length(da) == 2 && length(db) == 3) {
    # [M,K] x [B,K,N] -> [B,M,N]
    stopifnot(da[2] == db[2])
    bmm_rb(A, B, db[1], da[1], da[2], db[3])

  } else if (length(da) == 3 && length(db) == 3) {
    # [B,M,K] x [B,K,N] -> [B,M,N]
    stopifnot(da[1] == db[1], da[3] == db[2])
    bmm_bb(A, B, da[1], da[2], da[3], db[3])

  } else {
    stop("Invalid dimensions for %bmm%")
  }
}



## nullZ --------------------------------------------

# Null space basis from the reduced row echelon form, as null(A, "r") in MATLAB: each basis vector
# has a 1 in one non-pivot column and 0 in the others, integer if the echelon form is
.nullZ <- function(A, tol = sqrt(.Machine$double.eps)) {
  r <- .rref(A, tol)
  piv <- r$pivots
  free <- setdiff(seq_len(ncol(A)), piv)
  Z <- matrix(0, ncol(A), length(free))
  Z[cbind(free, seq_along(free))] <- 1
  if (length(piv)) Z[piv, ] <- -r$rref[seq_along(piv), free, drop = FALSE]
  Z
}


## rref ---------------------------------------------

# Reduced row echelon form by Gauss-Jordan elimination with partial pivoting: the form with entries
# below tol set to zero, and the pivot columns
.rref <- function(A, tol = sqrt(.Machine$double.eps)) {
  if (!is.numeric(A) || !is.matrix(A))
    stop("'A' must be a numeric matrix", call. = FALSE)
  m <- nrow(A)
  pivots <- integer(0)
  i <- 1L
  for (j in seq_len(ncol(A))) {
    if (i > m) break
    k <- i - 1L + which.max(abs(A[i:m, j]))
    if (abs(A[k, j]) <= tol) {
      A[i:m, j] <- 0
      next
    }
    pivots <- c(pivots, j)
    A[c(i, k), ] <- A[c(k, i), ]
    A[i, ] <- A[i, ] / A[i, j]
    others <- seq_len(m)[-i]
    A[others, ] <- A[others, ] - outer(A[others, j], A[i, ])
    i <- i + 1L
  }
  A[abs(A) <= tol] <- 0
  list(rref = A, pivots = pivots)
}




## submatrix --------------------------------------------

# M[rows, cols] that stays a matrix
.submatrix <- function(M, rows = seq_len(nrow(M)), cols = seq_len(ncol(M))) {
  M[rows, cols, drop = FALSE]
}




## .matchNum ---------------------------------------------

# Match with numeric tolerance; on a strictly increasing `y` by bisection, ties to the lower index
.matchNum <- function(x, y, tol = 1e-8) {
  n <- length(y)
  if (length(x) && n && !anyNA(x) && !anyNA(y) && !is.unsorted(y, strictly = TRUE)) {
    lo <- pmax(findInterval(x, y), 1L)
    hi <- pmin(lo + 1L, n)
    k <- ifelse(abs(y[hi] - x) < abs(y[lo] - x), hi, lo)
    k[abs(y[k] - x) > tol] <- NA_integer_
    return(k)
  }
  sapply(x, function(xi) {
    d <- abs(y - xi)
    if (min(d) > tol) return(NA_integer_)
    which.min(d)
  })
}
