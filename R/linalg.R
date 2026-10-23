#' Batched Matrix Multiplication
#'
#' Multiplies matrices batch by batch through BLAS. The batch index is the
#' first dimension of a 3D array.
#'
#' @section Supported Shapes:
#' \describe{
#'   \item{`[B,M,K] x [K,N] -> [B,M,N]`}{left operand batched}
#'   \item{`[M,K] x [B,K,N] -> [B,M,N]`}{right operand batched}
#'   \item{`[B,M,K] x [B,K,N] -> [B,M,N]`}{both operands batched}
#' }
#'
#' @param A Numeric matrix or 3D array.
#' @param B Numeric matrix or 3D array.
#' @return Numeric 3D array with dimensions `c(B, M, N)`.
#' @export
#' @examples
#' A <- array(rnorm(2 * 3 * 4), c(2, 3, 4))
#' M <- matrix(rnorm(4 * 5), 4, 5)
#' dim(A %bmm% M)
#' all.equal((A %bmm% M)[1, , ], A[1, , ] %*% M)
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

#' Null Space Basis from the Reduced Row Echelon Form
#'
#' Computes a basis of the null space of `A` from its reduced row echelon
#' form, like `null(A, "r")` in MATLAB: each basis vector has a 1 in one
#' non-pivot column and 0 in the others. The basis is integer if the reduced
#' row echelon form is.
#'
#' @param A Numeric matrix.
#' @param tol Numeric tolerance below which a pivot counts as zero, passed to
#'   [rref()].
#' @return A matrix with `ncol(A)` rows and one column per basis vector.
#'
#' @seealso [rref()]
#' @examples
#' A <- rbind(c(1, -1, 0), c(0, 1, -1))
#' nullZ(A)
#' A %*% nullZ(A)
#'
#' @export
nullZ <- function(A, tol = sqrt(.Machine$double.eps)) {
  r <- rref(A, tol)
  piv <- r$pivots
  free <- setdiff(seq_len(ncol(A)), piv)
  Z <- matrix(0, ncol(A), length(free))
  Z[cbind(free, seq_along(free))] <- 1
  if (length(piv)) Z[piv, ] <- -r$rref[seq_along(piv), free, drop = FALSE]
  Z
}


## rref ---------------------------------------------

#' Reduced Row Echelon Form of a Matrix
#'
#' Computes the reduced row echelon form of a numeric matrix by Gauss-Jordan
#' elimination with partial pivoting, like the MATLAB function `rref`.
#'
#' @param A Numeric matrix.
#' @param tol Numeric tolerance below which a pivot counts as zero. Defaults to
#'   `sqrt(.Machine$double.eps)`.
#'
#' @return A list with the elements `rref`, the reduced row echelon form of
#'   `A` with entries below `tol` set to zero, and `pivots`, the integer indices
#'   of the pivot columns.
#'
#' @seealso [nullZ()]
#' @examples
#' A <- rbind(c(1, 2, 3), c(2, 4, 7), c(1, 2, 4))
#' rref(A)
#'
#' @export
rref <- function(A, tol = sqrt(.Machine$double.eps)) {
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

#' Submatrix That Stays a Matrix
#'
#' @param M A matrix.
#' @param rows Index vector of rows. Defaults to all.
#' @param cols Index vector of columns. Defaults to all.
#' @return The matrix `M[rows, cols, drop = FALSE]`, with its dimnames.
#' @examples
#' M <- matrix(1:6, 2, 3, dimnames = list(c("a", "b"), c("x", "y", "z")))
#' submatrix(M, rows = 1)
#' @export
submatrix <- function(M, rows = 1:nrow(M), cols = 1:ncol(M)) {
  M[rows, cols, drop = FALSE]
}




## .matchNum ---------------------------------------------

# Match with numeric tolerance 
.matchNum <- function(x, y, tol = 1e-8) {
  sapply(x, function(xi) {
    d <- abs(y - xi)
    if (min(d) > tol) return(NA_integer_)
    which.min(d)
  })
}
