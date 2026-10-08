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



## nullZ (moved from symbolics.R) --------------------------------------------

#' Null Space Basis from the Reduced Row Echelon Form
#'
#' Computes a basis of the null space of `A` from its reduced row echelon
#' form, like `null(A, "r")` in MATLAB: each basis vector has a 1 in one
#' non-pivot column and 0 in the others. The basis is integer if the reduced
#' row echelon form is.
#'
#' @param A Numeric matrix.
#' @param tol Tolerance for pivots. Not used: [rref()] is called with its
#'   default tolerance.
#' @return A matrix with `ncol(A)` rows and one column per basis vector.
#'
#' @author Malenka Mader, \email{Malenka.Mader@@fdm.uni-freiburg.de}
#' @seealso [rref()]
#' @examples
#' A <- rbind(c(1, -1, 0), c(0, 1, -1))
#' nullZ(A)
#' A %*% nullZ(A)
#'
#' @export
nullZ <- function(A, tol=sqrt(.Machine$double.eps)) {
  
  ret <- rref(A) # compute reduced row echelon form of A
  ret[[1]] -> R # matrix A in rref 
  ret[[2]] -> pivcol #columns in which a pivot was found
  
  n <- ncol(A) # number of columns of A
  r <- length(pivcol) # rank of reduced row echelon form
  nopiv <- 1:n
  nopiv <- nopiv[-pivcol]  # columns in which no pivot was found
  
  Z <- mat <- matrix(0, nrow = n, ncol = n-r) # matrix containing the vectors spanning the null space
  if ( n>r ) {
    Z[nopiv,] <- diag(1, n-r, n-r)
    if ( r>0 ) {
      Z[pivcol,] <- -R[1:r,nopiv]
    }
  }
  return (Z) 
  
}




## rref (moved from symbolics.R) ---------------------------------------------

#' Reduced Row Echelon Form of a Matrix
#'
#' Computes the reduced row echelon form of a numeric matrix by Gauss-Jordan
#' elimination with partial pivoting, like the MATLAB function `rref`.
#'
#' @param A Numeric matrix.
#' @param tol Numeric tolerance below which a pivot counts as zero. Defaults to
#'   `sqrt(.Machine$double.eps)`.
#' @param verbose Logical. Not used.
#' @param fractions Logical. Not used.
#'
#' @return An unnamed list with two elements: the reduced row echelon form of
#'   `A`, and the indices of the pivot columns as a one-row matrix (`NULL` if
#'   there is none).
#'
#' @author Malenka Mader, \email{Malenka.Mader@@fdm.uni-freiburg.de}. The
#'   signature and the argument check follow the function `rref()` by John Fox.
#' @references John Fox, `rref()`, R-help mailing list, posted by Scott Hyde,
#'   2007-09-01,
#'   <https://stat.ethz.ch/pipermail/r-help/2007-September/139923.html>.
#' @seealso [nullZ()]
#' @examples
#' A <- rbind(c(1, 2, 3), c(2, 4, 7), c(1, 2, 4))
#' rref(A)
#'
#' @export
rref <- function(A, tol=sqrt(.Machine$double.eps), verbose=FALSE, fractions=FALSE){
  ## Signature and argument check after John Fox
  if ((!is.matrix(A)) || (!is.numeric(A)))
    stop("argument must be a numeric matrix")
  m <- nrow(A)
  n <- ncol(A)
  
  i <- 1 # row index
  j <- 1 # column index
  pivcol <- c() # vector of columns in which nozero pivots are found
  while ((i <= m) & (j <= n)){
    # find pivot in column j
    which <- which.max(abs(A[i:m,j])) # column in which pivot is
    k <- i+which-1 #row index, in which pivot is
    pivot <- A[k, j] # pivot of column j
    
    if ( abs(pivot) <= tol ) {
      A[i:m,j] =matrix(0,m-i+1,1) # column is negligible, zero it out
      j <- j+1
    } else {
      # remember column index
      pivcol <- cbind(pivcol,j)
      
      # swap i-th and k-th column
      A[cbind(i,k),j:n] = A[cbind(k, i),j:n];
      
      # divide pivot row by pivot element.
      A[i,j:n] = A[i,j:n]/A[i,j];
      
      # subtract multiples of pivot row from all other rows.
      otherRows <- 1:m
      otherRows <- otherRows[-i]
      for (u in otherRows) {
        A[u,j:n] = A[u,j:n] - A[u,j]*A[i,j:n];
      }
      i = i + 1;
      j = j + 1;
    }
  }
  return (list(A,pivcol))
}




## submatrix (moved from tools.R) --------------------------------------------

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




## .matchNum (moved from data.R) ---------------------------------------------

# Match with numeric tolerance 
.matchNum <- function(x, y, tol = 1e-8) {
  sapply(x, function(xi) {
    d <- abs(y - xi)
    if (min(d) > tol) return(NA_integer_)
    which.min(d)
  })
}
