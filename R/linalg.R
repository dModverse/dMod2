#' Batched Matrix Multiplication
#'
#' Efficient batched matrix multiplication with BLAS-backend.
#' Batch index is the FIRST dimension.
#'
#' @section Supported Contractions:
#' \describe{
#'   \item{\code{[B,M,K] x [K,N] -> [B,M,N]}}{Left batched}
#'   \item{\code{[M,K] x [B,K,N] -> [B,M,N]}}{Right batched}
#'   \item{\code{[B,M,K] x [B,K,N] -> [B,M,N]}}{Both batched}
#' }
#'
#' @param A Numeric matrix or 3D array
#' @param B Numeric matrix or 3D array
#' @return Numeric 3D array with dim \code{[B, M, N]}
#' @export
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

#' Find integer-null space of matrix A 
#  this function is written along the lines of matlab's function null(A,'r'), where 'r' specifies that integer solutions are returned instead of orthonormal vectors
#' 
#' @param A matrix for which the null space is searched
#' @param tol tolerance to find pivots in rref-function below
#' @return null space of A with only integers in it
#' 
#' @author Malenka Mader, \email{Malenka.Mader@@fdm.uni-freiburg.de}
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

#' Transform a matrix into reduced row echelon form
#'
#' This function computes the reduced row echelon form (RREF) of a numeric matrix.
#' It is written along the lines of the MATLAB function \code{rref}.
#'
#' @param A Numeric matrix for which the reduced row echelon form is computed.
#' @param tol Numeric tolerance used to identify pivot elements. Defaults to \code{sqrt(.Machine$double.eps)}.
#' @param verbose Logical. If \code{TRUE}, prints detailed information during computation.
#' @param fractions Logical. Currently not used.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{[[1]]}: The reduced row echelon form of \code{A}.
#'   \item \code{[[2]]}: The indices of the columns in which pivots were found.
#' }
#'
#' @author Malenka Mader, \email{Malenka.Mader@@fdm.uni-freiburg.de}
#'
#' @export
rref <- function(A, tol=sqrt(.Machine$double.eps), verbose=FALSE, fractions=FALSE){
  ## Written by John Fox
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

#' Submatrix of a matrix returning ALWAYS a matrix
#' 
#' @param M matrix
#' @param rows Index vector
#' @param cols Index vector
#' @return The matrix `M[rows, cols]`, keeping/adjusting attributes like ncol nrow and dimnames.
#' @export
submatrix <- function(M, rows = 1:nrow(M), cols = 1:ncol(M)) {
  
 M[rows, cols, drop = FALSE] 
  
  # myrows <- (structure(1:nrow(M), names = rownames(M)))[rows]
  # mycols <- (structure(1:ncol(M), names = colnames(M)))[cols]
  # 
  # if(any(is.na(myrows)) | any(is.na(mycols))) stop("subscript out of bounds")
  # 
  # matrix(M[myrows, mycols], 
  #        nrow = length(myrows), ncol = length(mycols), 
  #        dimnames = list(rownames(M)[myrows], colnames(M)[mycols]))

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
