\dontrun{
  ## ===============================================================
  ## Example: parvec class and derivative handling
  ## ===============================================================
  
  ## --- 1. Basic creation ------------------------------------------------------
  
  # create a simple parameter vector
  p <- as.parvec(c(a = 1, b = 2, c = 3))
  print(p)
  
  # extract the Jacobian (identity by default)
  getDerivs(p)
  
  # extract the Hessian (zero tensor by default)
  getDerivs2(p)
  
  
  ## --- 2. Custom derivatives --------------------------------------------------
  
  # create a Jacobian (deriv) and Hessian (deriv2) manually
  J <- matrix(c(1, 0, 0,
                0, 2, 0,
                0, 0, 3),
              nrow = 3, byrow = TRUE,
              dimnames = list(c("a", "b", "c"), c("x", "y", "z")))
  
  H <- array(0, dim = c(3, 3, 3),
             dimnames = list(c("a", "b", "c"),
                             c("x", "y", "z"),
                             c("x", "y", "z")))
  H["a","x","x"] <- 1
  H["b","y","y"] <- 2
  H["c","z","z"] <- 3
  
  # combine into a parvec
  p2 <- as.parvec(c(a = 10, b = 20, c = 30), deriv = J, deriv2 = H)
  
  # print summary
  print(p2)
  
  # extract first and second derivatives
  getDerivs(p2)
  getDerivs2(p2)
  
  
  ## --- 3. Subsetting ----------------------------------------------------------
  
  # select subset of parameters
  p_sub <- p2[c("a", "b")]
  print(p_sub)
  
  # confirm that Jacobian and Hessian were subset correctly
  getDerivs(p_sub)
  getDerivs2(p_sub)
  
  
  ## --- 4. Concatenation -------------------------------------------------------
  
  # combine multiple parvecs
  p3 <- as.parvec(c(d = 4, e = 5))
  p_combined <- c(p2, p3)
  
  print(p_combined)
  dim(getDerivs(p_combined))
  dim(getDerivs2(p_combined))
  
  
  ## --- 5. Consistency check ---------------------------------------------------
  
  stopifnot(
    all.equal(names(p_combined), c("a", "b", "c", "d", "e")),
    all(dim(getDerivs(p_combined))[1:2] == length(p_combined))
  )
  
  cat("\nAll checks passed.\n")
  
  ## --- END -------------------------------------------------------------------
}