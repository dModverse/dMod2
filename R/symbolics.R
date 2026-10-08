#' Forcing Functions as Strings
#'
#' Returns the expression of a common input function of `time`:
#' \describe{
#'   \item{`"Gauss"`}{Gaussian peak with parameters `scale`, `mu`, `tau`.}
#'   \item{`"Fermi"`}{decreasing step with parameters `scale`, `mu`, `tau`.}
#'   \item{`"1-Fermi"`}{increasing step with parameters `scale`, `mu`, `tau`.}
#'   \item{`"MM"`}{saturating increase with parameters `slope`, `vmax`.}
#'   \item{`"Signal"`}{rise and decay with parameters `max1`, `max2`, `tau1`,
#'     `tau2`.}
#'   \item{`"Dose"`}{smooth box of total amount `Dose` from `Tlag + Tinit`
#'     to `Tlag + Tinit + Tduration`.}
#' }
#'
#' @param type The function, one of the names above. Defaults to `"Gauss"`.
#' @param parameters Named character or numeric vector. Replaces the
#'   parameters named by it with its values. Defaults to `NULL`.
#' @return A character string with the expression.
#' @examples
#' forcingsSymb("Gauss")
#' forcingsSymb("Fermi", parameters = c(scale = 1, tau = "0.5"))
#' @export
forcingsSymb <- function(type =c("Gauss", "Fermi", "1-Fermi", "MM", "Signal", "Dose"), parameters = NULL) {
  
  type <- match.arg(type)
  
  # INPUT1 (differentiable box)
  fn1 <- "(.5*(1.-tanh(.5*k*(time-T1))))" # T1 = start, T2 = end, k/4 = +-steepness in T1 and T2
  fn2 <- "(.5*(1.+tanh(.5*k*(time-T2))))" # T1 = start, T2 = end, k/4 = +-steepness in T1 and T2
  integral <- "(Tduration/(exp(20)-1))" # 100 = k*Tduration
  
  k <- "(20/Tduration)"
  T1 <- "(Tlag+Tinit)"
  T2 <- "(Tlag+Tinit+Tduration)"
  
  INPUT <- paste0("Dose*(", fn1, ")*(", fn2, ")/", integral)
  INPUT <- replaceSymbols(c("k", "T1", "T2"), c(k, T1, T2), INPUT)
  INPUT <- paste0("(", INPUT, ")")
  
  
  
  fun <- switch(type,
                "Gauss"   = "(scale*exp(-(time-mu)^2/(2*tau^2))/(tau*2.506628))",
                "Fermi"   = "(scale/(exp((time-mu)/tau)+1))",
                "1-Fermi" = "(scale*exp((time-mu)/tau)/(exp((time-mu)/tau)+1))",
                "MM"      = "(slope*time/(1 + slope*time/vmax))",
                "Signal"  = "(max1*max2*(1-exp(-time/tau1))*exp(-time*tau2))",
                "Dose"   = INPUT
  )
  
  if(!is.null(parameters)) {
    fun <- replaceSymbols(names(parameters), parameters, fun)
  }
  
  return(fun)
  
}


#' Get coefficients from a character
#'
#' @param char character, e.g. "2*x + y"
#' @param symbol single character, e.g. "x" or "y"
#' @return numeric vector with the coefficients
#' @noRd
getCoefficients <- function(char, symbol) {
  
  pdata <- getParseData(parse(text = char, keep.source = TRUE))
  pdata <- pdata[pdata$terminal == TRUE, ] #  subset(pdata, terminal == TRUE)
  symbolPos <- which(pdata$text == symbol)
  coefficients <- rep(1, length(symbolPos))
  
  hasCoefficient <- rep(FALSE, length(symbolPos))
  hasCoefficient[symbolPos > 1] <- (pdata$text[symbolPos[symbolPos > 1] - 1] == "*")
  coefficients[hasCoefficient] <- pdata$text[symbolPos[hasCoefficient]-2]
  
  return(as.numeric(coefficients))
  
  
  
  
  
}


#' Substitute Earlier Entries into Later Ones
#'
#' Replaces, for each entry, its name in all later entries by the entry in
#' parentheses, so that no entry depends on an earlier name, e.g. for steady
#' state expressions.
#'
#' @param variables Named character vector.
#' @return Named character vector of the same length as `variables`.
#' @examples resolveRecurrence(c(A = "k1*B/k2", C = "A*k3+k4", D="A*C*k5"))
#' @export
resolveRecurrence <- function (variables) {
  if(length(variables) > 1) {
    for (i in 1:(length(variables) - 1)) {
      newvariables <- c(variables[1:i], 
                        unlist(replaceSymbols(names(variables)[i],
                                              paste("(", variables[i], ")", sep = ""), 
                                              variables[(i + 1):length(variables)])))
      names(newvariables) <- names(variables)
      variables <- newvariables
    }
  }
  
  return(variables)
}



## getElements (moved from toolsMarcus.R) ------------------------------------

#' Get Symbols and Numeric Constants from a Character
#'
#' @param char Character vector of expressions. Entries `"0"` are skipped.
#' @param exclude Character vector of symbols to drop from the result.
#'   Defaults to `NULL`.
#' @return Character vector of the symbols and numeric constants in order of
#'   appearance, with repetitions; `NULL` if `char` is `NULL`.
#' @seealso [getSymbols()]
#' @export
#' 
#' @examples getElements(c("A*AB+B^2"))
#' 
getElements <- function (char, exclude = NULL) 
{
  if (is.null(char)) 
    return(NULL)
  char <- char[char != "0"]
  out <- parse(text = char, keep.source = TRUE)
  out <- utils::getParseData(out)
  names <- out$text[out$token == "SYMBOL" | out$token == "NUM_CONST"]
  if (!is.null(exclude)) 
    names <- names[!names %in% exclude]
  return(names)
}




## blockdiagSymb (moved from tools.R) ----------------------------------------

#' Embed Two Matrices into One Block Diagonal Matrix
#'
#' @param M Matrix, or `NULL`.
#' @param N Matrix, or `NULL`.
#' @return Matrix with `M` as upper left and `N` as lower right block and 0
#'   elsewhere, with the dimnames of both. If one argument is `NULL`, the
#'   other; if both are, `NULL`.
#' @examples
#' M <- matrix(1:9, 3, 3, dimnames = list(letters[1:3], letters[1:3]))
#' N <- matrix(1:4, 2, 2, dimnames = list(LETTERS[1:2], LETTERS[1:2]))
#' blockdiagSymb(M, N)
#' @export
blockdiagSymb <- function(M, N) {
  
  red <- sapply(list(M, N), is.null)
  if(all(red)) {
    return()
  } else if(red[1]) {
    return(N)
  } else if(red[2]) {
    return(M)
  }
  
  A <- matrix(0, ncol=dim(N)[2], nrow=dim(M)[1])
  B <- matrix(0, ncol=dim(M)[2], nrow=dim(N)[1])
  result <- rbind(cbind(M, A), cbind(B, N))
  colnames(result) <- c(colnames(M), colnames(N))
  rownames(result) <- c(rownames(M), rownames(N))
  
  return(result)
  
}
