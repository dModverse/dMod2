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
  pdata <- pdata[pdata$terminal == TRUE, ]
  symbolPos <- which(pdata$text == symbol)
  coefficients <- rep(1, length(symbolPos))
  
  hasCoefficient <- rep(FALSE, length(symbolPos))
  hasCoefficient[symbolPos > 1] <- (pdata$text[symbolPos[symbolPos > 1] - 1] == "*")
  coefficients[hasCoefficient] <- pdata$text[symbolPos[hasCoefficient]-2]
  
  return(as.numeric(coefficients))
  
  
  
  
  
}


#' Symbols in Expressions
#'
#' Lists the variables that occur in character expressions, each once, in
#' order of first appearance. Function names, argument names and `pi`, the
#' number, are not variables.
#'
#' @param char Character vector of expressions. Entries `"0"` are skipped.
#' @param exclude Character vector of symbols to drop. Default `NULL`.
#' @return Character vector of the symbols, `NULL` if `char` is `NULL`.
#' @seealso [replaceSymbols()], [getElements()]
#' @examples
#' getSymbols(c("k1*A/(Km + A)", "exp(-k2*time)"))
#' getSymbols("k1*A + k2*B", exclude = "A")
#' @export
getSymbols <- function(char, exclude = NULL) {
  if (is.null(char)) return(NULL)
  char <- char[char != "0"]
  pd <- utils::getParseData(parse(text = char, keep.source = TRUE))
  out <- unique(pd$text[pd$token == "SYMBOL"])
  out[!out %in% c("pi", exclude)]
}


#' Replace Symbols in Expressions
#'
#' Substitutes symbols in character expressions. A replacement that is not a
#' single symbol, a non-negative number, a function call or already bracketed
#' is put in parentheses, so `x/a` with `a = "b+c"` becomes `x/(b+c)`. Only
#' variables are replaced, not function names, argument names or list
#' elements; the rest of the text is left as it is.
#'
#' @param what Character vector, the symbols to replace.
#' @param by Character or numeric vector, the replacements, recycled to the
#'   length of `what`.
#' @param x Character vector of expressions, one per element. Names and
#'   attributes are kept; `NA` and `"0"` are returned unchanged.
#' @return `x` with the symbols replaced.
#' @examples
#' replaceSymbols("a", "b+c", c(r1 = "x/a", r2 = "exp(a)*a^2"))
#' replaceSymbols(c("a", "b"), c("b", "a"), "a - b")
#' @export
replaceSymbols <- function(what, by, x) {
  if (!length(x) || !length(what)) return(x)
  by <- rep(as.character(by), length.out = length(what))
  names(by) <- what

  todo <- which(!is.na(x) & x != "0")
  if (!length(todo)) return(x)

  bracketed <- vapply(by, .bracketCompound, character(1))
  src <- gsub("\t", " ", unclass(x)[todo])
  if (any(grepl("\n", src, fixed = TRUE)))
    stop("replaceSymbols(): expressions must be single lines.", call. = FALSE)
  pd <- utils::getParseData(parse(text = src, keep.source = TRUE))
  pd <- pd[pd$terminal, c("line1", "col1", "col2", "token", "text")]
  pd <- pd[order(pd$line1, pd$col1), ]
  # the name after `$` or `@` is an element, not a variable
  before <- c("", pd$token[-nrow(pd)])
  next_  <- c(pd$token[-1], "")
  isHit  <- pd$token == "SYMBOL" & pd$text %in% what & !before %in% c("'$'", "'@'")
  # a symbol that is the whole expression or alone in parentheses needs no bracket
  single   <- pd$line1 %in% which(tabulate(pd$line1, length(src)) == 1L)
  pd$alone <- single | (before == "'('" & next_ == "')'")
  hit <- pd[isHit, ]

  out <- src
  for (ln in unique(hit$line1)) {
    h <- hit[hit$line1 == ln, ]
    s <- src[ln]
    for (k in rev(seq_len(nrow(h))))
      s <- paste0(substr(s, 1, h$col1[k] - 1),
                  if (h$alone[k]) by[[h$text[k]]] else bracketed[[h$text[k]]],
                  substr(s, h$col2[k] + 1, nchar(s)))
    out[ln] <- s
  }
  x[todo] <- out
  x
}

# A replacement that binds as one operand stays as it is, anything else is
# bracketed.
.bracketCompound <- function(expr) {
  e <- tryCatch(parse(text = expr, keep.source = FALSE), error = function(err) NULL)
  if (length(e) != 1L) return(expr)
  e <- e[[1]]
  if (is.name(e) || (is.numeric(e) && e >= 0)) return(expr)
  if (is.call(e) && is.name(e[[1]])) {
    fn <- as.character(e[[1]])
    if (fn == "(" || grepl("^[A-Za-z.][A-Za-z0-9._]*$", fn)) return(expr)
  }
  paste0("(", expr, ")")
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



## getElements ------------------------------------

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




## blockdiagSymb ----------------------------------------

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
