# Stoichiometric coefficients of one side of a reaction, a sum of terms
# `species` or `<numeric expression>*species`, named by species in order of appearance
.sideCoefficients <- function(side) {
  side <- trimws(side)
  if (!nzchar(side) || side == "0") return(numeric(0))
  terms <- function(e) {
    if (is.call(e) && identical(e[[1]], quote(`+`)) && length(e) == 3L)
      return(c(terms(e[[2]]), terms(e[[3]])))
    list(e)
  }
  out <- numeric(0)
  for (t in terms(str2lang(side))) {
    coef <- 1
    if (is.call(t) && identical(t[[1]], quote(`*`))) {
      coef <- t[[2]]
      t <- t[[3]]
    }
    if (!is.name(t) || length(all.vars(coef)) || identical(t, quote(pi)))
      stop("'", side, "' is not a sum of terms 'coefficient*species'.", call. = FALSE)
    nm <- as.character(t)
    out[nm] <- sum(out[nm], eval(coef, baseenv()), na.rm = TRUE)
  }
  out
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
#' @seealso [replaceSymbols()]
#' @examples
#' getSymbols(c("k1*A/(Km + A)", "exp(-k2*time)"))
#' getSymbols("k1*A + k2*B", exclude = "A")
#' @export
getSymbols <- function(char, exclude = NULL) {
  if (is.null(char)) return(NULL)
  pd <- utils::getParseData(parse(text = char[!char %in% "0"], keep.source = TRUE))
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







