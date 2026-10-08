# setdiff/intersect for character vectors that are already unique. base's
# versions coerce with as.vector() and run a unique() pass; on parameter-name
# vectors touched once per condition that is most of their cost.
.setdiffU <- function(a, b) a[match(a, b, 0L) == 0L]
.intersectU <- function(a, b) a[match(a, b, 0L) > 0L]

# Moves arguments given under a deprecated name, `renames` maps old to new,
# onto the current name in the calling frame with one warning each; both
# names together are an error. Returns `dots` without the old names.
.renameArgs <- function(dots, renames, label, strict = FALSE,
                        env = parent.frame()) {
  old <- intersect(names(dots), names(renames))
  for (o in old) {
    new <- renames[[o]]
    if (!eval(call("missing", as.name(new)), env))
      stop(sprintf("%s: give '%s' only; '%s' is its deprecated name.",
                   label, new, o), call. = FALSE)
    warning(sprintf("%s: '%s' is deprecated, use '%s'.", label, o, new),
            call. = FALSE)
    assign(new, dots[[o]], envir = env)
  }
  if (length(old)) dots <- dots[!names(dots) %in% old]
  if (strict && length(dots)) {
    nm <- names(dots)
    if (is.null(nm)) nm <- rep("", length(dots))
    stop(sprintf("%s: unused argument(s) %s.", label,
                 paste(ifelse(nzchar(nm), nm, "<unnamed>"), collapse = ", ")),
         call. = FALSE)
  }
  dots
}

# Arguments that no longer have an effect: each warns once, anything else left
# in `dots` is an error.
.droppedArgs <- function(dots, dropped, label) {
  nm <- names(dots)
  if (is.null(nm)) nm <- rep("", length(dots))
  for (d in intersect(nm, dropped))
    warning(sprintf("%s: '%s' is deprecated and ignored.", label, d),
            call. = FALSE)
  .renameArgs(dots[!nm %in% dropped], character(0), label, strict = TRUE)
  invisible(NULL)
}

## utils.R, general-purpose utility functions

#' Compare Two Objects and Return Differences
#'
#' Works on two objects or on a list of objects. A list is compared entry by
#' entry with a reference entry, including the attributes `"equations"`,
#' `"parameters"`, `"events"` and `"forcings"` of the entries.
#'
#' @param vec1 An [eqnvec], character vector or data frame, or a list of such
#'   objects or of objects with the attributes above.
#' @param vec2 Same as `vec1`. Not used if `vec1` is a list.
#' @param reference Integer, the reference entry of a list. Defaults to 1.
#' @param ... Arguments passed to the methods.
#' @return A data frame of the entries that differ, are missing in `vec2`
#'   or are additional in `vec2`, or `NULL` if there are none. For a list, a
#'   list of those per compared object and attribute.
#' 
#' @export
#' @examples
#' ## Compare equation vectors
#' eq1 <- eqnvec(a = "-k1*a + k2*b", b = "k2*a - k2*b")
#' eq2 <- eqnvec(a = "-k1*a", b = "k2*a - k2*b", c = "k2*b")
#' compare(eq1, eq2)
#' 
#' ## Compare character vectors
#' c1 <- c("a", "b")
#' c2 <- c("b", "c")
#' compare(c1, c2)
#' 
#' ## Compare data.frames
#' d1 <- data.frame(var = "a", time = 1, value = 1:3, method = "replace")
#' d2 <- data.frame(var = "a", time = 1, value = 2:4, method = "replace")
#' compare(d1, d2)
#' 
#' ## Compare structures like prediction functions
#' fn1 <- function(x) x^2
#' attr(fn1, "equations") <- eq1
#' attr(fn1, "parameters") <- c1
#' attr(fn1, "events") <- d1
#' 
#' fn2 <- function(x) x^3
#' attr(fn2, "equations") <- eq2
#' attr(fn2, "parameters") <- c2
#' attr(fn2, "events") <- d2
#' 
#' mylist <- list(f1 = fn1, f2 = fn2)
#' compare(mylist)

compare <- function(vec1, ...) {
  UseMethod("compare", vec1)
}

#' @export
#' @rdname compare
compare.list <- function(vec1, vec2 = NULL, reference = 1, ...) {
  
  index <- (1:length(vec1))[-reference]
  diffable.attributes <- c("equations", "parameters", "forcings", "events")
  
  
  out.total <- lapply(index, function(i) {
    
    # Compare objects if possible
    vec1.inner <- vec1[[reference]]
    vec2.inner <- vec1[[i]]
    out1 <- NULL
    if(any(class(vec1.inner) %in% c("eqnvec", "data.frame"))) {
      out1 <- list(compare(vec1.inner, vec2.inner))
      names(out1) <- "object"
    }
      
    # Compare comparable attributes of the object if available
    out2 <- NULL
    attributes1 <- attributes(vec1.inner)[diffable.attributes]
    attributes2 <- attributes(vec2.inner)[diffable.attributes]
    slots <- names(attributes1)[!is.na(names(attributes1))]
    out2 <- lapply(slots, function(n) {
      compare(attributes1[[n]], attributes2[[n]])
    })
    names(out2) <- slots
    
    c(out1, out2)
    
  })
  names(out.total) <- names(vec1)[index]
  
  ## Do resorting of the list
  innernames <- names(out.total[[1]])
  out.total <- lapply(innernames, function(n) {
    out <- lapply(out.total, function(out) out[[n]])
    out[!sapply(out, is.null)]
  })
  names(out.total) <- innernames
  
  
  return(out.total)
  
  
  
}

#' @export
#' @rdname compare
compare.character <- function(vec1, vec2 = NULL, ...) {
  missing <- setdiff(vec1, vec2)
  additional <- setdiff(vec2, vec1)
  
  out <- do.call(rbind, 
          list(different = NULL, 
               missing = data.frame(name = missing), 
               additional = data.frame(name = additional)
          )
  )
  
  if(nrow(out) == 0) out <- NULL
  return(out)
  
  
  
}

#' @export
#' @rdname compare
compare.eqnvec <- function(vec1, vec2 = NULL, ...) {

  names1 <- names(vec1)
  names2 <- names(vec2)
  
  missing <- setdiff(names1, names2)
  additional <- setdiff(names2, names1)
  joint <- intersect(names1, names2)
  
  # Compare joint equations
  v1 <- format(vec1)
  v2 <- format(vec2)
  not.coincide <- which(as.character(v1[joint]) != as.character(v2[joint]))
  
  different <- data.frame(name = names(v2[not.coincide]), equation = as.character(v2[not.coincide]))
  missing <- data.frame(name = names(v2[missing]), equation = as.character(v2[missing]))
  additional <- data.frame(name = names(v2[additional]), equation = as.character(v2[additional]))
  
  out <- do.call(rbind, list(different = different, missing = missing, additional = additional))
  if(nrow(out) == 0) out <- NULL
  return(out)
  
  
}

#' @export
#' @rdname compare
compare.data.frame <- function(vec1, vec2 = NULL, ...) {
  
  additional <- !duplicated(rbind(vec1, vec2))[-(1:nrow(vec1))]
  missing <- !duplicated(rbind(vec2, vec1))[-(1:nrow(vec2))]
  
  out <- do.call(rbind, list(different = character(0), missing = vec1[missing, ], additional = vec2[additional, ]))
  if(nrow(out) == 0) out <- NULL
  return(out)
  
  
}



#' Alternative version of expand.grid
#' @param seq1 Vector, numeric or character
#' @param seq2 Vector, numeric or character
#' @return Matrix of combinations of elements of `seq1` and `seq2`
#' @noRd
expand.grid.alt <- function(seq1, seq2) {
  cbind(Var1=rep.int(seq1, length(seq2)), Var2=rep(seq2, each=length(seq1)))
}




#' List dMod Objects in an Environment
#'
#' @param classlist Character vector of classes. Defaults to `"odemodel"`,
#'   `"parfn"`, `"prdfn"`, `"obsfn"`, `"objfn"` and `"datalist"`.
#' @param envir Environment to search. Defaults to the global environment.
#' @return Named character vector of the names of the objects of each class;
#'   the names of the vector are the classes, numbered when a class has more
#'   than one object.
#' @examples
#' env <- new.env()
#' assign("data", as.datalist(data.frame(name = "A", time = 0, value = 1,
#'                                       sigma = 1, condition = "C1")),
#'        envir = env)
#' lsdMod(envir = env)
#'
#' @export
lsdMod <- function(classlist = c("odemodel", "parfn", "prdfn", "obsfn", "objfn", "datalist"), envir = .GlobalEnv){
  glist <- as.list(envir)
  out <- list()
  for (a in classlist) {
    flist <- which(sapply(glist, function(f) any(class(f) == a)))
    out[[a]] <- names(glist[flist])
    #cat(a,": ")
    #cat(paste(out[[a]], collapse = ", "),"\n")
  }
  
  unlist(out)
  
  
  
}




#' Select Attributes
#'
#' Keeps or removes attributes of an object.
#'
#' @param x An object.
#' @param which Character vector of attribute names. `NULL` (default) stands for
#'   `"class"`, `"dim"`, `"dimnames"`, `"names"`, `"col.names"` and
#'   `"row.names"`.
#' @param keep Logical. `TRUE` (default) keeps the attributes in `which` and
#'   removes all others; `FALSE` removes those in `which`.
#' @param ... `atr` is deprecated, use `which`.
#'
#' @return `x` with the selected attributes.
#' @examples
#' x <- structure(1:3, names = c("a", "b", "c"), note = "temporary")
#' attrs(x)
#' attrs(x, "note", keep = FALSE)
#'   
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#' @author Mirjam Fehling-Kaschek, \email{mirjam.fehling@@physik.uni-freiburg.de}
#'   
#' @export
attrs <- function(x, which = NULL, keep = TRUE, ...) {

  .renameArgs(list(...), c(atr = "which"), "attrs", strict = TRUE)

  if (is.null(which)) {
    which <- c("class", "dim", "dimnames", "names", "col.names", "row.names")
  }
  
  xattr <- names(attributes(x))
  if (keep == TRUE) {
    attributes(x)[!xattr %in% which] <- NULL
  } else {
    attributes(x)[xattr %in% which] <- NULL
  }
  
  return(x)
}





# Match a `derivMode` argument against what an entry point supports. Mirrors
# cppDE's helper of the same shape: "forward" and "reverse" are separate build
# products and combine.
.matchDerivMode <- function(x, choices) {
  if ("symbolic" %in% x)
    stop('derivMode = "symbolic" is gone: every derivative is generated by ',
         'AD at code-generation time. Use "forward".', call. = FALSE)
  x <- unique(match.arg(x, choices, several.ok = TRUE))
  if (!length(x))
    stop("'derivMode' must name at least one of ",
         paste0('"', choices, '"', collapse = ", "), call. = FALSE)
  choices[choices %in% x]
}
