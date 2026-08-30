#' Time-course data for the JAK-STAT cell signaling pathway
#'
#' Phosphorylated Epo receptor (pEpoR), phosphorylated STAT in the
#' cytoplasm (tpSTAT) and total STAT (tSTAT) in the cytoplasmhave been 
#' measured at times 0, ..., 60.
#'
#' @name jakstat
#' @docType data
#' @keywords data
NULL


#' Time-course data for the Bile-Acid demonstration model
#'
#' @name badata
#' @docType data
#' @keywords data
NULL


#' Time-course data for STAT5 dimerisation in BaF3-EpoR cells
#'
#' Relative amounts of phosphorylated STAT5A (`pSTAT5A_rel`), phosphorylated
#' STAT5B (`pSTAT5B_rel`) and of STAT5A within the total STAT5 pool
#' (`rSTAT5A_rel`), measured at 16 time points after stimulation with Epo.
#' `sigma` is `NA` throughout: the measurement uncertainties are estimated
#' along with the model parameters.
#'
#' @source Boehm ME, Adlung L, Schilling M, Roth S, Klingmueller U, Lehmann WD
#' (2014), Identification of isoform-specific dynamics in
#' phosphorylation-dependent STAT5 dimerization by quantitative mass
#' spectrometry and mathematical modeling. J Proteome Res 13(12):5685-5694.
#' Values as distributed in the PEtab benchmark collection, model
#' `Boehm_JProteomeRes2014`.
#' @name boehm
#' @docType data
#' @keywords data
NULL



## combine (moved from tools.R) ----------------------------------------------

#' Combine several data.frames by rowbind
#' 
#' @param ... data.frames or matrices with not necessarily overlapping colnames
#' @details This function is useful when separating models into independent csv model files,
#' e.g.~a receptor model and several downstream pathways. Then, the models can be recombined 
#' into one model by `combine()`.
#' 
#' @return A `data.frame`
#' @export
#' @examples
#' data1 <- data.frame(Description = "reaction 1", Rate = "k1*A", A = -1, B = 1)
#' data2 <- data.frame(Description = "reaction 2", Rate = "k2*B", B = -1, C = 1)
#' combine(data1, data2)
#' @export
combine <- function(...) {
  
  # List of input data.frames
  mylist <- list(...)
  # Remove empty slots
  is.empty <- sapply(mylist, is.null)
  mylist <- mylist[!is.empty]
  
  mynames <- unique(unlist(lapply(mylist, function(S) colnames(S))))
  
  mylist <- lapply(mylist, function(l) {
    
    if(is.data.frame(l)) {
      i <- sapply(l, is.factor)
      l[i] <- lapply(l[i], as.character)
      present.list <- as.list(l)
      missing.names <- setdiff(mynames, names(present.list))
      missing.list <- structure(as.list(rep(NA, length(missing.names))), names = missing.names)
      combined.data <- do.call(function(...) cbind.data.frame(..., stringsAsFactors = FALSE), c(present.list, missing.list))
      rownames(combined.data) <- rownames(l)
    }
    if(is.matrix(l)) {
      present.matrix <- as.matrix(l)
      missing.names <- setdiff(mynames, colnames(present.matrix))
      missing.matrix <- matrix(0, nrow = nrow(present.matrix), ncol = length(missing.names), 
                             dimnames = list(NULL, missing.names))
      combined.data <- submatrix(cbind(present.matrix, missing.matrix), cols = mynames)
      rownames(combined.data) <- rownames(l)
    }
    
    return(combined.data)
  })
  
  out <- do.call(rbind, mylist)
  
  return(out)
  
  
}




## wide2long (moved from tools.R) --------------------------------------------

#' Translate wide output format (e.g., from ODE solver) into long format
#'
#' Converts simulation output in wide format into a tidy long format suitable for
#' plotting or further analysis (e.g., with \pkg{ggplot2}). The function assumes
#' that the first column of \code{out} represents a time-like variable and the
#' remaining columns contain values.
#'
#' @param out A \code{data.frame}, \code{matrix}, or a \code{list} of matrices in wide format.
#' @param keep Integer vector specifying the column indices to keep (default is \code{1}).
#' @param na.rm Logical. If \code{TRUE}, missing values are removed in the long-format output.
#'
#' @details
#' If \code{out} is a list, the list names are added as an additional column named
#' \code{"condition"}. This is particularly useful for plotting results from multiple
#' simulation conditions with \pkg{ggplot2}.
#'
#' @return A \code{data.frame} in long format with the following columns:
#' \itemize{
#'   \item \code{"time"} -- values from \code{out[, 1]}.
#'   \item \code{"name"} -- column names from \code{out[, -1]}.
#'   \item \code{"value"} -- corresponding numeric values.
#'   \item \code{"condition"} -- if \code{out} was a list, contains the list names.
#' }
#'
#' @export
wide2long <- function(out, keep = 1, na.rm = FALSE) {
  
  UseMethod("wide2long", out)
  
  
}

#' @rdname wide2long
#' @export
wide2long.data.frame <- function(out, keep = 1, na.rm = FALSE) {
  
  wide2long.matrix(out, keep = keep, na.rm = na.rm)
  
}

#' @rdname wide2long
#' @export
wide2long.matrix <- function(out, keep = 1, na.rm = FALSE) {
  
  timenames <- colnames(out)[keep]
  allnames <- colnames(out)[-keep]
  if (any(duplicated(allnames))) warning("Found duplicated colnames in out. Duplicates were removed.")
  times <- out[,keep]
  ntimes <- nrow(out)
  values <- unlist(out[,allnames])
  outlong <- data.frame(times, 
                        name = factor(rep(allnames, each = ntimes), levels = allnames), 
                        value = as.numeric(values))
  colnames(outlong)[1:length(keep)] <- timenames
  
  if (na.rm) outlong <- outlong[!is.na(outlong$value),]
  
  return(outlong)
  
}

#' @rdname wide2long
#' @export
wide2long.list <- function(out, keep = 1, na.rm = FALSE) {
  
  conditions <- names(out)
  
  outlong <- do.call(rbind, lapply(1:max(c(length(conditions), 1)), function(cond) {
    
    cbind(wide2long.matrix(out[[cond]]), condition = conditions[cond])
    
  }))
  
  
  
  return(outlong)
  
}




## long2wide (moved from tools.R) --------------------------------------------

#' Translate long to wide format (inverse of wide2long.matrix) 
#' 
#' @param out data.frame in long format 
#' @return data.frame in wide format 
#' @export
long2wide <- function(out) {
  
  timename <- colnames(out)[1]
  times <- unique(out[,1])
  allnames <- unique(as.character(out[,2]))
  M <- matrix(out[,3], nrow=length(times), ncol=length(allnames))
  M <- cbind(times, M)
  colnames(M) <- c(timename, allnames)
  
  return(M)
  
}




## lbind (moved from tools.R) ------------------------------------------------

#' Bind named list of data.frames into one data.frame
#' 
#' @param mylist A named list of data.frame. The data.frames are expected to have the same structure.
#' @details Each data.frame ist augented by a "condition" column containing the name attributed of
#' the list entry. Subsequently, the augmented data.frames are bound together by `rbind`.
#' @return data.frame with the originial columns augmented by a "condition" column.
#' @export
lbind <- function(mylist) {
  
  conditions <- names(mylist)
  #numconditions <- suppressWarnings(as.numeric(conditions))
  #
  # if(!any(is.na(numconditions))) 
  #   numconditions <- as.numeric(numconditions) 
  # else 
  numconditions <- conditions

  
  outlong <- do.call(rbind, lapply(1:length(conditions), function(cond) {
    
    myout <- mylist[[cond]]
    if (nrow(myout) > 0)
      myout[["condition"]] <- numconditions[cond]
    else
      myout[["condition"]] <- character(0)
    
    return(myout)
    
  }))
  
  return(outlong)
  
}
