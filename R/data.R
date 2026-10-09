#' Time-Course Data for the JAK-STAT Signaling Pathway
#'
#' Phosphorylated Epo receptor (`pEpoR`), total phosphorylated STAT in the
#' cytoplasm (`tpSTAT`) and total STAT in the cytoplasm (`tSTAT`), measured
#' between minute 0 and 60 after stimulation with Epo.
#'
#' @format A data frame with 47 rows and the columns `time` (minutes),
#'   `name`, `value`, `sigma` and `condition` (`"Epo"`).
#' @source Swameye I, Mueller TG, Timmer J, Sandra O, Klingmueller U (2003),
#'   Identification of nucleocytoplasmic cycling as a remote sensor in cellular
#'   signaling by databased modeling. PNAS 100(3):1028-1033.
#' @seealso [as.datalist()]
#' @examples
#' data(jakstat)
#' plotData(as.datalist(jakstat))
#' @name jakstat
#' @docType data
#' @keywords data
NULL


#' Time-Course Data for the Bile Acid Demonstration Model
#'
#' Bile acid amounts in the buffer (`buffer`) and in the cells (`cellular`)
#' between 0.1 and 41 time units, in two conditions, `"closed"` and `"open"`.
#'
#' @format A data frame with 32 rows and the columns `name`, `time`, `value`,
#'   `sigma`, `lloq` (lower limit of quantification, `-Inf` throughout) and
#'   `condition`.
#' @seealso `system.file("examples", "BA_transport.R", package = "dMod2")`
#'   builds and fits the model.
#' @examples
#' data(badata)
#' plotData(as.datalist(badata))
#' @name badata
#' @docType data
#' @keywords data
NULL


#' Time-Course Data for JAK2-STAT5 Signaling in CFU-E Cells
#'
#' 541 measurements from thirteen experiments on Epo stimulated erythroid
#' progenitor cells: phosphorylated JAK2 and Epo receptor, total and
#' phosphorylated STAT5, the feedback proteins CIS, SOCS3 and SHP1, their
#' transcripts, and dose responses at four times.
#'
#' @format A data frame with 541 rows and 11 columns:
#' \describe{
#'   \item{`time`}{minutes after stimulation.}
#'   \item{`name`}{observable, 20 levels.}
#'   \item{`value`}{signal divided by the maximum of its observable.}
#'   \item{`sigma`}{`NA` throughout; the errors are estimated in the fit.}
#'   \item{`condition`}{one of 36 experimental conditions.}
#'   \item{`experiment`}{one of 13 experiments; selects the scale and offset
#'     parameters.}
#'   \item{`epo_level`}{Epo dose.}
#'   \item{`ActD`}{1 where actinomycin D was added, else 0. The PEtab copy
#'     of the problem codes it the other way round.}
#'   \item{`CISoe`, `SOCS3oe`, `SHP1oe`}{1 where CIS, SOCS3 or SHP1 is
#'     overexpressed, else 0.}
#' }
#'
#' @source Bachmann J, Raue A, Schilling M, Boehm ME, Kreutz C, Kaschek D,
#' Busch H, Gretz N, Lehmann WD, Timmer J, Klingmueller U (2011),
#' Division of labor by dual feedback regulators controls JAK2/STAT5
#' signaling over broad ligand range. Mol Syst Biol 7:516. Values as
#' distributed in the PEtab benchmark collection, model
#' `Bachmann_MSB2011`, which is installed with the package, see
#' `system.file("extdata", "petab_bachmann", package = "dMod2")`.
#' @seealso `system.file("examples", "example_BachmannMSB2011.R", package =
#'   "dMod2")` builds the model.
#' @examples
#' data(bachmann)
#' head(bachmann)
#' @name bachmann
#' @docType data
#' @keywords data
NULL


#' Time-Course Data for STAT5 Dimerization in BaF3-EpoR Cells
#'
#' Relative amounts of phosphorylated STAT5A (`pSTAT5A_rel`), phosphorylated
#' STAT5B (`pSTAT5B_rel`) and of STAT5A within the total STAT5 pool
#' (`rSTAT5A_rel`), measured at 16 time points after stimulation with Epo.
#'
#' @format A data frame with 48 rows and the columns `time` (minutes),
#'   `name`, `value`, `sigma` (`NA` throughout; the errors are estimated in
#'   the fit) and `condition` (`"Boehm2014"`).
#' @source Boehm ME, Adlung L, Schilling M, Roth S, Klingmueller U, Lehmann WD
#' (2014), Identification of isoform-specific dynamics in
#' phosphorylation-dependent STAT5 dimerization by quantitative mass
#' spectrometry and mathematical modeling. J Proteome Res 13(12):5685-5694.
#' Values as distributed in the PEtab benchmark collection, model
#' `Boehm_JProteomeRes2014`.
#' @seealso `system.file("examples", "example_Boehm_JProteomeRes2014.R",
#'   package = "dMod2")` builds the model.
#' @examples
#' data(boehm)
#' plotData(as.datalist(boehm))
#' @name boehm
#' @docType data
#' @keywords data
NULL



## combine ----------------------------------------------

#' Combine Data Frames or Matrices by Rows
#'
#' Binds data frames or matrices by rows after adding the columns each one
#' lacks: `NA` for data frames, 0 for matrices.
#'
#' @param ... Data frames or matrices, with column names that need not
#'   overlap. `NULL` entries are dropped.
#'
#' @return A data frame or matrix with the union of the column names.
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




## wide2long --------------------------------------------

#' Translate Wide Format into Long Format
#'
#' Converts output in wide format, e.g. of a prediction, into long format with
#' one row per value.
#'
#' @param out A data frame, matrix, or list of those, in wide format.
#' @param keep Integer vector, the columns kept as identifiers. Defaults to 1,
#'   the time column.
#' @param na.rm Logical. `TRUE` removes rows with missing values. Defaults to
#'   `FALSE`.
#'
#' @return A data frame with the kept columns and the columns
#' \describe{
#'   \item{`name`}{the column names of the other columns of `out`.}
#'   \item{`value`}{their values.}
#'   \item{`condition`}{for a list, its names, or the list positions if it
#'     has none.}
#' }
#' @seealso [long2wide()]
#' @examples
#' out <- cbind(time = 0:2, A = c(1, 0.5, 0.25), B = c(0, 0.5, 0.75))
#' wide2long(out)
#' wide2long(list(C1 = out, C2 = out))
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
  
  # An unnamed list still needs a condition column; without one `cbind` would
  # drop it and every caller reading `$condition` would come up empty.
  conditions <- names(out)
  if (is.null(conditions)) conditions <- as.character(seq_along(out))
  
  outlong <- do.call(rbind, lapply(seq_along(out), function(cond) {
    
    cbind(wide2long.matrix(out[[cond]], keep = keep, na.rm = na.rm),
          condition = conditions[cond])
    
  }))
  
  
  
  return(outlong)
  
}




## long2wide --------------------------------------------

#' Translate Long Format into Wide Format
#'
#' Inverse of [wide2long()] for a single condition.
#'
#' @param out Data frame in long format with time, name and value as its
#'   first three columns, each name having a value at every time.
#' @return A matrix in wide format, one column per name after the time column.
#' @seealso [wide2long()]
#' @examples
#' out <- cbind(time = 0:2, A = c(1, 0.5, 0.25), B = c(0, 0.5, 0.75))
#' long2wide(wide2long(out))
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




## lbind ------------------------------------------------

#' Bind a Named List of Data Frames into One Data Frame
#'
#' Adds to each data frame a column `condition` with its name in the list and
#' binds them by rows.
#'
#' @param mylist A named list of data frames with the same columns.
#' @return A data frame with the original columns and `condition`.
#' @examples
#' lbind(list(C1 = data.frame(x = 1:2), C2 = data.frame(x = 3)))
#' @export
lbind <- function(mylist) {
  
  conditions <- names(mylist)
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
