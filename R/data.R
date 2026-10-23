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
#' @seealso `demo("BA_transport", package = "dMod2")` builds and fits the
#'   model.
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
#' @seealso `demo("BachmannMSB2011", package = "dMod2")` builds the model.
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
#' @seealso `demo("Boehm_JProteomeRes2014", package = "dMod2")` builds the
#'   model.
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

  mylist <- Filter(Negate(is.null), list(...))
  if (!length(mylist)) return(NULL)
  mynames <- unique(unlist(lapply(mylist, colnames)))

  if (all(vapply(mylist, is.matrix, TRUE))) {
    out <- do.call(rbind, lapply(mylist, function(l) {
      full <- matrix(0, nrow(l), length(mynames), dimnames = list(rownames(l), mynames))
      full[, colnames(l)] <- l
      full
    }))
    return(out)
  }

  mylist <- lapply(mylist, function(l) {
    l <- as.data.frame(l, stringsAsFactors = FALSE)
    l[] <- lapply(l, function(column) if (is.factor(column)) as.character(column) else column)
    l
  })
  out <- as.data.frame(data.table::rbindlist(mylist, use.names = TRUE, fill = TRUE))[mynames]
  # Row names as rbind() gives them: automatic ones are numbered, the others
  # made unique.
  if (!all(vapply(mylist, function(l) .row_names_info(l) < 0, TRUE))) {
    rn <- unlist(lapply(mylist, function(l)
      if (.row_names_info(l) < 0) rep(NA_character_, nrow(l)) else rownames(l)))
    rn[is.na(rn)] <- which(is.na(rn))
    rownames(out) <- make.unique(rn, sep = "")
  }
  out
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
  if (anyDuplicated(allnames))
    stop("Duplicated column names: ", paste(unique(allnames[duplicated(allnames)]), collapse = ", "),
         call. = FALSE)
  times <- out[, keep, drop = FALSE]
  ntimes <- nrow(out)
  values <- unlist(out[,allnames])
  outlong <- data.frame(times,
                        name = factor(rep(allnames, each = ntimes), levels = allnames),
                        value = as.numeric(values), row.names = NULL)
  colnames(outlong)[1:length(keep)] <- timenames
  
  if (na.rm) outlong <- outlong[!is.na(outlong$value),]
  
  return(outlong)
  
}

#' @rdname wide2long
#' @export
wide2long.list <- function(out, keep = 1, na.rm = FALSE) {
  # An unnamed list still gets a condition column, which callers read.
  conditions <- names(out)
  if (is.null(conditions)) conditions <- as.character(seq_along(out))
  pieces <- lapply(out, wide2long.matrix, keep = keep, na.rm = na.rm)
  outlong <- as.data.frame(data.table::rbindlist(pieces))
  outlong$condition <- rep(conditions, vapply(pieces, nrow, 1L))
  outlong
}








## lbind ------------------------------------------------

#' Bind a Named List of Data Frames into One Data Frame
#'
#' Adds to each data frame a column `condition` with its name in the list and
#' binds them by rows.
#'
#' @param mylist A named list of data frames; a column missing in one of them
#'   is `NA` there.
#' @return A data frame with the original columns and `condition`.
#' @examples
#' lbind(list(C1 = data.frame(x = 1:2), C2 = data.frame(x = 3)))
#' @export
lbind <- function(mylist) {
  if (!length(mylist)) return(NULL)
  conditions <- names(mylist)
  if (is.null(conditions)) conditions <- as.character(seq_along(mylist))
  out <- as.data.frame(data.table::rbindlist(mylist, use.names = TRUE, fill = TRUE))
  out$condition <- rep(conditions, vapply(mylist, nrow, 1L))
  out
}


#' Time-Course Data for IFN-alpha Signalling in Huh7.5 Cells
#'
#' 1918 measurements on human hepatoma cells stimulated with IFN-alpha at 0 h,
#' at 24 h or at both: phosphorylated and total STAT1 and STAT2 and IRF9 in the
#' cytoplasm and the nucleus, USP18, seven target gene transcripts and the
#' molecules per cell of STAT1, STAT2, IRF9 and USP18. Two experiments
#' perturb the USP18 feedback, by knock-down and by overexpression.
#'
#' @format A data frame with 1918 rows and 8 columns:
#' \describe{
#'   \item{`time`}{hours after the start of the experiment, when the first
#'     dose and the doxycycline are given.}
#'   \item{`name`}{observable, 22 levels. Immunoblots (`*_Cyt`, `*_Nuc`) and
#'     qPCR (`*_mRNA`, fold change against untreated cells) are on the
#'     natural log, molecules per cell (`mc_*`) on the linear scale.}
#'   \item{`value`}{the measurement, blots scaled across replicates.}
#'   \item{`sigma`}{standard deviation from the combined scaling and error
#'     model.}
#'   \item{`condition`}{one of 25 conditions, `<experiment>_<ifn1>_<ifn2>`,
#'     for example `"WT_10_5000"` for 10 U/ml at 0 h and 5000 U/ml at 24 h.}
#'   \item{`experiment`}{`"WT"`, `"USP18kd"` (USP18 knock-down by siRNA) or
#'     `"USP18oe"` (USP18 overexpression, induced by doxycycline from 0 h on).}
#'   \item{`ifn1`, `ifn2`}{IFN-alpha dose at 0 h and at 24 h in U/ml.}
#' }
#'
#' @source Kok F, Rosenblatt M, Teusel M, et al. (2020), Disentangling
#' molecular mechanisms regulating sensitization of interferon alpha signal
#' transduction. Mol Syst Biol 16:e8955, \doi{10.15252/msb.20198955}. Values
#' as fitted for the Huh7.5 core model of the publication.
#' @seealso `demo("KokMSB2020", package = "dMod2")` builds the model.
#' @examples
#' data(kok)
#' head(kok)
#' @name kok
#' @docType data
#' @keywords data
NULL
