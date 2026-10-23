
## Class "eqnlist" and its constructor ------------------------------------------


#' Generate or Coerce to an Equation List
#' @description An eqnlist stores an ODE as a stoichiometric matrix, rate
#' expressions, state names and compartment information. `eqnlist()` builds
#' one from these parts, `as.eqnlist()` from a data frame of reactions.
#' @param ... not used.
#' @details The data frame for `as.eqnlist()` contains the columns
#' "Description" (character), "Rate" (character) and one column per state,
#' the columns of the stoichiometric matrix. Layout attributes of the data
#' frame, as set by [as.data.frame.eqnlist()], are used unless `volumes` is
#' given.
#' @rdname eqnlist
#' @export
as.eqnlist <- function(data, volumes, ...) {
  UseMethod("as.eqnlist", data)
}

#' @export
#' @param data data frame with the columns "Description", "Rate" and one
#' column per state.
#' @rdname eqnlist
as.eqnlist.data.frame <- function(data, volumes = NULL, compartments = NULL, compartmentOf = NULL,
                                   reactionCompartment = NULL, amountStates = NULL, ...) {
  description <- as.character(data$Description)
  rates <- as.character(data$Rate)
  states <- setdiff(colnames(data), c("Description", "Rate"))
  smatrix <- as.matrix(data[, states]); colnames(smatrix) <- states

  # An explicit `volumes` states the layout and is not overruled by a layout
  # the data.frame brings along as attributes.
  inherit_layout <- is.null(volumes)
  if (is.null(volumes))             volumes             <- attr(data, "volumes")
  if (is.null(compartments) && inherit_layout)  compartments  <- attr(data, "compartments")
  if (is.null(compartmentOf) && inherit_layout) compartmentOf <- attr(data, "compartmentOf")
  if (is.null(reactionCompartment)) reactionCompartment <- attr(data, "reactionCompartment")
  if (is.null(amountStates))        amountStates        <- attr(data, "amountStates")

  eqnlist(smatrix, states, rates, volumes, description,
          compartments = compartments, compartmentOf = compartmentOf,
          reactionCompartment = reactionCompartment, amountStates = amountStates)

}


#' @export
#' @rdname eqnlist
#' @param x object to test.
is.eqnlist <- function(x) {

  required <- c("smatrix", "states", "rates", "volumes", "description",
                "compartments", "compartmentOf", "reactionCompartment")

  # A reaction-less list is valid and may already have a compartment layout.
  if (is.null(x[["smatrix"]])) {
    empty <- !length(x$states) && !length(x$rates) && !length(x$description) &&
      is.null(x$volumes) && is.null(x$reactionCompartment)
    hasComp <- !is.null(x$compartments)
    if (hasComp != !is.null(x$compartmentOf)) return(FALSE)
    return(empty && (!hasComp || all(x$compartmentOf %in% names(x$compartments))))
  }

  refs <- !is.null(x$compartments) && !is.null(x$compartmentOf) &&
          all(x$states %in% names(x$compartmentOf)) &&
          all(x$compartmentOf %in% names(x$compartments))

  # `$volumes` is derived, so an in-place edit of `$compartments` leaves it stale.
  volumes <- !refs ||
    identical(unname(x$volumes),
              unname(.derivedVolumes(x$compartments, x$compartmentOf, x$states)))

  reactions <- is.null(x$reactionCompartment) ||
    (length(x$reactionCompartment) == length(x$rates) &&
     all(is.na(x$reactionCompartment) |
         x$reactionCompartment %in% names(x$compartments)))

  amounts <- is.null(x$amountStates) || all(x$amountStates %in% x$states)

  # is.matrix() first: dim() on a non-matrix yields logical(0) under `&&`.
  inherits(x, "eqnlist") &&
    all(required %in% names(x)) &&
    is.matrix(x$smatrix) &&
    identical(as.character(colnames(x$smatrix)), as.character(x$states)) &&
    nrow(x$smatrix) == length(x$rates) &&
    ncol(x$smatrix) == length(x$states) &&
    refs && volumes && reactions && amounts
}


## Class "eqnlist" and its methods ------------------------------------------

#' Determine Conserved Quantities from the Stoichiometric Matrix
#'
#' Each vector of the null space basis of `S` ([nullZ()]) gives one linear
#' combination of the states that the reactions leave unchanged.
#'
#' @param S numeric matrix, one row per reaction and one named column per
#'   state. `NA` entries count as 0.
#' @param weight One of `"none"` (default) or `"volume"`. When `"volume"`, the
#'   columns of `S` are multiplied by their compartment volume before the kernel
#'   is computed, so the returned quantities are conserved in *amount* rather
#'   than concentration. Requires `volumes` to be supplied as numeric.
#' @param volumes Optional named numeric vector of volume values keyed by state,
#'   aligned with `colnames(S)`. Only consulted when `weight = "volume"`.
#' @return `NULL` if nothing is conserved, otherwise a data frame with one
#'   character column of expressions such as `"1*A+1*pA"`, one row per
#'   conserved quantity, and their number in the attribute `n`.
#'
#' @example inst/examples/equations.R
#' @export
conservedQuantities <- function(S, weight = c("none", "volume"), volumes = NULL) {
  weight <- match.arg(weight)
  S[is.na(S)] <- 0
  if (weight == "volume") S <- .volumeWeighted(S, volumes)

  basis <- nullZ(S)
  n <- ncol(basis)
  if (n == 0L) return(NULL)
  states <- colnames(S)
  if (is.null(states)) stop("Columns of stoichiometric matrix not named.")

  expr <- vapply(seq_len(n), function(k) {
    w <- basis[, k]
    used <- w != 0
    paste(w[used], states[used], sep = "*", collapse = "+")
  }, character(1))
  # gsub, not sub: a basis vector may have several negative weights.
  out <- data.frame(gsub("+-", "-", expr, fixed = TRUE))
  names(out) <- paste("Conserved quantities:", n)
  attr(out, "n") <- n
  out
}

# Scale the columns of S by the numeric volume of their state.
.volumeWeighted <- function(S, volumes) {
  if (is.null(volumes))
    stop("`weight = \"volume\"` requires a named numeric `volumes` argument.")
  if (is.null(colnames(S)))
    stop("`S` must have column names when `weight = \"volume\"`.")
  missing_vol <- setdiff(colnames(S), names(volumes))
  if (length(missing_vol) > 0L)
    stop("`volumes` missing entries for: ", paste(missing_vol, collapse = ", "))
  v_num <- suppressWarnings(as.numeric(volumes[colnames(S)]))
  if (anyNA(v_num))
    stop("`weight = \"volume\"` requires all volumes to be numeric; got symbolic expression(s).")
  sweep(S, 2, v_num, "*")
}


## ---- Totals (CQ basis) management on eqnlist ------------------------------

## Append `_2`, `_3`, ... until `base` is unique against `taken`.
#' @keywords internal
.uniqueName <- function(base, taken) {
  if (!(base %in% taken)) return(base)
  i <- 2L
  while (paste0(base, "_", i) %in% taken) i <- i + 1L
  paste0(base, "_", i)
}

## Longest common substring across all strings (empty if none).
#' @keywords internal
.lcs <- function(strings) {
  if (!length(strings)) return("")
  if (length(strings) == 1L) return(strings[1L])
  lcs2 <- function(a, b) {
    na <- nchar(a); nb <- nchar(b)
    if (na == 0L || nb == 0L) return("")
    best <- ""
    for (i in seq_len(na)) for (j in seq_len(nb)) {
      k <- 0L
      while (i + k <= na && j + k <= nb &&
             substr(a, i + k, i + k) == substr(b, j + k, j + k)) k <- k + 1L
      if (k > nchar(best)) best <- substr(a, i, i + k - 1L)
    }
    best
  }
  Reduce(lcs2, strings)
}

## `totalXxx` from LCS of CQ species (require >= 2 chars after trimming
## underscores), else `total_<index>`. Disambiguated against `used`/`parameters`.
#' @keywords internal
.smartTotalName <- function(species, used, parameters, index) {
  fallback <- .uniqueName(paste0("total_", index), c(used, parameters))
  if (length(species) < 2L) return(fallback)
  lcs <- gsub("^_+|_+$", "", .lcs(species))
  if (nchar(lcs) < 2L) return(fallback)
  .uniqueName(paste0("total", lcs), c(used, parameters))
}

## Linear coefficients of `expr` over `states` by 1-hot evaluation, checked
## for linearity at the all-twos point (value must equal 2 * sum(coefs)).
## NULL on unknown symbols or nonlinearity.
#' @keywords internal
.linearCoefs <- function(expr, states, tol = 1e-8) {
  syms <- getSymbols(expr)
  if (length(setdiff(syms, states))) return(NULL)
  parsed <- parse(text = expr)
  e0 <- setNames(as.list(rep(0, length(states))), states)
  base <- tryCatch(eval(parsed, envir = e0), error = function(e) NA_real_)
  if (!is.finite(base)) return(NULL)
  v <- setNames(numeric(length(states)), states)
  for (s in intersect(syms, states)) {
    e1 <- e0; e1[[s]] <- 1
    v[s] <- eval(parsed, envir = e1) - base
  }
  e2 <- setNames(as.list(rep(2, length(states))), states)
  check <- tryCatch(eval(parsed, envir = e2), error = function(e) NA_real_)
  if (!is.finite(check) || abs(check - base - 2 * sum(v)) > tol)
    return(NULL)
  v
}

#' Validate Conservation Expressions Against an Eqnlist
#'
#' Each expression must be linear in `eqnlist$states` and a true conservation
#' quantity (its coefficient vector lies in the left null space of the
#' stoichiometric matrix). The full basis must have rank equal to
#' `nrow(conservedQuantities(eqnlist$smatrix))`.
#'
#' @param totals Named list of conservation expressions.
#' @param eqnlist The [eqnlist] to validate against.
#' @param tol Numerical tolerance for the null-space check.
#' @return `TRUE` on success; otherwise `stop()`s with a descriptive error.
#' @keywords internal
.validateTotals <- function(totals, eqnlist, tol = 1e-8) {
  if (!is.list(totals) || is.null(names(totals)) || any(!nzchar(names(totals))))
    stop("`totals` must be a fully named list of conservation expressions.",
         call. = FALSE)
  states <- eqnlist$states
  S <- eqnlist$smatrix; S[is.na(S)] <- 0

  V <- matrix(0, length(totals), length(states),
              dimnames = list(names(totals), states))
  for (k in seq_along(totals)) {
    v <- .linearCoefs(totals[[k]], states)
    if (is.null(v))
      stop("totals[['", names(totals)[k], "']] is not linear in the model states.",
           call. = FALSE)
    if (max(abs(S %*% v)) > tol)
      stop("totals[['", names(totals)[k], "']] is not a conservation quantity ",
           "(S %*% v is not zero).", call. = FALSE)
    V[k, ] <- v
  }

  auto <- conservedQuantities(S)
  n_required <- if (is.null(auto)) 0L else nrow(auto)
  if (length(totals) != n_required)
    stop("Expected ", n_required, " conservation quantities, got ",
         length(totals), ".", call. = FALSE)
  if (n_required > 0L && qr(V)$rank < n_required)
    stop("Supplied `totals` are linearly dependent; rank ", qr(V)$rank,
         " < ", n_required, ".", call. = FALSE)
  TRUE
}

#' Conserved Quantities of an Equation List
#'
#' Returns the user-defined totals set by [customTotals()], otherwise the
#' basis found by [conservedQuantities()]. Totals without a user name are
#' named after their species.
#'
#' @param eqnlist object of class [eqnlist].
#' @return Named list of conservation expressions (character). Empty list when
#'   the stoichiometric matrix admits no conserved quantities.
#' @seealso [customTotals()], [conservedQuantities()]
#' @export
getTotals <- function(eqnlist) {
  if (!is.null(eqnlist$totals)) return(eqnlist$totals)
  S <- eqnlist$smatrix
  if (is.null(S) || ncol(S) == 0L) return(list())
  cq <- conservedQuantities(S)
  if (is.null(cq) || nrow(cq) == 0L) return(list())
  out <- list()
  for (i in seq_len(nrow(cq))) {
    expr <- as.character(cq[i, 1])
    sp   <- intersect(getSymbols(expr), eqnlist$states)
    nm   <- .smartTotalName(sp, names(out), eqnlist$states, i)
    out[[nm]] <- expr
  }
  out
}

#' Set or Reset User-Defined Conserved Quantities
#'
#' Attaches a named list of conservation expressions to the eqnlist. They
#' replace the detected basis, also as parameters of [Pimpl()]. Each
#' expression must be linear in the states and conserved by the reactions, and
#' together they must span as many quantities as [conservedQuantities()] finds.
#'
#' @param eqnlist object of class [eqnlist].
#' @param totals named list of expressions, or `NULL` or `list()` to return to
#'   the detected basis.
#' @return The eqnlist with `$totals` updated.
#' @seealso [getTotals()]
#' @export
customTotals <- function(eqnlist, totals) {
  if (is.null(totals) || (is.list(totals) && length(totals) == 0L)) {
    eqnlist$totals <- NULL
    return(eqnlist)
  }
  .validateTotals(totals, eqnlist)
  attr(totals, "custom") <- TRUE
  eqnlist$totals <- totals
  eqnlist
}


#' Generate a Table of Reactions from an Equation List
#' 
#' @param eqnlist object of class [eqnlist]
#' @return `data.frame` with educts, products, rate and description. The first
#' column is a check if the reactions comply with reaction kinetics.
#' 
#' @example inst/examples/equations.R
#' @export
getReactions <- function(eqnlist) {

  S <- eqnlist$smatrix
  if (is.null(S)) return()
  rates <- eqnlist$rates

  side <- function(v) paste(ifelse(v == 1, names(v), paste0(v, "*", names(v))), collapse = " + ")
  reactions <- apply(S, 1, function(v) {
    v <- v[!is.na(v)]
    c(side(-v[v < 0]), side(v[v > 0]))
  })

  # "!" marks a reaction whose rate does not contain all its educts.
  check <- vapply(seq_along(rates), function(i) {
    educts <- eqnlist$states[which(S[i, ] < 0)]
    if (all(educts %in% getSymbols(rates[i]))) " " else "!"
  }, "")

  out <- data.frame(check, reactions[1, ], "->", reactions[2, ], rates, eqnlist$description,
                    stringsAsFactors = FALSE)
  colnames(out) <- c("Check", "Educt", "->", "Product", "Rate", "Description")
  rownames(out) <- seq_len(nrow(out))
  out
}


#' Add Reaction to an Equation List
#'
#' @param eqnlist equation list, see [eqnlist]
#' @param from character with the left hand side of the reaction, e.g. "2*A + B"
#' @param to character with the right hand side of the reaction, e.g. "C + 2*D"
#' @param rate character. The rate associated with the reaction. The name is employed as a description
#' of the reaction.
#' @param description Optional description instead of `names(rate)`.
#' @param compartment character, compartment ID for the states this reaction
#' introduces. Default `"defaultComp"`, created with volume `"1"` if new.
#' States that already have a compartment keep it; see [assignCompartment()].
#' @param rateCompartment compartment ID of the frame in which `rate` is a
#' concentration rate. Default `NA` lets [getFluxes()] infer the frame from the
#' educts. If the educts span several compartments and `compartment` is given
#' explicitly, `compartment` is the frame.
#' @return An object of class [eqnlist].
#' @examples
#' f <- eqnlist()
#' f <- addReaction(f, "2*A+B", "C + 2*D", "k1*B*A^2")
#' f <- addReaction(f, "C + A", "B + A", "k2*C*A")
#'
#' @example inst/examples/equations.R
#' @export
#' @rdname addReaction
addReaction <- function(eqnlist, from, to, rate, description = names(rate),
                         compartment = "defaultComp", rateCompartment = NA_character_) {


  if (missing(eqnlist)) eqnlist <- eqnlist()

  volumes <- eqnlist$volumes
  compartments_in <- eqnlist$compartments
  compartmentOf_in <- eqnlist$compartmentOf
  reactionCompartment_in <- eqnlist$reactionCompartment

  # Analyze the reaction character expressions
  eductCoef <- .sideCoefficients(from)
  educts <- names(eductCoef)
  productCoef <- .sideCoefficients(to)
  products <- names(productCoef)
  if (!length(educts) && !length(products))
    stop("A reaction needs at least one species.", call. = FALSE)
  reserved <- intersect(c(educts, products), c("Rate", "Description"))
  if (length(reserved))
    stop("'", reserved[1], "' is reserved and cannot name a state.", call. = FALSE)


  # States introduced by this reaction
  states <- unique(c(educts, products))

  # Description
  if(is.null(description)) description <- ""

  # Stoichiometric matrix
  smatrix <- matrix(NA, nrow = 1, ncol=length(states)); colnames(smatrix) <- states
  if(length(educts)>0) smatrix[,educts] <- -eductCoef
  if(length(products)>0) {
    filled <- !is.na(smatrix[,products])
    smatrix[,products[filled]] <- smatrix[,products[filled]] + productCoef[filled]
    smatrix[,products[!filled]] <- productCoef[!filled]
  }


  smatrix[smatrix == "0"] <- NA


  # data.frame
  mydata <- cbind(data.frame(Description = description, Rate = as.character(rate)), as.data.frame(smatrix))
  row.names(mydata) <- NULL


  if(!is.null(eqnlist)) {
    mydata0 <- as.data.frame(eqnlist)
    mydata <- combine(mydata0, mydata)
  }

  # Extend compartment assignment for brand-new states with the `compartment` arg.
  new_states <- setdiff(states, names(compartmentOf_in))
  compartments_out <- compartments_in
  compartmentOf_out <- compartmentOf_in
  if (length(new_states) > 0L) {
    if (is.null(compartments_out)) compartments_out <- list()
    if (is.null(compartmentOf_out)) compartmentOf_out <- character(0)
    if (!compartment %in% names(compartments_out)) {
      compartments_out[[compartment]] <- list(volume = "1", rule = NULL)
    }
    compartmentOf_out <- c(compartmentOf_out,
                           setNames(rep(compartment, length(new_states)), new_states))
  }

  # Educts in different compartments need a frame; an explicit `compartment` is one.
  if (is.na(rateCompartment) && !missing(compartment)) {
    educt_comps <- unique(unname(compartmentOf_out[educts]))
    if (length(educt_comps) > 1L) rateCompartment <- compartment
  }

  # Extend reactionCompartment with the value for this new reaction. When the
  # input list has no annotations (NULL), pad with NA for the existing rates
  # so the final vector lines up with the combined data.frame rows.
  existing_n <- length(eqnlist$rates)
  if (is.null(reactionCompartment_in)) reactionCompartment_in <- rep(NA_character_, existing_n)
  reactionCompartment_out <- c(reactionCompartment_in, as.character(rateCompartment))
  if (all(is.na(reactionCompartment_out))) reactionCompartment_out <- NULL

  new_el <- as.eqnlist(mydata, volumes = volumes,
                       compartments = compartments_out,
                       compartmentOf = compartmentOf_out,
                       reactionCompartment = reactionCompartment_out,
                       amountStates = eqnlist$amountStates)

  ## Preserve user-customized totals if they survive the structural change.
  old_totals <- eqnlist$totals
  if (!is.null(old_totals) && isTRUE(attr(old_totals, "custom"))) {
    new_el <- tryCatch(customTotals(new_el, unclass(old_totals)),
                       error = function(e) {
                         warning("customTotals invalidated by addReaction(): ",
                                 conditionMessage(e), ". Resetting to auto.",
                                 call. = FALSE)
                         new_el
                       })
  }
  new_el
}


# Rebuild an eqnlist around a new layout, so `$volumes` is re-derived.
.withCompartments <- function(eqnlist, compartments, compartmentOf) {

  if (is.null(eqnlist$smatrix))
    return(eqnlist(compartments = compartments, compartmentOf = compartmentOf))

  out <- eqnlist(smatrix = eqnlist$smatrix, states = eqnlist$states,
                 rates = eqnlist$rates, description = eqnlist$description,
                 compartments = compartments, compartmentOf = compartmentOf,
                 reactionCompartment = eqnlist$reactionCompartment,
                 amountStates = eqnlist$amountStates)

  old_totals <- eqnlist$totals
  if (!is.null(old_totals) && isTRUE(attr(old_totals, "custom")))
    out <- tryCatch(customTotals(out, unclass(old_totals)), error = function(e) out)

  out
}


#' Assign States to Compartments
#'
#' @description Declares the compartment of a state, before or after the state
#' is introduced by [addReaction()].
#'
#' @param eqnlist object of class [eqnlist]
#' @param ... named arguments `state = "compartment"`, or a single named
#' character vector of the same shape. Compartments that do not exist yet are
#' created with volume `"1"`.
#' @param volume Optional volume expression for the target compartment. Only
#' allowed when `...` names a single compartment; use [setCompartmentVolume()]
#' for several.
#' @param rule optional volume rule (`dV/dt`) for the target compartment, see
#' [setCompartmentVolume()]. Default `NULL`.
#' @return An object of class [eqnlist].
#' @seealso [setCompartmentVolume()], [addReaction()]
#' @examples
#' f <- eqnlist() |>
#'   assignCompartment(TGFb = "extraCell", volume = "V_ext") |>
#'   addReaction("R + TGFb", "R_TGFb", "k_on*R*TGFb", compartment = "Cell")
#' f$compartmentOf
#' @export
assignCompartment <- function(eqnlist, ..., volume = NULL, rule = NULL) {

  assignment <- unlist(list(...))
  if (length(assignment) == 0L) return(eqnlist)
  if (is.null(names(assignment)) || any(!nzchar(names(assignment))))
    stop("`...` must be named: assignCompartment(eqnlist, state = \"compartment\").")
  assignment <- setNames(as.character(assignment), names(assignment))

  compartments <- eqnlist$compartments
  if (is.null(compartments)) compartments <- list()
  targets <- unique(unname(assignment))
  for (cid in setdiff(targets, names(compartments)))
    compartments[[cid]] <- list(volume = "1", rule = NULL)

  if (!is.null(volume) || !is.null(rule)) {
    if (length(targets) != 1L)
      stop("`volume`/`rule` apply to a single compartment; `...` names ",
           paste(targets, collapse = ", "), ". Use setCompartmentVolume().")
    entry <- compartments[[targets]]
    if (!is.null(volume)) entry$volume <- as.character(volume)
    if (!is.null(rule))   entry$rule   <- as.character(rule)
    compartments[[targets]] <- entry
  }

  compartmentOf <- eqnlist$compartmentOf
  if (is.null(compartmentOf)) compartmentOf <- character(0)
  compartmentOf[names(assignment)] <- unname(assignment)

  .withCompartments(eqnlist, compartments, compartmentOf)

}


#' Set Compartment Volumes
#'
#' @description Sets the volume expression, and optionally the volume rule, of
#' one or more compartments and updates the per-state volumes.
#'
#' @param eqnlist object of class [eqnlist]
#' @param ... named arguments `compartment = "volume expression"`, or a single
#' named character vector of the same shape. Compartments that do not exist yet
#' are created.
#' @param rules Optional named character vector `compartment = "dV/dt"`, adding
#' a volume rule. `NA` removes an existing rule.
#' @return An object of class [eqnlist].
#' @seealso [assignCompartment()]
#' @examples
#' f <- eqnlist() |>
#'   addReaction("A", "B", "k*A", compartment = "cyt") |>
#'   setCompartmentVolume(cyt = "V_cyt")
#' f$volumes
#' @export
setCompartmentVolume <- function(eqnlist, ..., rules = NULL) {

  vols <- unlist(list(...))
  if (length(vols) > 0L && (is.null(names(vols)) || any(!nzchar(names(vols)))))
    stop("`...` must be named: setCompartmentVolume(eqnlist, compartment = \"V\").")
  if (length(vols) == 0L && is.null(rules)) return(eqnlist)

  compartments <- eqnlist$compartments
  if (is.null(compartments)) compartments <- list()

  entry_of <- function(cid) {
    if (cid %in% names(compartments)) compartments[[cid]] else list(volume = "1", rule = NULL)
  }
  for (cid in names(vols)) {
    entry <- entry_of(cid)
    entry$volume <- as.character(vols[[cid]])
    compartments[[cid]] <- entry
  }
  for (cid in names(rules)) {
    entry <- entry_of(cid)
    entry$rule <- if (is.na(rules[[cid]])) NULL else as.character(rules[[cid]])
    compartments[[cid]] <- entry
  }

  compartmentOf <- eqnlist$compartmentOf
  if (is.null(compartmentOf)) compartmentOf <- character(0)

  .withCompartments(eqnlist, compartments, compartmentOf)

}


# Reference compartment per reaction: the frame in which its rate is a
# concentration-rate. Shared with .volumeScaledReactions() so the ODE and the
# csv the steady-state backend reads cannot drift apart.
.refCompartments <- function(SMatrix, compOf, reactionCompartment, description) {

  # A model without reactions has no stoichiometry to read frames off.
  if (is.null(SMatrix)) return(character(0))

  vref_cid <- rep(NA_character_, nrow(SMatrix))
  for (i in seq_len(nrow(SMatrix))) {
    if (!is.null(reactionCompartment) && !is.na(reactionCompartment[i])) {
      vref_cid[i] <- reactionCompartment[i]
      next
    }
    row_i <- SMatrix[i, ]
    educt_idx <- which(!is.na(row_i) & row_i < 0)
    product_idx <- which(!is.na(row_i) & row_i > 0)
    cand <- if (length(educt_idx) > 0) unique(compOf[educt_idx])
            else if (length(product_idx) > 0) unique(compOf[product_idx])
            else character(0)
    if (length(cand) == 1L) {
      vref_cid[i] <- cand
    } else if (length(cand) > 1L) {
      stop(sprintf(
        "Reaction %d (\"%s\") spans compartments (%s). Pass `reactionCompartment` to name the frame in which the rate is a concentration-rate.",
        i, description[i], paste(cand, collapse = ", ")))
    } else {
      # Net-zero reaction ("A -> A"): no flux to scale, so no frame to name.
      vref_cid[i] <- NA_character_
    }
  }
  vref_cid
}


# Volumes of the reference compartments returned by .refCompartments(), NA-safe.
.refVolumes <- function(vref_cid, compartments) {
  vapply(vref_cid, function(cid) if (is.na(cid)) NA_character_ else compartments[[cid]]$volume,
         character(1))
}


# Reactions for the steady-state backend, with getFluxes()' V_ref / V_X scaling
# split into one row per distinct ratio of a reaction. Attribute "volumes": V_X
# per state, "1" for amount states; rows of one reaction share rate * V_X.
.volumeScaledReactions <- function(eqnlist) {

  data <- as.data.frame(eqnlist)
  SMatrix <- eqnlist$smatrix
  variables <- eqnlist$states
  compartments <- eqnlist$compartments
  compartmentOf <- eqnlist$compartmentOf
  if (is.null(SMatrix) || is.null(compartments) || is.null(compartmentOf))
    return(data)

  # A volume rule adds a dilution term -[X]*V'/V, which is not a reaction flux.
  ruled <- names(compartments)[vapply(compartments,
    function(cmp) !is.null(cmp$rule) && nzchar(cmp$rule), logical(1))]
  if (length(ruled))
    stop("Compartment(s) ", paste(ruled, collapse = ", "), " have a volume ",
         "rule, whose dilution term -[X]*V'/V is not a reaction flux.",
         call. = FALSE)

  compOf <- compartmentOf[variables]
  vol <- vapply(compOf, function(cid) compartments[[cid]]$volume, character(1))
  isAmount <- variables %in% eqnlist$amountStates
  vol[isAmount] <- "1"
  vref_cid <- .refCompartments(SMatrix, compOf, eqnlist$reactionCompartment,
                               eqnlist$description)
  vref_vol <- .refVolumes(vref_cid, compartments)

  # as.data.frame.eqnlist(): Description, Rate, then states in order.
  stateCols <- seq_along(variables) + 2L

  rows <- lapply(seq_len(nrow(SMatrix)), function(i) {
    touched <- which(!is.na(SMatrix[i, ]))
    ratio <- ifelse(!is.na(vref_cid[i]) & vref_cid[i] == compOf[touched] & !isAmount[touched], "",
                    paste0("*(", vref_vol[i], "/", vol[touched], ")"))
    lapply(unique(ratio), function(rr) {
      keep <- touched[ratio == rr]
      row <- data[i, , drop = FALSE]
      if (nzchar(rr)) row$Rate <- paste0("(", data$Rate[i], ")", rr)
      row[1, stateCols] <- NA
      row[1, stateCols[keep]] <- as.numeric(SMatrix[i, keep])
      row
    })
  })

  out <- do.call(rbind, unlist(rows, recursive = FALSE))
  rownames(out) <- NULL
  attr(out, "volumes") <- vol
  out
}


#' Generate List of Fluxes from Equation List
#'
#' @param eqnlist object of class [eqnlist].
#' @param type `"conc"` (default) or `"amount"` for fluxes in units of
#' concentration or of amount.
#' @return Named list with one named character vector per state, its in- and
#' out-fluxes.
#' @example inst/examples/equations.R
#' @export
getFluxes <- function(eqnlist, type = c("conc", "amount")) {

  type <- match.arg(type)

  description <- eqnlist$description
  rate <- eqnlist$rates
  variables <- eqnlist$states
  SMatrix <- eqnlist$smatrix
  compartments <- eqnlist$compartments
  compartmentOf <- eqnlist$compartmentOf
  reactionCompartment <- eqnlist$reactionCompartment

  if (is.null(SMatrix)) return()

  # Fallback when compartment info was not populated: treat every state as
  # living in an implicit "defaultComp" compartment with volume "1".
  if (is.null(compartments) || is.null(compartmentOf)) {
    compOf <- setNames(rep("defaultComp", length(variables)), variables)
    compartments <- list(defaultComp = list(volume = "1", rule = NULL))
  } else {
    compOf <- compartmentOf[variables]
  }
  volumes <- vapply(compOf, function(cid) compartments[[cid]]$volume, character(1))
  names(volumes) <- variables
  # Amount states are not divided by a volume: their V_X is 1.
  isAmount <- setNames(variables %in% eqnlist$amountStates, variables)
  volumes[isAmount] <- "1"

  # Reference compartment V_ref per reaction: reactionCompartment[i] if non-NA,
  # else the unique educt compartment, else the unique product compartment for
  # pure synthesis. Educts in several compartments without annotation error.
  vref_cid <- .refCompartments(SMatrix, compOf, reactionCompartment, description)
  vref_vol <- .refVolumes(vref_cid, compartments)

  # generate equation expressions
  terme <- lapply(seq_along(variables), function(j) {
    v <- SMatrix[,j]
    nonZeros <- which(!is.na(v))
    var.description <- description[nonZeros]
    positives <- which(v > 0)
    destin_cid <- compOf[[j]]
    destin_vol <- volumes[[j]]
    destin_amount <- isAmount[[j]]

    # Uniform flux formula: flux_X = stoich_X * rate * V_ref / V_X for every
    # species in every reaction.
    switch(type,
           conc = {
             volumes.ratios <- if (destin_amount) paste0("*(", vref_vol, ")")
                               else paste0("*(", vref_vol, "/", destin_vol, ")")
             volumes.ratios[is.na(vref_cid) | (!destin_amount & vref_cid == destin_cid)] <- ""
           },
           amount = {
             volumes.ratios <- paste0("*(", vref_vol, ")")
             volumes.ratios[is.na(vref_cid)] <- ""
           }
    )

    numberchar <- as.character(v)
    if (nonZeros[1] %in% positives) {
      numberchar[positives] <- paste(c("", rep("+", length(positives)-1)), numberchar[positives], sep = "")
    } else {
      numberchar[positives] <- paste("+", numberchar[positives], sep = "")
    }
    var.flux <- if (length(nonZeros) == 0L) character(0)
                else paste0(numberchar[nonZeros], "*(", rate[nonZeros], ")", volumes.ratios[nonZeros])
    names(var.flux) <- var.description

    # A volume rule adds -[X]*(dV/dt)/V to d[X]/dt. Amounts have no such term.
    r <- if (type == "conc" && !destin_amount) compartments[[destin_cid]]$rule else NULL
    if (!is.null(r) && nzchar(r)) {
      dilution <- paste0("-(", variables[j], ")*(", r, ")/(", destin_vol, ")")
      names(dilution) <- paste0("dilution_", destin_cid)
      var.flux <- c(var.flux, dilution)
    }

    return(var.flux)
  })

  fluxes <- terme
  names(fluxes) <- variables

  return(fluxes)


}






#' Coerce Equation List into a Data Frame
#' 
#' @param x object of class [eqnlist].
#' @param ... not used.
#' @return A `data.frame` with the columns "Description" (character),
#' "Rate" (character) and one column per state, the columns of the
#' stoichiometric matrix. The compartment layout is kept as attributes and
#' read back by [as.eqnlist()]. `NULL` for an empty equation list.
#' @export
as.data.frame.eqnlist <- function(x, ...) {

  eqnlist <- x

  if(is.null(eqnlist$smatrix)) return()

  data <- data.frame(Description = eqnlist$description,
                     Rate = eqnlist$rate,
                     eqnlist$smatrix,
                     stringsAsFactors = FALSE)

  attr(data, "volumes") <- eqnlist$volumes
  attr(data, "compartments") <- eqnlist$compartments
  attr(data, "compartmentOf") <- eqnlist$compartmentOf
  attr(data, "reactionCompartment") <- eqnlist$reactionCompartment
  attr(data, "amountStates") <- eqnlist$amountStates

  return(data)
}

#' Write Equation List into a CSV File
#'
#' Compartments and reaction frames are not written; a warning is issued if
#' the equation list has any.
#'
#' @param eqnlist object of class [eqnlist].
#' @param ... arguments going to [utils::write.csv()]. `row.names` defaults to
#' `FALSE` and `na` to `""`.
#' @return `NULL`, invisibly. Called for the file it writes.
#' @seealso [as.data.frame.eqnlist()]
#'
#' @export
#' @importFrom utils file.edit getParseData install.packages installed.packages read.csv str tail write.csv
write.eqnlist <- function(eqnlist, ...) {
  
  
  arglist <- list(...)
  argnames <- names(arglist)
  if (!"row.names" %in% argnames) arglist$row.names <- FALSE
  if (!"na" %in% argnames) arglist$na <- ""
  
  if (!length(eqnlist$rates)) stop("The eqnlist has no reaction to write.", call. = FALSE)
  arglist$x <- as.data.frame(eqnlist)

  # The csv has no place for the layout; don't drop it silently.
  if (!.trivialCompartments(eqnlist$compartments) || !is.null(eqnlist$reactionCompartment))
    warning("write.eqnlist(): the csv format has no place for compartments; ",
            "volumes and reaction frames are lost on read-back.", call. = FALSE)

  do.call(write.csv, arglist)
  
}


#' Subset of an Equation List
#' 
#' @param x object of class [eqnlist].
#' @param ... logical expression for subsetting.
#' @details The expression in `...` can refer to `Educt`, `Product`, `Rate`
#' and `Description`. Within it, `%in%` searches the educts and products of
#' each reaction (see examples).
#'
#' @return Object of class [eqnlist].
#' @seealso [addReaction()], [as.eqnlist()]
#' @examples
#' reactions <- data.frame(Description = c("Activation", "Deactivation"), 
#'                         Rate = c("act*A", "deact*pA"), A=c(-1,1), pA=c(1, -1) )
#' f <- as.eqnlist(reactions)
#' subset(f, "A" %in% Educt)
#' subset(f, "pA" %in% Product)
#' subset(f, grepl("act", Rate))
#' @export subset.eqnlist
#' @export
subset.eqnlist <- function(x, ...) {
  
  eqnlist <- x
  callerEnv <- parent.frame()
  
  # Do selection on data.frame
  data <- getReactions(eqnlist)
  if (is.null(data)) return(eqnlist())
  
  data.list <- list(Educt = lapply(data$Educt, getSymbols), 
                    Product = lapply(data$Product, getSymbols),
                    Rate = data$Rate,
                    Description = data$Description,
                    Check = data$Check)
  
  # The condition resolves against the reaction table, then the calling frame.
  # `%in%` is overloaded for Educt and Product, lists of symbol vectors, and is
  # the base operator otherwise.
  env <- list2env(data.list, envir = new.env(parent = callerEnv))
  env$"%in%" <- function(x, table) {
    if (!is.list(table)) return(base::`%in%`(x, table))
    vapply(table, function(el) any(base::`%in%`(x, el)), logical(1))
  }
  condition <- eval(substitute(alist(...)))[[1L]]
  select <- which(eval(condition, env))
  if (length(select) == 0) return(eqnlist())
  
  # Translate subsetting on eqnlist entries
  # smatrix
  smatrix <- submatrix(eqnlist$smatrix, rows = select)
  empty <- sapply(1:ncol(smatrix), function(i) all(is.na(smatrix[, i])))
  smatrix <- submatrix(smatrix, cols = !empty)
  
  # states and rates
  states <- colnames(smatrix)
  rates <- eqnlist$rates[select]

  # volumes (derived view; filter to surviving states)
  volumes <- eqnlist$volumes
  if(!is.null(volumes)) volumes <- volumes[intersect(names(volumes),  states)]

  # description
  description <- eqnlist$description[select]

  reactionCompartment <- if (!is.null(eqnlist$reactionCompartment)) eqnlist$reactionCompartment[select] else NULL
  if (!is.null(reactionCompartment) && all(is.na(reactionCompartment))) reactionCompartment <- NULL

  # Restrict to surviving states; a compartment named as a reaction frame counts
  # as referenced too. `%in%` is locally shadowed above, so use base::`%in%`.
  compartmentOf <- eqnlist$compartmentOf
  compartments <- eqnlist$compartments
  if (!is.null(compartmentOf)) {
    compartmentOf <- compartmentOf[intersect(names(compartmentOf), states)]
    if (!is.null(compartments)) {
      used_cids <- unique(c(compartmentOf, stats::na.omit(reactionCompartment)))
      compartments <- compartments[base::`%in%`(names(compartments), used_cids)]
    }
  }

  new_el <- eqnlist(smatrix, states, rates, volumes, description,
                    compartments = compartments, compartmentOf = compartmentOf,
                    reactionCompartment = reactionCompartment,
                    amountStates = intersect(eqnlist$amountStates, states))

  old_totals <- eqnlist$totals
  if (!is.null(old_totals) && isTRUE(attr(old_totals, "custom"))) {
    new_el <- tryCatch(customTotals(new_el, unclass(old_totals)),
                       error = function(e) {
                         warning("customTotals invalidated by subset(): ",
                                 conditionMessage(e), ". Resetting to auto.",
                                 call. = FALSE)
                         new_el
                       })
  }
  new_el
}


#' Print or Pander Equation List
#'
#' @param x object of class [eqnlist].
#' @param pander logical, use pander for output (R Markdown). Default `FALSE`.
#' @param ... not used.
#' @return `x`, invisibly.
#' 
#' @export
print.eqnlist <- function(x, pander = FALSE, ...) {

  reactions <- getReactions(x)

  if (pander) {
    if (!requireNamespace("pander", quietly = TRUE)) {
      print(reactions)
      return(invisible(x))
    }
    .panderTableOptions()
    reactions$Check <- NULL
    reactions$Rate <- unname(format.eqnvec(as.character(reactions$Rate)))
    pander::pander(reactions)
    return(invisible(x))
  }

  totals <- getTotals(x)
  if (length(totals)) {
    custom <- isTRUE(attr(x$totals, "custom"))
    cat("Conserved quantities", if (custom) " (custom)", ":\n", sep = "")
    cat(sprintf("  %s = %s\n", names(totals), unlist(totals, use.names = FALSE)), sep = "")
  }
  layout <- .format_compartments(x$compartments, x$compartmentOf, x$states, x$amountStates)
  if (length(layout)) cat("", layout, sep = "\n")
  cat("\n")
  print(reactions)
  invisible(x)
}


# Left-aligned pander tables that are never split.
.panderTableOptions <- function() {
  pander::panderOptions("table.alignment.default", "left")
  pander::panderOptions("table.split.table", Inf)
  pander::panderOptions("table.split.cells", Inf)
}


# Internal: is this the implicit single unit-volume compartment?
.trivialCompartments <- function(compartments) {
  if (is.null(compartments)) return(TRUE)
  if (length(compartments) != 1L) return(FALSE)
  only <- compartments[[1L]]
  identical(only$volume, "1") && is.null(only$rule)
}


# Internal: render a compact compartment summary for print.eqnlist.
# Returns character(0) when the model has exactly one default-volume
# compartment and nothing else to say, so prints stay short.
.format_compartments <- function(compartments, compartmentOf, states = names(compartmentOf),
                                 amountStates = NULL) {
  if (is.null(compartments) || is.null(compartmentOf)) return(character(0))
  declared <- setdiff(names(compartmentOf), states)
  if (.trivialCompartments(compartments) && length(declared) == 0L &&
      length(amountStates) == 0L)
    return(character(0))

  header <- "Compartments:"
  comp_entries <- vapply(names(compartments), function(cid) {
    entry <- compartments[[cid]]
    rule_txt <- if (!is.null(entry$rule) && nzchar(entry$rule)) paste0(", rule=", entry$rule) else ""
    sprintf("  %s (V=%s%s)", cid, entry$volume, rule_txt)
  }, character(1))

  by_comp <- split(states, unname(compartmentOf[states]))
  assign_lines <- sprintf("  %s: %s", names(by_comp),
                          vapply(by_comp, function(sts) paste(sts, collapse = ", "), character(1)))
  out <- c(header, comp_entries, "States by compartment:", assign_lines)

  if (length(declared) > 0L) {
    by_decl <- split(declared, unname(compartmentOf[declared]))
    out <- c(out, "Declared, not used by any reaction yet:",
             sprintf("  %s: %s", names(by_decl),
                     vapply(by_decl, function(sts) paste(sts, collapse = ", "), character(1))))
  }

  if (length(amountStates) > 0L)
    out <- c(out, paste("States in substance units:", paste(amountStates, collapse = ", ")))

  out
}



## Class "eqnvec" and its constructors --------------------------------------------



#' Coerce to an Equation Vector
#'
#' @param x object of class `character` or `eqnlist`.
#' @param ... for the `eqnlist` method, arguments going to [getFluxes()],
#' e.g. `type`; otherwise not used.
#' @details For an `eqnlist`, the right-hand sides of the ODE are assembled
#' from [getFluxes()].
#' @return Object of class [eqnvec].
#' @export
as.eqnvec <- function(x, ...) {
  UseMethod("as.eqnvec", x)
}

#' @param names character, the left-hand sides of the equations. Default
#' `NULL` uses `names(x)`.
#' @rdname as.eqnvec
#' @export
as.eqnvec.character <- function(x = NULL, names = NULL, ...) {
  
  equations <- x
  
  if (is.null(equations)) return(NULL)
  
  if (is.null(names)) names <- names(equations)
  if (is.null(names)) stop("equations need names")
  if (length(names) != length(equations)) stop("Length of names and equations do not coincide")
  if (any(is.na(names) | !nzchar(names))) stop("every equation needs a name")
  if (anyDuplicated(names)) stop("equation names are duplicated: ",
                                 paste(unique(names[duplicated(names)]), collapse = ", "))
  tryCatch(parse(text = equations), error = function(e)
    stop("equations cannot be parsed: ", conditionMessage(e), call. = FALSE))
  
  out <- structure(equations, names = names)
  class(out) <- c("eqnvec", "character")
  
  return(out)
  
}



#' @rdname as.eqnvec
#' @export
as.eqnvec.eqnlist <- function(x, ...) {
  
  eqnlist <- x
  
  terme <- getFluxes(eqnlist, ...)
  if(is.null(terme)) return()
  terme <- lapply(terme, function(t) if (length(t) == 0L) "0" else paste(t, collapse=" "))
  
  
  terme <- do.call(c, terme)
  
  as.eqnvec(terme, names(terme))
  
}

#' @export
c.eqnlist <- function(...) {

  inputs <- list(...)
  inputs <- inputs[!vapply(inputs, is.null, logical(1))]
  if (length(inputs) == 0L) return(eqnlist())
  hasReactions <- !vapply(inputs, function(x) is.null(x$smatrix), logical(1))

  # Merge compartments with conflict detection
  all_compartments <- list()
  for (el in inputs) {
    if (is.null(el$compartments)) next
    for (cid in names(el$compartments)) {
      new_entry <- el$compartments[[cid]]
      if (cid %in% names(all_compartments)) {
        old_entry <- all_compartments[[cid]]
        if (!identical(old_entry$volume, new_entry$volume)) {
          stop(sprintf("Compartment conflict: '%s' has volume '%s' in one eqnlist and '%s' in another.",
                       cid, old_entry$volume, new_entry$volume))
        }
        if (!identical(old_entry$rule, new_entry$rule)) {
          stop(sprintf("Compartment conflict: '%s' has different `rule` in the input eqnlists.", cid))
        }
      } else {
        all_compartments[[cid]] <- new_entry
      }
    }
  }

  all_compartmentOf <- character(0)
  for (el in inputs) {
    if (is.null(el$compartmentOf)) next
    for (st in names(el$compartmentOf)) {
      cid <- unname(el$compartmentOf[[st]])
      old <- if (st %in% names(all_compartmentOf)) unname(all_compartmentOf[[st]])
      # an explicit compartment wins over the default one
      if (is.null(old) || identical(old, "defaultComp")) {
        all_compartmentOf[st] <- cid
      } else if (!identical(old, cid) && !identical(cid, "defaultComp")) {
        stop(sprintf("State '%s' assigned to different compartments across input eqnlists.", st))
      }
    }
  }

  if (length(all_compartments) == 0L) all_compartments <- NULL
  if (length(all_compartmentOf) == 0L) all_compartmentOf <- NULL

  # Concatenate reactionCompartment annotations. If any input has them, the result needs
  # to produce a vector of length nrow(combined). Missing entries become NA.
  any_rc <- any(vapply(inputs, function(el) !is.null(el$reactionCompartment), logical(1)))
  if (any_rc) {
    all_rc <- unlist(lapply(inputs, function(el) {
      if (is.null(el$reactionCompartment)) rep(NA_character_, length(el$rates))
      else el$reactionCompartment
    }))
    if (all(is.na(all_rc))) all_rc <- NULL
  } else {
    all_rc <- NULL
  }

  out <- if (!any(hasReactions)) {
    eqnlist(compartments = all_compartments, compartmentOf = all_compartmentOf)
  } else {
    as.eqnlist(Reduce(combine, lapply(inputs[hasReactions], as.data.frame)),
               compartments = all_compartments, compartmentOf = all_compartmentOf,
               reactionCompartment = all_rc,
               amountStates = unique(unlist(lapply(inputs, function(el) el$amountStates))))
  }

  # Custom totals survive, completed by the detected totals of the other inputs,
  # as long as they form a basis of the result.
  custom <- vapply(inputs, function(el) isTRUE(attr(el$totals, "custom")), logical(1))
  if (!any(custom)) return(out)
  totals <- unlist(lapply(inputs[custom], function(el) unclass(el$totals)))
  auto <- unlist(lapply(inputs[!custom & hasReactions], function(el) unclass(getTotals(el))))
  if (length(auto))
    totals <- c(totals, setNames(auto, make.unique(c(names(totals), rep("total", length(auto))),
                                                   sep = "_")[-seq_along(totals)]))
  tryCatch(customTotals(out, as.list(totals)), error = function(e) {
    warning("customTotals invalidated by c(): ", conditionMessage(e),
            ". Resetting to auto.", call. = FALSE)
    out
  })

}


#' @export
#' @param x object of any class
#' @rdname eqnvec
is.eqnvec <- function(x) inherits(x, "eqnvec") && length(names(x)) == length(x)


## Class "eqnvec" and its methods --------------------------------------------





#' Format Equation Vector with Spaces
#' 
#' @param x object of class [eqnvec], or a named parsable character vector.
#' @param ... not used.
#' @return Named character vector with spaces around `+`, `-`, `*` and `/`.
#' @export format.eqnvec
#' @export
format.eqnvec <- function(x, ...) {
  
  # pads the operator tokens only, so a number like 1e-3 stays whole
  vapply(x, function(eqn) {
    tokens <- getParseData(parse(text = eqn, keep.source = TRUE))
    tokens <- tokens[tokens$terminal, ]
    text <- ifelse(tokens$token == "'^'", "^", tokens$text)
    op <- tokens$token %in% c("'+'", "'-'", "'*'", "'/'")
    text[op] <- paste0(" ", text[op], " ")
    paste(text, collapse = "")
  }, character(1))
}

#' Print Equation Vector
#'
#' @param x object of class [eqnvec].
#' @param width numeric, width of the print-out. Default 140.
#' @param pander logical, use pander for output (R Markdown). Default `FALSE`.
#' @param ... not used.
#' @return `x`, invisibly.
#' 
#' @import stringr
#' @export
print.eqnvec <- function(x, width = 140, pander = FALSE, ...) {

  if (pander) {
    tab <- data.frame(unname(format.eqnvec(unclass(x))), row.names = names(x))
    names(tab) <- ""
    if (!requireNamespace("pander", quietly = TRUE)) {
      print(tab)
      return(invisible(x))
    }
    .panderTableOptions()
    pander::pander(tab)
    return(invisible(x))
  }

  # Columns: index, right-aligned name, " <- ", right-hand side. Spaces are
  # dropped from the equations; continuation lines start under the first one.
  lhs <- names(x)
  ord <- order(lhs)
  nameWidth <- max(nchar(lhs), 5L)
  indent <- 3L + 1L + nameWidth + 4L
  rhsWidth <- max(width, nameWidth + 10) - indent
  rhs <- stringr::str_wrap(gsub(" ", "", unclass(x)[ord], fixed = TRUE),
                           width = rhsWidth, exdent = indent)
  idx <- stringr::str_pad(c("Idx", seq_along(ord)), width = 3L)
  key <- stringr::str_pad(c("Inner", lhs[ord]), width = nameWidth)
  cat(paste0(idx, " ", key, " <- ", c("Outer", rhs)), sep = "\n")
  invisible(x)
}



#' Summary of an Equation Vector
#'
#' Prints, for each equation, the symbols it depends on.
#'
#' @param object object of class [eqnvec].
#' @param ... not used.
#' @return `object`, invisibly.
#' 
#' @export
summary.eqnvec <- function(object, ...) {
  symbols <- vapply(object, function(eqn) paste(getSymbols(eqn), collapse = ", "),
                    character(1))
  cat(paste0(names(object), " = f( ", symbols, ")"), sep = "\n")
  invisible(object)
}


#' @export
c.eqnvec <- function(...) {
 
  joined <- unlist(lapply(list(...), unclass))
  if (anyDuplicated(names(joined))) stop("Names must be unique")
  as.eqnvec(joined)
}

#' @export
"[.eqnvec" <- function(x, ...) {
  out <- unclass(x)[...]
  class(out) <- c("eqnvec", "character")
  return(out)
}









## eqnvec / eqnlist constructors ----------------------------------------

## Equation classes -------------------------------------------------------

#' Generate Equation Vector Object
#'
#' @description The eqnvec object stores explicit algebraic equations, like the
#' right-hand sides of an ODE, observation functions or parameter transformations
#' as named character vectors.
#' @param ... mathematical expressions as characters to be coerced,
#' the right-hand sides of the equations
#' @return Object of class `eqnvec`, a named character vector.
#' @example inst/examples/eqnvec.R
#' @seealso [eqnlist]
#' @export
eqnvec <- function(...) {
  x <- c(...)
  if (!length(x)) return(NULL)
  as.eqnvec(.fillNames(vapply(x, as.character, ""), "eqn"))
}

#' @rdname eqnlist
#' @export
#' @param smatrix numeric stoichiometric matrix, one row per reaction and one
#'   column per state. Default `NULL` gives an empty equation list.
#' @param states character, the state names. Default `colnames(smatrix)`.
#' @param rates character, the rate expressions, one per reaction.
#' @param volumes named character of state volumes. Without `compartments`
#'   and `compartmentOf`, each distinct expression becomes a compartment with
#'   ID `c1`, `c2`, ...
#' @param description character, one description per reaction. Default
#'   `NULL` numbers the reactions.
#' @param compartments named list keyed by compartment ID; each entry is a
#'   volume expression (character) or a list with fields `volume` and `rule`
#'   (`dV/dt`, see [setCompartmentVolume()]).
#' @param compartmentOf named character vector mapping states to compartment IDs.
#'   States not listed default to compartment `"defaultComp"` with volume `"1"`.
#'   May also name states that do not exist yet; see [assignCompartment()].
#' @param reactionCompartment character of length `nrow(smatrix)`, the
#'   compartment in whose frame each rate is given; `NA` infers it from the
#'   educts. Required for reactions whose educts span several compartments.
#'   Default `NULL`.
#' @param amountStates character, the states given in amounts rather than
#'   concentrations (SBML's `hasOnlySubstanceUnits`). Their fluxes are not
#'   divided by a compartment volume. Default `NULL`.
#' @param totals named list of user-defined conservation expressions, see
#'   [customTotals()]. Default `NULL` uses the detected basis.
#' @return Object of class `eqnlist`, a list with the entries above and the
#'   per-state `volumes`. `is.eqnlist()` returns a logical.
#' @seealso [addReaction()], [as.eqnlist()], [assignCompartment()],
#'   [setCompartmentVolume()]
#' @example inst/examples/eqnlist.R
eqnlist <- function(smatrix = NULL, states = colnames(smatrix), rates = NULL,
                    volumes = NULL, description = NULL,
                    compartments = NULL, compartmentOf = NULL,
                    reactionCompartment = NULL, amountStates = NULL, totals = NULL) {

  # A model without species has no stoichiometry;
  # canonicalise it to NULL so is.eqnlist() recognises the empty list.
  if (length(states) == 0L && length(rates) == 0L) smatrix <- NULL

  # One column per state, one row per rate; "not involved" is NA, as in addReaction().
  if (!is.null(smatrix)) {
    smatrix <- as.matrix(smatrix)
    if (ncol(smatrix) != length(states))
      stop("Number of states does not coincide with number of columns of stoichiometric matrix")
    if (nrow(smatrix) != length(rates))
      stop("Number of rates does not coincide with number of rows of stoichiometric matrix")
    colnames(smatrix) <- states
    smatrix[!is.na(smatrix) & smatrix == 0] <- NA
    if (is.null(description)) description <- seq_len(nrow(smatrix))
  }

  resolved <- .resolve_compartments(as.character(states), compartments, compartmentOf, volumes)

  # Reaction-compartment annotation: optional per-reaction V_ref. NA means
  # "infer from educts" in getFluxes(). Default NULL (no annotation anywhere).
  if (!is.null(rates) && !is.null(reactionCompartment)) {
    reactionCompartment <- as.character(reactionCompartment)
    if (length(reactionCompartment) != length(rates))
      stop("`reactionCompartment` must match the number of rates (", length(rates), ").")
    bad <- setdiff(stats::na.omit(reactionCompartment), names(resolved$compartments))
    if (length(bad) > 0L)
      stop("`reactionCompartment` references undefined compartments: ", paste(bad, collapse = ", "))
  }

  if (!is.null(amountStates)) {
    amountStates <- as.character(amountStates)
    bad <- setdiff(amountStates, as.character(states))
    if (length(bad) > 0L)
      stop("`amountStates` names states that do not exist: ", paste(bad, collapse = ", "))
    if (length(amountStates) == 0L) amountStates <- NULL
  }

  out <- list(smatrix = smatrix,
              states = as.character(states),
              rates = as.character(rates),
              volumes = resolved$volumes,
              description = as.character(description),
              compartments = resolved$compartments,
              compartmentOf = resolved$compartmentOf,
              reactionCompartment = reactionCompartment,
              amountStates = amountStates,
              totals = totals)
  class(out) <- c("eqnlist", "list")

  return(out)
}


# The per-state `$volumes` view of the layout. Declared-but-absent states have none.
.derivedVolumes <- function(compartments, compartmentOf, states) {
  if (length(states) == 0L) return(NULL)
  setNames(vapply(compartmentOf[states], function(id) compartments[[id]]$volume, character(1)),
           states)
}


# Resolve compartment inputs into canonical (compartments, compartmentOf, volumes);
# all NULL when nothing is known at all. `compartmentOf` may name states that do
# not exist yet; those declarations are kept behind the present states.
.resolve_compartments <- function(states, compartments, compartmentOf, volumes) {

  if (length(states) == 0L && length(compartmentOf) == 0L && length(compartments) == 0L) {
    return(list(compartments = NULL, compartmentOf = NULL, volumes = NULL))
  }

  normalize_entry <- function(e) {
    if (is.character(e) && length(e) == 1L) return(list(volume = unname(e), rule = NULL))
    if (is.list(e) && !is.null(e$volume)) {
      return(list(volume = as.character(e$volume), rule = e$rule))
    }
    stop("Each compartment entry must be a character volume expression or a list with `$volume`.")
  }

  given_volumes <- volumes

  if (!is.null(compartmentOf) || (!is.null(compartments) && is.null(volumes))) {
    compartments <- lapply(compartments, normalize_entry)
    if (is.null(compartmentOf)) compartmentOf <- setNames(character(0), character(0))
    compartmentOf <- setNames(as.character(compartmentOf), names(compartmentOf))
    if (length(compartmentOf) && (is.null(names(compartmentOf)) || any(!nzchar(names(compartmentOf)))))
      stop("`compartmentOf` must be a fully named character vector (names = state IDs).")
    missing_states <- setdiff(states, names(compartmentOf))
    if (length(missing_states) > 0L) {
      if (!"defaultComp" %in% names(compartments))
        compartments[["defaultComp"]] <- list(volume = "1", rule = NULL)
      compartmentOf <- c(compartmentOf, setNames(rep("defaultComp", length(missing_states)), missing_states))
    }
    bad <- setdiff(compartmentOf, names(compartments))
    if (length(bad) > 0L)
      stop("compartmentOf references undefined compartments: ", paste(unique(bad), collapse = ", "))
  } else if (!is.null(volumes)) {
    v <- as.character(volumes)
    nms <- names(volumes)
    if (is.null(nms) || any(!nzchar(nms)))
      stop("`volumes` must be a fully named character vector.")
    names(v) <- nms
    # First appearance, not split(): generated IDs must not depend on the locale.
    compartments <- list()
    compartmentOf <- character(0)
    for (expr in unique(v)) {
      cid <- paste0("c", length(compartments) + 1L)
      compartments[[cid]] <- list(volume = expr, rule = NULL)
      compartmentOf[names(v)[v == expr]] <- cid
    }
    missing_states <- setdiff(states, names(compartmentOf))
    if (length(missing_states) > 0L) {
      compartments[["defaultComp"]] <- list(volume = "1", rule = NULL)
      compartmentOf <- c(compartmentOf, setNames(rep("defaultComp", length(missing_states)), missing_states))
    }
  } else {
    compartments <- list(defaultComp = list(volume = "1", rule = NULL))
    compartmentOf <- setNames(rep("defaultComp", length(states)), states)
  }

  # States first (in state order), declarations for absent states behind them.
  compartmentOf <- compartmentOf[c(states, setdiff(names(compartmentOf), states))]

  vols <- .derivedVolumes(compartments, compartmentOf, states)

  # The explicit layout wins over `volumes`; disagreement means a stale `$volumes`.
  if (!is.null(given_volumes) && !is.null(vols)) {
    common <- intersect(names(given_volumes), names(vols))
    clash <- common[as.character(given_volumes[common]) != vols[common]]
    if (length(clash) > 0L)
      warning("`volumes` disagrees with the compartment layout for ",
              paste(clash, collapse = ", "), "; the layout wins. ",
              "Use setCompartmentVolume() to change a compartment volume.",
              call. = FALSE)
  }

  list(compartments = compartments, compartmentOf = compartmentOf, volumes = vols)
}





#' Transform an ODE to log10 Coordinates
#'
#' Takes `f` as the right-hand side of \eqn{\dot x = f(x)} and returns the
#' same system in \eqn{x_{l10} = \log_{10} x},
#' \deqn{\dot x_{l10} = \frac{f(x)}{\log(10)\, x}, \qquad x = 10^{x_{l10}}.}
#' Read results back with `10^x_l10`.
#'
#' @param f `eqnvec` or named character vector, the right-hand sides.
#' @param suffix character appended to each state name. Default `"_l10"`.
#' @param simplify logical, cancel the division symbolically with the Python
#'   module `sympy`. Default `TRUE`. With `FALSE`, or without `sympy`, the
#'   plain quotient is returned.
#' @return An `eqnvec` over the renamed states.
#' @seealso [odemodel()]
#' @export
#' @examplesIf reticulate::py_module_available("sympy")
#' log10Transform(eqnvec(x = "-k*x"))
#' log10Transform(eqnvec(A = "-k1*A + k2*B", B = "k1*A - k2*B"))
log10Transform <- function(f, suffix = "_l10", simplify = TRUE) {

  f <- as.eqnvec(f)
  states <- names(f)
  renamed <- paste0(states, suffix)

  spy <- NULL
  if (simplify)
    spy <- tryCatch(reticulate::import("sympy", convert = TRUE),
                    error = function(e) NULL)

  rhs <- vapply(states, function(k) {
    quotient <- paste0("(", f[[k]], ")/(", k, ")")
    if (!is.null(spy)) {
      # cancel removes the state where it divides out, expand then distributes
      # the remainder over the division. Every name stays a plain symbol, also
      # one sympy knows as a constant such as E or I.
      syms <- getSymbols(quotient)
      cancelled <- tryCatch(
        gsub("\\*\\*", "^", as.character(
          spy$expand(spy$cancel(spy$sympify(gsub("\\^", "**", quotient),
                                            locals = setNames(lapply(syms, spy$Symbol), syms)))))),
        error = function(e) NULL)
      if (!is.null(cancelled)) quotient <- cancelled
    }
    quotient <- replaceSymbols(states, paste0("(10^(", renamed, "))"), quotient)
    # a bare symbol needs no parentheses; anything else does
    if (grepl("^-?[A-Za-z0-9_.]+$", quotient)) paste0(quotient, "/log(10)")
    else paste0("(", quotient, ")/log(10)")
  }, character(1))

  as.eqnvec(setNames(rhs, renamed))

}
