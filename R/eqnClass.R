
## Class "eqnlist" and its constructor ------------------------------------------


#' Coerce to an equation list
#' @description Translates a reaction network, e.g. defined by a data.frame, into an equation list object.
#' @param ... additional arguments to be passed to or from methods.
#' @details If `data` is a `data.frame`, it must contain columns "Description" (character),
#' "Rate" (character), and one column per ODE state with the state names.
#' The state columns correspond to the stoichiometric matrix.
#' @return Object of class [eqnlist]
#' @rdname eqnlist
#' @export
as.eqnlist <- function(data, volumes, ...) {
  UseMethod("as.eqnlist", data)
}

#' @export
#' @param data data.frame with columns Description, Rate, and one colum for each state
#' reflecting the stoichiometric matrix
#' @rdname eqnlist
as.eqnlist.data.frame <- function(data, volumes = NULL, compartments = NULL, compartmentOf = NULL,
                                   reactionCompartment = NULL, amountStates = NULL, ...) {
  description <- as.character(data$Description)
  rates <- as.character(data$Rate)
  states <- setdiff(colnames(data), c("Description", "Rate"))
  smatrix <- as.matrix(data[, states]); colnames(smatrix) <- states

  # An explicit `volumes` is the legacy way of stating the layout: it must not
  # be overruled by a layout the data.frame merely carries along as attributes.
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
#' @param x object of class `eqnlist`
is.eqnlist <- function(x) {

  required <- c("smatrix", "states", "rates", "volumes", "description",
                "compartments", "compartmentOf", "reactionCompartment")

  # A reaction-less list is valid and may already carry a compartment layout.
  if (is.null(x$smatrix))
    return(length(x$states) == 0 &&
           length(x$rates) == 0 &&
           is.null(x$volumes) &&
           length(x$description) == 0 &&
           is.null(x$reactionCompartment) &&
           ((is.null(x$compartments) && is.null(x$compartmentOf)) ||
            (!is.null(x$compartments) && !is.null(x$compartmentOf) &&
             all(x$compartmentOf %in% names(x$compartments)))))

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
    all(names(x$smatrix) == names(x$states)) &&
    nrow(x$smatrix) == length(x$rates) &&
    ncol(x$smatrix) == length(x$states) &&
    refs && volumes && reactions && amounts
}


## Class "eqnlist" and its methods ------------------------------------------

#' Determine conserved quantites by finding the kernel of the stoichiometric
#' matrix
#'
#' @param S Stoichiometric matrix
#' @param weight One of `"none"` (default) or `"volume"`. When `"volume"`, the
#'   columns of `S` are multiplied by their compartment volume before the kernel
#'   is computed, so the returned quantities are conserved in *amount* rather
#'   than concentration. Requires `volumes` to be supplied as numeric.
#' @param volumes Optional named numeric vector of volume values keyed by state,
#'   aligned with `colnames(S)`. Only consulted when `weight = "volume"`.
#' @return Data frame with conserved quantities carrying an attribute with the
#'   number of conserved quantities.
#' @author Malenka Mader, \email{Malenka.Mader@@fdm.uni-freiburg.de}
#'
#' @example inst/examples/equations.R
#' @export
conservedQuantities <- function(S, weight = c("none", "volume"), volumes = NULL) {
  weight <- match.arg(weight)
  # Get kernel of S
  S[is.na(S)] <- 0
  if (weight == "volume") {
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
    S <- sweep(S, 2, v_num, "*")
  }
  v <- nullZ(S)
  n_cq <-  ncol(v)
  
  # Iterate over conserved quantities, removes 0s, etc.
  if (n_cq > 0) {
    if (is.null(colnames(S))) stop("Columns of stoichiometric matrix not named.") else variables <- colnames(S)
    cq <- matrix(nrow = ncol(v), ncol = 1)
    for (iCol in 1:ncol(v)) {
      is.zero <- v[, iCol] == 0
      # gsub, not sub: a kernel vector may carry more than one negative weight.
      cq[iCol, 1] <- gsub("+-", "-", paste0(v[!is.zero, iCol], "*", variables[!is.zero], collapse = "+"), fixed = TRUE)
    }
    
    colnames(cq) <- paste0("Conserved quantities: ", n_cq)
    cq <- as.data.frame(cq)
    attr(x = cq, which = "n") <- n_cq
    
  } else {
    cq <- c()
  }
  
  return(cq)
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

## Extract the linear coefficient vector of `expr` over `states` via 1-hot
## evaluation, then verify linearity by re-evaluating at the all-twos point
## (linear sum must equal 2 * sum(coefs); any product term breaks this).
## Returns NULL on unknown symbols or detected nonlinearity.
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

#' Validate a list of conservation expressions against an eqnlist
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

#' Conservation-quantity basis of an `eqnlist`
#'
#' Returns the conservation expressions associated with the model. If the
#' eqnlist carries user-defined `$totals` (set via [customTotals]) those are
#' returned verbatim; otherwise the auto-detected basis from
#' [conservedQuantities()] is rendered with smart `totalXxx` names from the
#' longest common substring of each CQ's species.
#'
#' @param eqnlist An [eqnlist].
#' @return Named list mapping `total_name` -> conservation expression
#'   (character). Empty list when the smatrix admits no CQs.
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

#' Set or reset user-defined conservation-quantity totals
#'
#' Attaches a named list of conservation expressions to the eqnlist; these
#' override auto-detection and flow into [Pimpl] / [Pequil] as the new
#' parameter basis. Each expression is validated against the stoichiometric
#' matrix (must lie in the left null space of `S`) and the basis as a whole
#' must have the same rank as `conservedQuantities(S)`. Pass `NULL` or
#' `list()` to reset to auto-detection.
#'
#' @param eqnlist An [eqnlist].
#' @param totals Named list of expressions, or `NULL` / `list()` to reset.
#' @return The eqnlist with `$totals` updated.
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


#' Generate a table of reactions (data.frame) from an equation list
#' 
#' @param eqnlist object of class [eqnlist]
#' @return `data.frame` with educts, products, rate and description. The first
#' column is a check if the reactions comply with reaction kinetics.
#' 
#' @example inst/examples/equations.R
#' @export
getReactions <- function(eqnlist) {
  
  # Extract information from eqnlist
  S <- eqnlist$smatrix
  rates <- eqnlist$rates
  description <- eqnlist$description
  variables <- eqnlist$states
  
  # Determine lhs and rhs of reactions
  if(is.null(S)) return()
  
  reactions <- apply(S, 1, function(v) {
    
    numbers <- v[which(!is.na(v))]
    educts <- -numbers[numbers < 0]
    products <- numbers[numbers > 0]
    educts <- paste(paste(educts, names(educts), sep = "*"), collapse=" + ")
    products <- paste(paste(products, names(products), sep = "*"), collapse=" + ")
    educts <- gsub("1*", "", educts, fixed = TRUE)
    products <- gsub("1*", "", products, fixed = TRUE)
    
    reaction <- paste(educts, "->", products)
    return(c(educts, products))
    
  })
  educts <- reactions[1,]
  products <- reactions[2,]
  
  # Check for consistency
  exclMarks.logical <- unlist(lapply(1:length(rates), function(i) {
    
    myrate <- rates[i]
    parsedRate <- getParseData(parse(text=myrate, keep.source = TRUE))
    symbols <- parsedRate$text[parsedRate$token=="SYMBOL"]
    
    educts <- variables[which(S[i,]<0)]
    
    !all(unlist(lapply(educts, function(e) any(e==symbols))))
    
  }))
  exclMarks <- rep(" ", ncol(reactions))
  exclMarks[exclMarks.logical] <- "!"
  

  # Generate data.frame  
  out <- data.frame(exclMarks, educts, "->", products, rates, description, stringsAsFactors = FALSE)
  colnames(out) <- c("Check", "Educt",  "->",  "Product", "Rate", "Description")
  rownames(out) <- 1:nrow(out)
  
  return(out)
  
}


#' Add reaction to reaction table
#'
#' @param eqnlist equation list, see [eqnlist]
#' @param from character with the left hand side of the reaction, e.g. "2*A + B"
#' @param to character with the right hand side of the reaction, e.g. "C + 2*D"
#' @param rate character. The rate associated with the reaction. The name is employed as a description
#' of the reaction.
#' @param description Optional description instead of `names(rate)`.
#' @param compartment Character, compartment ID for the states this reaction
#' *introduces*, and the frame the reaction is written in. Defaults to
#' `"defaultComp"`; created with volume `"1"` if new. States that already have a
#' compartment keep it -- use [assignCompartment()] to place a species that does
#' not belong to the compartment of the reaction first mentioning it.
#' @param rateCompartment Optional compartment ID naming the frame in which `rate`
#' is a concentration-rate. Needed when educts span multiple compartments (e.g.
#' membrane binding `L_ext + R_cyt -> Complex`); leave as `NA` (the default) to
#' let [getFluxes()] infer the frame from the educts. When the educts do span
#' compartments and `compartment` was given explicitly, that compartment is
#' taken as the frame.
#' @return An object of class [eqnlist].
#' @examples
#' f <- eqnlist()
#' f <- addReaction(f, "2*A+B", "C + 2*D", "k1*B*A^2")
#' f <- addReaction(f, "C + A", "B + A", "k2*C*A")
#'
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
  educts <- getSymbols(from)
  eductCoef <- 0
  if(length(educts) > 0) eductCoef <- sapply(educts, function(e) sum(getCoefficients(from, e)))
  products <- getSymbols(to)
  productCoef <- 0
  if(length(products) > 0) productCoef <- sapply(products, function(p) sum(getCoefficients(to, p)))


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


#' Assign states to compartments
#'
#' @description Declares which compartment a state lives in, independently of
#' the reactions that use it. [addReaction()] only ever assigns states it
#' introduces, so without an explicit declaration a species inherits the
#' compartment of whichever reaction happens to mention it first -- which makes
#' the model depend on the order in which it is written. `assignCompartment()`
#' removes that dependency: it works before the state exists (the declaration is
#' remembered and applied when the reaction arrives) as well as afterwards.
#'
#' @param eqnlist object of class [eqnlist]
#' @param ... named arguments `state = "compartment"`, or a single named
#' character vector of the same shape. Compartments that do not exist yet are
#' created with volume `"1"`.
#' @param volume Optional volume expression for the target compartment. Only
#' allowed when `...` names a single compartment; use [setCompartmentVolume()]
#' for several.
#' @param rule Optional volume rule (`dV/dt`) for the target compartment, see
#' [eqnlist].
#' @return An object of class [eqnlist].
#' @seealso [setCompartmentVolume()], [addReaction()]
#' @examples
#' # TGFb belongs to the extracellular space although the first reaction that
#' # mentions it is a cell-surface binding step.
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


#' Set compartment volumes
#'
#' @description Changes the volume expression (and optionally the volume rule)
#' of one or more compartments. Use this rather than assigning into
#' `eqnlist$compartments` directly: the per-state `$volumes` view is derived
#' from the layout and has to be recomputed, which in-place assignment does not
#' do.
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


# Reactions as a data.frame for the steady-state backend, which sees only rates
# and stoichiometry while getFluxes() scales every flux by V_ref / V_X. One csv
# row carries one rate, so a reaction touching states at different ratios is
# split into one row per ratio.
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
    stop("Compartment(s) ", paste(ruled, collapse = ", "), " carry a volume ",
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
  out
}


#' Generate list of fluxes from equation list
#' 
#' @param eqnlist object of class [eqnlist].
#' @param type "conc." or "amount" for fluxes in units of concentrations or
#' number of molecules. 
#' @return list of named characters, the in- and out-fluxes for each state.
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

  # Resolve per-reaction reference compartment V_ref (concentration-rate frame).
  # Priority: (1) user-supplied reactionCompartment[i] if non-NA, (2) unique
  # educt compartment, (3) unique product compartment for pure synthesis.
  # When educts span multiple compartments and no annotation is given, we
  # error with a clear message pointing the user at `reactionCompartment`.
  vref_cid <- .refCompartments(SMatrix, compOf, reactionCompartment, description)
  vref_vol <- .refVolumes(vref_cid, compartments)

  # generate equation expressions
  terme <- lapply(1:length(variables), function(j) {
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



#' Symbolic time derivative of equation vector given an equation list
#' 
#' The time evolution of the internal states is defined in the equation list.
#' Time derivatives of observation functions are expressed in terms of the
#' rates of the internal states.
#' 
#' @param observable named character vector or object of type [eqnvec]
#' @param eqnlist equation list
#' @details Observables are translated into an ODE
#' @return An object of class [eqnvec]
#' @example inst/examples/equations.R
#' @export
dot <- function(observable, eqnlist) {
  
 
  # Analyze the observable character expression
  symbols <- getSymbols(observable)
  states <- intersect(symbols, eqnlist$states)
  derivatives <- lapply(observable, function(obs) {
    out <- lapply(as.list(states), function(x) paste(deparse(D(parse(text=obs), x), width.cutoff = 500),collapse=""))
    names(out) <- states
    return(out)
  })
  
  # Generate equations from eqnist
  f <- as.eqnvec(eqnlist)
  
  newodes <- sapply(derivatives, function(der) {
    
    prodSymb(matrix(der, nrow = 1), matrix(f[names(der)], ncol = 1))
    
#     
#     out <- sapply(names(der), function(n) {
#       d <- der[n]
#       
#       if (d != "0") {
#         prodSymb(matrix(d, nrow = 1), matrix(f[names(d)], ncol = 1))
#       } else  {
#         return("0")
#       }
#         
#       
#       
#       #paste( paste("(", d, ")", sep="") , paste("(", f[names(d)], ")",sep=""), sep="*") else return("0")
#     })
#     out <- paste(out, collapse = "+")
#     
#     return(out)
    
  })
  
  as.eqnvec(newodes)
}



#' Coerce equation list into a data frame
#' 
#' @param x object of class [eqnlist]
#' @param ... other arguments
#' @return a `data.frame` with columns "Description" (character), 
#' "Rate" (character), and one column per ODE state with the state names. 
#' The state columns correspond to the stoichiometric matrix.
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

#' Write equation list into a csv file
#' 
#' @param eqnlist object of class [eqnlist]
#' @param ... Arguments going to [write.table][utils::write.table]
#' 
#' @export
#' @importFrom utils file.edit getParseData install.packages installed.packages read.csv str tail write.csv
write.eqnlist <- function(eqnlist, ...) {
  
  
  arglist <- list(...)
  argnames <- names(arglist)
  if (!"row.names" %in% argnames) arglist$row.names <- FALSE
  if (!"na" %in% argnames) arglist$na <- ""
  
  arglist$x <- as.data.frame(eqnlist)

  # The csv has no place for the layout; don't drop it silently.
  if (!.trivialCompartments(eqnlist$compartments) || !is.null(eqnlist$reactionCompartment))
    warning("write.eqnlist(): the csv format has no place for compartments; ",
            "volumes and reaction frames are lost on read-back.", call. = FALSE)

  do.call(write.csv, arglist)
  
}


#' subset of an equation list
#' 
#' @param x the equation list
#' @param ... logical expression for subsetting
#' @details The argument `...` can contain "Educt", "Product", "Rate" and "Description".
#' The "%in%" operator is modified to allow searches in Educt and Product (see examples).
#' 
#' @return An object of class [eqnlist]
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
  
  # Do selection on data.frame
  data <- getReactions(eqnlist)
  if(is.null(data)) return()
  
  data.list <- list(Educt = lapply(data$Educt, getSymbols), 
                    Product = lapply(data$Product, getSymbols),
                    Rate = data$Rate,
                    Description = data$Description,
                    Check = data$Check)
  
  "%in%" <- function(x, table) sapply(table, function(mytable) any(x == mytable))
  select <- which(eval(substitute(...), data.list))
  if (length(select) == 0) return(NULL)
  
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


#' Print or pander equation list
#' 
#' @param x object of class [eqnlist]
#' @param pander logical, use pander for output (used with R markdown)
#' @param ... additional arguments
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#' @author Daniel Kaschek, \email{daniel.kaschek@@physik.uni-freiburg.de}
#' 
#' @export
print.eqnlist <- function(x, pander = FALSE, ...) {

  eqnlist <- x

  totals <- getTotals(eqnlist)
  r <- getReactions(eqnlist)
  comp_lines <- .format_compartments(eqnlist$compartments, eqnlist$compartmentOf,
                                     eqnlist$states, eqnlist$amountStates)

  if (!pander) {
    if (length(totals)) {
      tag <- if (isTRUE(attr(eqnlist$totals, "custom"))) " (custom)" else ""
      cat("Conserved quantities", tag, ":\n", sep = "")
      for (nm in names(totals))
        cat("  ", nm, " = ", totals[[nm]], "\n", sep = "")
    }
    if (length(comp_lines) > 0L) {
      cat("\n")
      cat(comp_lines, sep = "\n")
    }
    cat("\n")
    print(r)
  } else if (requireNamespace("pander", quietly = TRUE)) {
    pander::panderOptions("table.alignment.default", "left")
    pander::panderOptions("table.split.table", Inf)
    pander::panderOptions("table.split.cells", Inf)
    exclude <- "Check"
    r <- r[, setdiff(colnames(r), exclude)]
    r$Rate <- paste0(format.eqnvec(as.character(r$Rate)))
    pander::pander(r)
  } else {
    print(r)
  }
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



#' Coerce to an equation vector
#' 
#' @param x object of class `character` or `eqnlist`
#' @param ... arguments going to the corresponding methods
#' @details If `x` is of class `eqnlist`, [getFluxes] is called and coerced
#' into a vector of equations.
#' @return object of class [eqnvec].
#' @export
as.eqnvec <- function(x, ...) {
  UseMethod("as.eqnvec", x)
}

#' Generate equation vector object
#'
#' @param names character, the left-hand sides of the equation
#' @rdname as.eqnvec
#' @export
as.eqnvec.character <- function(x = NULL, names = NULL, ...) {
  
  equations <- x
  
  if (is.null(equations)) return(NULL)
  
  if (is.null(names)) names <- names(equations)
  if (is.null(names)) stop("equations need names")
  if (length(names) != length(equations)) stop("Length of names and equations do not coincide")
  try.parse <- try(parse(text = equations), silent = TRUE)
  if (inherits(try.parse, "try-error")) stop("equations cannot be parsed: ", try.parse)
  
  out <- structure(equations, names = names)
  class(out) <- c("eqnvec", "character")
  
  return(out)
  
}



#' Transform equation list into vector of equations
#' 
#' @description An equation list stores an ODE in a list format. The function
#' translates this list into the right-hand sides of the ODE.
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
  inputs <- inputs[!vapply(inputs, function(x) is.null(x) || is.null(x$smatrix), logical(1))]
  if (length(inputs) == 0L) return(eqnlist())

  # Merge stoichiometry / rates / description via the data.frame path
  out <- Reduce(combine, lapply(inputs, as.data.frame))

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
      if (st %in% names(all_compartmentOf)) {
        if (!identical(unname(all_compartmentOf[[st]]), cid)) {
          stop(sprintf("State '%s' assigned to different compartments across input eqnlists.", st))
        }
      } else {
        all_compartmentOf[st] <- cid
      }
    }
  }

  if (length(all_compartments) == 0L) all_compartments <- NULL
  if (length(all_compartmentOf) == 0L) all_compartmentOf <- NULL

  # Concatenate reactionCompartment annotations. If any input has them, we need
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

  as.eqnlist(out, compartments = all_compartments, compartmentOf = all_compartmentOf,
             reactionCompartment = all_rc,
             amountStates = unique(unlist(lapply(inputs, function(el) el$amountStates))))

}


#' @export
#' @param x obect of any class
#' @rdname eqnvec
is.eqnvec <- function(x) {
  if (inherits(x, "eqnvec") &&
      length(x) == length(names(x))
  )
    return(TRUE)
  
  else
    return(FALSE)
}


## Class "eqnvec" and its methods --------------------------------------------





#' Encode equation vector in format with sufficient spaces
#' 
#' @param x object of class [eqnvec]. Alternatively, a named parsable character vector.
#' @param ... additional arguments
#' @return named character
#' @export format.eqnvec
#' @export
format.eqnvec <- function(x, ...) {
  
  eqnvec <- x
  
  eqns <- sapply(eqnvec, function(eqn) {
    parser.out <- getParseData(parse(text = eqn, keep.source = TRUE))
    parser.out <- subset(parser.out, terminal == TRUE)
    # parser.out$text[parser.out$text == "*"] <- "*" (avoid non-ASCII characters for CRAN)
    out <- paste(parser.out$text, collapse = "")
    return(out)
  })
  
  patterns <- c("+", "-", "*", "/")
  for (p in patterns) eqns <- gsub(p, paste0(" ", p, " "), eqns, fixed = TRUE)
  
  return(eqns)
    
  
}

#' Print equation vector
#' 
#' @param x object of class [eqnvec].
#' @param width numeric, width of the print-out
#' @param pander logical, use pander for output (used with R markdown)
#' @param ... not used right now
#' 
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#' 
#' @import stringr
#' @export
print.eqnvec <- function(x, width = 140, pander = FALSE, ...) {
  
  eqnvec <- x

  # Stuff to print
  m_odr <- "Idx"
  m_rel <- " <- "
  m_sep <- " "
  m_species <- names(eqnvec)
  
  # Width of stuff to print
  m_odrWidth <- max(3, nchar(m_odr))
  m_speciesWidth <- max(nchar(m_species), nchar("outer"))
  m_lineWidth <- max(width, m_speciesWidth + 10)
  m_relWidth <- nchar(m_rel)
  m_sepWidth <- nchar(m_sep)
  
  # Compound widths
  m_frontWidth <- m_odrWidth + m_speciesWidth + m_relWidth + m_sepWidth
  m_eqnWidth <- m_lineWidth - m_frontWidth
  
  # Order of states for alphabetical for print out. The Idx column is the
  # sequential print order (1, 2, 3, ...), not the pre-sort position.
  m_eqnOrder <- order(m_species)

  # Iterate over species
  m_msgEqn <- do.call(c, mapply(function(eqn, spec, odr) {
    return(paste0(
      str_pad(string = odr, side = "left", width = m_odrWidth),
      m_sep,
      str_pad(string = spec, side = "left", width = m_speciesWidth),
      m_rel,
      str_wrap(string = gsub(x = eqn, pattern = " ", replacement = "", fixed = TRUE),
               width = m_eqnWidth, exdent = m_frontWidth)
    ))
  }, eqn = eqnvec[m_eqnOrder], spec = m_species[m_eqnOrder],
     odr = seq_along(m_eqnOrder), SIMPLIFY = FALSE))
  
  # Print to command line or to pander
  if (!pander) {
    cat(paste0(str_pad(string = m_odr, side = "left", width = m_odrWidth),
               m_sep,
               str_pad(string = "Inner", side = "left", width = m_speciesWidth),
               m_rel,
               "Outer\n"))
    cat(m_msgEqn, sep = "\n")
  } else if (requireNamespace("pander", quietly = TRUE)) {
    pander::panderOptions("table.alignment.default", "left")
    pander::panderOptions("table.split.table", Inf)
    pander::panderOptions("table.split.cells", Inf)
    out <- as.data.frame(unclass(eqnvec), stringsAsFactors = FALSE)
    colnames(out) <- "" #  as.character(substitute(eqnvec))
    out[, 1] <- format.eqnvec(out[, 1])
    pander::pander(out)

  } else {
    out <- as.data.frame(unclass(eqnvec), stringsAsFactors = FALSE)
    colnames(out) <- ""
    out[, 1] <- format.eqnvec(out[, 1])
    print(out)
  }
  

}



#' Summary of an equation vector
#' 
#' @param object of class [eqnvec].
#' @param ... additional arguments
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
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
 
  out <- lapply(list(...), unclass)
  out <- do.call(c, out)
  if (any(duplicated(names(out)))) {
    stop("Names must be unique")
  }
  
  as.eqnvec(out)
}

#' @export
"[.eqnvec" <- function(x, ...) {
  out <- unclass(x)[...]
  class(out) <- c("eqnvec", "character")
  return(out)
}



#' Identify linear variables in an equation vector using sympy
#'
#' @param eqnvec An object of class `eqnvec`, representing a set of equations.
#' @details This function calls Python's `sympy` library via `reticulate` to symbolically analyze equations and determine if variables appear linearly in all equations.
#'
#' @return A character vector of variables that occur linearly in all equations.
#'
#' @examples
#' eqnvec <- as.eqnvec(
#'   c("-k1*A", "k1*A - k2*B", "-k3*B*C/(Km+C) + k4*pC", "k3*B*C/(Km+C) - k4*pC"),
#'   names = c("A", "B", "C", "pC")
#' )
#' getLinVars(eqnvec)
#'
#' @export
getLinVars <- function(eqnvec) {
  if (!inherits(eqnvec, "eqnvec")) {
    stop("Input 'eqnvec' must be of class 'eqnvec'.")
  }
  .require_ns("reticulate", "getLinVars()")

  sympy <- reticulate::import("sympy")
  sympy_zero <- sympy$Integer(0)
  
  variables <- names(eqnvec)
  sympy_vars <- lapply(variables, sympy$Symbol)
  sympy_eqns <- lapply(as.character(eqnvec), sympy$simplify)
  
  is_linear_in_eq <- function(eqn, var) {
    first_derivative <- sympy$diff(eqn, var)
    second_derivative <- sympy$diff(first_derivative, var)
    is_second_derivative_zero <- sympy$simplify(second_derivative) == sympy_zero
    is_first_derivative_nonzero <- sympy$simplify(first_derivative) != sympy_zero
    is_second_derivative_zero && is_first_derivative_nonzero
  }
  
  linear_vars <- sapply(seq_along(sympy_vars), function(i) {
    var <- sympy_vars[[i]]
    all(sapply(sympy_eqns, function(eqn) is_linear_in_eq(eqn, var)))
  })
  
  variables[linear_vars]
}






## eqnvec / eqnlist constructors (moved from classes.R) ----------------------------------------

## Equation classes -------------------------------------------------------

#' Generate equation vector object
#'
#' @description The eqnvec object stores explicit algebraic equations, like the
#' right-hand sides of an ODE, observation functions or parameter transformations
#' as named character vectors.
#' @param ... mathematical expressions as characters to be coerced,
#' the right-hand sides of the equations
#' @return object of class `eqnvec`, basically a named character.
#' @example inst/examples/eqnvec.R
#' @seealso [eqnlist]
#' @export
eqnvec <- function(...) {

  mylist <- list(...)
  if (length(mylist) > 0) {
    mynames <- paste0("eqn", 1:length(mylist))
    is.available <- !is.null(names(mylist))
    mynames[is.available] <- names(mylist)[is.available]

    names(mylist) <- mynames
    out <- unlist(mylist)

    return(as.eqnvec(out))

  } else {

    return(NULL)

  }

}

#' Generate eqnlist object
#'
#' @description The eqnlist object stores an ODE as a stoichiometric matrix,
#' rate expressions, state names, and compartment information.
#' @export
#' @param smatrix Numeric stoichiometric matrix; one row per reaction, one
#'   column per state.
#' @param states Character vector of state names.
#' @param rates Character vector of rate expressions.
#' @param volumes Named character of state volumes (kept for back-compat; when
#'   supplied without `compartments`/`compartmentOf`, distinct expressions are
#'   auto-assigned IDs `c1`, `c2`, ...).
#' @param description Character vector describing each reaction.
#' @param compartments Named list keyed by compartment ID; each entry is a
#'   volume expression (character) or a list with fields `volume` and `rule`
#'   (`rule` reserved for future dynamic-volume support, must be `NULL`).
#' @param compartmentOf Named character vector mapping state → compartment ID.
#'   States not listed default to compartment `"defaultComp"` with volume `"1"`.
#'   May also name states that do not exist yet; see [assignCompartment()].
#' @param reactionCompartment Optional character vector of length
#'   `nrow(smatrix)`. Per-reaction reference compartment ID; use `NA` to infer
#'   from educts. Required when educts span multiple compartments.
#' @param amountStates Optional character vector of states that carry substance
#'   units (amounts) rather than concentrations, SBML's `hasOnlySubstanceUnits`.
#'   Their fluxes are not divided by a compartment volume.
#' @param totals Optional named list of user-defined conservation-quantity
#'   expressions, as produced by [customTotals]. `NULL` leaves the basis to be
#'   auto-detected from the stoichiometric matrix.
#' @return An object of class `eqnlist`, basically a list.
#' @example inst/examples/eqnlist.R
eqnlist <- function(smatrix = NULL, states = colnames(smatrix), rates = NULL,
                    volumes = NULL, description = NULL,
                    compartments = NULL, compartmentOf = NULL,
                    reactionCompartment = NULL, amountStates = NULL, totals = NULL) {

  # Dimension checks and preparations for non-empty argument list.
  if (all(!is.null(c(smatrix, states, rates)))) {
    #Dimension checks
    d1 <- dim(smatrix)
    l2 <- length(states)
    l3 <- length(rates)
    if (l2 != d1[2]) stop("Number of states does not coincide with number of columns of stoichiometric matrix")
    if (l3 != d1[1]) stop("Number of rates does not coincide with number of rows of stoichiometric matrix")

    # Prepare variables
    smatrix <- as.matrix(smatrix)
    colnames(smatrix) <- states
    # "not involved" is NA, not 0, as in addReaction().
    smatrix[!is.na(smatrix) & smatrix == 0] <- NA
    if (is.null(description)) {
      description <- 1:nrow(smatrix)
    }
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

  if (!is.null(compartments) && !is.null(compartmentOf)) {
    compartments <- lapply(compartments, normalize_entry)
    compartmentOf <- setNames(as.character(compartmentOf), names(compartmentOf))
    if (is.null(names(compartmentOf)) || any(!nzchar(names(compartmentOf))))
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





#' Transform an ODE to log10 coordinates
#'
#' Reads `f` as the right-hand side of \eqn{\dot x = f(x)} and returns the same
#' system in \eqn{x_{l10} = \log_{10} x}. By the chain rule
#'
#' \deqn{\dot x_{l10} = \frac{f(x)}{\log(10)\, x}, \qquad x = 10^{x_{l10}}}
#'
#' so each right-hand side is divided by its own state and every remaining state
#' replaced by `10^` of the transformed one. The division is cancelled
#' symbolically first, which is what makes the result readable: `x = "-k*x"`
#' becomes `x_l10 = "-k/log(10)"`, not a ratio that happens to simplify at
#' runtime.
#'
#' Integrating in these coordinates cannot produce a negative state, since a power
#' of ten is positive for every real exponent -- positivity is the geometry of the
#' chart rather than a constraint checked afterwards. Where a trajectory would have
#' crossed \eqn{x = 0}, the transformed variable escapes to `-Inf` and the solver
#' stops there instead.
#'
#' The back-substitution is written `(10^(x_l10))`, not `exp10(x_l10)`, although
#' the latter reads better and works inside a [P()] transformation: `exp10` is a
#' C99 function with no entry in R's derivatives table, so a model built on it has
#' no symbolic Jacobian and [Xs()] cannot generate sensitivities. The parentheses
#' are load-bearing -- R's `^` is right-associative, so a bare `10^x_l10`
#' substituted into `x^2` would mean `10^(x_l10^2)`.
#'
#' Read results back in R with `10^x_l10`.
#'
#' @param f `eqnvec` (or named character vector) holding the right-hand sides.
#' @param suffix appended to each state name.
#' @param simplify cancel the division symbolically with sympy. `FALSE`, or an
#'   unavailable sympy, leaves the plain quotient, which is equivalent but wordy.
#' @return An `eqnvec` over the renamed states.
#' @seealso [odemodel()], which compiles the result.
#' @export
#' @examples
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
      # what is left over the division: (-A*k1 + B*k2)/A becomes -k1 + B*k2/A
      cancelled <- tryCatch(
        gsub("\\*\\*", "^", as.character(
          spy$expand(spy$cancel(spy$sympify(gsub("\\^", "**", quotient)))))),
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
