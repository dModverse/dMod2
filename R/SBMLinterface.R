#' Import an SBML model
#'
#' Reads an SBML Level 3 file via a Python helper (`inst/code/sbmlImport.py`)
#' that uses `python-libsbml`. The Python environment is provisioned
#' automatically by `reticulate` on first use, no manual venv setup is
#' required. Users who want to point dMod at an existing interpreter
#' (e.g. a hand-managed venv or conda env) can set the
#' `DMOD_LIBSBML_PYTHON` environment variable to its absolute path; that
#' bypasses reticulate entirely.
#'
#' Kinetic laws coming out of SBML are **extensive** (amount/time) by SBML
#' convention, whereas dMod stores rates in **concentration-style**. On import,
#' each kinetic law `K` is divided by the volume of its home compartment (the
#' compartment shared by the educts, or the product compartment for pure
#' synthesis) to produce `rate_dMod = K / V`. Combined with the volume-ratio
#' factors emitted by [getFluxes()], this preserves the SBML semantics:
#' `K_SBML = rate_dMod * V`.
#'
#' @param modelpath Path to the sbml file
#'
#' @return list of eqnlist, parameters and inits
#' @export
#' @importFrom stringr str_replace_all
importSbml <- function(modelpath) {

  .require_ns("reticulate", "SBML import")
  .require_ns("rjson", "SBML import")
  importscript <- system.file("code/sbmlImport.py", package = "dMod2")
  tmpfile_json <- tempfile()
  modelpath <- normalizePath(modelpath, mustWork = TRUE)

  venv_python <- .dmod_libsbml_python()

  status <- system2(venv_python,
                    args = c(shQuote(importscript),
                             shQuote(modelpath),
                             shQuote(tmpfile_json)))
  if (status != 0L || !file.exists(tmpfile_json))
    stop("SBML import failed (exit ", status, "). ",
         "Check that the libsbml virtualenv has python-libsbml installed.")
  json_content <- rjson::fromJSON(file = tmpfile_json)
  json_content <- .renameNonSyntactic(json_content)

  # `S` from the python side is a list-of-lists with shape
  # `[n_species][n_reactions]`. rjson::fromJSON collapses fully-nested
  # arrays whose inner length is 1 into a bare vector, which breaks
  # `cbind`'s list-of-columns contract. The reshape goes against the known
  # state/reaction counts so the resulting matrix is always
  # `[n_reactions rows x n_states cols]`, matching `eqnlist`'s convention.
  S_raw <- json_content[["S"]]
  n_states <- length(json_content[["stateNames"]])
  n_rxns   <- length(json_content[["v"]])
  # A model whose dynamics live entirely in RateRules has no reactions; the
  # empty JSON payload would otherwise cbind into a list-mode matrix.
  S <- if (n_rxns == 0L || n_states == 0L) NULL
       else if (is.list(S_raw)) do.call(cbind, S_raw)
       else matrix(unlist(S_raw), nrow = n_rxns, ncol = n_states)
  if (!is.null(S)) {
    storage.mode(S) <- "double"
    S[S == 0] <- NA
  }

  # libsbml L3 emits natural log as `ln(x)`, but both R's stats::D() and the
  # C math library expect `log(x)` for natural log. Apply this normalisation
  # to every formula channel, rates, initial-value expressions, and
  # AssignmentRule RHSs, at a single point. The leading boundary
  # `(^|[^A-Za-z0-9_.])` prevents matching `eln`, `arcln`, etc.
  .normalise_formula <- function(s) {
    if (length(s) == 0L) return(s)
    s <- stringr::str_replace_all(s, "\\*\\*", "^")
    s <- stringr::str_replace_all(s, "(^|[^A-Za-z0-9_.])ln\\(", "\\1log(")
    s
  }

  v <- .normalise_formula(json_content[["v"]])

  states <- json_content[["stateNames"]]

  # Build compartment records + state-to-compartment map from the libsbml payload.
  compartments <- list()
  compartmentOf <- character(0)
  comp_json <- json_content[["compartments"]]
  spc_json  <- json_content[["speciesCompartments"]]
  if (!is.null(comp_json) && length(comp_json) > 0L) {
    for (c in comp_json) {
      # Compartments with size = 1 (and no rule) carry no symbolic content --
      # storing them as the literal "1" keeps the compartment ID out of the
      # kinetic laws, which is what dMod's roundtrip expects when the source
      # eqnlist had volume "1". Otherwise use the SBML compartment ID as the
      # volume symbol so the trafo can override it.
      trivial <- !is.null(c$size) && is.numeric(c$size) && isTRUE(c$size == 1)
      # An InitialAssignment on the compartment is its symbolic volume.
      volume <- if (!is.null(c$sizeAssignment)) .normalise_formula(c$sizeAssignment)
                else if (trivial) "1" else c$id
      compartments[[c$id]] <- list(volume = volume, rule = NULL)
    }
    if (!is.null(spc_json) && length(spc_json) > 0L) {
      compartmentOf <- unlist(spc_json)
      compartmentOf <- compartmentOf[intersect(names(compartmentOf), states)]
    }
  }
  if (length(compartmentOf) == 0L) {
    compartments <- NULL
    compartmentOf <- NULL
  }

  # Normalize each kinetic law: rate_dMod = K / V_ref. Any compartment works as
  # V_ref, getFluxes() multiplies it back in, as long as the choice is
  # recorded, which is what `reactionCompartment` is for. The educt compartment
  # is the natural one; pure-synthesis reactions fall back to the product side.
  reactionCompartment <- NULL
  if (!is.null(compartmentOf) && !is.null(S)) {
    reactionCompartment <- rep(NA_character_, length(v))
    for (i in seq_along(v)) {
      row_i <- S[i, ]
      educt_idx <- which(!is.na(row_i) & row_i < 0)
      product_idx <- which(!is.na(row_i) & row_i > 0)
      home_cids <- if (length(educt_idx) > 0)   unique(compartmentOf[states[educt_idx]])
                   else if (length(product_idx) > 0) unique(compartmentOf[states[product_idx]])
                   else character(0)
      if (length(home_cids) == 0L) next
      if (length(home_cids) > 1L) reactionCompartment[i] <- home_cids[1]
      home_vol <- compartments[[home_cids[1]]]$volume
      if (!identical(home_vol, "1"))
        v[i] <- paste0("(", v[i], ")/(", home_vol, ")")
    }
    if (all(is.na(reactionCompartment))) reactionCompartment <- NULL
  }

  # SBML species with hasOnlySubstanceUnits carry amounts, not concentrations.
  amountStates <- intersect(unlist(json_content[["amountSpecies"]]), states)

  pars <- setNames(json_content[["p"]], json_content[["parameterNames"]])
  x0 <- setNames(.normalise_formula(json_content[["x0"]]),
                 json_content[["stateNames"]])

  # Inline AssignmentRules from the SBML model. Each rule `lhs := rhs` becomes
  # an algebraic substitution applied to all rates and species initials. PEtab
  # benchmark models (e.g. Boehm_JProteomeRes2014) use rules to encode
  # time-varying inputs like `BaF3_Epo := 1.25e-7 * exp(-k * time)`.
  # Iterate to a fixed point so chained rules resolve. After inlining, the
  # LHS symbols are no longer free parameters and are dropped from `pars`.
  rules <- json_content[["assignmentRules"]]
  if (length(rules)) {
    rule_lhs <- names(rules)
    rule_rhs <- .normalise_formula(unlist(rules, use.names = FALSE))
    # Wrap each RHS in parens so substitution into a sub-expression keeps
    # operator precedence intact (e.g. `1.25e-7 * exp(...)` inside `a * lhs`).
    rule_rhs <- paste0("(", rule_rhs, ")")
    max_iter <- length(rule_lhs) + 1L
    for (it in seq_len(max_iter)) {
      new_v  <- replaceSymbols(rule_lhs, rule_rhs, v)
      new_x0 <- replaceSymbols(rule_lhs, rule_rhs, x0)
      if (identical(new_v, v) && identical(new_x0, x0)) break
      v <- new_v; x0 <- new_x0
    }
    pars <- pars[setdiff(names(pars), rule_lhs)]
  }

  # --- rate rules ---
  # `<rateRule variable="X">` defines dX/dt = rhs. Per SBML spec, X cannot
  # also be produced/consumed by reactions, so the new column is independent
  # of existing kinetic laws. dC/dt is intensive already, so the rate is
  # appended AFTER the volume-division loop above (no /V wrap). RateRules on
  # non-species (parameter / compartment) are skipped with a warning.
  rate_rules <- json_content[["rateRules"]]
  if (length(rate_rules)) {
    rr_lhs <- names(rate_rules)
    rr_rhs <- .normalise_formula(unlist(rate_rules, use.names = FALSE))
    # Inline assignment-rule LHSs into the rate-rule RHSs too, so a RateRule
    # that references a rule-defined symbol doesn't carry it as a free var.
    if (length(rules)) {
      for (it in seq_len(max_iter)) {
        new_rr <- replaceSymbols(rule_lhs, rule_rhs, rr_rhs)
        if (identical(new_rr, rr_rhs)) break
        rr_rhs <- new_rr
      }
    }
    # A RateRule may target a non-constant parameter or compartment. That is a
    # dynamic quantity, so promote it to a state: its SBML value becomes the
    # initial condition and the rule becomes its rate. States without a
    # compartment land in the unit-volume default, which leaves dX/dt = rhs.
    promote <- setdiff(rr_lhs, states)
    if (length(promote)) {
      par_ia <- json_content[["parameterAssignments"]]
      init <- vapply(promote, function(nm) {
        ia <- par_ia[[nm]]
        if (!is.null(ia) && nzchar(ia)) return(.normalise_formula(ia))
        v0 <- suppressWarnings(as.numeric(pars[nm]))
        if (length(v0) != 1L || is.na(v0)) "0" else format(v0, digits = 17)
      }, character(1))
      if (!is.null(S)) S <- cbind(S, matrix(NA_real_, nrow(S), length(promote)))
      states <- c(states, promote)
      x0     <- c(x0, setNames(unname(init), promote))
      pars   <- pars[setdiff(names(pars), promote)]
    }

    for (k in seq_along(rr_lhs)) {
      var <- rr_lhs[k]; rhs <- rr_rhs[k]
      # Append a virtual reaction: stoichiometry +1 on `var`, 0 elsewhere;
      # rate string = rhs. After eqnlist construction this contributes
      # +rhs to the RHS row of `var`, which is exactly the rate rule.
      # An amount species has V_X = 1, so divide its volume out to leave dn/dt = rhs.
      if (var %in% amountStates && !is.null(compartmentOf)) {
        var_vol <- compartments[[unname(compartmentOf[var])]]$volume
        if (!identical(var_vol, "1")) rhs <- paste0("(", rhs, ")/(", var_vol, ")")
      }
      new_col <- rep(NA_real_, length(states))
      new_col[match(var, states)] <- 1
      S <- if (is.null(S)) matrix(new_col, nrow = 1L)
           else rbind(S, new_col)
      v <- c(v, rhs)
    }
  }

  # Virtual rate-rule reactions have no annotation; their frame is inferred.
  if (!is.null(reactionCompartment) && length(v) > length(reactionCompartment))
    reactionCompartment <- c(reactionCompartment,
                             rep(NA_character_, length(v) - length(reactionCompartment)))

  reactions <- eqnlist(smatrix = S, states = states, rates = v,
                       compartments = compartments, compartmentOf = compartmentOf,
                       reactionCompartment = reactionCompartment,
                       amountStates = amountStates)

  observables <- json_content[["observables"]]
  # An observable may read a rule-defined symbol; the rule LHS is no longer a
  # parameter after inlining, so the formula has to be inlined as well.
  if (length(rules) && length(observables)) {
    obs_chr <- unlist(observables, use.names = FALSE)
    for (it in seq_len(length(rule_lhs) + 1L)) {
      new_obs <- replaceSymbols(rule_lhs, rule_rhs, obs_chr)
      if (identical(new_obs, obs_chr)) break
      obs_chr <- new_obs
    }
    observables <- setNames(as.list(obs_chr), names(observables))
  }

  # --- events ---
  # SBML <event> -> dMod eventlist (one row per <eventAssignment>). Triggers
  # of the form `time >=/== T` (numeric or symbolic T) populate `time`; other
  # triggers fall back to a root expression of the form `lhs - (rhs)` so the
  # ODE solver can detect the zero crossing. Method is always "replace"
  # (SBML eventAssignments are assignment-style by spec).
  events_json <- json_content[["events"]]
  events_df <- NULL
  events_src <- NULL
  if (length(events_json)) {
    rows <- list()
    # The source rows carry only what the SBML declared. The derived
    # compartment rescales are regenerated on the next import, so exporting
    # them would apply them twice.
    rows_src <- list()
    cmp_pat <- "^\\s*(.+?)\\s*(>=|>|<=|<|==|!=)\\s*(.+?)\\s*$"
    for (ev in events_json) {
      tt <- ev[["triggerTime"]]
      tf <- ev[["triggerFormula"]]
      tt_num <- suppressWarnings(as.numeric(tt))
      time_val <- if (length(tt) && !is.null(tt) && !is.na(tt_num)) tt_num
                  else if (length(tt) && nzchar(tt)) as.character(tt)
                  else NA
      root_val <- NA_character_
      if (is.na(time_val) && length(tf) && nzchar(tf)) {
        m <- regmatches(tf, regexec(cmp_pat, tf))[[1L]]
        root_val <- if (length(m) == 4L)
                      paste0("(", m[2L], ") - (", m[4L], ")")
                    else
                      tf
      }
      for (a in ev$assignments) {
        # An event assignment to a compartment changes its size while species
        # amounts stay put, so concentrations scale inversely. Emit the
        # rescale first, while the old volume is still in effect.
        new_val <- .normalise_formula(a$formula)
        if (!is.null(compartmentOf) && a$variable %in% names(compartments)) {
          old_vol <- compartments[[a$variable]]$volume
          inside  <- setdiff(names(compartmentOf)[compartmentOf == a$variable],
                             amountStates)
          for (sp in intersect(inside, states))
            rows[[length(rows) + 1L]] <- data.frame(
              var    = sp,
              time   = time_val,
              value  = paste0(sp, " * (", old_vol, ") / (", new_val, ")"),
              root   = root_val,
              method = "replace",
              stringsAsFactors = FALSE)
        }
        src_row <- data.frame(
          var    = a$variable,
          time   = time_val,
          value  = new_val,
          root   = root_val,
          method = "replace",
          stringsAsFactors = FALSE
        )
        rows[[length(rows) + 1L]] <- src_row
        rows_src[[length(rows_src) + 1L]] <- src_row
      }
    }
    if (length(rows)) {
      events_df <- as.eventlist(do.call(rbind, rows))
      events_src <- as.eventlist(do.call(rbind, rows_src))
    }
  }

  out <- list(reactions = reactions, pars = pars, inits = x0,
              observables = observables, assignmentRules = rules,
              events = events_df,
              eventsSource = events_src,
              renamed = attr(json_content, "renamed") %||% character(0))
  out
}


# SBML ids may start with an underscore or a digit, which R cannot parse, and
# every symbolic step downstream goes through `parse()`. The rename happens on
# the raw JSON, before anything reads a formula.
# Reserved in C++, so a state or compartment carrying one of these names cannot
# reach the code generator unrenamed.
.CPP_KEYWORDS <- c(
  "alignas", "alignof", "and", "asm", "auto", "bool", "break", "case", "catch",
  "char", "class", "compl", "concept", "const", "consteval", "constexpr",
  "continue", "decltype", "default", "delete", "do", "double", "else", "enum",
  "explicit", "export", "extern", "false", "float", "for", "friend", "goto",
  "if", "inline", "int", "long", "mutable", "namespace", "new", "noexcept",
  "not", "nullptr", "operator", "or", "private", "protected", "public",
  "register", "requires", "return", "short", "signed", "sizeof", "static",
  "struct", "switch", "template", "this", "throw", "true", "try", "typedef",
  "typeid", "typename", "union", "unsigned", "using", "virtual", "void",
  "volatile", "wchar_t", "while", "xor", "std")


# Reserved in R, so a bare identifier of this name cannot be parsed.
.R_RESERVED <- c("if", "else", "repeat", "while", "function", "for", "next",
                 "break", "TRUE", "FALSE", "NULL", "Inf", "NaN", "NA",
                 "NA_integer_", "NA_real_", "NA_character_", "in")


# An identifier both R and SBML accept: an SBML SId allows letters, digits and
# underscores only, so `make.names()` is no good here, it produces dots.
.safeName <- function(x) {
  y <- gsub("[^A-Za-z0-9_]", "_", x)
  lead <- !grepl("^[A-Za-z]", y)
  y[lead] <- paste0("X", y[lead])
  y
}


# Internal: apply a rename map to expressions, textually. A non-syntactic id
# cannot go through `replaceSymbols()`, which parses, so the substitution has
# to work on the string with identifier boundaries of its own.
.renameIds <- function(x, map) {
  if (!length(map) || !length(x)) return(x)
  old <- names(map)
  for (i in seq_along(old))
    x <- gsub(paste0("(?<![A-Za-z0-9_.])", old[i], "(?![A-Za-z0-9_.])"),
              map[[i]], x, perl = TRUE)
  x
}

.renameNonSyntactic <- function(js) {
  ids <- unique(c(
    unlist(js[c("stateNames", "parameterNames")], use.names = FALSE),
    vapply(js$compartments %||% list(), function(z) z$id %||% "", ""),
    names(js$assignmentRules), names(js$rateRules), names(js$observables)))
  # `NA` is a legal SBML id and models use it for Avogadro's number. R parses
  # it as the logical constant, so it never shows up as a symbol and reaches
  # the code generator verbatim; it has to be renamed like any reserved word.
  ids <- unique(ids[!is.na(ids) & nzchar(ids)])
  bad <- ids[.safeName(ids) != ids | ids %in% .R_RESERVED]

  # A species or compartment named like a C++ keyword collides with the
  # generated sources. PEtab parameters are left alone: their names are part of
  # the problem description and must survive an export unchanged.
  own <- unique(c(unlist(js$stateNames, use.names = FALSE),
                  vapply(js$compartments %||% list(),
                         function(z) z$id %||% "", "")))
  bad <- unique(c(bad, intersect(own, .CPP_KEYWORDS)))
  if (!length(bad)) return(js)
  new <- .safeName(bad)
  # A name that is only reserved, in R or in C++, stays valid but has to
  # differ; `.safeName()` would leave it untouched.
  res <- new %in% c(.R_RESERVED, .CPP_KEYWORDS)
  new[res] <- paste0(new[res], "_")
  clash <- new %in% setdiff(ids, bad)
  new[clash] <- paste0(new[clash], "_")
  warning("SBML id(s) renamed to stay parsable in R: ",
          paste(bad, "->", new, collapse = ", "), call. = FALSE)

  # textual, on every character leaf and on every name in the tree
  fix <- function(x) .renameIds(x, setNames(new, bad))
  walk <- function(z) {
    if (is.character(z)) z <- fix(z)
    if (!is.null(names(z))) names(z) <- fix(names(z))
    if (is.list(z)) z[] <- lapply(z, walk)
    z
  }
  out <- walk(js)
  # The PEtab tables outside the model refer to these ids too and have to be
  # renamed with the same map.
  attr(out, "renamed") <- setNames(new, bad)
  out
}


#' Export an eqnlist to an SBML Level 3 file
#'
#' Serialises an [eqnlist] plus parameter values and initial concentrations to
#' SBML Level 3 Version 2. Each reaction's kinetic law is emitted as
#' `V_ref * rate_dMod` to restore SBML's extensive-flux convention
#' (`K_SBML = rate_dMod * V_ref`). Symbolic volumes are written as an
#' `<initialAssignment>` on the compartment. Requires a Python environment with `libsbml`
#' installed; the default location matches the one used by [importSbml()].
#'
#' @param eqnlist Object of class [eqnlist] to export.
#' @param parameters Named numeric vector of parameter values (including any
#'   compartment-size parameters referenced in `eqnlist$compartments`). An entry
#'   whose name is a compartment or a species is that element's size or initial
#'   value, not a separate `<parameter>`. Pass `NULL` to write parameters
#'   without values.
#' @param inits Named numeric *or* character vector of initial values keyed
#'   by state name. Numeric entries (or character entries that parse as
#'   numeric) are written as `initialConcentration`, `initialAmount` for
#'   states listed in `eqnlist$amountStates`; non-numeric character
#'   entries are emitted as `<initialAssignment>` formulas and let the SBML
#'   simulator resolve the expression against `parameters` at sim time.
#'   Missing states default to 0.
#' @param filepath Path to the SBML output file.
#' @param modelID SBML model identifier. Defaults to `"dMod_export"`.
#' @param events Object of class [eventlist] to emit as `<listOfEvents>`, or
#'   `NULL`. Rows sharing a trigger become one `<event>`, since SBML applies an
#'   event's assignments together. A compartment an event assigns to is written
#'   as non-constant.
#' @return `filepath`, invisibly.
#' @export
exportSbml <- function(eqnlist, parameters = NULL, inits = NULL, filepath,
                         modelID = "dMod_export", events = NULL) {

  .require_ns("reticulate", "SBML export")
  .require_ns("rjson", "SBML export")
  stopifnot(is.eqnlist(eqnlist))
  # A model without species needs no compartment layout; it exports as a
  # parameter-only SBML document.
  if (length(eqnlist$states) > 0L &&
      (is.null(eqnlist$compartments) || is.null(eqnlist$compartmentOf)))
    stop("`eqnlist` must have populated compartments/compartmentOf. Use the updated constructor.")

  # Compartments: numeric size when the volume expression parses as numeric;
  # otherwise 1.0 and rely on the homonymous parameter for the symbolic value.
  ev_targets <- if (is.null(events)) character(0)
                else unique(as.character(events$var))
  comp_list <- lapply(names(eqnlist$compartments), function(cid) {
    entry <- eqnlist$compartments[[cid]]
    vol <- entry$volume
    size <- suppressWarnings(as.numeric(vol))
    # A compartment whose volume expression is its own symbol carries the
    # value in the parameter vector; write that as the size instead of a
    # self-referential assignment.
    if (is.na(size) && identical(vol, cid) && cid %in% names(parameters))
      size <- suppressWarnings(as.numeric(parameters[[cid]]))
    out <- list(id = cid, size = if (is.na(size)) 1.0 else size,
                spatialDimensions = 3L,
                constant = !(cid %in% ev_targets))
    if (is.na(size)) out$sizeAssignment <- vol
    out
  })

  species_list <- lapply(eqnlist$states, function(st) {
    raw <- if (!is.null(inits) && st %in% names(inits)) inits[[st]] else 0
    num <- suppressWarnings(as.numeric(raw))
    amount <- st %in% eqnlist$amountStates
    base <- list(id = st, compartment = unname(eqnlist$compartmentOf[[st]]),
                 hasOnlySubstanceUnits = amount)
    if (!is.na(num)) {
      if (amount) c(base, list(initialAmount = num)) else c(base, list(initialConcentration = num))
    } else {
      # symbolic: emit as <initialAssignment>; the formula may reference any
      # parameter declared on the SBML side (incl. compartment-volume IDs).
      c(base, list(initialAssignment = as.character(raw)))
    }
  })

  # A name shared with a compartment or a species is already declared by that
  # element; re-declaring it as a parameter would duplicate the SBML id.
  declared <- c(names(eqnlist$compartments), eqnlist$states)
  param_list <- list()
  if (!is.null(parameters)) {
    nms  <- setdiff(names(parameters), declared)
    vals <- suppressWarnings(as.numeric(unlist(parameters[nms], use.names = FALSE)))
    # SBML needs a number on every <parameter>. A missing one would reach
    # libsbml as a null and fail there without naming the culprit.
    bad <- nms[!is.finite(vals)]
    if (length(bad))
      stop("exportSbml: no finite value for parameter(s) ",
           paste(bad, collapse = ", "), ".", call. = FALSE)
    for (i in seq_along(nms))
      param_list[[length(param_list) + 1L]] <- list(id = nms[i], value = vals[i])
  }

  # Each row of the stoichiometric matrix becomes one reaction. The kinetic
  # law is `rate * V_ref`, the α-bridge identity in the export direction.
  smatrix <- eqnlist$smatrix
  vref_vol <- .refVolumes(.refCompartments(smatrix, eqnlist$compartmentOf[eqnlist$states],
                                           eqnlist$reactionCompartment, eqnlist$description),
                          eqnlist$compartments)
  # unname(): a named index vector would make rjson emit an object, not an array.
  rxn_rows <- if (is.null(smatrix)) integer(0)
              else unname(which(apply(!is.na(smatrix), 1L, any)))
  rxn_list <- lapply(rxn_rows, function(i) {
    row_i <- smatrix[i, ]
    # `which()` on a named vector preserves names, which would propagate
    # through lapply() into a *named* list, rjson then serialises it as
    # a JSON object, breaking the array-of-dicts contract dmodToSbml.py
    # expects. unname() the indices.
    educt_idx <- unname(which(!is.na(row_i) & row_i < 0))
    product_idx <- unname(which(!is.na(row_i) & row_i > 0))

    educts <- lapply(educt_idx, function(j)
      list(species = eqnlist$states[j], stoich = as.numeric(abs(row_i[j]))))
    products <- lapply(product_idx, function(j)
      list(species = eqnlist$states[j], stoich = as.numeric(row_i[j])))

    home_vol <- vref_vol[i]
    kinetic_law <- paste0("(", home_vol, ") * (", eqnlist$rates[i], ")")

    list(id = paste0("r", i),
         reactants = educts,
         products = products,
         kineticLaw = kinetic_law)
  })

  # Rows sharing a trigger belong to one <event>: SBML applies an event's
  # assignments together, using the values from before the event.
  event_list <- list()
  if (!is.null(events) && nrow(events)) {
    ev <- as.data.frame(events, stringsAsFactors = FALSE)
    trig <- ifelse(is.na(ev$time), as.character(ev$root),
                   paste("time >=", ev$time))
    for (t in unique(trig)) {
      rows <- ev[trig == t, , drop = FALSE]
      event_list[[length(event_list) + 1L]] <- list(
        trigger = t,
        assignments = lapply(seq_len(nrow(rows)), function(i)
          list(variable = rows$var[i], formula = as.character(rows$value[i]))))
    }
  }

  spec <- list(modelId = modelID,
               compartments = comp_list,
               species = species_list,
               parameters = param_list,
               reactions = rxn_list,
               events = event_list,
               outfile = normalizePath(filepath, mustWork = FALSE))

  spec_json <- tempfile(fileext = ".json")
  writeLines(rjson::toJSON(spec), spec_json)

  script <- system.file("code/dmodToSbml.py", package = "dMod2")
  venv_python <- .dmod_libsbml_python()
  err <- system2(venv_python, args = c(shQuote(script), shQuote(spec_json)),
                 stdout = TRUE, stderr = TRUE)
  status <- attr(err, "status")
  if (!is.null(status) && status != 0L)
    stop("SBML export failed (exit ", status, "): ",
         paste(utils::tail(err, 4), collapse = " | "), call. = FALSE)

  invisible(filepath)
}


.dmod_libsbml_python <- function() {
  # Explicit override wins (existing user envs, conda, CI with prebuilt
  # interpreters, ...). Skip reticulate provisioning entirely.
  override <- Sys.getenv("DMOD_LIBSBML_PYTHON", unset = "")
  python <- if (nzchar(override)) {
    if (!file.exists(override))
      stop("DMOD_LIBSBML_PYTHON=", override, " does not exist.")
    override
  } else {
    # `python-libsbml` was declared via reticulate::py_require() in
    # .onLoad(). py_exe() materialises the managed env (downloads Python +
    # installs the requirement on first call) and returns the interpreter
    # path used by the system2() calls below.
    tryCatch(reticulate::py_exe(), error = function(e) {
      stop("Could not provision a Python with python-libsbml via ",
           "reticulate (", conditionMessage(e), "). Set ",
           "DMOD_LIBSBML_PYTHON to point at a Python interpreter that ",
           "has python-libsbml installed.")
    })
  }

  # Probe `import libsbml` once per session. For the override path this
  # catches a wrong interpreter early; for the reticulate path it is a
  # cheap sanity check that the requirement actually resolved. Cached via
  # an env var so the ~30 ms python spawn does not repeat across
  # importSbml / exportSbml calls.
  if (!identical(Sys.getenv("DMOD_LIBSBML_OK", unset = ""), "1")) {
    status <- suppressWarnings(
      system2(python, args = c("-c", shQuote("import libsbml")),
              stdout = FALSE, stderr = FALSE))
    if (status != 0L)
      stop("Python at ", python, " could not `import libsbml` ",
           "(status ", status, "). ",
           if (nzchar(override))
             "Install python-libsbml into that env."
           else
             "reticulate did not provision python-libsbml as expected.")
    Sys.setenv(DMOD_LIBSBML_OK = "1")
  }
  python
}
