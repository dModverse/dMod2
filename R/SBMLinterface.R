#' Import an SBML Model
#'
#' Reads an SBML Level 3 file. Requires \pkg{rjson} and \pkg{reticulate}, which
#' provides a Python with `python-libsbml`; set the environment variable
#' `DMOD_LIBSBML_PYTHON` to the path of a Python interpreter to use that one
#' instead. Kinetic laws are divided by the volume of their compartment.
#'
#' @param modelpath Path to the SBML file.
#' @param keep Parameter ids that keep their own symbol although the SBML sets
#'   them by an `<initialAssignment>`; a constant assignment becomes their
#'   default value. Any other initial assignment is substituted into the rates,
#'   initial values and observables. Default `NULL`.
#'
#' @return A list with components
#'   \describe{
#'     \item{`reactions`}{The model as an [eqnlist].}
#'     \item{`pars`}{Named numeric, the parameter values.}
#'     \item{`inits`}{Named character, the initial value of each state.}
#'     \item{`observables`}{Named list of observable formulas, possibly
#'       empty.}
#'     \item{`assignmentRules`}{Named list of the formulas of the assignment
#'       rules and of the parameter initial assignments that `keep` does not
#'       name, already substituted into rates, initial values and observables.}
#'     \item{`events`}{An [eventlist] or `NULL`, including the rescaling of
#'       states when an event changes a compartment volume.}
#'     \item{`eventsSource`}{An [eventlist] or `NULL`, the event assignments
#'       of the SBML file only.}
#'     \item{`compartmentInits`}{Named character, the size at the start of each
#'       compartment an event assigns to; such a compartment keeps the id as volume.}
#'     \item{`renamed`}{Named character, the new ids of SBML ids that are no
#'       valid R or C++ names, named by the old ids.}
#'   }
#' @seealso [exportSbml()], [importPEtab()]
#' @examplesIf requireNamespace("reticulate", quietly = TRUE) && requireNamespace("rjson", quietly = TRUE) && reticulate::py_module_available("libsbml")
#' f <- addReaction(eqnlist(), from = "A", to = "B", rate = "k1*A")
#' f <- addReaction(f, from = "B", to = "", rate = "k2*B")
#' file <- file.path(tempdir(), "model.xml")
#' exportSbml(f, parameters = c(k1 = 0.5, k2 = 0.1), inits = c(A = 1, B = 0),
#'            filepath = file)
#' model <- importSbml(file)
#' model$reactions
#' model$pars
#' @export
#' @importFrom stringr str_replace_all
importSbml <- function(modelpath, keep = NULL) {

  .require_ns("reticulate", "SBML import")
  .require_ns("rjson", "SBML import")
  sbml <- .renameNonSyntactic(.sbmlPayload(normalizePath(modelpath, mustWork = TRUE)))

  # rjson collapses inner arrays of length 1 into a bare vector, so `S` is
  # reshaped from the known counts into `[n_reactions x n_states]`.
  S_raw <- sbml[["S"]]
  n_states <- length(sbml[["stateNames"]])
  n_rxns   <- length(sbml[["v"]])
  # A model whose dynamics live entirely in RateRules has no reactions; the
  # empty JSON payload would otherwise cbind into a list-mode matrix.
  S <- if (n_rxns == 0L || n_states == 0L) NULL
       else if (is.list(S_raw)) do.call(cbind, S_raw)
       else matrix(unlist(S_raw), nrow = n_rxns, ncol = n_states)
  if (!is.null(S)) {
    storage.mode(S) <- "double"
    S[S == 0] <- NA
  }

  # libsbml writes natural log as `ln(x)`, R and C expect `log(x)`. Every
  # formula passes through here; the leading boundary leaves `eln` intact.
  .normalizeFormula <- function(s) {
    if (length(s) == 0L) return(s)
    s <- stringr::str_replace_all(s, "\\*\\*", "^")
    s <- stringr::str_replace_all(s, "(^|[^A-Za-z0-9_.])ln\\(", "\\1log(")
    s
  }

  v <- .normalizeFormula(sbml[["v"]])

  states <- sbml[["stateNames"]]

  # Build compartment records + state-to-compartment map from the libsbml payload.
  compartments <- list()
  compartmentOf <- character(0)
  comp_json <- sbml[["compartments"]]
  spc_json  <- sbml[["speciesCompartments"]]
  # A compartment an event assigns to keeps the id as volume, so the kinetic laws follow the event;
  # the size at the start is in `compartmentInits`.
  ev_targets <- unlist(lapply(sbml[["events"]], function(ev)
    vapply(ev$assignments, function(a) a$variable, "")))
  compartmentInits <- character(0)
  if (!is.null(comp_json) && length(comp_json) > 0L) {
    for (c in comp_json) {
      # A unit-size compartment becomes the literal "1", keeping the compartment id out of
      # the kinetic laws; otherwise the id is a volume symbol the trafo can set.
      trivial <- !is.null(c$size) && is.numeric(c$size) && isTRUE(c$size == 1)
      if (c$id %in% ev_targets) {
        compartmentInits[[c$id]] <- if (!is.null(c$sizeAssignment)) .normalizeFormula(c$sizeAssignment)
                                    else if (is.numeric(c$size)) format(c$size, digits = 17) else "1"
        compartments[[c$id]] <- list(volume = c$id, rule = NULL)
        next
      }
      # An InitialAssignment on the compartment is the symbolic volume.
      volume <- if (!is.null(c$sizeAssignment)) .normalizeFormula(c$sizeAssignment)
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

  # rate_dMod = K / V_ref. Any compartment works as V_ref once recorded in
  # `reactionCompartment`; the educt side is preferred, the product side is the
  # fallback for pure synthesis.
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

  # SBML species with hasOnlySubstanceUnits hold amounts, not concentrations.
  amountStates <- intersect(unlist(sbml[["amountSpecies"]]), states)

  pars <- setNames(sbml[["p"]], sbml[["parameterNames"]])
  x0 <- setNames(.normalizeFormula(sbml[["x0"]]),
                 sbml[["stateNames"]])

  # AssignmentRules `lhs := rhs` are inlined into rates and initials up to a
  # fixed point, so chained rules resolve; their LHS are dropped from `pars`.
  rules <- sbml[["assignmentRules"]]
  # An initialAssignment on a constant parameter is substituted like a rule,
  # with species standing for their initial values. A parameter driven by a
  # RateRule is promoted to a state below instead.
  par_ia <- sbml[["parameterAssignments"]]
  for (nm in intersect(keep, names(par_ia))) {
    v0 <- tryCatch(eval(str2lang(.normalizeFormula(par_ia[[nm]])), .formulaEnv()),
                   error = function(e) NULL)
    if (is.numeric(v0) && length(v0) == 1L && nm %in% names(pars)) pars[[nm]] <- v0
  }
  ia_lhs <- setdiff(names(par_ia),
                    c(names(rules), names(sbml[["rateRules"]]), keep))
  if (length(ia_lhs)) {
    ia_rhs <- .normalizeFormula(unlist(par_ia[ia_lhs], use.names = FALSE))
    if (length(states))
      ia_rhs <- replaceSymbols(states, paste0("(", x0[states], ")"), ia_rhs)
    rules <- c(rules, setNames(as.list(ia_rhs), ia_lhs))
  }
  if (length(rules)) {
    rule_lhs <- names(rules)
    rule_rhs <- .normalizeFormula(unlist(rules, use.names = FALSE))
    # Parenthesized so substitution into a sub-expression keeps precedence.
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

  # ---- rate rules ----
  # A RateRule defines dX/dt = rhs and SBML forbids reactions on X. It is added
  # after the volume division above, since dC/dt is already intensive.
  rate_rules <- sbml[["rateRules"]]
  if (length(rate_rules)) {
    rr_lhs <- names(rate_rules)
    rr_rhs <- .normalizeFormula(unlist(rate_rules, use.names = FALSE))
    # Inline assignment-rule LHSs into the rate-rule RHSs too, so a RateRule
    # that references a rule-defined symbol doesn't keep it as a free var.
    if (length(rules)) {
      for (it in seq_len(max_iter)) {
        new_rr <- replaceSymbols(rule_lhs, rule_rhs, rr_rhs)
        if (identical(new_rr, rr_rhs)) break
        rr_rhs <- new_rr
      }
    }
    # A RateRule on a non-constant parameter or compartment promotes it to a
    # state with the SBML value as initial condition. Without a compartment it
    # gets unit volume, so dX/dt = rhs.
    promote <- setdiff(rr_lhs, states)
    if (length(promote)) {
      par_ia <- sbml[["parameterAssignments"]]
      init <- vapply(promote, function(nm) {
        ia <- par_ia[[nm]]
        if (!is.null(ia) && nzchar(ia)) return(.normalizeFormula(ia))
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
      # A virtual reaction with stoichiometry +1 on `var` and rate rhs. An amount
      # species has V_X = 1, so the volume is divided out to leave dn/dt = rhs.
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

  observables <- sbml[["observables"]]
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

  # ---- events ----
  # One "replace" row per eventAssignment. A `time >= T` trigger sets `time`,
  # any other trigger becomes the root `lhs - (rhs)`, which cppDE fires at crossings in both
  # directions; a pair `lhs >= rhs`, `lhs <= rhs` with the same assignments is one root event.
  events_json <- sbml[["events"]]
  events_df <- NULL
  events_src <- NULL
  if (length(events_json)) {
    rows <- list()
    # The source rows contain only what the SBML declared. The derived
    # compartment rescales are regenerated on the next import, so exporting
    # them would apply them twice.
    rows_src <- list()
    down <- down_src <- logical(0)
    cmp_pat <- "^\\s*(.+?)\\s*(>=|>|<=|<|==|!=)\\s*(.+?)\\s*$"
    for (ev in events_json) {
      tt <- ev[["triggerTime"]]
      tf <- ev[["triggerFormula"]]
      tt_num <- suppressWarnings(as.numeric(tt))
      time_val <- if (length(tt) && !is.null(tt) && !is.na(tt_num)) tt_num
                  else if (length(tt) && nzchar(tt)) as.character(tt)
                  else NA
      root_val <- NA_character_
      falling <- FALSE
      if (is.na(time_val) && length(tf) && nzchar(tf)) {
        m <- regmatches(tf, regexec(cmp_pat, tf))[[1L]]
        falling <- length(m) == 4L && m[3L] %in% c("<=", "<")
        root_val <- if (length(m) == 4L)
                      paste0("(", m[2L], ") - (", m[4L], ")")
                    else
                      tf
      }
      for (a in ev$assignments) {
        # An event assignment to a compartment changes the compartment size while species
        # amounts stay put, so concentrations scale inversely. Emit the
        # rescale first, while the old volume is still in effect.
        new_val <- .normalizeFormula(a$formula)
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
          down <- c(down, falling)
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
        down <- c(down, falling)
        down_src <- c(down_src, falling)
      }
    }
    # a falling row whose rising twin exists is the second half of one root event
    paired <- function(df, falling) {
      key <- do.call(paste, c(df[c("var", "value", "root")], sep = "\r"))
      df[!(falling & key %in% key[!falling]), , drop = FALSE]
    }
    if (length(rows)) {
      events_df <- as.eventlist(paired(do.call(rbind, rows), down))
      events_src <- as.eventlist(paired(do.call(rbind, rows_src), down_src))
    }
  }

  out <- list(reactions = reactions, pars = pars, inits = x0,
              observables = observables, assignmentRules = rules,
              events = events_df,
              eventsSource = events_src,
              compartmentInits = compartmentInits,
              renamed = attr(sbml, "renamed") %||% character(0))
  out
}


# Reserved in C++, so a state or compartment with one of these names cannot
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


# Formulas of SBML and PEtab files are evaluated here: an environment on emptyenv() with the
# arithmetic, the comparisons and the functions of the expression language of cppDE, nothing else.
.formulaEnv <- local({
  env <- NULL
  function() {
    if (!is.null(env)) return(env)
    e <- new.env(parent = emptyenv())
    base <- c("(", "+", "-", "*", "/", "^", "<", "<=", ">", ">=", "==", "!=", "&&", "||", "!",
              "&", "|", "exp", "log", "log2", "log10", "sqrt", "sin", "cos", "tan", "asin", "acos",
              "atan", "atan2", "sinh", "cosh", "tanh", "asinh", "acosh", "atanh", "abs", "min", "max",
              "floor", "ceiling", "round", "trunc", "sign", "ifelse", "pi")
    for (f in base) assign(f, get(f, envir = baseenv()), envir = e)
    fns <- list(
      ln = log, exp2 = function(x) 2^x, exp10 = function(x) 10^x, cbrt = function(x) sign(x) * abs(x)^(1/3),
      pow = function(x, y) x^y, root = function(x, n) x^(1/n), ceil = ceiling,
      cot = function(x) 1/tan(x), sec = function(x) 1/cos(x), csc = function(x) 1/sin(x),
      acot = function(x) atan(1/x), asec = function(x) acos(1/x), acsc = function(x) asin(1/x),
      coth = function(x) 1/tanh(x), sech = function(x) 1/cosh(x), csch = function(x) 1/sinh(x),
      acoth = function(x) atanh(1/x), asech = function(x) acosh(1/x), acsch = function(x) asinh(1/x),
      erf = function(x) 2 * stats::pnorm(x * sqrt(2)) - 1, erfc = function(x) 2 * stats::pnorm(-x * sqrt(2)),
      Heaviside = function(x) as.numeric(x > 0),
      piecewise = function(...) {
        a <- list(...); n <- length(a)
        for (i in seq_len(n %/% 2L)) if (isTRUE(as.logical(a[[2L * i]]))) return(a[[2L * i - 1L]])
        if (n %% 2L) a[[n]] else 0
      })
    list2env(fns, envir = e)
    lockEnvironment(e, bindings = TRUE)
    env <<- e
    e
  }
})


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
# to work on the string with separate identifier boundaries.
.renameIds <- function(x, map) {
  if (!length(map) || !length(x)) return(x)
  old <- names(map)
  for (i in seq_along(old))
    x <- gsub(paste0("(?<![A-Za-z0-9_.])", old[i], "(?![A-Za-z0-9_.])"),
              map[[i]], x, perl = TRUE)
  x
}

# The model as the list that code/sbmlImport.py writes to JSON under libsbml's
# Python. rjson reads arrays of length one as bare values.
.sbmlPayload <- function(path) {
  out <- tempfile(fileext = ".json")
  on.exit(unlink(out))
  script <- system.file("code/sbmlImport.py", package = "dMod2")
  status <- system2(.dmod_libsbml_python(), shQuote(c(script, path, out)))
  if (status != 0L || !file.exists(out))
    stop("SBML import failed (exit ", status, "). ",
         "Check that the libsbml virtualenv has python-libsbml installed.")
  rjson::fromJSON(file = out)
}


# SBML ids may start with an underscore or a digit, which R cannot parse. They
# are renamed on the raw JSON, before anything reads a formula.
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


#' Export an Equation List to SBML
#'
#' Writes an [eqnlist] with parameter values and initial values to an SBML
#' Level 3 Version 2 file. Kinetic laws are multiplied by the volume of their
#' compartment. Symbolic volumes are written as an `<initialAssignment>` of the
#' compartment. Requires the same Python environment as [importSbml()].
#'
#' @param eqnlist Object of class [eqnlist] to export.
#' @param parameters Named numeric vector of parameter values (including any
#'   compartment-size parameters referenced in `eqnlist$compartments`). An entry
#'   named like a compartment whose volume is the compartment symbol is the size
#'   of that compartment; entries named like other compartments or like species
#'   are not written. Default `NULL`: no `<parameter>` is written.
#' @param inits Named numeric or character vector of initial values keyed
#'   by state name. Numbers, also as character, are written as
#'   `initialConcentration`, or `initialAmount` for the states in
#'   `eqnlist$amountStates`; other expressions as `<initialAssignment>`
#'   formulas. Missing states default to 0. Default `NULL`.
#' @param filepath Path of the SBML file to write.
#' @param modelID SBML model identifier. Defaults to `"dMod_export"`.
#' @param events Object of class [eventlist] to emit as `<listOfEvents>`, or
#'   `NULL`. Rows sharing a trigger become one `<event>`. A root event becomes two events with
#'   the triggers `root >= 0` and `root <= 0`, initially true, so that the SBML model fires at
#'   crossings in both directions and not at the start, as cppDE does. A compartment an event
#'   assigns to is written as non-constant. Default `NULL`.
#' @return `filepath`, invisibly.
#' @seealso [importSbml()], [exportPEtab()]
#' @inherit importSbml examples
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
    stop("exportSbml: `eqnlist` has no compartments; assign them with assignCompartment().",
         call. = FALSE)

  # Compartments: numeric size when the volume expression parses as numeric;
  # otherwise 1.0 and rely on the homonymous parameter for the symbolic value.
  ev_targets <- if (is.null(events)) character(0)
                else unique(as.character(events$var))
  comp_list <- lapply(names(eqnlist$compartments), function(cid) {
    entry <- eqnlist$compartments[[cid]]
    vol <- entry$volume
    size <- suppressWarnings(as.numeric(vol))
    # A compartment whose volume expression is the compartment symbol has the
    # value in the parameter vector; write that as the size instead of a
    # self-referential assignment.
    if (is.na(size) && identical(vol, cid)) {
      if (!cid %in% names(parameters))
        stop("exportSbml: no size for compartment `", cid, "`; pass it in ",
             "`parameters`.", call. = FALSE)
      size <- suppressWarnings(as.numeric(parameters[[cid]]))
    }
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
  # law is `rate * V_ref`, the inverse of the division by the volume on import.
  smatrix <- eqnlist$smatrix
  vref_vol <- .refVolumes(.refCompartments(smatrix, eqnlist$compartmentOf[eqnlist$states],
                                           eqnlist$reactionCompartment, eqnlist$description),
                          eqnlist$compartments)
  # unname(): a named index vector would make rjson emit an object, not an array.
  rxn_rows <- if (is.null(smatrix)) integer(0)
              else unname(which(apply(!is.na(smatrix), 1L, any)))
  rxn_list <- lapply(rxn_rows, function(i) {
    row_i <- smatrix[i, ]
    # Unnamed indices, so rjson writes arrays rather than objects.
    educt_idx <- unname(which(!is.na(row_i) & row_i < 0))
    product_idx <- unname(which(!is.na(row_i) & row_i > 0))

    educts <- lapply(educt_idx, function(j)
      list(species = eqnlist$states[j], stoich = as.numeric(abs(row_i[j]))))
    products <- lapply(product_idx, function(j)
      list(species = eqnlist$states[j], stoich = as.numeric(row_i[j])))

    home_vol <- vref_vol[i]
    kinetic_law <- paste0("(", home_vol, ") * (", eqnlist$rates[i], ")")
    # SBML requires every species in the kinetic law to be a reactant, a
    # product or a declared modifier (a rate rule imported as a reaction
    # names the other states only in the rate).
    modifiers <- setdiff(intersect(getSymbols(eqnlist$rates[i]), eqnlist$states),
                         eqnlist$states[c(educt_idx, product_idx)])

    list(id = paste0("r", i),
         reactants = educts,
         products = products,
         modifiers = as.list(modifiers),
         kineticLaw = kinetic_law)
  })

  # Rows sharing a trigger belong to one <event>: SBML applies an event's
  # assignments together, using the values from before the event.
  event_list <- list()
  if (!is.null(events) && nrow(events)) {
    ev <- as.data.frame(events, stringsAsFactors = FALSE)
    # an SBML trigger fires when it turns true: a root that cppDE fires on crossings in both
    # directions becomes two triggers, true beforehand so that neither fires at the start
    isRoot <- is.na(ev$time)
    key <- ifelse(isRoot, as.character(ev$root), paste("time >=", ev$time))
    # SBML assigns a value: an added or multiplied value becomes the formula on the old one
    formula <- function(var, value, method)
      switch(if (is.null(method) || is.na(method)) "replace" else as.character(method),
             replace = as.character(value),
             add = sprintf("%s + (%s)", var, value), multiply = sprintf("%s * (%s)", var, value),
             stop("exportSbml: event method ", method, " is not replace, add or multiply.", call. = FALSE))
    for (k in unique(key)) {
      rows <- ev[key == k, , drop = FALSE]
      assignments <- lapply(seq_len(nrow(rows)), function(i)
        list(variable = rows$var[i], formula = formula(rows$var[i], rows$value[i], rows$method[i])))
      triggers <- if (isRoot[match(k, key)]) paste0("(", k, c(") >= 0", ") <= 0")) else k
      for (t in triggers)
        event_list[[length(event_list) + 1L]] <- list(trigger = t, initialValue = !identical(t, k),
                                                      assignments = assignments)
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
  # An explicit interpreter bypasses reticulate provisioning.
  override <- Sys.getenv("DMOD_LIBSBML_PYTHON", unset = "")
  python <- if (nzchar(override)) {
    if (!file.exists(override))
      stop("DMOD_LIBSBML_PYTHON=", override, " does not exist.")
    override
  } else {
    # `python-libsbml` was declared via reticulate::py_require() in .onLoad(),
    # but py_exe() reports a path only once Python is initialized. Before that
    # it falls back to discovery, misses the managed env and returns "".
    tryCatch({
      if (!reticulate::py_available(initialize = TRUE))
        stop("reticulate could not start Python")
      exe <- reticulate::py_exe()
      if (!nzchar(exe)) stop("reticulate reported no interpreter path")
      exe
    }, error = function(e) {
      stop("Could not provision a Python with python-libsbml via ",
           "reticulate (", conditionMessage(e), "). Set ",
           "DMOD_LIBSBML_PYTHON to point at a Python interpreter that ",
           "has python-libsbml installed.")
    })
  }

  # `import libsbml` is probed once per session, cached in an env var, to catch
  # a wrong interpreter or an unresolved requirement early.
  if (!identical(Sys.getenv("DMOD_LIBSBML_OK", unset = ""), "1")) {
    probe <- function() suppressWarnings(
      system2(python, args = c("-c", shQuote("import libsbml")),
              stdout = FALSE, stderr = FALSE))
    status <- probe()
    # an existing ~/.virtualenvs/r-reticulate ignores py_require(); install there
    if (status != 0L && !nzchar(override)) {
      message("Installing python-libsbml into ", python, ".")
      try(reticulate::py_install("python-libsbml", pip = TRUE), silent = TRUE)
      status <- probe()
    }
    if (status != 0L)
      stop("Python at ", python, " could not `import libsbml` ",
           "(status ", status, "). ",
           if (nzchar(override))
             "Install python-libsbml into that env."
           else
             "Install it with `pip install python-libsbml`.")
    Sys.setenv(DMOD_LIBSBML_OK = "1")
  }
  python
}
