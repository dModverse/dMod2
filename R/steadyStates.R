#' Analytical Steady States
#'
#' Computes symbolic steady-state expressions for use in a parameter
#' transformation, following the AlyssaPetit method \[1\], \[2\]. Each balance
#' equation is solved for a state or a rate constant such that the result is
#' a ratio of sums of positive terms. A state whose equation was spent on a
#' rate constant stays a free parameter. The solver is a Python module called
#' through \pkg{reticulate}; it needs \pkg{numpy}, \pkg{sympy} and, from
#' version `"1.2"`, \pkg{scipy}.
#'
#' @param model An [eqnlist], or the path of a csv file describing the model.
#' @param file Character, path of the RDS file the result is written to with
#'   `saveRDS()`. For an `eqnlist` the model is also written to
#'   `<file>_model.csv` for the solver. Default `NULL`: for an `eqnlist` this
#'   is `"reactions_for_Alyssa"`, so both files land in the working directory;
#'   for a csv model nothing is written.
#' @param rates Not used.
#' @param forcings Character vector, names of forcings. These states are held
#'   at zero and treated as exogenous. Default `NULL`.
#' @param givenCQs Unnamed character vector of conserved quantities, either as
#'   `c("A + pA = totA", "B + pB = totB")` or as `c("A + pA", "B + pB")`.
#'   Default `NULL`: a basis is derived, or taken from [customTotals()] when the
#'   `eqnlist` has them. In `"1.4"` each given quantity keeps one of its states
#'   free, the first where possible.
#' @param neglect Character vector, states and rate parameters the solver must
#'   not solve for. Default `NULL`.
#' @param sparsifyLevel Numeric, upper bound on the length of the linear
#'   combinations used to simplify the stoichiometric matrix. Default `NULL`:
#'   `2` for versions `"1.0"` and `"1.1"`, unused from `"1.2"`.
#' @param outputFormat `"R"` (default) for dMod, or `"M"` for d2d \[3\].
#' @param testSteady How the solution is verified: `"fast"` (default), a
#'   probabilistic check modulo a prime; `"exact"`, symbolic substitution;
#'   `"skip"`.
#' @param walltime Integer, time budget of the solver in seconds, `0`
#'   (default) for none.
#' @param simplify `TRUE` (default) simplifies every expression once with
#'   sympy, `FALSE` skips it, `"full"` simplifies further at a higher cost.
#' @param solveQuadratic Logical, whether a balance quadratic in its own state
#'   may be solved by its positive root instead of for a rate constant, which
#'   introduces `sqrt()` terms. Default `FALSE`.
#' @param positive Positivity assumption for root and pivot selection: `TRUE`
#'   (default) all symbols are positive, `FALSE` none, or a character vector of
#'   the positive symbols. With `TRUE` no result contains a subtraction.
#' @param branches Logical. If `TRUE`, a quadratic with two positive roots is
#'   returned with a selector symbol `branch_<state>` taking `-1` or `+1`.
#'   Requires `solveQuadratic = TRUE`. Default `FALSE`.
#' @param priority Character vector of state and rate-parameter names, most
#'   preferred first: named states are solved earlier, named rate parameters
#'   are preferred as pivots. Default `NULL`.
#' @param resolve Logical, whether the result is passed through
#'   [resolveRecurrence()] so that no equation refers to another. Default
#'   `TRUE`.
#' @param verbose `TRUE` (default) reports progress in a few lines and the
#'   result, `FALSE` only the result, `"full"` every step of the solver.
#' @param version Solver version, one of `"1.0"`, `"1.1"`, `"1.2"`, `"1.3"`
#'   (default) and `"1.4"`. `"1.4"` solves every unknown from a combination of
#'   balances, so every result is a ratio of positive sums; it needs
#'   `positive = TRUE` and `outputFormat = "R"` and treats compartment volumes
#'   as fixed.
#'
#' @details
#' Arguments by version:
#' \tabular{ll}{
#'   `"1.0"` \tab `sparsifyLevel`; the solution is always tested.\cr
#'   `"1.1"` \tab adds `testSteady` (`"fast"` runs as `"exact"`).\cr
#'   `"1.2"`, `"1.3"` \tab add `walltime`, `simplify`, `solveQuadratic`,
#'     `positive`, `branches`, `priority`; `sparsifyLevel` is unused.\cr
#'   `"1.3"`, `"1.4"` \tab honour `verbose` themselves; older versions are
#'     silenced with `verbose = FALSE`.\cr
#'   `"1.4"` \tab ignores `branches`; `testSteady = "exact"` runs as
#'     `"fast"`.\cr
#' }
#'
#' @return Named character vector of steady-state equations in dMod format, or
#'   `0` if no solution was found. An entry whose value is its own name is a
#'   free parameter.
#'
#' @references \[1\] <https://pmc.ncbi.nlm.nih.gov/articles/PMC4863410/>
#' @references \[2\] <https://github.com/marcusrosenblatt/AlyssaPetit>
#' @references \[3\] <https://github.com/Data2Dynamics/d2d>
#'
#' @author Marcus Rosenblatt, \email{marcus.rosenblatt@@fdm.uni-freiburg.de}
#'
#' @seealso [repar()], [P()]
#'
#' @export
#' @importFrom utils write.table
#' @examplesIf requireNamespace("reticulate", quietly = TRUE) && reticulate::py_module_available("sympy") && reticulate::py_module_available("scipy")
#' reactions <- eqnlist()
#' reactions <- addReaction(reactions, "Tca_buffer", "Tca_cyto",
#'                          "import_Tca*Tca_buffer", "Basolateral uptake")
#' reactions <- addReaction(reactions, "Tca_cyto", "Tca_buffer",
#'                          "export_Tca_baso*Tca_cyto", "Basolateral efflux")
#' reactions <- addReaction(reactions, "Tca_cyto", "Tca_canalicular",
#'                          "export_Tca_cana*Tca_cyto", "Canalicular efflux")
#' reactions <- addReaction(reactions, "Tca_canalicular", "Tca_buffer",
#'                          "transport_Tca*Tca_canalicular", "Transport bile")
#'
#' steadies <- steadyStates(reactions, file = file.path(tempdir(), "steady"),
#'                          verbose = FALSE)
#' steadies
#'
#' # Parameter transformation that puts the model in steady state
#' parameters <- getParameters(reactions)
#' trafo <- repar("x ~ y", eqnvec(setNames(parameters, parameters)),
#'                x = names(steadies), y = steadies)
#' trafo
steadyStates <- function(model, file = NULL, rates = NULL, forcings = NULL,
                         givenCQs = NULL, neglect = NULL, sparsifyLevel = NULL,
                         outputFormat = "R", testSteady = c("fast", "exact", "skip"),
                         walltime = 0L, simplify = TRUE, solveQuadratic = FALSE,
                         positive = TRUE, branches = FALSE, priority = NULL,
                         version = "1.3", resolve = TRUE, verbose = TRUE) {

  .require_ns("reticulate", "steadyStates()")
  # Validate version and verification mode
  version <- match.arg(version, choices = c("1.0", "1.1", "1.2", "1.3", "1.4"))
  testSteady <- match.arg(testSteady)
  if (!(isTRUE(verbose) || isFALSE(verbose) || identical(verbose, "full")))
    stop("verbose must be TRUE, FALSE or \"full\".")

  # Forward customTotals() as givenCQs (the CSV drops totals metadata, so the
  # backend would otherwise re-derive an arbitrary conserved-quantity basis).
  if (is.null(givenCQs) && inherits(model, "eqnlist") &&
      !is.null(model$totals) && isTRUE(attr(model$totals, "custom"))) {
    givenCQs <- unname(paste(unlist(model$totals), names(model$totals),
                             sep = " = "))
  }

  # Default sparsifyLevel depends on version: v1.0/v1.1 still use sparsify
  # (default 2), v1.2+ ignores it (default 0, avoids the info print).
  if (is.null(sparsifyLevel)) {
    sparsifyLevel <- if (version %in% c("1.2", "1.3", "1.4")) 0 else 2
  }

  # Check if model is an equation list
  if (inherits(model, "eqnlist")) {
    # 1.4 would take a volume factor for a rate constant: volumes are not unknowns
    if (version == "1.4" && length(model$volumes))
      neglect <- union(neglect, getSymbols(unique(unlist(model$volumes))))
    if (is.null(file)) file <- "reactions_for_Alyssa"
    # Not write.eqnlist(): the backend never sees the volumes, so the
    # V_ref / V_X factors getFluxes() applies have to be folded in first.
    tab <- .volumeScaledReactions(model)
    # State volumes let the backend test for sink clusters in amounts, across
    # the rows a reaction is split into
    volumes <- attr(tab, "volumes")
    # Names sympy resolves to its own objects (Ci, Si, E, S, Q, gamma, ...) reach the
    # backend under an alias: its implicit sympify() of strings would otherwise build
    # Symbol*Ci and fail. Mapped back on the result below.
    symAlias <- .ssSympyAliases(unique(c(names(tab)[-(1:2)], getSymbols(tab$Rate),
                                         forcings, neglect, priority,
                                         if (is.character(positive)) positive,
                                         getSymbols(givenCQs), getSymbols(volumes))))
    if (length(symAlias)) {
      ren <- function(x) .ssRename(x, symAlias)
      tab$Rate <- ren(tab$Rate)
      hit <- names(tab) %in% names(symAlias)
      names(tab)[hit] <- symAlias[names(tab)[hit]]
      volumes <- setNames(ren(volumes), ren(names(volumes)))
      forcings <- ren(forcings); neglect <- ren(neglect); priority <- ren(priority)
      if (is.character(positive)) positive <- ren(positive)
      if (length(givenCQs)) givenCQs <- ren(givenCQs)
    }
    utils::write.csv(tab, file = paste0(file, "_model.csv"),
                     row.names = FALSE, na = "")
    model <- paste0(file, "_model.csv")
  } else symAlias <- volumes <- character(0)
  if (!is.null(givenCQs) && length(names(givenCQs)) > 0)
    stop("givenCQs must not have names. Please unname() them.")

  # Ensure Python dependencies are available
  reticulate::py_require("numpy")
  reticulate::py_require("sympy")
  # v1.2+ uses scipy.optimize.linprog for structural sink-cluster detection
  # (states whose combined mass leaks monotonically and must therefore be 0).
  # Older versions don't need scipy.
  if (version %in% c("1.2", "1.3", "1.4")) reticulate::py_require("scipy")

  pymodule <- paste0("AlyssaPetit_ver", gsub("\\.", "_", version))
  ap <- reticulate::import_from_path(pymodule,
                                     path = system.file("code", package = "dMod2"))
  # 1.3 honours verbose itself; the older backends are silenced from outside
  alyssa <- function(...) {
    if (version %in% c("1.3", "1.4") || !isFALSE(verbose)) return(ap$Alyssa(...))
    res <- NULL
    reticulate::py_capture_output(res <- ap$Alyssa(...))
    res
  }

  # Version-specific Python signatures:
  #   v1.0: Alyssa(filename, injections, givenCQs, neglect, sparsifyLevel, outputFormat)
  #   v1.1: Alyssa(filename, injections, givenCQs, neglect, sparsifyLevel, outputFormat, testSteady)
  #   v1.2/v1.3 (same signature): Alyssa(filename, injections, givenCQs, neglect, sparsifyLevel, outputFormat, testSteady, walltime, simplify, solveQuadratic, positive, branches, priority)
  #        -- v1.2 additionally runs structural sink-cluster detection a priori,
  #          and (when `solveQuadratic=TRUE`) attempts a closed-form quadratic
  #          state-side solve before resorting to flux-parameter pivots.
  #        -- v1.3 records solutions lazily (one textual resolution at output
  #          time) and guards direct solves against lock deadlocks.
  if (version == "1.0") {
    if (testSteady == "skip")
      message("Note: version 1.0 does not support testSteady='skip', test will always run.")
    if (isTRUE(solveQuadratic))
      message("Note: version 1.0 does not support solveQuadratic=TRUE, ignored.")
    if (length(priority) > 0)
      message("Note: version 1.0 does not support 'priority', ignored.")
    m_ss <- alyssa(model, as.list(forcings), as.list(givenCQs),
                      as.list(neglect), sparsifyLevel, outputFormat)

  } else if (version == "1.1") {
    if (isTRUE(solveQuadratic))
      message("Note: version 1.1 does not support solveQuadratic=TRUE, ignored.")
    if (length(priority) > 0)
      message("Note: version 1.1 does not support 'priority', ignored.")
    if (testSteady == "fast")
      message("Note: version 1.1 has no 'fast' test, using 'exact' instead.")
    # v1.1 backend takes the legacy "T"/"F" tokens.
    m_ss <- alyssa(model, as.list(forcings), as.list(givenCQs),
                      as.list(neglect), sparsifyLevel, outputFormat,
                      if (testSteady == "skip") "F" else "T")

  } else {
    # v1.2 / v1.3 / v1.4 (shared signature)
    # 1.4 reports what it ignores after its own summary, in one line
    ignored <- character(0)
    if (version == "1.4") {
      ignored <- c(if (isTRUE(branches)) "branches",
                   if (testSteady == "exact") "testSteady = \"exact\" (runs \"fast\")")
    }
    # simplify can be TRUE / FALSE / "full" -- pass through untouched so the
    # Python side sees either a Python bool or the literal string "full".
    if (is.character(simplify)) {
      simplify <- match.arg(tolower(simplify), choices = "full")
      simplify_arg <- simplify
    } else {
      simplify_arg <- as.logical(simplify)
    }
    # bool -> Python bool, character vector -> Python list of symbol names.
    positive_arg <- if (is.character(positive)) as.list(positive) else as.logical(positive)
    args <- list(model,
                 injections     = as.list(forcings),
                 givenCQs       = as.list(givenCQs),
                 neglect        = as.list(neglect),
                 sparsifyLevel  = as.integer(sparsifyLevel),
                 outputFormat   = outputFormat,
                 testSteady     = testSteady,
                 walltime       = as.integer(walltime),
                 simplify       = simplify_arg,
                 solveQuadratic = as.logical(solveQuadratic),
                 positive       = positive_arg,
                 branches       = as.logical(branches),
                 priority       = as.list(as.character(priority)))
    if (version %in% c("1.3", "1.4")) {
      args$verbose <- verbose
      args$volumes <- as.list(volumes)
    }
    m_ss <- do.call(alyssa, args)
    if (length(ignored))
      message("  note: version 1.4 ignores ", paste(ignored, collapse = ", "))
  }

  if (is.null(m_ss) || identical(m_ss, 0L)) return(0)

  # All versions return a list of "lhs=rhs" strings.
  # Parse into named character vector (dMod format).
  m_ssChar <- do.call(c, lapply(strsplit(m_ss, "="), function(eq) {
    out <- trimws(eq[2])
    names(out) <- trimws(eq[1])
    return(out)
  }))

  if (length(m_ssChar) == 0) return(0)
  if (length(symAlias)) {
    back <- setNames(names(symAlias), symAlias)
    m_ssChar <- setNames(.ssRename(m_ssChar, back),
                         ifelse(names(m_ssChar) %in% names(back), back[names(m_ssChar)],
                                names(m_ssChar)))
  }

  # Versions 1.2+ resolve on the backend side; this covers the older ones. Only
  # rewrite when there is something to resolve, resolveRecurrence() reformats
  # every equation it touches, which would break printed == returned.
  if (resolve && outputFormat == "R" && length(m_ssChar) > 1) {
    # Entries kept as their own name are free parameters, not definitions.
    passthrough <- names(m_ssChar)[trimws(m_ssChar) == names(m_ssChar)]
    recurrent <- function(x) vapply(seq_along(x), function(i) {
      any(setdiff(getSymbols(x[i]), passthrough) %in% names(x)[-i])
    }, logical(1))
    if (any(recurrent(m_ssChar))) {
      m_ssChar <- resolveRecurrence(m_ssChar)
      left <- recurrent(m_ssChar)
      if (any(left))
        warning("Steady-state equations still reference each other after ",
                "resolveRecurrence(): ", paste(names(m_ssChar)[left], collapse = ", "),
                ". The returned vector cannot be used as a parameter ",
                "transformation as is.")
    }
  }

  # Write steady states to disk
  if (!is.null(file) && is.character(file))
    saveRDS(object = m_ssChar, file = file)

  return(m_ssChar)
}


# Model names that sympy's sympify() turns into its own objects (functions, constants,
# the S registry): each gets an alias that is a plain symbol and clashes with no other
# name. Named character vector original -> alias; empty when nothing clashes.
.ssSympyAliases <- function(names) {
  names <- unique(names[nzchar(names)])
  if (!length(names)) return(character(0))
  spy <- tryCatch(reticulate::import("sympy", convert = FALSE), error = function(e) NULL)
  if (is.null(spy)) return(character(0))
  bad <- names[vapply(names, function(nm) {
    ok <- tryCatch(reticulate::py_to_r(spy$sympify(nm)$is_Symbol),
                   error = function(e) FALSE)
    !isTRUE(ok) || !identical(reticulate::py_to_r(spy$sympify(nm)$name), nm)
  }, logical(1))]
  if (!length(bad)) return(character(0))
  alias <- paste0(bad, "_dModSym")
  while (any(alias %in% names)) alias <- paste0(alias, "_")
  setNames(alias, bad)
}

# whole-name replacement old -> new in plain strings (also "a + b = tot")
.ssRename <- function(x, map) {
  if (!length(x) || !length(map)) return(x)
  for (nm in names(map))
    x <- gsub(paste0("(?<![A-Za-z0-9_.])", gsub(".", "\\.", nm, fixed = TRUE),
                     "(?![A-Za-z0-9_.])"), map[[nm]], x,
              perl = TRUE)
  x
}
