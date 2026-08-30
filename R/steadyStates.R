#' Calculate analytical steady states
#'
#' Computes symbolic steady-state expressions tailored to parameter estimation,
#' following the AlyssaPetit method \[1\], \[2\]. The solver itself is a Python
#' module (Python 3.x), called through reticulate.
#'
#' Rather than solving every balance equation for its own state, the solver
#' spends each equation on whichever state or rate constant keeps the result a
#' ratio of sums of positive terms. A state whose equation went to a rate
#' constant is absent from the result and stays a free parameter of the
#' transformation.
#'
#' @param model An `eqnlist`, or the name of a csv file describing the model.
#' @param file Character, path the result is written to with `saveRDS()`, and
#'   base name of the intermediate csv handed to the backend. For an `eqnlist`
#'   model it defaults to `"reactions_for_Alyssa"`.
#' @param rates Unused, retained for backward compatibility.
#' @param forcings Character vector, names of the forcings. These states are
#'   held at zero and treated as exogenous.
#' @param givenCQs Unnamed character vector of conserved quantities, either as
#'   `c("A + pA = totA", "B + pB = totB")` or as `c("A + pA", "B + pB")`. `NULL`
#'   (default) derives a basis automatically.
#' @param neglect Character vector, states and rate parameters the solver must
#'   not resolve. A neglected state stays a free parameter of the transformation;
#'   a neglected rate parameter is never used as a pivot.
#' @param sparsifyLevel Numeric, upper bound on the length of the linear
#'   combinations used to simplify the stoichiometric matrix. Versions `"1.0"`
#'   and `"1.1"` only.
#' @param outputFormat Character, `"R"` (default) for dMod-compatible output, or
#'   `"M"` for d2d \[3\].
#' @param testSteady Character, how the solution is verified. One of `"fast"`
#'   (default; probabilistic Schwartz-Zippel check over GF(p), with negligible
#'   error probability), `"exact"` (symbolic substitution, slow on large `sqrt`
#'   solutions) or `"skip"`. `"fast"` requires version `"1.2"` or later and
#'   falls back to `"exact"` on `"1.1"`; version `"1.0"` always tests.
#' @param walltime Integer, wall-clock budget in seconds for the solver, `0`
#'   (default) for unlimited. Version `"1.2"` and later.
#' @param simplify Final-simplification mode. `TRUE` (default) applies
#'   `sympy.simplify` once per expression, `FALSE` skips it and returns bulkier
#'   output, `"full"` adds a `cancel`/`posify`/`factor` pipeline that is slower
#'   but more compact. Version `"1.2"` and later.
#' @param solveQuadratic Logical, whether a cycle whose final equation is
#'   quadratic in its own state may be closed by the positive root of
#'   `a*X^2 + b*X + c = 0` instead of by a rate-parameter pivot. This keeps the
#'   pivoted rate constants out of the result, at the price of `sqrt(...)` terms
#'   that some workflows cannot consume in a parameter transformation. Default
#'   `FALSE`. Version `"1.2"` and later.
#' @param positive Positivity assumption used for root and pivot selection.
#'   `TRUE` (default) treats all parameters, initial values and totals as
#'   positive, `FALSE` assumes nothing, and a character vector names the symbols
#'   to treat as positive. Version `"1.2"` and later.
#'
#'   With `TRUE` the result is manifestly non-negative: no expression contains a
#'   subtraction, so no parameter choice can drive a state negative. Where
#'   solving for a state would yield a difference, its equation is spent on one
#'   of its own rate constants instead. That solve is linear and therefore has a
#'   unique root, which is why a pivot, not root selection, is the remedy. If no
#'   pivot succeeds either, a diagnosis is printed and `0` returned.
#' @param branches Logical, whether a quadratic with two positive roots
#'   (bistability) is emitted with a selector symbol `branch_<state>` taking
#'   `-1` or `+1` to recover both steady states. Requires
#'   `solveQuadratic = TRUE`. `FALSE` (default) emits only the provably unique
#'   positive root and pivots ambiguous cases. Version `"1.2"` and later.
#' @param priority Character vector of state and rate-parameter names,
#'   most-preferred first, biasing the resolution order: a named state is
#'   resolved earlier, a named rate parameter is preferred as pivot. Structural
#'   constraints take precedence, so this is a preference, not a guarantee.
#'   Unmatched names are reported and ignored. Version `"1.2"` and later.
#' @param resolve Logical, whether the equations are passed through
#'   [resolveRecurrence()] so that none refers to another, as a dMod parameter
#'   transformation substitutes all entries at once. `FALSE` keeps the compact
#'   recurrent form. Versions `"1.2"` and later already resolve in the backend,
#'   so this matters for `"1.0"` and `"1.1"`.
#' @param version Character, backend version. One of `"1.0"` (original), `"1.1"`
#'   (adds `testSteady`), `"1.2"` (sink-cluster detection, `walltime`,
#'   priority-table cycle breaking, `simplify` toggle, optional quadratic
#'   state-side solve), or `"1.3"` (default; same interface as `"1.2"`, but
#'   solutions are recorded lazily and resolved once at output time, with a
#'   lock guard replacing most rollbacks -- typically orders of magnitude
#'   faster on feedback-heavy networks).
#'
#' @return Named character vector of steady-state equations in dMod format, or
#'   `0` if no solution was found. An entry whose value is its own name denotes a
#'   free parameter.
#'
#' @references \[1\] <https://pmc.ncbi.nlm.nih.gov/articles/PMC4863410/>
#' @references \[2\] <https://github.com/marcusrosenblatt/AlyssaPetit>
#' @references \[3\] <https://github.com/Data2Dynamics/d2d>
#'
#' @author Marcus Rosenblatt, \email{marcus.rosenblatt@@fdm.uni-freiburg.de}
#'
#' @export
#' @importFrom utils write.table
#' @example inst/examples/steadystates.R
steadyStates <- function(model, file = NULL, rates = NULL, forcings = NULL,
                         givenCQs = NULL, neglect = NULL, sparsifyLevel = NULL,
                         outputFormat = "R", testSteady = c("fast", "exact", "skip"),
                         walltime = 0L, simplify = TRUE, solveQuadratic = FALSE,
                         positive = TRUE, branches = FALSE, priority = NULL,
                         version = "1.3", resolve = TRUE) {

  .require_ns("reticulate", "steadyStates()")
  # Validate version and verification mode
  version <- match.arg(version, choices = c("1.0", "1.1", "1.2", "1.3"))
  testSteady <- match.arg(testSteady)

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
    sparsifyLevel <- if (version %in% c("1.2", "1.3")) 0 else 2
  }

  # Check if model is an equation list
  if (inherits(model, "eqnlist")) {
    if (is.null(file)) file <- "reactions_for_Alyssa"
    # Not write.eqnlist(): the backend never sees the volumes, so the
    # V_ref / V_X factors getFluxes() applies have to be folded in first.
    utils::write.csv(.volumeScaledReactions(model),
                     file = paste0(file, "_model.csv"),
                     row.names = FALSE, na = "")
    model <- paste0(file, "_model.csv")
  }
  if (!is.null(givenCQs) && length(names(givenCQs)) > 0)
    stop("givenCQs must not have names. Please unname() them.")

  # Ensure Python dependencies are available
  reticulate::py_require("numpy")
  reticulate::py_require("sympy")
  # v1.2+ uses scipy.optimize.linprog for structural sink-cluster detection
  # (states whose combined mass leaks monotonically and must therefore be 0).
  # Older versions don't need scipy.
  if (version %in% c("1.2", "1.3")) reticulate::py_require("scipy")

  pymodule <- paste0("AlyssaPetit_ver", gsub("\\.", "_", version))
  ap <- reticulate::import_from_path(pymodule,
                                     path = system.file("code", package = "dMod2"))

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
    m_ss <- ap$Alyssa(model, as.list(forcings), as.list(givenCQs),
                      as.list(neglect), sparsifyLevel, outputFormat)

  } else if (version == "1.1") {
    if (isTRUE(solveQuadratic))
      message("Note: version 1.1 does not support solveQuadratic=TRUE, ignored.")
    if (length(priority) > 0)
      message("Note: version 1.1 does not support 'priority', ignored.")
    if (testSteady == "fast")
      message("Note: version 1.1 has no 'fast' test, using 'exact' instead.")
    # v1.1 backend takes the legacy "T"/"F" tokens.
    m_ss <- ap$Alyssa(model, as.list(forcings), as.list(givenCQs),
                      as.list(neglect), sparsifyLevel, outputFormat,
                      if (testSteady == "skip") "F" else "T")

  } else {
    # v1.2 / v1.3 (shared signature)
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
    m_ss <- ap$Alyssa(model,
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

  # Versions 1.2+ resolve on the backend side; this covers the older ones. Only
  # rewrite when there is something to resolve -- resolveRecurrence() reformats
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
