#' Analytical Steady States
#'
#' Computes symbolic steady-state expressions for use in a parameter
#' transformation, following the AlyssaPetit method \[1\], \[2\]. Each balance
#' equation is solved for a state or a rate constant such that the result is
#' a ratio of sums of positive terms. A state whose equation was spent on a
#' rate constant stays a free parameter. Every unknown is solved from a
#' combination of balances, so every result is a ratio of positive sums, and
#' compartment volumes are fixed. The solver is a Python module called through
#' \pkg{reticulate}; it needs \pkg{numpy}, \pkg{sympy} and \pkg{scipy}.
#'
#' @param model An [eqnlist], or the path of a csv file describing the model.
#' @param file Character, path of the RDS file the result is written to with
#'   `saveRDS()`. For an `eqnlist` the model is also written to
#'   `<file>_model.csv` for the solver. Default `NULL` writes nothing; the
#'   solver then reads the model from a temporary file.
#' @param forcings Character vector, names of forcings. These states are held
#'   at zero and treated as exogenous. Default `NULL`.
#' @param givenCQs Unnamed character vector of conserved quantities, either as
#'   `c("A + pA = totA", "B + pB = totB")` or as `c("A + pA", "B + pB")`.
#'   Default `NULL`: a basis is derived, or taken from [customTotals()] when the
#'   `eqnlist` has them. Each given quantity keeps one of its states free, the
#'   first where possible.
#' @param neglect Character vector, states and rate parameters the solver must
#'   not solve for. Default `NULL`.
#' @param testSteady How the solution is verified: `"fast"` (default), a
#'   probabilistic check modulo a prime, or `"skip"`.
#' @param walltime Integer, time budget of the solver in seconds, `0`
#'   (default) for none.
#' @param simplify `TRUE` (default) simplifies every expression once with
#'   sympy, `FALSE` skips it, `"full"` simplifies further at a higher cost.
#' @param solveQuadratic Logical, whether a balance quadratic in its own state
#'   may be solved by its positive root instead of for a rate constant, which
#'   introduces `sqrt()` terms. Default `FALSE`.
#' @param priority Character vector of state and rate-parameter names, most
#'   preferred first: named states are solved earlier, named rate parameters
#'   are preferred as pivots. Default `NULL`.
#' @param resolve Logical, whether the result is passed through
#'   [resolveRecurrence()] so that no equation refers to another. Default
#'   `TRUE`.
#' @param verbose `TRUE` reports progress in a few lines and the result,
#'   `FALSE` (default) only the result, `"full"` every step of the solver.
#'
#' @return Named character vector of steady-state equations in dMod format, or
#'   `0` if no solution was found. An entry whose value is its own name is a
#'   free parameter.
#'
#' @references \[1\] Rosenblatt, M., Timmer, J. and Kaschek, D. (2016).
#'   Customized steady-state constraints for parameter estimation in
#'   non-linear ordinary differential equation models.
#'   \emph{Frontiers in Cell and Developmental Biology} 4, 41.
#'   \doi{10.3389/fcell.2016.00041}
#' @references \[2\] <https://github.com/marcusrosenblatt/AlyssaPetit>
#'
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
#' steadies <- steadyStates(reactions)
#' steadies
#'
#' # Parameter transformation that puts the model in steady state
#' parameters <- getParameters(reactions)
#' trafo <- repar(eqnvec(setNames(parameters, parameters)), "x ~ y",
#'                x = names(steadies), y = steadies)
#' trafo
steadyStates <- function(model, file = NULL, forcings = NULL,
                         givenCQs = NULL, neglect = NULL, testSteady = c("fast", "skip"),
                         walltime = 0L, simplify = TRUE, solveQuadratic = FALSE,
                         priority = NULL, resolve = TRUE, verbose = FALSE) {

  .require_ns("reticulate", "steadyStates()")
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

  # Check if model is an equation list
  if (inherits(model, "eqnlist")) {
    # the solver would take a volume factor for a rate constant: volumes are not unknowns
    if (length(model$volumes))
      neglect <- union(neglect, getSymbols(unique(unlist(model$volumes))))
    csvBase <- if (is.null(file)) tempfile("steadyStates") else file
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
                                         getSymbols(givenCQs), getSymbols(volumes))))
    if (length(symAlias)) {
      ren <- function(x) .ssRename(x, symAlias)
      tab$Rate <- ren(tab$Rate)
      hit <- names(tab) %in% names(symAlias)
      names(tab)[hit] <- symAlias[names(tab)[hit]]
      volumes <- setNames(ren(volumes), ren(names(volumes)))
      forcings <- ren(forcings); neglect <- ren(neglect); priority <- ren(priority)
      if (length(givenCQs)) givenCQs <- ren(givenCQs)
    }
    model <- paste0(csvBase, "_model.csv")
    utils::write.csv(tab, file = model, row.names = FALSE, na = "")
    if (is.null(file)) on.exit(unlink(model), add = TRUE)
  } else symAlias <- volumes <- character(0)
  if (length(names(givenCQs)))
    stop("givenCQs are equations \"lhs = total\" without names, drop the names with unname().")

  reticulate::py_require(c("numpy", "sympy", "scipy"))
  ap <- reticulate::import_from_path("AlyssaPetit",
                                     path = system.file("code", package = "dMod2"))
  # simplify is TRUE, FALSE or "full", passed on as a Python bool or the string
  if (is.character(simplify)) simplify <- match.arg(tolower(simplify), choices = "full")
  else simplify <- as.logical(simplify)
  m_ss <- ap$Alyssa(model,
                    injections     = as.list(forcings),
                    givenCQs       = as.list(givenCQs),
                    neglect        = as.list(neglect),
                    testSteady     = testSteady,
                    walltime       = as.integer(walltime),
                    simplify       = simplify,
                    solveQuadratic = as.logical(solveQuadratic),
                    priority       = as.list(as.character(priority)),
                    verbose        = verbose,
                    volumes        = as.list(volumes))

  if (is.null(m_ss) || identical(m_ss, 0L)) return(0)

  # The solver answers with "lhs=rhs" strings, the result is rhs named by lhs.
  sides <- strsplit(as.character(unlist(m_ss)), "=", fixed = TRUE)
  steady <- setNames(trimws(vapply(sides, `[`, "", 2L)), trimws(vapply(sides, `[`, "", 1L)))

  if (length(steady) == 0) return(0)
  if (length(symAlias)) {
    back <- setNames(names(symAlias), symAlias)
    steady <- setNames(.ssRename(steady, back),
                       ifelse(names(steady) %in% names(back), back[names(steady)],
                              names(steady)))
  }

  # The solver substitutes its results into each other; this catches what is left.
  # Only rewrite when there is something to resolve, resolveRecurrence() reformats
  # every equation it touches, which would break printed == returned.
  if (resolve && length(steady) > 1) {
    # Entries kept as their own name are free parameters, not definitions.
    passthrough <- names(steady)[trimws(steady) == names(steady)]
    recurrent <- function(x) vapply(seq_along(x), function(i) {
      any(setdiff(getSymbols(x[i]), passthrough) %in% names(x)[-i])
    }, logical(1))
    if (any(recurrent(steady))) {
      steady <- resolveRecurrence(steady)
      left <- recurrent(steady)
      if (any(left))
        warning("Steady-state equations still reference each other after ",
                "resolveRecurrence(): ", paste(names(steady)[left], collapse = ", "),
                ". The returned vector cannot be used as a parameter ",
                "transformation as is.")
    }
  }

  # Write steady states to disk
  if (!is.null(file) && is.character(file))
    saveRDS(object = steady, file = file)

  steady
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
