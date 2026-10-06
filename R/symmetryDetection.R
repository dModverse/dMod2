#' Structural non-identifiabilities of an ODE model
#'
#' @description Finds the directions in parameters and initial values along which
#'   the observations of a model do not change: the nullspace of the observability
#'   matrix over finite fields, with certified rank and exact closed forms. The model
#'   `f`, the observables `g` and an optional `trafo` are equations in any form
#'   [as.eqnvec()] accepts. The computation runs in the Python package `symident`
#'   through `reticulate`, installed on first use from
#'   <https://github.com/dModverse/symident>.
#'
#'   Right-hand sides, observables and initial values may contain `exp()`,
#'   `exp10()`, `b^x`, hyperbolic and trigonometric functions, free exponents `x^n`,
#'   and `log()`, fractional powers and `abs()` of `positive` coordinates. `time` is
#'   the model time.
#'
#' @param f Right-hand sides: an [eqnlist], an [eqnvec] or a named character vector.
#'   `reduceCQ` and `equilibrate` use the reactions of an [eqnlist].
#' @param g Observables: one [eqnvec] or named character vector for all conditions,
#'   or a list with one per row of `conditions`. An observable `a*log(h) + c` with a
#'   number `a` is analysed through `h`.
#' @param trafo Parameter transformation: one [eqnvec] for all conditions or a list
#'   with one per condition. An entry named like a parameter is substituted into `f`,
#'   `g`, event values and event times; an entry named like a state is its initial
#'   value. On the right-hand side a state name is the initial value of that state.
#'   A parameter that enters only through `k = "exp(lk)"`, `"exp10(lk)"` or
#'   `"10^lk"` is analysed in `b^lk` and reported in `lk`.
#' @param parameters Character vector of additional symbols treated as parameters.
#' @param fixed Character vector of known symbols: a fixed parameter is a known
#'   constant, a fixed state has no unknown initial value.
#' @param gaugePreference `FALSE` (default) reconstructs every direction. Otherwise
#'   one coordinate per scaling is fixed first and only the general directions are
#'   reconstructed in that gauge. `NULL` chooses the gauge that leaves the general
#'   directions on the fewest coordinates; a character vector ranks the coordinates
#'   to fix (`*` as wildcard). The gauge is returned as `$gauge`.
#' @param forcings Character vector of input states, 0 at the start and at rest.
#' @param events An [eventlist] with event times as numbers or expressions in the
#'   parameters, in chronological order. Root events are not supported. An event
#'   value named like a column of `conditions` is read from it.
#' @param conditions Data frame with one row per condition and columns named by
#'   model symbols or event values. A numeric cell fixes the symbol, a character cell
#'   renames it. A column named like a dynamic state sets its initial value. A
#'   direction is reported only if it is non-identifiable in every condition.
#' @param equilibrate Logical. Start at a steady state of `f` with the inputs at 0;
#'   the earliest events apply on top. Coupled steady states may need msolve, which
#'   `symident.install_msolve()` builds or `SYMIDENT_MSOLVE` points to.
#' @param reduceCQ Logical, [eqnlist] only. Report the freedom of a conserved moiety
#'   on its total (`TRUE`) or on the initial value of one species (`FALSE`).
#' @param freeInitial States, at most one per conserved quantity, carrying the free
#'   resting value of their moiety under `equilibrate = TRUE, reduceCQ = FALSE`.
#' @param reconstruct Logical. Return general directions as exact rational functions
#'   instead of their support.
#' @param positive Coordinates known to be positive: `TRUE` (default, all), `FALSE`
#'   (none) or a character vector.
#' @param verify Logical. Check that the rank has saturated where the Lie order is
#'   not certified; the result is in `$info$verification`.
#' @param cores Number of threads and worker processes.
#' @param control A [reconstControl()] list.
#' @param scalingsOnly Logical. Only the scaling symmetries, from an exact integer
#'   kernel; `identifiable`, `rank` and `dim` are `NA`.
#' @param symEngine `"modular"` (default, finite fields) or `"symbolic"` (sympy, small
#'   models without `equilibrate` and later events).
#' @param verbose Logical. Print the result.
#'
#' @return An object of class `symmetrydetection`:
#'   \describe{
#'     \item{`identifiable`}{`TRUE` or `FALSE`, `NA` with `scalingsOnly`.}
#'     \item{`rank`, `dim`}{rank of the observability matrix and number of
#'       coordinates.}
#'     \item{`symmetries`}{the directions, each a generator
#'       \eqn{X = \sum_i \eta_i \partial_{z_i}}: `generator` (\eqn{\eta} by
#'       coordinate), `weights` (a scaling's integer weights), `type` (`"scaling"` or
#'       `"general"`), `degree`, `support`, `explicit`, `certified`, `reason`,
#'       `route`, `display`, `completeGenerator` and `factor`.}
#'     \item{`info`}{`engine`, the Lie order and its certification
#'       (`lieOrderUsed`, `lieOrderDriver`, `liePlateau`, `lieCertified`),
#'       `rankProven`, `rankProof`, `gapOrderUsed`, `conditions`, `segments`,
#'       `coordinates`, `settings`, `elapsed` and `verification`.}
#'     \item{`gauge`}{the coordinates fixed under `gaugePreference`.}
#'     \item{`call`}{the matched call.}
#'   }
#'   `print()` shows the verdict and the generators, `summary()` also the
#'   computation. Both take `verbose` and `width`.
#'
#' @details The analysis starts at the earliest event; later events split the time
#'   line into segments. Exponentials and trigonometric functions of states enter as
#'   auxiliary states; directions are reported in the original functions. Switches
#'   of symident are read from the environment variables `SYMIDENT_*` and from the
#'   options `dMod.sym.*`.
#'
#' @example inst/examples/symmetryDetection.R
#' @export
symmetryDetection <- function(f = NULL, g = NULL, trafo = NULL, parameters = NULL,
                              fixed = NULL, gaugePreference = FALSE, forcings = NULL,
                              events = NULL, conditions = NULL, equilibrate = FALSE,
                              reduceCQ = FALSE, freeInitial = NULL, reconstruct = FALSE,
                              positive = TRUE, verify = TRUE, cores = 1,
                              control = reconstControl(), scalingsOnly = FALSE,
                              symEngine = c("modular", "symbolic"), verbose = TRUE) {
  symEngine <- match.arg(symEngine)
  if (is.null(f))
    stop("Provide the model right-hand sides via `f` ",
         "(eqnlist, eqnvec or named character vector).", call. = FALSE)
  spec <- .symSpec(f, g, trafo, conditions, events)
  syms <- list(parameters = parameters, fixed = fixed, forcings = forcings,
               free_initial = freeInitial)
  for (k in names(syms))
    if (!is.null(syms[[k]])) spec[[k]] <- as.list(as.character(syms[[k]]))
  spec$gauge_preference <- if (is.null(gaugePreference)) "AUTO"
    else if (isFALSE(gaugePreference)) FALSE else as.list(as.character(gaugePreference))
  spec$equilibrate <- isTRUE(equilibrate)
  spec$reduce_cq <- isTRUE(reduceCQ)
  spec$reconstruct <- isTRUE(reconstruct)
  spec$verify <- isTRUE(verify)
  spec$scalings_only <- isTRUE(scalingsOnly)
  spec$positive <- if (is.logical(positive)) isTRUE(positive) else as.list(as.character(positive))
  spec$sym_engine <- symEngine
  spec$cores <- as.integer(max(1L, cores))
  spec$control <- setNames(lapply(unclass(control), function(x)
    if (is.numeric(x) && is.infinite(x)) "Inf" else x), .symSnake(names(control)))
  res <- .symFromPy(.symCall("detect_json", spec), match.call())
  if (isTRUE(verbose)) print(res)
  invisible(res)
}

#' Settings of the observability engine
#'
#' Saturation and closed-form reconstruction settings for
#' `symmetryDetection(control = reconstControl())`. Raise the caps to reconstruct
#' wide or high-degree directions, at the cost of more samples.
#'
#' @param relevanceCap Maximum number of coordinates in one entry of a direction
#'   for the dense fit; wider entries use sparse interpolation.
#' @param relevanceCapDir Maximum number of coordinates in one direction; a wider
#'   one is reported by its support.
#' @param relevanceCapSparse Maximum number of coordinates in one entry for the
#'   sparse (Ben-Or/Tiwari) fit; a wider entry is reported by its support.
#' @param degreeCap Total degree bound of the dense rational fit.
#' @param sampleSlack Number of samples beyond the minimum of the fit.
#' @param probeRetries Number of retries of the relevance probe when a
#'   perturbation changes the pivots.
#' @param laurentDegNum,laurentDegDen Numerator degree and monomial denominator
#'   degree bounds of the sparse Laurent fit.
#' @param laurentCandCap Maximum number of candidate monomials of the Laurent and
#'   general sparse fits.
#' @param termCap Maximum number of terms of a sparse entry.
#' @param generalDegNum,generalDegDen Numerator and denominator degree bounds of the
#'   general sparse rational fit.
#' @param gapOrderCap Maximum order of the power series in the time between
#'   events.
#' @param minsupportCandCap Maximum number of column subsets searched for
#'   directions with small support.
#' @param perprimeCap Maximum number of samples per prime for the reconstruction
#'   under `equilibrate = TRUE`.
#' @param perprimeMinPrimes Minimum number of primes with samples for a
#'   reconstruction under `equilibrate = TRUE`.
#' @param timeout Time limit in seconds for the reconstruction. Directions not
#'   finished in time are reported by their support. `Inf` (default) sets no limit.
#' @return A `reconstControl` list.
#' @seealso [symmetryDetection()]
#' @export
reconstControl <- function(relevanceCap = 6L, relevanceCapDir = 24L, relevanceCapSparse = 30L,
                           degreeCap = 4L, sampleSlack = 5L, probeRetries = 8L,
                           laurentDegNum = 4L, laurentDegDen = 2L, laurentCandCap = 200000L,
                           termCap = 60L, generalDegNum = 4L, generalDegDen = 3L,
                           gapOrderCap = 8L, minsupportCandCap = 20000L, perprimeCap = 120L,
                           perprimeMinPrimes = 3L, timeout = Inf) {
  stopifnot(relevanceCap >= 0L, relevanceCapDir >= 1L,
            relevanceCapSparse >= relevanceCap, degreeCap >= 0L,
            sampleSlack >= 0L, probeRetries >= 1L, termCap >= 1L,
            laurentDegNum >= 1L, laurentDegDen >= 0L, laurentCandCap >= 1L,
            generalDegNum >= 1L, generalDegDen >= 1L, gapOrderCap >= 0L,
            minsupportCandCap >= 1L, perprimeCap >= 1L, perprimeMinPrimes >= 2L,
            is.numeric(timeout), timeout > 0)
  structure(list(relevanceCap = as.integer(relevanceCap),
                 relevanceCapDir = as.integer(relevanceCapDir),
                 relevanceCapSparse = as.integer(relevanceCapSparse),
                 degreeCap = as.integer(degreeCap), sampleSlack = as.integer(sampleSlack),
                 probeRetries = as.integer(probeRetries),
                 laurentDegNum = as.integer(laurentDegNum),
                 laurentDegDen = as.integer(laurentDegDen),
                 laurentCandCap = as.integer(laurentCandCap), termCap = as.integer(termCap),
                 generalDegNum = as.integer(generalDegNum),
                 generalDegDen = as.integer(generalDegDen),
                 gapOrderCap = as.integer(gapOrderCap),
                 minsupportCandCap = as.integer(minsupportCandCap),
                 perprimeCap = as.integer(perprimeCap),
                 perprimeMinPrimes = as.integer(perprimeMinPrimes), timeout = timeout),
            class = c("reconstcontrol", "list"))
}


# ---- the interface to the Python package symident --------------------------------

# Where symident is installed from until it is on PyPI
.symidentSource <- "git+https://github.com/dModverse/symident"

# A module of symident. A missing symident is requested into reticulate's
# ephemeral environment; under RETICULATE_PYTHON it has to be installed there.
.symident <- function(module = "rjson") {
  .require_ns("reticulate", "symmetryDetection()")
  .require_ns("jsonlite", "symmetryDetection()")
  if (!reticulate::py_module_available("symident"))
    suppressWarnings(try(reticulate::py_require(paste("symident @", .symidentSource)),
                         silent = TRUE))
  tryCatch(reticulate::import(paste0("symident.", module)),
           error = function(e)
    stop("symmetryDetection() needs the Python package symident. Install it into ",
         "the Python environment of reticulate with `pip install ", .symidentSource,
         "` (", conditionMessage(e), ").", call. = FALSE))
}

# camelCase to snake_case, the argument names of symident
.symSnake <- function(x) tolower(gsub("([a-z0-9])([A-Z])", "\\1_\\2", x))

# snake_case to camelCase, the field names of the R objects
.symCamel <- function(x) gsub("_([a-z0-9])", "\\U\\1", x, perl = TRUE)

# field names to camelCase, recursively; the fields in `data` are named by model
# symbols and keep their names
.symCamelFields <- function(x, data = c("generator", "complete_generator", "weights",
                                        "display", "factor", "pins", "survivor_meaning",
                                        "carrier_domain", "trafo", "totals")) {
  if (!is.list(x)) return(x)
  nm <- names(x)
  for (i in seq_along(x))
    if (is.null(nm) || !nm[i] %in% data) x[i] <- list(.symCamelFields(x[[i]], data))
  if (!is.null(nm)) names(x) <- .symCamel(nm)
  x
}

# one call into symident: the arguments as JSON in, the result as JSON out, with its
# warnings raised and the JSON kept as attribute "json". The SYMIDENT_* environment
# and the dMod.sym.* options travel with it, since Python reads its environment once
.symCall <- function(fun, ...) {
  si <- .symident()
  env <- Sys.getenv()
  env <- as.list(env[grepl("^SYMIDENT_", names(env)) & nzchar(env)])
  opt <- options()
  opt <- opt[grepl("^dMod\\.sym\\.", names(opt))]
  names(opt) <- .symSnake(sub("^dMod\\.sym\\.", "", names(opt)))
  args <- lapply(list(...), function(a) if (is.character(a) && length(a) == 1L) a
                 else as.character(jsonlite::toJSON(a, auto_unbox = TRUE, null = "null",
                                                    digits = NA, force = TRUE)))
  out <- do.call(si[[fun]], c(args, list(
    env = as.character(jsonlite::toJSON(env, auto_unbox = TRUE)),
    options = as.character(jsonlite::toJSON(opt, auto_unbox = TRUE, null = "null")))))
  res <- jsonlite::fromJSON(out, simplifyVector = FALSE)
  for (w in unlist(res$warnings)) warning(w, call. = FALSE)
  res$warnings <- NULL
  structure(res, json = out)
}

# the result as symident wrote it, for its own report and for the reduction
.symRaw <- function(x) {
  raw <- attr(x, "symident")
  if (is.null(raw))
    stop("the object does not carry the symident result it was built from.", call. = FALSE)
  raw
}

# the report lines of symident for a result, written to the console
.symReport <- function(fun, x, ...) {
  json <- reticulate::import("json", convert = FALSE)$loads(.symRaw(x))
  writeLines(.symident("report")[[fun]](json, ...))
}

# the model as the JSON specification symident reads
.symSpec <- function(f, g, trafo, conditions, events) {
  eqn <- function(x) if (is.null(x)) NULL else
    as.list(setNames(as.character(x), names(x)))
  out <- list()
  if (inherits(f, "eqnlist")) {
    S <- f$smatrix; S[is.na(S)] <- 0
    out$reactions <- list(species = as.list(colnames(f$smatrix)),
                          smatrix = lapply(seq_len(nrow(S)), function(i) I(as.numeric(S[i, ]))),
                          rates = as.list(as.character(f$rates)))
    tt <- getTotals(f)
    if (length(tt)) out$totals <- lapply(tt, as.character)
  }
  out$f <- eqn(as.eqnvec(f))
  out$g <- if (is.list(g) && !inherits(g, "eqnvec"))
    lapply(unname(g), function(x) eqn(as.eqnvec(x)))
    else if (!is.null(g)) eqn(as.eqnvec(g))
  if (!is.null(trafo))
    out$trafo <- if (is.list(trafo) && !inherits(trafo, "eqnvec"))
      lapply(unname(trafo), function(x) if (is.null(x)) NULL else eqn(as.eqnvec(x)))
      else eqn(as.eqnvec(trafo))
  if (!is.null(conditions)) {
    cd <- as.data.frame(conditions, stringsAsFactors = FALSE)
    out$conditions <- list(rows = as.list(rownames(cd)), cols = lapply(as.list(cd), function(v)
      as.list(if (is.factor(v)) as.character(v) else v)))
  }
  if (!is.null(events) && nrow(as.data.frame(events))) {
    ed <- as.data.frame(events, stringsAsFactors = FALSE)
    out$events <- lapply(seq_len(nrow(ed)), function(i) list(
      var = as.character(ed$var[i]), time = as.character(ed$time[i]),
      value = as.character(ed$value[i]), method = as.character(ed$method[i]),
      root = if (is.null(ed$root) || is.na(ed$root[i])) NULL else as.character(ed$root[i])))
  }
  out
}

# the result of symident as a symmetrydetection object, the raw result attached
.symFromPy <- function(res, call) {
  raw <- attr(res, "json")
  res <- .symCamelFields(unclass(res))
  attr(res, "json") <- NULL
  chr <- function(x) if (is.null(x)) NULL else as.character(unlist(x))
  vec <- function(x) if (is.null(x)) NULL else
    as.eqnvec(setNames(as.character(unlist(x)), names(x)))
  res$symmetries <- lapply(res$symmetries, function(d) {
    d$generator <- vec(d$generator)
    d$completeGenerator <- vec(d$completeGenerator)
    d$support <- chr(d$support)
    if (!is.null(d$display)) d$display <- lapply(d$display, as.character)
    structure(d, class = "symmetrygenerator")
  })
  info <- res$info
  info$kernel <- NULL
  if (!is.null(info$lieOrderDriver)) info$lieOrderDriver <- as.integer(info$lieOrderDriver) + 1L
  if (length(info$lieUncertified)) info$lieUncertified <- as.integer(unlist(info$lieUncertified)) + 1L
  info$lieBlockOrders <- if (length(info$lieBlockOrders)) as.integer(unlist(info$lieBlockOrders))
  info$coordinates <- chr(info$coordinates)
  info$modelExprs <- chr(info$modelExprs)
  if (!is.null(info$gaugeSuggestion)) {
    gs <- info$gaugeSuggestion
    info$gaugeSuggestion <- list(gauge = chr(gs$gauge), sizes = as.integer(unlist(gs$sizes)),
                                 plain = chr(gs$plain), plainSizes = as.integer(unlist(gs$plainSizes)))
  }
  if (!is.null(info$settings$positive) && is.list(info$settings$positive))
    info$settings$positive <- chr(info$settings$positive)
  res$info <- info
  if (!is.null(res$gauge)) res$gauge <- chr(res$gauge)
  for (k in c("rank", "dim")) res[[k]] <- if (is.null(res[[k]])) NA_integer_ else as.integer(res[[k]])
  if (is.null(res$identifiable)) res$identifiable <- NA
  res$call <- call
  structure(res, class = "symmetrydetection", symident = raw)
}

#' @export
print.symmetrydetection <- function(x, verbose = FALSE, width = getOption("width"), ...) {
  .symReport("result_lines", x, verbose = isTRUE(verbose), width = as.integer(width))
  invisible(x)
}

#' @export
summary.symmetrydetection <- function(object, verbose = FALSE, width = getOption("width"), ...)
  structure(list(object = object, verbose = isTRUE(verbose), width = as.integer(width)),
            class = "summary.symmetrydetection")

#' @export
print.summary.symmetrydetection <- function(x, ...) {
  .symReport("summary_lines", x$object, verbose = x$verbose, width = x$width)
  invisible(x)
}
