#' Structural Non-Identifiabilities of an ODE Model
#'
#' @description Finds the directions in parameters and initial values along which
#'   the observations of a model do not change: the nullspace of the observability
#'   matrix over finite fields, with exact closed forms and a level that says whether
#'   each one is proven. The model `f`, the observables `g` and an optional `trafo` are
#'   equations in any form [as.eqnvec()] accepts. The computation runs in the Python
#'   package `symident` through `reticulate`, installed from PyPI on first use.
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
#' @param regimes For a model with switches, the regimes an experiment visits: a
#'   character vector of comparisons such as `"x > c"`, each selecting the regimes on
#'   that side of a switch, or a list of such vectors whose comparisons hold together;
#'   the first entry is the regime at the start. A list named by the rows of
#'   `conditions` gives them per condition. Without it every regime counts; a switch on
#'   parameters alone needs its side named here.
#' @param equilibrate Logical. Start at a steady state of `f` with the inputs at 0;
#'   the earliest events apply on top. Coupled steady states may need msolve; when
#'   it is missing, an interactive session asks whether to build it into the
#'   symident cache (on Windows with Rtools). Defaults to `FALSE`.
#' @param reduceCQ Logical, [eqnlist] only. Report the freedom of a conserved moiety
#'   on its total (`TRUE`) or on the initial value of one species (`FALSE`,
#'   default).
#' @param freeInitial States, at most one per conserved quantity, holding the free
#'   resting value of their moiety under `equilibrate = TRUE, reduceCQ = FALSE`.
#' @param reconstruct Logical. Return general directions as exact rational functions
#'   instead of their support. Defaults to `FALSE`.
#' @param positive Coordinates known to be positive: `TRUE` (default, all), `FALSE`
#'   (none) or a character vector.
#' @param verify Logical. Check that the rank has saturated where the Lie order is
#'   not proven sufficient; the result is in `$info$verification`. Defaults to `TRUE`.
#' @param cores Number of threads and worker processes. Defaults to 1.
#' @param control A [reconstControl()] list.
#' @param scalingsOnly Logical. Only the scaling symmetries, from an exact integer
#'   kernel; `identifiable`, `rank` and `dim` are `NA`. Defaults to `FALSE`.
#' @param symEngine `"modular"` (default, finite fields) or `"symbolic"` (sympy, small
#'   models without `equilibrate` and later events).
#' @param verbose Logical. Print the result. Defaults to `FALSE`.
#'
#' @return An object of class `symmetrydetection`:
#'   \describe{
#'     \item{`identifiable`}{`TRUE` or `FALSE`, `NA` with `scalingsOnly`; `FALSE` holds as far as
#'       `info$rankStatus` says.}
#'     \item{`rank`, `dim`}{rank of the observability matrix and number of
#'       coordinates.}
#'     \item{`symmetries`}{the directions, each a generator
#'       \eqn{X = \sum_i \xi_i \partial_{\theta_i}}: `generator` (\eqn{\xi} by
#'       coordinate), `weights` (a scaling's integer weights), `type` (`"scaling"` or
#'       `"general"`), `degree`, `support`, `explicit`, `level`, `stateComponents`,
#'       `reason`, `route`, `display`, `completeGenerator` and `factor`. `level` is
#'       `"exact"` for a direction proven over the rationals, `"verified"` for one
#'       checked at random points modulo a large prime at the Lie and gap orders used
#'       and `"open"` without closed form; `stateComponents` holds, for an exact
#'       direction, the components on the states of the field that proves it.}
#'     \item{`info`}{`engine`, the Lie order and whether it is proven sufficient
#'       (`lieOrderUsed`, `lieOrderDriver`, `liePlateau`, `lieCertified`),
#'       the gap order of the events after the start and whether it is proven
#'       sufficient (`gapOrderUsed`, `gapCertified`, `NULL` without such events, and
#'       `gapCertificate`, why not), `rankProven`, `rankProof`, `rankFields` (a rank
#'       the invariant fields establish, checked at random points modulo a prime),
#'       `rankFieldsBound` (the probability that a wrong field passes), `rankStatus`
#'       (`status` `"proven"`, `"verified"` or `"lower bound"`, `why` and
#'       `probability`), `proof` (why each `"verified"` direction is not proven),
#'       `bound` (the probability that a wrong closed form passes the check, only where
#'       both orders are proven sufficient or the rank is settled), `limits` (what
#'       stopped part of the analysis), `errors` (exceptions turned into failed
#'       routes), `branches` (resting states with several roots modulo a prime, of
#'       which the first is taken), `conditions`, `segments`, `coordinates`,
#'       `settings`, `elapsed` and `verification`.}
#'     \item{`gauge`}{the coordinates fixed under `gaugePreference`.}
#'     \item{`call`}{the matched call.}
#'   }
#'   `print()` shows the verdict and the generators, `summary()` also the
#'   computation. Both take `verbose` and `width`.
#'
#' @details The analysis starts at the earliest event; later events split the time
#'   line into segments. Directions are reported in the original coordinates.
#'
#'   Every environment variable `SYMIDENT_<NAME>` and every option
#'   `dMod.sym.<name>` set in R is passed to symident for the call, the option
#'   name converted to snake case. The environment variables a user may need
#'   are `SYMIDENT_MSOLVE`, the path of the msolve executable (unset by
#'   default, msolve is then searched in the symident cache and on the
#'   `PATH`), and `SYMIDENT_CACHE`, the directory of that cache (by default
#'   `~/.cache/symident` on Linux). The other switches are diagnostic, see the
#'   symident documentation at <https://dmodverse.github.io/symident/>.
#'
#' @seealso [symmetryReduction()], [reconstControl()]
#' @examplesIf requireNamespace("reticulate", quietly = TRUE) && reticulate::py_module_available("symident")
#' f <- eqnvec(m = "ktx - dm*m", p = "ktl*m - dp*p")
#' g <- eqnvec(y = "p")
#' res <- symmetryDetection(f, g, reconstruct = TRUE)
#' res
#' summary(res)
#' res$identifiable
#' @export
symmetryDetection <- function(f, g, trafo = NULL, parameters = NULL,
                              fixed = NULL, gaugePreference = FALSE, forcings = NULL,
                              events = NULL, conditions = NULL, regimes = NULL,
                              equilibrate = FALSE, reduceCQ = FALSE, freeInitial = NULL,
                              reconstruct = FALSE,
                              positive = TRUE, verify = TRUE, cores = 1,
                              control = reconstControl(), scalingsOnly = FALSE,
                              symEngine = c("modular", "symbolic"), verbose = FALSE) {
  symEngine <- match.arg(symEngine)
  if (missing(f) || is.null(f))
    stop("symmetryDetection: give the right-hand sides `f` ",
         "(eqnlist, eqnvec or named character vector).", call. = FALSE)
  if (missing(g))
    stop("symmetryDetection: give the observables `g` ",
         "(eqnvec, named character vector or a list of those).", call. = FALSE)
  spec <- .symSpec(f, g, trafo, conditions, events)
  syms <- list(parameters = parameters, fixed = fixed, forcings = forcings,
               free_initial = freeInitial)
  for (k in names(syms))
    if (!is.null(syms[[k]])) spec[[k]] <- as.list(as.character(syms[[k]]))
  if (!is.null(regimes)) spec$regimes <- .symRegimes(regimes)
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
    if (is.numeric(x) && is.infinite(x)) "Inf" else x), .symControlSnake(names(control)))
  res <- .symFromPy(.symCall("detect_json", spec), match.call())
  if (isTRUE(verbose)) print(res)
  invisible(res)
}

#' Settings of the Observability Engine
#'
#' Saturation and closed-form reconstruction settings for
#' `symmetryDetection(control = reconstControl())`. Raise the caps to reconstruct
#' wide or high-degree directions, at the cost of more samples.
#'
#' @param relevanceCap Maximum number of coordinates in one entry of a direction
#'   for the dense fit; wider entries use sparse interpolation. Defaults to 6.
#' @param relevanceCapDir Maximum number of coordinates in one direction; a wider
#'   one is reported by its support. Defaults to 24.
#' @param relevanceCapSparse Maximum number of coordinates in one entry for the
#'   sparse (Ben-Or/Tiwari) fit; a wider entry is reported by its support.
#'   Defaults to 30.
#' @param degreeCap Total degree bound of the dense rational fit. Defaults to 4.
#' @param sampleSlack Number of samples beyond the minimum of the fit. Defaults
#'   to 5.
#' @param probeRetries Number of retries of the relevance probe when a
#'   perturbation changes the pivots. Defaults to 8.
#' @param laurentDegNum,laurentDegDen Numerator degree and monomial denominator
#'   degree bounds of the sparse Laurent fit. Default to 4 and 2.
#' @param laurentCandCap Maximum number of candidate monomials of the Laurent and
#'   general sparse fits. Defaults to 200000.
#' @param termCap Maximum number of terms of a sparse entry. Defaults to 60.
#' @param generalDegNum,generalDegDen Numerator and denominator degree bounds of the
#'   general sparse rational fit. Default to 4 and 3.
#' @param gapOrderCap Maximum order of the power series in the time between
#'   events. Defaults to 8.
#' @param minSupportCandCap Maximum number of column subsets searched for
#'   directions with small support. Defaults to 20000.
#' @param perPrimeCap Maximum number of samples per prime for the reconstruction
#'   under `equilibrate = TRUE`. Defaults to 120.
#' @param perPrimeMinPrimes Minimum number of primes with samples for a
#'   reconstruction under `equilibrate = TRUE`. Defaults to 3.
#' @param homotopy Whether directions still open under `equilibrate = TRUE` are
#'   fitted along lines in parameter space as a last route. Defaults to `TRUE`.
#' @param timeout Time limit in seconds for the reconstruction. Directions not
#'   finished in time are reported by their support. `Inf` (default) sets no limit.
#' @param homotopyDegree Total degree bound of the entries the homotopy route fits.
#' @param degreeCapMax Degree up to which an entry past `degreeCap` is probed again,
#'   the bound doubled each time.
#' @param exactProof Whether a closed form that is not its own lift is proven through
#'   the fields of each regime, one detection per regime.
#' @param exactWorkCap Work bound of that proof: regimes times coordinates times tape
#'   instructions; above it the proof is skipped and the limit noted.
#' @param seed Seed of the random choices.
#' @details `homotopyDegree`, `degreeCapMax`, `exactProof`, `exactWorkCap` and `seed`
#'   default to `NULL`, which keeps the value of symident (4, 8, `TRUE`, 500000 and
#'   0).
#' @return A list of class `reconstcontrol` with the settings.
#' @seealso [symmetryDetection()]
#' @export
#' @examples
#' reconstControl(degreeCap = 6L, timeout = 60)
reconstControl <- function(relevanceCap = 6L, relevanceCapDir = 24L, relevanceCapSparse = 30L,
                           degreeCap = 4L, sampleSlack = 5L, probeRetries = 8L,
                           laurentDegNum = 4L, laurentDegDen = 2L, laurentCandCap = 200000L,
                           termCap = 60L, generalDegNum = 4L, generalDegDen = 3L,
                           gapOrderCap = 8L, minSupportCandCap = 20000L, perPrimeCap = 120L,
                           perPrimeMinPrimes = 3L, homotopy = TRUE, timeout = Inf,
                           homotopyDegree = NULL, degreeCapMax = NULL, exactProof = NULL,
                           exactWorkCap = NULL, seed = NULL) {
  stopifnot(relevanceCap >= 0L, relevanceCapDir >= 1L,
            relevanceCapSparse >= relevanceCap, degreeCap >= 0L,
            sampleSlack >= 0L, probeRetries >= 1L, termCap >= 1L,
            laurentDegNum >= 1L, laurentDegDen >= 0L, laurentCandCap >= 1L,
            generalDegNum >= 1L, generalDegDen >= 1L, gapOrderCap >= 0L,
            minSupportCandCap >= 1L, perPrimeCap >= 1L, perPrimeMinPrimes >= 2L,
            isTRUE(homotopy) || isFALSE(homotopy), is.numeric(timeout), timeout > 0,
            is.null(homotopyDegree) || homotopyDegree >= 1L,
            is.null(degreeCapMax) || degreeCapMax >= 0L,
            is.null(exactProof) || isTRUE(exactProof) || isFALSE(exactProof),
            is.null(exactWorkCap) || exactWorkCap > 0, is.null(seed) || is.numeric(seed))
  int <- function(x) if (is.null(x)) NULL else as.integer(x)
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
                 minSupportCandCap = as.integer(minSupportCandCap),
                 perPrimeCap = as.integer(perPrimeCap),
                 perPrimeMinPrimes = as.integer(perPrimeMinPrimes), homotopy = homotopy,
                 timeout = timeout, homotopyDegree = int(homotopyDegree),
                 degreeCapMax = int(degreeCapMax), exactProof = exactProof,
                 exactWorkCap = if (is.null(exactWorkCap)) NULL else as.numeric(exactWorkCap),
                 seed = int(seed)),
            class = c("reconstcontrol", "list"))
}


# ---- the interface to the Python package symident --------------------------------

.symidentMin <- "0.10.1"
.symidentReq <- paste0("symident>=", .symidentMin)

# whether symident is installed in version .symidentMin or later
.symidentCurrent <- function() {
  v <- tryCatch(reticulate::import("importlib.metadata")$version("symident"),
                error = function(e) NA_character_)
  !is.na(v) && utils::compareVersion(v, .symidentMin) >= 0
}

# A module of symident, installed on first use.
.symident <- function(module = "rjson") {
  .require_ns("reticulate", "symmetryDetection()")
  .require_ns("jsonlite", "symmetryDetection()")
  if (!.symidentCurrent())
    suppressWarnings(try(reticulate::py_require(.symidentReq), silent = TRUE))
  if (!.symidentCurrent()) {
    message("Installing symident into ", reticulate::py_exe(), ".")
    try(reticulate::py_install(.symidentReq, pip = TRUE), silent = TRUE)
    reticulate::import("importlib")$invalidate_caches()
  }
  if (!.symidentCurrent())
    stop("symmetryDetection() needs the Python package symident (>= ", .symidentMin, "). ",
         "Install it into ", reticulate::py_exe(), " with `pip install -U symident`.", call. = FALSE)
  reticulate::import(paste0("symident.", module))
}

# camelCase to snake_case, the argument names of symident
.symSnake <- function(x) tolower(gsub("([a-z0-9])([A-Z])", "\\1_\\2", x))

# The names of reconstControl() in symident, where three differ from the snake case
.symControlSnake <- function(x) {
  py <- c(minSupportCandCap = "minsupport_cand_cap", perPrimeCap = "perprime_cap",
          perPrimeMinPrimes = "perprime_min_primes")
  ifelse(x %in% names(py), py[x], .symSnake(x))
}

# snake_case to camelCase, the field names of the R objects
.symCamel <- function(x) gsub("_([a-z0-9])", "\\U\\1", x, perl = TRUE)

# field names to camelCase, recursively; the fields in `data` are named by model
# symbols and keep their names
.symCamelFields <- function(x, data = c("generator", "complete_generator", "weights",
                                        "display", "factor", "pins", "survivor_meaning",
                                        "carrier_domain", "trafo", "totals",
                                        "state_components")) {
  if (!is.list(x)) return(x)
  nm <- names(x)
  for (i in seq_along(x))
    if (is.null(nm) || !nm[i] %in% data) x[i] <- list(.symCamelFields(x[[i]], data))
  if (!is.null(nm)) names(x) <- .symCamel(nm)
  x
}

# regimes as symident reads them: each entry one comparison or a list of comparisons, a
# list named by condition holding the entries of each condition
.symRegimes <- function(x) {
  entries <- function(r) lapply(if (is.list(r)) r else as.list(r), function(e)
    if (length(e) == 1L) as.character(e) else as.list(as.character(e)))
  if (is.list(x) && !is.null(names(x))) lapply(x, entries) else entries(x)
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
  args <- c(args, list(
    env = as.character(jsonlite::toJSON(env, auto_unbox = TRUE)),
    options = as.character(jsonlite::toJSON(opt, auto_unbox = TRUE, null = "null"))))
  out <- tryCatch(do.call(si[[fun]], args),
                  symident.sym.msolve.MsolveMissingError = function(e) {
                    if (!.offerMsolve())
                      stop("this steady state needs msolve: confirm its installation in an ",
                           "interactive session, or set SYMIDENT_MSOLVE to an msolve ",
                           "executable.", call. = FALSE)
                    do.call(si[[fun]], args)
                  })
  res <- jsonlite::fromJSON(out, simplifyVector = FALSE)
  for (w in unlist(res$warnings)) warning(w, call. = FALSE)
  res$warnings <- NULL
  structure(res, json = out)
}

# Asks whether to install msolve into the symident cache and does so; FALSE when the
# session is not interactive or the answer is no.
.offerMsolve <- function() {
  if (!interactive()) return(FALSE)
  inst <- .symident("install")
  q <- paste0("This steady state needs msolve, which is not installed. Download and build ",
              "msolve ", inst$MSOLVE_VERSION, " into ", inst$cache_dir(), "? (yes/no) ")
  repeat {
    a <- tolower(trimws(readline(q)))
    if (a %in% c("n", "no")) return(FALSE)
    if (a %in% c("y", "yes")) break
  }
  path <- if (.Platform$OS.type == "windows") .buildMsolveRtools(inst)
          else inst$install_msolve(ask = FALSE)
  message("msolve installed in ", path)
  TRUE
}

# Builds msolve with the compiler of R and the static GMP, MPFR and FLINT of Rtools
# into the layout the symident cache expects. symident runs msolve on one thread, so
# OpenMP is a stub.
.buildMsolveRtools <- function(inst) {
  config <- function(x) system2(file.path(R.home("bin"), "R"), c("CMD", "config", x),
                                stdout = TRUE)
  cc <- config("CC")
  libs <- file.path(config("LOCAL_SOFT"), "lib", c("libflint.a", "libmpfr.a", "libgmp.a"))
  if (!nzchar(Sys.which(cc)) || !all(file.exists(libs)))
    stop("building msolve on Windows needs Rtools.", call. = FALSE)
  version <- inst$MSOLVE_VERSION
  work <- tempfile("msolve")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  tarball <- file.path(work, "msolve.tar.gz")
  utils::download.file(sprintf(
    "https://github.com/algebraic-solving/msolve/archive/refs/tags/v%s.tar.gz", version),
    tarball, mode = "wb", quiet = TRUE)
  utils::untar(tarball, exdir = work)
  src <- file.path(work, paste0("msolve-", version))
  writeLines(c(sprintf('#define VERSION "%s"', version), "#define HAVE_FLINT_NMOD_H 1",
               "#define HAVE_TIMESPEC_GET 1"), file.path(src, "config.h"))
  dir.create(file.path(src, "ompstub"))
  writeLines(c(
    "#include <time.h>",
    "static inline int omp_get_thread_num(void) { return 0; }",
    "static inline int omp_get_max_threads(void) { return 1; }",
    "static inline void omp_set_num_threads(int n) { (void)n; }",
    "static inline double omp_get_wtime(void) {",
    "  struct timespec ts; timespec_get(&ts, TIME_UTC); return ts.tv_sec + 1e-9 * ts.tv_nsec;",
    "}"), file.path(src, "ompstub", "omp.h"))
  run <- function(args) {
    out <- suppressWarnings(system2(cc, args, stdout = TRUE, stderr = TRUE))
    if (!is.null(attr(out, "status")))
      stop("the msolve build failed:\n", paste(utils::tail(out, 20), collapse = "\n"),
           call. = FALSE)
  }
  # one translation unit per library; main.c includes libmsolve.c
  units <- c("usolve/usolve.c", "fglm/fglm_core.c", "neogb/gb.c", "msolve/main.c")
  objs <- file.path(work, sub("\\.c$", ".o", basename(units)))
  inc <- paste0("-I", shQuote(c(src, file.path(src, c("src", "ompstub"))), type = "cmd"))
  message("Compiling msolve ", version, ", this takes a minute.")
  for (i in seq_along(units))
    run(c("-std=gnu17", "-O2", "-DHAVE_CONFIG_H", inc, "-w", "-c",
          shQuote(file.path(src, "src", units[i]), type = "cmd"),
          "-o", shQuote(objs[i], type = "cmd")))
  prefix <- file.path(inst$cache_dir(), paste0("deps-msolve-", version))
  marker <- file.path(prefix, ".msolve-complete")
  unlink(marker)
  dir.create(file.path(prefix, "bin"), recursive = TRUE, showWarnings = FALSE)
  exe <- file.path(prefix, "bin", "msolve.exe")
  run(c("-static", "-s", "-o", shQuote(exe, type = "cmd"), shQuote(rev(objs), type = "cmd"),
        shQuote(libs, type = "cmd"), "-lm"))
  inst$msolve_check$cache_clear()
  if (is.null(inst$msolve_check(exe)))
    stop("the msolve build produced no working executable.", call. = FALSE)
  file.create(marker)
  normalizePath(exe)
}

# the result as symident wrote it, for its own report and for the reduction
.symRaw <- function(x) {
  raw <- attr(x, "symident")
  if (is.null(raw))
    stop("the object does not hold the symident result it was built from.", call. = FALSE)
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
    d$stateComponents <- vec(d$stateComponents)
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
