#' Structural non-identifiabilities of an ODE model
#'
#' @description Finds the directions in parameters and initial values along which
#'   the observations of a model do not change. The model `f`, the observables `g`
#'   and an optional `trafo` are equations in any form [as.eqnvec()] accepts. The
#'   computation runs in the Python module `symmetryDetection` through `reticulate`.
#'   `method` selects the engine:
#'
#'   * `"observability"` (default): the nullspace of the observability matrix. The
#'     only engine that is exhaustive and supports `equilibrate`, event times in
#'     parameters and closed-form reconstruction.
#'   * `"polynomial"`: the polynomial Lie-symmetry ansatz of Merkt et al. (2015),
#'     with generators and finite transformations. The cost grows with
#'     `polynomialControl(pMax =)`. Ignores `events` and `conditions`.
#'   * `"scaling"`: scaling symmetries only, from an exact integer kernel. Ignores
#'     `equilibrate`.
#'
#'   All engines accept `exp()`, `exp10()`, `b^x`, `sinh()`, `cosh()` and `tanh()`
#'   of states and parameters. `log()`, fractional powers and `abs()` of
#'   coordinates declared `positive` are analysed in \eqn{L_v = \log v} and
#'   reported in \eqn{v}. The symbol `time` is the model time, known at the start
#'   of the analysis. Otherwise `"observability"` with `symEngine = "modular"`
#'   requires rational expressions, up to free power exponents `x^n` and
#'   observables `a*log(h) + c` with a number `a`. The methods are described in
#'   `vignette("Symmetries")`.
#'
#' @param f Right-hand sides: an [eqnlist], an [eqnvec] or a named character vector.
#'   `reduceCQ` needs an [eqnlist].
#' @param g Observables: one [eqnvec] or named character vector for all conditions,
#'   or, for `"observability"` and `"scaling"`, a list with one per row of
#'   `conditions`. Use a list when the conditions measure different observables.
#'   Without `conditions` and a `trafo` list, its length sets the number of
#'   conditions. An observable `a*log(h) + c` with a number `a` is analysed through
#'   `h`.
#' @param trafo Parameter transformation: one [eqnvec] for all conditions or, for
#'   `"observability"`, a list with one per condition. An entry named like a
#'   parameter is substituted into `f`, `g`, event values and event times; an entry
#'   named like a state is its initial value. On the right-hand side a state name is
#'   the initial value of that state, so a steady state from [steadyStates()] can be
#'   passed as it is. A parameter that enters only through `k = "exp(lk)"`,
#'   `"exp10(lk)"` or `"10^lk"` is analysed in `b^lk` and reported in `lk`.
#' @param method `"observability"` (default), `"polynomial"` or `"scaling"`.
#' @param parameters Character vector of additional symbols treated as parameters.
#' @param forcings Character vector of input states. For `"observability"` an input
#'   starts at 0 and is held at 0 in the steady state; for `"polynomial"` and
#'   `"scaling"` it does not transform.
#' @param events An [eventlist], for `"observability"` and `"scaling"`. An event
#'   time is a number or, for `"observability"` with `symEngine = "modular"`, an
#'   expression in the parameters; list the events in chronological order. Root
#'   events are not supported. The analysis starts at the earliest numeric event
#'   time, or at 0 without events. An event value named like a column of
#'   `conditions` is read from it. A later event on a state inside an exponential
#'   must replace the state or add to it.
#' @param conditions Data frame with one row per experimental condition and columns
#'   named by model symbols or event values. A numeric cell fixes the symbol, a
#'   character cell renames it. A column named like a dynamic state sets its
#'   initial value, one named like a constant state sets its value. A direction is
#'   reported only if it is non-identifiable in every condition.
#' @param fixed Character vector of known symbols. For `"observability"` a fixed
#'   parameter is a known constant and a fixed state has no unknown initial value;
#'   for `"polynomial"` a fixed symbol does not transform.
#' @param equilibrate Logical, `"observability"` only. Start at a steady state of
#'   `f` with the inputs at 0 instead of at free initial values. The earliest events
#'   apply on top; initial values in `trafo` are ignored. Not available for
#'   exponentials of states.
#' @param reduceCQ Logical, [eqnlist] only. `FALSE` (default) keeps every species
#'   and reports the freedom of a conserved moiety on the initial value of one
#'   species. `TRUE` eliminates one species per conserved quantity and reports the
#'   freedom on a total named after [getTotals()], which `fixed` or `trafo` can fix.
#'   Both give the same verdict. Set to `FALSE` with a warning when `trafo` gives the
#'   initial value of a moiety species.
#' @param freeInitial Character vector of states, at most one per conserved
#'   quantity, that carry the free resting value of their moiety under
#'   `equilibrate = TRUE, reduceCQ = FALSE`. Invalid choices are dropped; ignored
#'   with a warning otherwise.
#' @param positive Coordinates known to be positive, as in [symmetryReduction()]:
#'   `TRUE` (default) for all, `FALSE` for none, or a character vector. Sets the
#'   domain on which `completeGenerator` has a flow for every `s`.
#' @param reconstruct Logical, `"observability"` only. Return general directions as
#'   exact rational functions instead of their support. A direction that cannot be
#'   reconstructed or verified keeps `explicit = FALSE`.
#' @param verify Logical (default `TRUE`), `"observability"` only. Where the Lie
#'   order is not certified (`$info$lieCertified`), check that the rank has
#'   saturated and warn if not. `DMOD_SYM_VERIFY_MARGIN` (default 6) sets how far
#'   the check looks; the result is in `$info$verification`.
#' @param cores Number of threads for `"observability"`, shared between the
#'   steady-state solves and the kernel. `"polynomial"` and `"scaling"` run
#'   serially.
#' @param control A [reconstControl()] list for the `"observability"` engine.
#' @param polynomial A [polynomialControl()] list for the `"polynomial"` engine.
#' @param scaling A [scalingControl()] list for the `"scaling"` engine.
#' @param symEngine For `"observability"`: `"modular"` (default) computes over
#'   finite fields and scales to large models; `"symbolic"` is an exact sympy
#'   computation for small models, without `equilibrate` and later events.
#' @param verbose Logical (default `TRUE`). Print the result on return.
#'
#' @return An object of class `symmetrydetection`:
#'   \describe{
#'     \item{`method`}{the engine.}
#'     \item{`identifiable`}{`TRUE` or `FALSE` for `"observability"`; `NA` for
#'       `"scaling"` and `"polynomial"`, whose search is not exhaustive.}
#'     \item{`rank`, `dim`}{rank of the observability matrix and number of
#'       coordinates; `NA` for `"scaling"` and `"polynomial"`.}
#'     \item{`symmetries`}{the directions, each a generator
#'       \eqn{X = \sum_i \eta_i \partial_{z_i}} with `generator` (the components
#'       \eqn{\eta_i} by coordinate), `weights` (integer weights of a scaling, else
#'       `NULL`), `type` (`"scaling"`, removed by fixing one coordinate of the support, or
#'       `"general"`, removed by [symmetryReduction()]), `degree` (`-1` if not
#'       polynomial), `support`, `explicit`, `reason`, `certified`,
#'       `transformation` (`"polynomial"` only), `verified`, `display` (factored
#'       components for printing), `completeGenerator` (`factor * generator`, with
#'       a flow for every `s` on the domain set by `positive`) and `factor`.}
#'     \item{`info`}{`engine`, the Lie order and its certification
#'       (`lieOrderUsed`, `lieOrderDriver`, `lieBudget`, `liePlateau`,
#'       `lieCertified`), `gapOrderUsed`, `conditions`, `segments`, `coordinates`,
#'       `settings`, `elapsed` and `verification`.}
#'     \item{`call`}{the matched call.}
#'   }
#'   `print()` shows the verdict and the generators, `summary()` adds the
#'   computation and `summary(verbose = TRUE)` the settings. Both take `fixed`, a
#'   candidate set of fixed coordinates, and report whether it removes every
#'   scaling direction, and `width`.
#'
#' @details Events after the earliest one split the time line into segments. The
#'   analysis starts at the earliest event, so only an event placed there shows its
#'   transient. A dose on a species that `reduceCQ = TRUE` eliminates is not seen.
#'
#'   An exponential of a state enters the `"observability"` engine as an auxiliary
#'   state whose initial value is a further coordinate, tied to the others by its
#'   differential. The `"polynomial"` and `"scaling"` engines treat it as an
#'   independent variable. Both are exact, since exponentials of exponents that are
#'   linearly independent over the rationals are algebraically independent (Ax
#'   1971).
#'
#' @note The interface, defaults and output structure may still change.
#'
#' @references
#'   Merkt B, Timmer J, Kaschek D (2015). Higher-order Lie symmetries in
#'   identifiability and predictability analysis of dynamic models. Physical
#'   Review E 92, 012920. \doi{10.1103/PhysRevE.92.012920}
#'
#'   Ax J (1971). On Schanuel's conjectures. Annals of Mathematics 93, 252-268.
#'
#' @examples
#' \dontrun{
#' # A reversible reaction observed through alpha * A: the scale is free.
#' eq <- eqnlist() |>
#'   addReaction("A", "B", "k1 * A") |>
#'   addReaction("B", "A", "k2 * B")
#'
#' out <- symmetryDetection(eq, eqnvec(Aobs = "alpha * A"))
#' summary(out)
#' out <- symmetryDetection(eq, eqnvec(Aobs = "alpha * A"), method = "polynomial")
#' out <- symmetryDetection(eq, eqnvec(Aobs = "alpha * A"), method = "scaling")
#'
#' # A state-named trafo entry is an initial value.
#' out <- symmetryDetection(eqnvec(x = "b - a*x"), eqnvec(y = "s*x"),
#'                          trafo = eqnvec(x = "b/a"), reconstruct = TRUE)
#'
#' # Two switch values separate k1 and k2, one does not.
#' fu <- eqnvec(A = "-(k1 + u * k2) * A", u = "0")
#' events <- addEvent(eventlist(), var = "u", time = -1, value = "var_u",
#'                    method = "replace")
#' grid <- data.frame(var_u = c(0, 1), row.names = c("ctrl", "stim"))
#' out <- symmetryDetection(fu, eqnvec(y = "A"), events = events, conditions = grid)
#' out$identifiable
#'
#' # Observables measured in different conditions: one eqnvec per condition.
#' fd <- eqnvec(A = "-k1*A", B = "-k2*B")
#' out <- symmetryDetection(fd, list(eqnvec(yA = "sA*A"), eqnvec(yB = "sB*B")),
#'                          conditions = data.frame(row.names = c("elisa", "wb")))
#'
#' # A known dose makes a, b and s identifiable.
#' dose <- addEvent(eventlist(), var = "x", time = 0, value = "dose",
#'                  method = "replace")
#' out <- symmetryDetection(eqnvec(x = "b - a*x"), eqnvec(y = "s*x"), events = dose,
#'                          conditions = data.frame(dose = 2, row.names = "stim"))
#'
#' # A state inside exp(): A + c, kA + kd*c and kE*exp(-c) give the same y.
#' out <- symmetryDetection(eqnvec(A = "kA - kd*A", B = "kB - kE*exp(A)*B"),
#'                          eqnvec(y = "s*B"), reconstruct = TRUE)
#'
#' # Worked scripts
#' file.edit(system.file("examples", "symmetryDetection.R", package = "dMod2"))
#' file.edit(system.file("examples", "symmetryExponentials.R", package = "dMod2"))
#' }
#' @export
symmetryDetection <- function(f = NULL, g = NULL, trafo = NULL,
                              method = c("observability", "polynomial", "scaling"),
                              parameters = NULL, fixed = NULL, forcings = NULL,
                              events = NULL, conditions = NULL,
                              equilibrate = FALSE,
                              reduceCQ = FALSE, freeInitial = NULL,
                              reconstruct = FALSE, positive = TRUE,
                              verify = TRUE, cores = 1,
                              control = reconstControl(),
                              polynomial = polynomialControl(),
                              scaling = scalingControl(),
                              symEngine = c("modular", "symbolic"),
                              verbose = TRUE) {

  # ==== arguments, engine settings and the report delivery ==========================
  .require_ns("reticulate", "symmetryDetection()")
  .require_ns("lpSolve", "symmetryDetection()")
  symEngine <- match.arg(symEngine)
  method <- match.arg(method)
  equilibrate <- isTRUE(equilibrate)
  # for the summary header and the elapsed time
  .symCall <- match.call()
  .symT0 <- Sys.time()
  # settings for summary(), which shows those of the chosen method
  .symSettings <- list(positive = positive,
                        reduceCQ = isTRUE(reduceCQ), equilibrate = isTRUE(equilibrate),
                        reconstruct = isTRUE(reconstruct), verify = isTRUE(verify),
                        symEngine = symEngine, degreeCap = control$degreeCap,
                        certifyPoly = isTRUE(control$certifyPoly),
                        ansatz = polynomial$ansatz, pMax = polynomial$pMax,
                        polyBackend = if (is.null(polynomial$backend)) "symengine"
                                      else polynomial$backend)
  # set by the `trafo` block when a right-hand side names a state: the analysis runs
  # on a renamed initial-value coordinate and the report maps it back
  icFrom <- icTo <- character(0)
  logArg <- NULL            # set by the log chart below: log_v analysed, v reported
  chartCoords <- character(0)
  # every engine returns through here: map names back, finalise, print unless
  # verbose = FALSE, return invisibly
  deliver <- function(raw, method, coordinates = NULL) {
    raw <- .symLogArgBack(raw, logArg, sd)
    if (!is.null(logArg) && length(coordinates))
      coordinates <- replaceSymbols(logArg$L, logArg$v, as.character(coordinates))
    raw <- .symRenameResult(raw, icFrom, icTo)
    if (length(icFrom) && length(coordinates))
      coordinates <- replaceSymbols(icTo, icFrom, as.character(coordinates))
    res <- .symFinalize(raw, method, .symSettings, .symCall,
                         elapsed = as.numeric(Sys.time() - .symT0, units = "secs"),
                         coordinates = coordinates)
    if (isTRUE(verbose)) print(res)
    invisible(res)
  }

  # warn about arguments the chosen engine ignores
  supplied <- setdiff(names(match.call())[-1], "")
  # NULL counts as the default, so callers forwarding a fixed argument list stay quiet
  supplied <- supplied[vapply(supplied, function(a)
    !is.null(mget(a, envir = environment(), ifnotfound = list(NULL))[[1]]), logical(1))]
  applies <- switch(method,
    observability = c("events", "conditions", "equilibrate", "control"),
    polynomial = "polynomial",
    scaling = c("scaling", "events", "conditions"))
  methodSpecific <- c("events", "conditions", "equilibrate",
                      "control", "polynomial", "scaling")
  ignored <- intersect(supplied, setdiff(methodSpecific, applies))
  if (length(ignored))
    warning("symmetryDetection(): argument(s) ", paste(ignored, collapse = ", "),
            " do not apply to method = \"", method, "\" and are ignored.",
            call. = FALSE)

  # ==== the model: f, the observation g, and the parameter transformation ===========
  # model right-hand sides as an eqnvec; f is eqnlist, eqnvec or named character
  if (is.null(f))
    stop("Provide the model right-hand sides via `f` ",
         "(eqnlist, eqnvec or named character vector).")
  feqnlist <- if (inherits(f, "eqnlist")) f else NULL
  fdyn <- as.eqnvec(f)
  states <- names(fdyn)

  # `g` is one eqnvec for every condition or a list with one per condition (for
  # observables that differ between conditions). Held as the list `gset`.
  gPerCond <- !is.null(g) && is.list(g) && !inherits(g, "eqnvec")
  gset <- if (gPerCond) lapply(g, as.eqnvec)
          else list(if (is.null(g)) NULL else as.eqnvec(g))
  if (gPerCond && !length(gset))
    stop("symmetryDetection(): `g` is an empty list; give one eqnvec per condition.",
         call. = FALSE)
  # flat character/name vectors of the observation; not c() on the eqnvecs, whose
  # names may repeat across conditions
  gChar  <- function() unlist(lapply(gset, as.character), use.names = FALSE)
  gNames <- function() unique(unlist(lapply(gset, names), use.names = FALSE))
  gLines <- function(k = 1L) .symEqnLines(gset[[k]])
  # applies a substitution to the optional model pieces; `fdyn` is handled per pass
  substModel <- function(fn) {
    gset <<- lapply(gset, fn)
    initial <<- fn(initial)
    if (!is.null(condInitial))
      condInitial <<- lapply(condInitial, fn)
  }
  # time in f or g is a clock: time' = 1 from the start of the analysis, the
  # earliest numeric event time or 0
  if (!"time" %in% states &&
      "time" %in% getSymbols(c(as.character(fdyn), gChar()))) {
    evt <- if (!is.null(events) && nrow(as.data.frame(events)))
      suppressWarnings(as.numeric(as.character(as.data.frame(events)$time)))
    t0 <- if (length(evt) && any(!is.na(evt))) min(evt, na.rm = TRUE) else 0
    fdyn <- as.eqnvec(c(setNames(as.character(fdyn), states), time = "1"))
    states <- names(fdyn)
    if (method == "observability") {
      addClock <- function(tr) {
        tr <- if (is.null(tr)) character(0) else as.eqnvec(tr)
        as.eqnvec(c(setNames(as.character(tr), names(tr)),
                    time = format(t0, digits = 15)))
      }
      trafo <- if (!is.null(trafo) && is.list(trafo) && !inherits(trafo, "eqnvec"))
        lapply(trafo, addClock) else addClock(trafo)
    } else fixed <- unique(c(fixed, "time"))
  }
  # an observable named like a state would define that symbol twice in the model lines
  if (length(clashObs <- intersect(gNames(), states)))
    stop("symmetryDetection(): observable(s) ", paste(clashObs, collapse = ", "),
         " in `g` carry the name of a state in `f`. Rename the observable(s) ",
         "(e.g. ", clashObs[1], "_obs).", call. = FALSE)
  parameters <- parameters %||% character(0)

  # freeInitial: moiety species that keep a free resting initial value. Only meaningful
  # with equilibrate = TRUE, reduceCQ = FALSE; otherwise every species already has a
  # free initial value or the moiety is reduced to a `total`.
  freeInitial <- freeInitial %||% character(0)
  if (length(freeInitial)) {
    if (length(bad <- setdiff(freeInitial, states)))
      warning("symmetryDetection(): freeInitial names non-state(s) ",
              paste(bad, collapse = ", "), "; ignored.", call. = FALSE)
    freeInitial <- intersect(freeInitial, states)
    if (length(freeInitial) && !(method == "observability" &&
                                 equilibrate && !isTRUE(reduceCQ))) {
      warning("symmetryDetection(): freeInitial applies only to ",
              "method = \"observability\" with equilibrate = TRUE and ",
              "reduceCQ = FALSE (the held-variable moiety parameterisation); ignored.",
              call. = FALSE)
      freeInitial <- character(0)
    }
  }

  # `trafo`: one eqnvec or a list with one per condition. A parameter-named entry is a
  # substitution, a state-named entry an initial value. A single trafo is substituted
  # up front, a list per condition (composing with the condition grid).
  initial <- NULL
  condSubs <- NULL          # per-condition parameter substitutions (trafo list)
  condInitial <- NULL       # per-condition initial conditions (trafo list)
  trafoSyms <- character(0) # symbols the trafo introduces (substitution targets/params)
  trafoSubs <- NULL         # single-trafo parameter substitutions, kept for the CQ pass
  trafoHit <- character(0)  # substitution targets that actually occurred somewhere
  trafoList <- !is.null(trafo) && is.list(trafo) && !inherits(trafo, "eqnvec")
  if (trafoList) {
    trafos <- lapply(trafo, as.eqnvec)
    icSplit <- .symInitialValueSplit(trafos, states, taken = c(parameters, gNames()))
    trafos <- icSplit$trafos; icFrom <- icSplit$from; icTo <- icSplit$to
    condSubs <- lapply(trafos, function(tr) {
      se <- setdiff(names(tr), states)
      as.list(setNames(as.character(tr[se]), se))
    })
    condInitial <- lapply(trafos, function(tr) {
      ic <- intersect(names(tr), states)
      if (length(ic)) as.eqnvec(tr[ic]) else NULL
    })
    trafoSyms <- unique(unlist(lapply(trafos,
                        function(tr) getSymbols(as.character(tr)))))
    if (method == "polynomial")
      warning("symmetryDetection(): a per-condition `trafo` list is not yet ",
              "supported for method = \"polynomial\"; pass a single trafo.", call. = FALSE)
  } else if (!is.null(trafo)) {
    icSplit <- .symInitialValueSplit(list(as.eqnvec(trafo)), states,
                                     taken = c(parameters, gNames()))
    trafo <- icSplit$trafos[[1]]; icFrom <- icSplit$from; icTo <- icSplit$to
    icEntries  <- intersect(names(trafo), states)
    subEntries <- setdiff(names(trafo), states)
    subs <- trafo[subEntries]
    sub <- .symSubst(names(subs), subs)
    # Record which targets occur, to report entries that substitute into nothing.
    # Identity entries (`x = "x"`) are never reported; symbols in event values or the
    # condition grid count as present. The CQ reduction adds a second pass for totals.
    trafoSubs <- subs[trimws(as.character(subs)) != names(subs)]
    if (length(trafoSubs)) {
      # as.character() per piece: an initial condition shares its state's name, which
      # c.eqnvec() would reject as a duplicate
      knownSyms <- c(getSymbols(c(as.character(fdyn), gChar(),
                                  as.character(trafo[icEntries]))),
                     if (!is.null(events) && nrow(as.data.frame(events)))
                       getSymbols(c(as.character(as.data.frame(events)$value),
                                    as.character(as.data.frame(events)$time))),
                     if (!is.null(conditions))
                       c(names(as.data.frame(conditions)),
                         unlist(lapply(as.data.frame(conditions),
                                       function(cl) if (is.character(cl)) cl else NULL))))
      trafoHit <- intersect(names(trafoSubs), knownSyms)
    }
    fdyn <- sub(fdyn)
    substModel(sub)
    if (length(icEntries)) initial <- sub(trafo[icEntries])
    # event values and times are in the same parameters
    if (!is.null(events) && nrow(as.data.frame(events)) && length(trafoSubs)) {
      events$value <- replaceSymbols(names(trafoSubs), as.character(trafoSubs),
                                     as.character(events$value))
      events$time <- replaceSymbols(names(trafoSubs), as.character(trafoSubs),
                                    as.character(events$time))
    }
  }

  # grid columns and event values naming a renamed initial value follow the rename
  if (length(icFrom)) {
    if (!is.null(conditions)) {
      cdf <- as.data.frame(conditions, stringsAsFactors = FALSE)
      hit <- match(colnames(cdf), icFrom)
      if (any(!is.na(hit))) colnames(cdf)[!is.na(hit)] <- icTo[hit[!is.na(hit)]]
      conditions <- cdf
    }
    if (!is.null(events) && nrow(as.data.frame(events)))
      events$value <- replaceSymbols(icFrom, icTo, as.character(events$value))
  }

  # A grid cell naming a state is substituted into f, where the name is the running
  # state, not its initial value; warn. Cells of state-named columns are initial values.
  if (!is.null(conditions)) {
    cdf <- as.data.frame(conditions, stringsAsFactors = FALSE)
    keep <- setdiff(colnames(cdf), states)
    chr <- keep[vapply(cdf[keep], function(cl)
      is.character(cl) || is.factor(cl), logical(1))]
    hit <- if (!length(chr)) character(0) else
      intersect(getSymbols(as.character(unlist(lapply(cdf[chr], as.character)))), states)
    if (length(hit))
      warning("symmetryDetection(): condition-grid value(s) name the state(s) ",
              paste(hit, collapse = ", "), ". A grid value is substituted into `f`, ",
              "where that name is the running state and not an initial value. ",
              "Pass it through a per-condition `trafo` instead.", call. = FALSE)
  }

  # a per-condition `g` list sets the condition count and must agree with the other
  # per-condition inputs
  nCondObs <- if (gPerCond) length(gset) else 0L
  if (gPerCond) {
    if (!method %in% c("observability", "scaling"))
      stop("symmetryDetection(): a per-condition `g` list applies to ",
           "method = \"observability\" and method = \"scaling\" only; ",
           "method = \"", method, "\" has no conditions. Pass a single eqnvec.",
           call. = FALSE)
    nGridG <- if (is.null(conditions)) 0L else nrow(as.data.frame(conditions))
    if (nGridG && nCondObs != nGridG)
      stop("symmetryDetection(): the per-condition `g` list length (", nCondObs,
           ") must match the condition grid rows (", nGridG, ").", call. = FALSE)
    if (length(condSubs) && nCondObs != length(condSubs))
      stop("symmetryDetection(): the per-condition `g` list length (", nCondObs,
           ") must match the per-condition `trafo` list length (",
           length(condSubs), ").", call. = FALSE)
    # the polynomial certificate needs a single g
    if (isTRUE(control$certifyPoly)) {
      warning("symmetryDetection(): reconstControl(certifyPoly = TRUE) needs a ",
              "single `g` and is skipped for a per-condition `g` list.", call. = FALSE)
      control$certifyPoly <- FALSE
      .symSettings$certifyPoly <- FALSE
    }
  }

  # equilibrate solves the initial conditions from f = 0, so any given in `trafo`
  # (single or per-condition) are dropped with a warning
  icNames <- unique(c(names(initial),
                      unlist(lapply(condInitial, names))))
  if (equilibrate && length(icNames)) {
    warning("symmetryDetection(): equilibrate solves the steady state from f = 0; ",
            "the initial condition(s) for ", paste(icNames, collapse = ", "),
            " in `trafo` are ignored.", call. = FALSE)
    initial <- NULL
    condInitial <- if (!is.null(condInitial))
      lapply(condInitial, function(x) NULL) else NULL
  }

  # reduceCQ = TRUE would eliminate a moiety species together with its supplied initial
  # value and over-report identifiability; force FALSE (equilibrate dropped them above)
  if (isTRUE(reduceCQ) && !equilibrate && !is.null(feqnlist) && length(icNames)) {
    moietyStates <- getSymbols(as.character(getTotals(feqnlist)))
    clash <- intersect(icNames, moietyStates)
    if (length(clash)) {
      warning("symmetryDetection(): a `trafo` initial condition was supplied for ",
              paste(clash, collapse = ", "), ", which participate(s) in a conserved ",
              "quantity; reduceCQ = TRUE would eliminate a moiety species and discard ",
              "that steady-state relation. Forcing reduceCQ = FALSE so the supplied ",
              "steady state is used as given.", call. = FALSE)
      reduceCQ <- FALSE
      .symSettings$reduceCQ <- FALSE
    }
  }

  # ==== conserved-quantity reduction and the substitutions it enables ===============
  # conserved-quantity reduction (eqnlist input only)
  cqTotals <- character(0)  # `total` parameters the reduction introduces
  if (!is.null(feqnlist) && isTRUE(reduceCQ)) {
    totals <- getTotals(feqnlist)
    if (length(totals)) {
      # a moiety species under a free exponent must survive the reduction as a
      # bare symbol; keep it and eliminate another species of its total instead
      avoidCQ <- .symFreeExponentBases(c(as.character(fdyn), gChar()), names(fdyn))
      cq <- .detect_and_substitute_cq(totals, TRUE, fdyn, names(fdyn),
                                      parameters, expressInTotals = TRUE,
                                      avoid = avoidCQ)
      fdyn <- cq$f[setdiff(names(cq$f), cq$elim_states)]
      parameters <- cq$parameters
      states <- names(fdyn)
      if (length(cq$cq_info)) {
        cqTotals <- vapply(cq$cq_info, function(ci) ci$total_name, character(1))
        keys <- vapply(cq$cq_info, function(ci) ci$elim_state, character(1))
        vals <- vapply(cq$cq_info, function(ci) unname(ci$recon_expr), character(1))
        substModel(.symSubst(keys, vals))
      }
    }
  }

  # `trafo` entries naming a `total` apply only now that the totals exist; pinned
  # totals leave `parameters`. A per-condition list is substituted later anyway.
  if (length(trafoSubs) && length(cqTotals)) {
    hitCQ <- intersect(names(trafoSubs), cqTotals)
    if (length(hitCQ)) {
      subsCQ <- trafoSubs[hitCQ]
      subCQ <- .symSubst(names(subsCQ), subsCQ)
      fdyn <- subCQ(fdyn)
      substModel(subCQ)
      parameters <- setdiff(parameters, hitCQ)
      trafoHit <- union(trafoHit, hitCQ)
    }
  }

  # a substitution target that occurs nowhere would silently leave the parameter free
  if (length(trafoSubs)) {
    missTrafo <- setdiff(names(trafoSubs), trafoHit)
    if (length(missTrafo))
      warning("symmetryDetection(): `trafo` ",
              if (length(missTrafo) == 1L) "entry " else "entries ",
              paste(missTrafo, collapse = ", "),
              if (length(missTrafo) == 1L) " does" else " do",
              " not occur in the model and had no effect. ",
              "Note that the `total` parameters of a conserved quantity only ",
              "exist with reduceCQ = TRUE, and are named after getTotals(f).",
              call. = FALSE)
  }

  code_dir <- system.file("code", package = "dMod2")
  sysmod <- reticulate::import("sys", convert = TRUE)
  if (!(code_dir %in% sysmod$path)) sysmod$path <- c(code_dir, sysmod$path)
  sd <- reticulate::import("symmetryDetection", convert = TRUE)

  # abs() and sign() resolved by the declared signs; max, min and steps refused
  if (any(grepl("\\b(abs|sign|max|min|pmax|pmin|ifelse|Heaviside)\\s*\\(",
                c(as.character(fdyn), gChar())))) {
    toL <- function(e) if (is.null(e) || !length(e)) list() else
      as.list(setNames(as.character(e), names(e)))
    sc <- sd$signChart(toL(fdyn), lapply(gset, toL),
                       if (isTRUE(positive)) TRUE else
                         as.list(if (isFALSE(positive)) character(0) else positive))
    if (!is.null(sc$why)) stop("symmetryDetection(): ", sc$why, ".", call. = FALSE)
    fromS <- function(e) if (!length(e)) NULL else
      as.eqnvec(setNames(gsub("**", "^", unlist(e), fixed = TRUE), names(e)))
    fdyn <- fromS(sc$f)
    gset <- lapply(sc$g, fromS)
  }

  # log() or a fractional power of a positive v: analysed in log_v = log(v), where
  # only exponentials remain, and mapped back in deliver(). Not for the scaling
  # engine, which takes log() itself and does not see translations.
  if (!isFALSE(positive) && method != "scaling") {
    toList <- function(e) if (is.null(e) || !length(e)) list() else
      as.list(setNames(as.character(e), names(e)))
    # replace and multiply events enter the chart like initial values
    evIc <- if (!is.null(events) && nrow(as.data.frame(events))) {
      evdf0 <- as.data.frame(events, stringsAsFactors = FALSE)
      keep <- as.character(evdf0$method) %in% c("replace", "multiply")
      lapply(which(keep), function(i)
        setNames(list(as.character(evdf0$value[i])), as.character(evdf0$var[i])))
    }
    lc <- sd$logArgChart(toList(fdyn), lapply(gset, toList),
                         if (isTRUE(positive)) TRUE else as.list(positive),
                         as.list(unique(c(parameters, gNames(), icTo))),
                         unname(c(list(toList(initial)), lapply(condInitial, toList),
                                  evIc)))
    if (!is.null(lc$why))
      stop("symmetryDetection(): ", lc$why, ".", call. = FALSE)
    if (!is.null(lc)) {
      vs <- vapply(lc$map, function(m) m$v, "")
      Ls <- vapply(lc$map, function(m) m$L, "")
      clash <- intersect(vs, c(forcings,
        if (!is.null(conditions)) names(as.data.frame(conditions)),
        unlist(lapply(condSubs, function(cs) c(names(cs), getSymbols(unlist(cs)))))))
      if (length(clash))
        stop("symmetryDetection(): ", paste(clash, collapse = ", "), " appear(s) ",
             "inside log() or under a fractional power and is analysed as ",
             "exp(log_<name>); forcings and conditions on ",
             "it are not supported.", call. = FALSE)
      fromPy <- function(e) if (!length(e)) NULL else
        as.eqnvec(setNames(gsub("**", "^", unlist(e), fixed = TRUE), names(e)))
      fdyn <- fromPy(lc$f)
      states <- names(fdyn)
      gset <- lapply(lc$g, fromPy)
      # log(p) of a prime p in a known number: a known constant
      lognum <- unlist(lc$consts)
      if (!is.null(events) && nrow(as.data.frame(events))) {
        evdf <- as.data.frame(events, stringsAsFactors = FALSE)
        for (i in seq_len(nrow(evdf))) {
          ce <- sd$logArgEvent(as.character(evdf$var[i]), as.character(evdf$value[i]),
                               as.character(evdf$method[i]), lc$map)
          if (!is.null(ce$why))
            stop("symmetryDetection(): ", ce$why, ".", call. = FALSE)
          evdf$var[i] <- ce$var
          evdf$value[i] <- gsub("**", "^", ce$value, fixed = TRUE)
          evdf$method[i] <- ce$method
          lognum <- c(lognum, unlist(ce$consts))
        }
        evdf$time <- replaceSymbols(vs, paste0("exp(", Ls, ")"), as.character(evdf$time))
        events <- as.eventlist(evdf)
      }
      initial <- fromPy(lc$ic[[1]])
      if (!is.null(condInitial))
        condInitial <- lapply(lc$ic[-1], fromPy)
      # the symbols of a and b stay coordinates although the chart absorbs them
      ab <- getSymbols(unlist(lapply(lc$map, function(m) c(m$a, m$b))))
      chartCoords <- setdiff(ab, c(states, fixed, vs))
      parameters <- unique(c(replaceSymbols(vs, Ls, parameters), chartCoords))
      fixed <- if (length(fixed)) replaceSymbols(vs, Ls, fixed) else fixed
      if (length(lognum)) fixed <- unique(c(fixed, lognum))
      logArg <- list(v = vs, L = Ls, map = lc$map)
    }
  }

  # an event fires at a known time, a number or an expression in the parameters; a
  # time given in parameters needs the modular observability engine
  if (!is.null(events) && nrow(as.data.frame(events))) {
    evdf <- as.data.frame(events)
    if (!is.null(evdf$root) && any(!is.na(evdf$root)))
      stop("symmetryDetection(): root events are not supported; give the event a time.",
           call. = FALSE)
    et <- as.character(evdf$time)
    symTime <- is.na(suppressWarnings(as.numeric(et)))
    if (any(symTime) && !(method == "observability" && symEngine == "modular"))
      stop("symmetryDetection(): an event time given in parameters (",
           paste(unique(et[symTime]), collapse = ", "), ") needs method = ",
           "\"observability\" with symEngine = \"modular\".", call. = FALSE)
  }

  # ==== engine: observability (conditions/segments, symbolic or modular) ============
  if (method == "observability") {
    # condition/event resolution shared by the symbolic and modular engines. Grid
    # targets include symbols only in initial values, event values or a trafo list.
    extraSyms <- c(if (!is.null(initial)) getSymbols(as.character(as.eqnvec(initial))),
                   if (!is.null(events)) getSymbols(c(as.character(as.data.frame(events)$value),
                                                      as.character(as.data.frame(events)$time))),
                   if (length(condInitial)) unlist(lapply(condInitial,
                     function(x) if (is.null(x)) NULL else getSymbols(as.character(x)))),
                   trafoSyms)
    symbols <- unique(c(states, gNames(),
                        getSymbols(c(as.character(fdyn), gChar())), extraSyms))
    # states with a right-hand side that evaluates to 0 are constant in time and
    # substituted by their per-condition value (boolean switches, held inputs)
    isZeroRHS <- vapply(as.character(fdyn), function(r)
      isTRUE(suppressWarnings(tryCatch(eval(parse(text = r)) == 0,
                                       error = function(e) FALSE))), logical(1))
    constStates <- names(fdyn)[isZeroRHS]
    # states forced to zero at the resting state are held at zero in the f = 0
    # solve but stay dynamic states in the observability tape
    equilZeroStates <- if (equilibrate && !is.null(feqnlist))
      intersect(.equil_zero_states(feqnlist, forcings), states) else character(0)
    # equilibrate without reduceCQ: f = 0 lacks one equation per moiety, so one pivot
    # species per moiety keeps its resting value free. Pivots stay dynamic states and
    # exclude zero, forced and free-exponent species.
    heldStateParams <- character(0)   # named: pivot state -> initial-value parameter
    if (equilibrate && !isTRUE(reduceCQ) && !is.null(feqnlist)) {
      totalsFV <- getTotals(feqnlist)
      if (length(totalsFV)) {
        avoidFV <- unique(c(equilZeroStates, forcings,
          .symFreeExponentBases(c(as.character(fdyn), gChar()), names(fdyn))))
        decFV <- .cq_pivot_decomposition(totalsFV, states, parameters, avoid = avoidFV,
                                         prefer = freeInitial)
        piv <- intersect(decFV$pivots[!is.na(decFV$pivots)], states)
        piv <- setdiff(piv, c(equilZeroStates, forcings))
        if (length(piv)) {
          # held under the pivot's own name, so the moiety freedom is reported as its
          # initial value (dMod convention), not as a `total`
          heldStateParams <- setNames(piv, piv)
        }
      }
    }
    res <- .symResolveConditions(conditions, events, initial, symbols, states,
                                   constStates, forcings, equilibrate = equilibrate,
                                   condSubs = condSubs, condInitial = condInitial,
                                   nCondObs = nCondObs)
    # one observation per segment, expanded along chainOf; NULL for a shared `g`
    segObs <- if (gPerCond) lapply(res$chainOf, gLines) else NULL

    # codimension of the specialisation, the budget of .symSaturateCertify(): each pinned
    # coordinate of (x0, theta) counts once; with equilibrate every state counts
    modelSyms <- getSymbols(c(as.character(fdyn), gChar()))
    pinnedIC <- unique(c(if (equilibrate) states else icNames,
                         .symEventPinnedStates(events, conditions),
                         .symGridPinnedStates(conditions, states),
                         intersect(c(forcings, fixed, constStates), states)))
    gridPars <- if (is.null(conditions)) character(0) else {
      cdfC <- as.data.frame(conditions, stringsAsFactors = FALSE)
      cols <- setdiff(colnames(cdfC), states)
      intersect(cols[vapply(cdfC[cols], is.numeric, logical(1))], modelSyms)
    }
    pinnedPar <- unique(c(names(trafoSubs), unlist(lapply(condSubs, names)),
                          setdiff(fixed, states), gridPars))
    codimSpec <- length(pinnedIC) + length(setdiff(pinnedPar, ""))

    spy <- tryCatch(reticulate::import("sympy", convert = TRUE),
                    error = function(err) NULL)

    # symbolic engine: the same matrix reduced exactly with sympy up to the saturation
    # bound, independent of the modular kernel; small models, no gaps or equilibrate
    if (symEngine == "symbolic") {
      if (isTRUE(equilibrate))
        stop("symEngine = \"symbolic\" does not support equilibrate; supply the ",
             "steady state explicitly through `trafo` (e.g. from steadyStates()), ",
             "which is an exact substitution, or use symEngine = \"modular\".",
             call. = FALSE)
      if (res$nGaps > 0L)
        stop("symEngine = \"symbolic\" handles single-segment conditions only; ",
             "later events (gaps) need symEngine = \"modular\".", call. = FALSE)
      sr <- sd$observabilitySympyMulti(
        model = .symEqnLines(fdyn), observation = gLines(),
        conditionSubs = res$subs, conditionIC0 = res$ic0,
        conditionObs = segObs,
        fixed = if (length(fixed)) fixed else NULL,
        parameters = if (length(parameters)) parameters else NULL,
        inputs = if (length(forcings)) forcings else NULL)
      if (!isTRUE(sr$ok))
        stop("symEngine = \"symbolic\": ",
             if (!is.null(sr$why)) sr$why else "could not build the symbolic system.",
             call. = FALSE)
      sr$method <- "observability"; sr$engine <- "symbolic"
      sr$lieOrderUsed <- as.integer(sr$lieOrder)
      sr$conditions <- as.integer(res$nConditions)
      sr$segments <- as.integer(res$nConditions)   # single-segment: one per condition
      sr$gapOrderUsed <- 0L
      sr$nonIdentifiable <- .symRelabelDirections(sr$nonIdentifiable, sd)
      if (isTRUE(control$certifyPoly))
        sr$nonIdentifiable <- .symCertifyPoly(sr$nonIdentifiable,
          .symEqnLines(fdyn), gLines(), forcings, fixed, parameters, control, sd)
      return(deliver(sr, method))
    }
    # equilibrate uses the implicit system: states stay coordinates and the steady
    # state enters as the tangency constraint df.xi = 0, one resting state per condition
    useImplicit <- isTRUE(equilibrate)
    runObs <- function(ui) {
      multi <- sd$compileObservabilityTapeMulti(
        model = .symEqnLines(fdyn), observation = gLines(),
        conditionSubs = res$subs, conditionIC0 = res$ic0,
        conditionObs = segObs,
        fixed = if (length(fixed)) fixed else NULL,
        parameters = if (length(parameters)) parameters else NULL,
        equilibrate = equilibrate, segEquilibrate = as.list(res$segEquil),
        forcings = if (length(forcings)) forcings else NULL,
        conditionEvents = res$segEvents, conditionT0Events = res$events0,
        conditionTimes = if (res$nGaps > 0L) res$times else NULL,
        keepCoords = if (length(chartCoords)) as.list(chartCoords) else NULL,
        jointSteadyState = isTRUE(ui),
        jointFixedStates = if (isTRUE(ui) && length(equilZeroStates))
          equilZeroStates else NULL,
        heldStateParams = if (isTRUE(ui) && length(heldStateParams))
          as.list(heldStateParams) else NULL)
      if (!isTRUE(multi$ok))
        return(list(ok = FALSE, nonrational = multi$nonrational, why = multi$why))
      list(ok = TRUE, result = .symLogParamBack(.observability_analytic_multi(multi, spy = spy,
             closedForm = reconstruct, sd = sd, cores = cores,
             equilZeroStates = equilZeroStates, t0events = res$events0,
             nConditions = res$nConditions, chainOf = res$chainOf,
             nGaps = res$nGaps, implicitSteadyState = isTRUE(ui), control = control,
             verify = verify, codimSpec = codimSpec), multi$logParams, sd))
    }
    ro <- runObs(useImplicit)
    if (useImplicit && !isFALSE(ro$ok) && is.null(ro$result))
      stop("symmetryDetection(): the implicit steady-state path could not be ",
           "evaluated (e.g. a singular resting Jacobian from an unreduced conserved ",
           "moiety). Reduce conserved moieties (reduceCQ = TRUE) or supply an ",
           "explicit steady state through `trafo` (from steadyStates()).",
           call. = FALSE)
    if (isFALSE(ro$ok) && !is.null(ro$why))
      stop("symmetryDetection(): ", ro$why, ".", call. = FALSE)
    if (isFALSE(ro$ok))
      stop("method = \"observability\" requires right-hand sides, observables and ",
           "initial conditions built from +, -, *, /, integer powers, exp(), b^x, ",
           "hyperbolic functions and free power exponents x^n, log() and fractional ",
           "powers of coordinates declared `positive`, and observables ",
           "a*log(h) + offset with a number a; anything else is not rational.\n  ",
           paste(unlist(ro$nonrational), collapse = "\n  "),
           "\nUse symEngine = \"symbolic\" for other functions.", call. = FALSE)
    res <- ro$result
    if (is.list(res)) {
      res$nonIdentifiable <- .symRelabelDirections(res$nonIdentifiable, sd)
      if (isTRUE(control$certifyPoly))
        res$nonIdentifiable <- .symCertifyPoly(res$nonIdentifiable,
          .symEqnLines(fdyn), gLines(), forcings, fixed, parameters, control, sd)
    }
    if (isTRUE(verify) && is.list(res) && is.list(res$verification) &&
        isFALSE(res$verification$ok))
      warning("symmetryDetection(verify = TRUE): the Schwartz-Zippel saturation guard ",
              "found the rank still growing past the reported Lie order; the ",
              "directions may be over-reported; inspect $verification (",
              res$verification$reason, ").", call. = FALSE)
    return(deliver(res, method))
  }

  # ==== engine: scaling (exact integer kernel) ======================================
  # A known dose pins its state (weight 0); conditions intersect their scaling
  # lattices. A scaling leaves f = 0 invariant, so `equilibrate` does not apply.
  if (method == "scaling") {
    fixedScal <- unique(c(fixed, .symEventPinnedStates(events, conditions),
                          .symGridPinnedStates(conditions, states)))
    syms <- unique(c(states, gNames(), getSymbols(c(as.character(fdyn), gChar()))))
    multiCond <- (!is.null(conditions) && nrow(as.data.frame(conditions)) > 1L) ||
                 length(condSubs) > 1L || nCondObs > 1L
    reticulate::py_capture_output(
      res <- if (multiCond) {
        pc <- .symPercondLines(fdyn, gset, gPerCond, conditions, condSubs, syms)
        sd$scalingSymmetriesMulti(
          perCondModel = lapply(pc, `[[`, "f"),
          perCondObs   = lapply(pc, `[[`, "g"),
          inputs = if (length(forcings)) forcings else NULL,
          fixed  = if (length(fixedScal)) fixedScal else NULL, logs = TRUE)
      } else {
        sd$symmetryDetectiondMod(
          model = .symEqnLines(fdyn), observation = gLines(),
          inputs = if (length(forcings)) forcings else NULL,
          fixed = fixedScal, backend = scaling$backend,
          parameters = if (length(parameters)) parameters else NULL, method = "scaling")
      })
    return(deliver(res, "scaling"))
  }

  # ==== engine: polynomial (Lie-symmetry ansatz) ====================================
  # Merkt et al. (2015) on the trafo-substituted f and g; the Python tag is "liesym"
  fld <- function(nm, default) if (is.null(polynomial[[nm]])) default else polynomial[[nm]]
  # the Python report is captured; print()/summary() display the result
  reticulate::py_capture_output(
    res <- sd$symmetryDetectiondMod(
      model       = .symEqnLines(fdyn),
      observation = gLines(),
      ansatz      = fld("ansatz", "uni"),
      pMax        = as.integer(fld("pMax", 2L)),
      inputs      = if (length(forcings)) forcings else NULL,
      fixed       = fixed,
      allTrafos   = fld("allTrafos", FALSE),
      lieOrder    = as.integer(fld("lieOrder", 0L)),
      exact       = fld("exact", TRUE),
      verify      = fld("verify", TRUE),
      backend     = if (is.null(polynomial$backend)) "symengine" else polynomial$backend,
      parameters  = if (length(parameters)) parameters else NULL,
      method      = "liesym"
    ))
  polyCoords <- setdiff(unique(c(states,
                                  getSymbols(c(as.character(fdyn), gChar())),
                                  parameters)),
                         c(forcings, fixed))
  deliver(res, "polynomial", coordinates = polyCoords)
}


# ---- Schwartz-Zippel saturation guard (verify = TRUE) ---------------------------
# A rank can plateau and grow again, so a premature Lie stop over-reports
# non-identifiability. Extends the Lie order past NtUsed on the existing kernel and
# base point (only the jet grows). The rank is monotone in the order, so the far end
# decides; the orders in between are built only on failure, to locate the growth.
.symSzSaturationGuard <- function(kcall, point0Solved, NtUsed, reportedRank,
                                     margin = as.integer(Sys.getenv("DMOD_SYM_VERIFY_MARGIN", "6"))) {
  P <- .symPrimes[1]                          # saturation prime: point0Solved's solve is cached
  r0 <- tryCatch(kcall(point0Solved, P, as.integer(NtUsed)), error = function(e) NULL)
  if (is.null(r0) || !isTRUE(r0$ok))
    return(list(ok = NA, method = "saturation guard",
                reason = "base point not re-evaluable"))
  base <- as.integer(.symRankOf(r0)); maxR <- base; growAt <- NA_integer_
  margin <- max(1L, margin)
  rankAt <- function(k) {                     # extend the Lie order; only the jet grows
    rk <- tryCatch(kcall(point0Solved, P, as.integer(NtUsed + k)), error = function(e) NULL)
    if (is.null(rk) || !isTRUE(rk$ok)) NA_integer_ else as.integer(.symRankOf(rk))
  }
  top <- rankAt(margin)
  if (!is.na(top)) maxR <- max(maxR, top)
  if (is.na(top) || top > base) {             # only a failure pays for the window
    for (k in seq_len(margin)) {
      rk <- rankAt(k)
      if (is.na(rk)) next
      if (rk > maxR) maxR <- rk
      if (rk > base) { growAt <- as.integer(NtUsed + k); break }
    }
    if (is.na(growAt) && !is.na(top) && top > base)
      growAt <- as.integer(NtUsed + margin)
  }
  ok <- is.na(growAt)
  list(ok = ok, method = "saturation guard", lieOrderUsed = as.integer(NtUsed),
       ordersChecked = as.integer(NtUsed + margin),
       kernelRank = base, kernelRankExtended = maxR, growAt = growAt,
       reason = if (ok)
         sprintf("rank %d stable through Lie order %d (%d orders beyond the reported saturation)",
                 base, NtUsed + margin, margin)
       else sprintf("rank grows %d -> %d at Lie order %d (the reported Lie order was premature)",
                    base, maxR, growAt))}


# states whose initial value a grid pins to a number, so they cannot scale
.symGridPinnedStates <- function(conditions, states) {
  if (is.null(conditions) || !length(states)) return(character(0))
  grid <- as.data.frame(conditions, stringsAsFactors = FALSE)
  cand <- intersect(colnames(grid), states)
  cand[vapply(cand, function(cl) {
    v <- as.character(unlist(grid[[cl]]))
    length(v) > 0L && all(!is.na(suppressWarnings(as.numeric(v))))
  }, logical(1))]
}


# States pinned by a numeric replace/add dose in any condition cannot scale. A multiply
# dose or a symbolic value pins nothing; a grid-column value is resolved to its cells.
.symEventPinnedStates <- function(events, conditions) {
  if (is.null(events)) return(character(0))
  ev <- as.data.frame(events, stringsAsFactors = FALSE)
  if (!nrow(ev)) return(character(0))
  grid <- if (is.null(conditions)) NULL else as.data.frame(conditions, stringsAsFactors = FALSE)
  cols <- if (is.null(grid)) character(0) else colnames(grid)
  isKnown <- function(v) {
    v <- as.character(v)
    vals <- if (v %in% cols) as.character(unlist(grid[[v]])) else v
    length(vals) > 0L && all(!is.na(suppressWarnings(as.numeric(vals))))
  }
  pinned <- character(0)
  for (i in seq_len(nrow(ev)))
    if (as.character(ev$method[i]) %in% c("replace", "add") && isKnown(ev$value[i]))
      pinned <- c(pinned, as.character(ev$var[i]))
  unique(pinned)
}


# "<name> = <rhs>" model lines, the form every Python engine entry point reads.
.symEqnLines <- function(e) {
  if (is.null(e)) return(NULL)
  e <- as.eqnvec(e)
  if (!length(e)) return(NULL)
  as.character(vapply(seq_along(e), function(i) paste(names(e)[i], "=", e[i]),
                      character(1)))
}


# A state name on the right of a `trafo` entry is that state's initial value, but
# substituted into f and g it would read as the running state (a steadyStates() trafo
# would collapse f to 0). Each such symbol is renamed to a fresh initial-value
# coordinate in all right-hand sides, and a state without its own initial condition
# gets one tying it to that coordinate. Returns the trafos and the rename `from`, `to`.
.symInitialValueSplit <- function(trafos, states, taken = character(0)) {
  trafos <- lapply(trafos, as.eqnvec)
  none <- list(trafos = trafos, from = character(0), to = character(0))
  if (!length(trafos) || !length(states)) return(none)
  conflict <- unique(unlist(lapply(trafos, function(tr) {
    subEntries <- setdiff(names(tr), states)
    if (!length(subEntries)) return(character(0))
    intersect(getSymbols(as.character(tr[subEntries])), states)
  })))
  conflict <- intersect(states, conflict)          # model order, for a stable rename
  if (!length(conflict)) return(none)

  taken <- unique(c(taken, states, unlist(lapply(trafos, names)),
                    unlist(lapply(trafos, function(tr) getSymbols(as.character(tr))))))
  to <- character(0)
  for (X in conflict) {
    cand <- paste0(X, "_init")
    while (cand %in% taken) cand <- paste0(cand, "0")
    taken <- c(taken, cand); to <- c(to, cand)
  }

  out <- lapply(trafos, function(tr) {
    tr <- as.eqnvec(setNames(replaceSymbols(conflict, to, as.character(tr)), names(tr)))
    add <- setdiff(conflict, intersect(names(tr), states))
    if (length(add))
      tr <- as.eqnvec(setNames(c(as.character(tr), to[match(add, conflict)]),
                               c(names(tr), add)))
    tr
  })
  list(trafos = out, from = conflict, to = to)
}

# Map the initial-value coordinates of .symInitialValueSplit() back onto the names the
# user gave, in the coordinate list and in every direction the engine returned.
.symRenameResult <- function(raw, from, to) {
  if (!length(from) || is.null(raw) || !is.list(raw)) return(raw)
  ren <- function(x) if (!length(x)) x else replaceSymbols(to, from, as.character(x))
  renNames <- function(v) {
    i <- match(names(v), to)
    if (any(!is.na(i))) names(v)[!is.na(i)] <- from[i[!is.na(i)]]
    v
  }
  renVals <- function(v) {
    if (is.list(v)) { v[] <- lapply(v, function(e) if (is.character(e)) ren(e) else e); v }
    else if (is.character(v)) { v[] <- ren(v); v } else v
  }
  fixDir <- function(d) {
    if (!is.list(d)) return(d)
    for (fld in c("vector", "infinitesimals"))
      if (!is.null(d[[fld]])) d[[fld]] <- renVals(renNames(d[[fld]]))
    if (!is.null(d$support)) d$support <- ren(d$support)
    d
  }
  if (!is.null(raw$coordinates)) raw$coordinates <- ren(raw$coordinates)
  if (!is.null(raw$nonIdentifiable))
    raw$nonIdentifiable <- lapply(raw$nonIdentifiable, fixDir)
  else if (is.null(raw$rank)) raw <- lapply(raw, fixDir)   # polynomial: bare list
  raw
}


# substitution keys -> vals as an eqnvec-preserving function; NULL or empty input
# passes through, so it maps over the optional model pieces unguarded
.symSubst <- function(keys, vals) {
  keys <- as.character(keys); vals <- as.character(vals)
  function(e) {
    if (is.null(e) || !length(e) || !length(keys)) return(e)
    e <- as.eqnvec(e)
    setNames(replaceSymbols(keys, vals, e), names(e))
  }
}


# Per-condition (f, g) lines for the scaling engine: grid cells and the trafo list
# substitute parameters only, so all conditions share the states. `gset` has one
# entry per condition when `gPerCond`, else one for all.
.symPercondLines <- function(fdyn, gset, gPerCond, conditions, condSubs, symbols) {
  grid <- if (is.null(conditions)) NULL else as.data.frame(conditions, stringsAsFactors = FALSE)
  nGrid <- if (is.null(grid)) 0L else nrow(grid)
  K <- max(length(condSubs), nGrid, if (gPerCond) length(gset) else 0L, 1L)
  # states are never baked into f (see above); a grid column naming one fixes that
  # state's initial value, which the scaling engine takes through `fixed`
  subCols <- if (is.null(grid)) character(0)
             else setdiff(intersect(colnames(grid), symbols), names(fdyn))
  lapply(seq_len(K), function(k) {
    subs <- list()
    for (col in subCols) subs[[col]] <- as.character(grid[k, col])
    if (length(condSubs) && !is.null(condSubs[[k]]))
      for (nm in names(condSubs[[k]])) subs[[nm]] <- condSubs[[k]][[nm]]
    sub <- function(e) { e <- as.eqnvec(e)
      if (!length(subs)) e else
        setNames(replaceSymbols(names(subs), unlist(subs), e), names(e)) }
    list(f = .symEqnLines(sub(fdyn)),
         g = .symEqnLines(sub(gset[[if (gPerCond) k else 1L]])))
  })
}


# species that are the base of a free exponent (`C3` in `C3^nhill`); eliminating one in
# the CQ reduction would put a sum under the exponent, so they are no pivots
.symFreeExponentBases <- function(exprs, states) {
  if (!length(exprs) || !length(states)) return(character(0))
  hits <- regmatches(exprs, gregexpr(
    "[A-Za-z_][A-Za-z0-9_]*\\s*\\^\\s*(?![0-9])", exprs, perl = TRUE))
  bases <- sub("\\s*\\^.*$", "", unlist(hits, use.names = FALSE))
  intersect(unique(bases), states)
}


# `expr` identically zero once the `dead` symbols are set to zero, tested at
# random positive values for the surviving symbols
.expr_is_zero <- function(expr, dead, ntry = 4, tol = 1e-9) {
  free <- setdiff(getSymbols(expr), dead)
  for (i in seq_len(ntry)) {
    env <- as.list(setNames(rep(0, length(dead)), dead))
    if (length(free))
      env <- c(env, as.list(setNames(runif(length(free), 0.5, 1.5), free)))
    val <- tryCatch(eval(parse(text = expr), env), error = function(e) NA_real_)
    if (is.na(val) || abs(val) > tol) return(FALSE)
  }
  TRUE
}

# Sink cluster among the live reactions: species with a nonnegative combination whose
# net production is <= 0 in every reaction and < 0 overall, so it vanishes at steady
# state. One LP per candidate species on the stoichiometry (as in steadyStates).
.equil_sink_cluster <- function(M, eps = 1e-8, Mbig = 1e4) {
  nF <- nrow(M); nS <- ncol(M)
  if (nF == 0L || nS == 0L) return(integer(0))
  c_obj <- colSums(M); id <- diag(nS)
  for (i in seq_len(nS)) {
    lb <- rep(0, nS); ub <- rep(Mbig, nS); lb[i] <- 1; ub[i] <- 1
    res <- tryCatch(
      lpSolve::lp("min", c_obj, rbind(M, id, id),
                  c(rep("<=", nF), rep(">=", nS), rep("<=", nS)),
                  c(rep(0, nF), lb, ub)),
      error = function(e) NULL)
    if (!is.null(res) && res$status == 0 && res$objval < -eps)
      return(which(res$solution > eps))
  }
  integer(0)
}

# states forced to zero at steady state with the forcings at zero: a species dies when
# all its producing rates vanish, iterated to a fixpoint; a stall tries a sink cluster
.equil_zero_states <- function(eqnlist, forcings) {
  S <- eqnlist$smatrix
  if (is.null(S) || !length(S)) return(character(0))
  S[is.na(S)] <- 0
  species <- colnames(S)
  rates <- as.character(eqnlist$rates)
  seed <- intersect(forcings, c(species, getSymbols(paste(rates, collapse = "+"))))
  dead <- seed
  repeat {
    live <- !vapply(rates, .expr_is_zero, logical(1), dead = dead)
    newdead <- character(0)
    for (sp in setdiff(species, dead)) {
      prod <- which(S[, sp] > 0)
      if (length(prod) && !any(live[prod])) newdead <- c(newdead, sp)
    }
    if (!length(setdiff(newdead, dead))) {
      alive <- setdiff(species, dead)
      cl <- .equil_sink_cluster(S[live, alive, drop = FALSE])
      if (length(cl)) newdead <- alive[cl] else break
    }
    if (!length(setdiff(newdead, dead))) break
    dead <- union(dead, newdead)
  }
  setdiff(dead, forcings)
}


# Analytic observability path: the matrix over GF(p) from the C++ kernel
# (src/symmetry_kernel.cpp), certified rank and directions, exact scalings peeled off,
# the rest reconstructed as rational functions by interpolation per prime, Chinese
# remaindering and rational reconstruction. No floating point.

# four primes < 2^31 (their product < 2^124, within unsigned __int128)
.symPrimes <- c(2147483647, 2147483629, 2147483587, 2147483579)

# a fifth prime, disjoint from the reconstruction primes, used only to certify a
# reconstructed direction against the nullspace at a fresh evaluation
.symVerifyPrime <- 2147483563

# TRUE once the reconstruction is past reconstControl(timeout =); `deadline` is a
# POSIXct set at its start, or NULL
.symExpired <- function(ctrl) {
  d <- ctrl$deadline
  !is.null(d) && Sys.time() > d
}

#' Settings of the observability engine
#'
#' Saturation and closed-form reconstruction settings for
#' `symmetryDetection(method = "observability", control = reconstControl())`. Raise
#' the caps to reconstruct wide or high-degree directions, at the cost of more
#' samples.
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
#'   directions with small support; wider directions go to the general fit.
#' @param perprimeCap Maximum number of samples per prime for the reconstruction
#'   under `equilibrate = TRUE`. A direction that needs more at degree 2 is reported
#'   by its support.
#' @param perprimeMinPrimes Minimum number of primes with samples for a
#'   reconstruction under `equilibrate = TRUE`.
#' @param certifyPoly Logical. Check with the `"polynomial"` engine whether each
#'   polynomial direction of degree up to `degreeCap` is a Lie point symmetry, and
#'   set `$certified`. Costly; `FALSE` by default.
#' @param certifyPolyDeg Ansatz degree for `certifyPoly`; `NULL` (default) uses the
#'   degree of each direction.
#' @param timeout Time limit in seconds for the reconstruction. Directions not
#'   finished in time are reported by their support. `Inf` (default) sets no limit.
#' @return A `reconstControl` list.
#' @seealso [symmetryDetection()]
#' @export
reconstControl <- function(relevanceCap       = 6L,
                           relevanceCapDir    = 24L,
                           relevanceCapSparse = 30L,
                           degreeCap          = 4L,
                           sampleSlack        = 5L,
                           probeRetries       = 8L,
                           laurentDegNum      = 4L,
                           laurentDegDen      = 2L,
                           laurentCandCap     = 200000L,
                           termCap            = 60L,
                           generalDegNum      = 4L,
                           generalDegDen      = 3L,
                           gapOrderCap        = 8L,
                           minsupportCandCap  = 20000L,
                           perprimeCap        = 120L,
                           perprimeMinPrimes  = 3L,
                           certifyPoly        = FALSE,
                           certifyPolyDeg     = NULL,
                           timeout            = Inf) {
  stopifnot(relevanceCap >= 0L, relevanceCapDir >= 1L,
            relevanceCapSparse >= relevanceCap, degreeCap >= 0L,
            sampleSlack >= 0L, probeRetries >= 1L, termCap >= 1L,
            laurentDegNum >= 1L, laurentDegDen >= 0L, laurentCandCap >= 1L,
            generalDegNum >= 1L, generalDegDen >= 1L, gapOrderCap >= 0L,
            minsupportCandCap >= 1L, perprimeCap >= 1L, perprimeMinPrimes >= 2L,
            is.logical(certifyPoly), is.null(certifyPolyDeg) || certifyPolyDeg >= 1L,
            is.numeric(timeout), timeout > 0)
  structure(list(relevanceCap       = as.integer(relevanceCap),
                 relevanceCapDir    = as.integer(relevanceCapDir),
                 relevanceCapSparse = as.integer(relevanceCapSparse),
                 degreeCap          = as.integer(degreeCap),
                 sampleSlack        = as.integer(sampleSlack),
                 probeRetries       = as.integer(probeRetries),
                 laurentDegNum      = as.integer(laurentDegNum),
                 laurentDegDen      = as.integer(laurentDegDen),
                 laurentCandCap     = as.integer(laurentCandCap),
                 termCap            = as.integer(termCap),
                 generalDegNum      = as.integer(generalDegNum),
                 generalDegDen      = as.integer(generalDegDen),
                 gapOrderCap        = as.integer(gapOrderCap),
                 minsupportCandCap  = as.integer(minsupportCandCap),
                 perprimeCap        = as.integer(perprimeCap),
                 perprimeMinPrimes  = as.integer(perprimeMinPrimes),
                 certifyPoly        = isTRUE(certifyPoly),
                 certifyPolyDeg     = if (is.null(certifyPolyDeg)) NULL
                                      else as.integer(certifyPolyDeg),
                 timeout            = timeout),
            class = c("reconstcontrol", "list"))
}


#' Settings of the polynomial engine
#'
#' Settings for `symmetryDetection(method = "polynomial", polynomial =
#' polynomialControl())`.
#'
#' @param ansatz Infinitesimal ansatz: `"uni"`, `"par"` or `"multi"`.
#' @param pMax Maximum degree of the ansatz, at least 1.
#' @param lieOrder Also require invariance of the Lie derivatives \eqn{L^k g} for
#'   \eqn{k = 1, \dots,} `lieOrder`.
#' @param exact Logical. Exact modular linear algebra (`TRUE`) or floating point.
#' @param verify Logical. Verify each generator symbolically.
#' @param allTrafos Logical. Keep transformations that share a common parameter
#'   factor.
#' @param backend `"symengine"` (falls back to sympy) or `"sympy"`.
#' @return A `polynomialControl` list.
#' @seealso [symmetryDetection()], [scalingControl()]
#' @export
polynomialControl <- function(ansatz    = c("uni", "par", "multi"),
                              pMax      = 2L,
                              lieOrder  = 0L,
                              exact     = TRUE,
                              verify    = TRUE,
                              allTrafos = FALSE,
                              backend   = c("symengine", "sympy")) {
  ansatz  <- match.arg(ansatz)
  backend <- match.arg(backend)
  stopifnot(pMax >= 1L, lieOrder >= 0L,
            is.logical(exact), is.logical(verify), is.logical(allTrafos))
  structure(list(ansatz = ansatz, pMax = as.integer(pMax),
                 lieOrder = as.integer(lieOrder), exact = isTRUE(exact),
                 verify = isTRUE(verify), allTrafos = isTRUE(allTrafos),
                 backend = backend),
            class = c("polynomialcontrol", "list"))
}


#' Settings of the scaling engine
#'
#' Settings for `symmetryDetection(method = "scaling", scaling = scalingControl())`.
#'
#' @param backend `"symengine"` (falls back to sympy) or `"sympy"`.
#' @return A `scalingControl` list.
#' @seealso [symmetryDetection()], [polynomialControl()]
#' @export
scalingControl <- function(backend = c("symengine", "sympy")) {
  backend <- match.arg(backend)
  structure(list(backend = backend),
            class = c("scalingcontrol", "list"))
}


# Sort in the C locale: the coordinate order fixes monomial tables, evaluation points
# and the gauge search, so a locale-dependent order could change the result.
.symSort <- function(x) sort(x, method = "radix")


# ---- modular linear algebra over GF(p) -----------------------------------------------

.symSieve <- function(n) {
  limit <- 200L
  repeat {
    is_p <- rep(TRUE, limit)
    is_p[1] <- FALSE
    for (i in 2:floor(sqrt(limit))) if (is_p[i]) is_p[seq(i * i, limit, by = i)] <- FALSE
    pr <- which(is_p)
    if (length(pr) >= n) return(pr[seq_len(n)])
    limit <- limit * 2L
  }
}


# A stream of distinct primes used as generic evaluation coordinates; grows on
# demand so the interpolation never runs out of sample points.
.symPool <- function() {
  cache <- .symSieve(1000L)
  function(k) {
    if (max(k) > length(cache)) cache <<- .symSieve(max(2L * length(cache), max(k)))
    cache[k]
  }
}


.symNullResidues <- function(res, freeCol, p) {
  nz <- res$dim
  piv <- res$pivots
  v <- integer(nz)
  v[freeCol + 1L] <- 1L
  for (ri in seq_along(piv))
    v[piv[ri] + 1L] <- as.integer((p - res$R[ri, freeCol + 1L]) %% p)
  v
}


# TRUE iff vector x lies in the column span of the nz-row matrix M over GF(p).
.symInSpan <- function(M, x, nz, p)
  ncol(M) > 0L && !is.null(symSolveMod(matrix(as.integer(M), nz), as.integer(x), p))


# a*b mod p < 2^31 in doubles: b is split at 15 bits so every product stays below
# 2^47 (exact to 2^53). Vectorised; b is a scalar or of the length of a.
.symMulmod <- function(a, b, p) {
  a <- a %% p; b <- b %% p
  hi <- (a * (b %/% 32768)) %% p
  ((hi * 32768) %% p + a * (b %% 32768)) %% p
}
# inverse by Fermat, computed with the split multiply so the squarings stay exact
.symInvmod <- function(a, p) {
  r <- 1; b <- a %% p; e <- p - 2
  while (e > 0) {
    if (e %% 2 == 1) r <- .symMulmod(r, b, p)
    e <- e %/% 2
    if (e > 0) b <- .symMulmod(b, b, p)
  }
  r
}

# RREF over GF(p) via symRrefMod(): reduced rows, 0-based pivots, rank. The zero rows
# the kernel drops are restored; .symMinsupportGauge searches them.
.symRrefModp <- function(M, p) {
  nr <- nrow(M); nc <- ncol(M)
  if (!nr || !nc)
    return(list(R = matrix(0, nr, nc), piv = integer(0), rank = 0L))
  rr <- symRrefMod(matrix(as.numeric(M) %% p, nr, nc), p)
  R <- matrix(0, nr, nc)
  if (rr$rank > 0L) R[seq_len(rr$rank), ] <- rr$R
  list(R = R, piv = rr$piv, rank = rr$rank)
}

# RREF of B preferring the columns `physCols0` (0-based) as pivots, so residual
# directions anchor on physical coordinates. Rows in the original column order.
.symPhysRref <- function(B, physCols0, nz, p) {
  if (nrow(B) == 0L) return(list(R = B, piv = integer(0)))
  ord <- c(sort(as.integer(physCols0)),
           sort(setdiff(0:(nz - 1L), as.integer(physCols0))))   # physical columns first
  rr <- .symRrefModp(B[, ord + 1L, drop = FALSE], p)
  Rback <- matrix(0, nrow(rr$R), nz)
  Rback[, ord + 1L] <- rr$R
  list(R = Rback, piv = ord[rr$piv + 1L])                        # 0-based original pivots
}

# k-by-k matrix inverse over GF(p) by Gauss-Jordan; NULL if singular.
.symMatinvModp <- function(B, p) {
  k <- nrow(B)
  A <- cbind(matrix(as.numeric(B) %% p, k, k), diag(k))
  for (col in seq_len(k)) {
    nz <- which(A[col:k, col] %% p != 0)
    if (!length(nz)) return(NULL)
    pr <- col + nz[1] - 1L
    if (pr != col) { tmp <- A[col, ]; A[col, ] <- A[pr, ]; A[pr, ] <- tmp }
    A[col, ] <- .symMulmod(A[col, ], .symInvmod(A[col, col], p), p)
    for (r in seq_len(k)) if (r != col) {
      f <- A[r, col] %% p
      if (f != 0) A[r, ] <- (A[r, ] - .symMulmod(A[col, ], f, p)) %% p
    }
  }
  A[, (k + 1L):(2L * k), drop = FALSE]
}

# reduce rows R modulo the RREF of the scaling lattice S, memoised per prime
.symReduceModRows <- function(S) {
  cache <- new.env(parent = emptyenv())
  function(R, p) {
    if (nrow(S) == 0L) return(R %% p)
    key <- as.character(p)
    Sr <- cache[[key]]
    if (is.null(Sr)) { Sr <- .symRrefModp(S, p); cache[[key]] <- Sr }
    for (j in seq_along(Sr$piv)) {
      fac <- R[, Sr$piv[j] + 1L] %% p
      for (i in seq_len(nrow(R)))
        if (fac[i] != 0) R[i, ] <- (R[i, ] - .symMulmod(Sr$R[j, ], fac[i], p)) %% p
    }
    R
  }
}

# Shared tail of the canon/logcoord gauges: `rowsFn(rp, p, zvals)` maps a kernel
# result to the k reduced direction rows (or NULL). Pin the k pivot columns of the
# reference rows to the identity, so each direction's entry reads off directly.
.symGaugeFromRows <- function(residualFree, rowsFn, refRows, P, k, nz) {
  raw <- list(anchors = residualFree, residueFns = vector("list", k))
  if (is.null(refRows)) return(raw)
  rr <- .symRrefModp(refRows, P)
  if (rr$rank < k) return(raw)
  dp <- rr$piv
  makeFn <- function(i) function(rp, p, zvals = NULL) {
    R <- rowsFn(rp, p, zvals)
    if (is.null(R)) return(NULL)
    Binv <- .symMatinvModp(R[, dp + 1L, drop = FALSE], p)
    if (is.null(Binv)) return(NULL)
    w <- numeric(nz)
    for (l in seq_len(k)) w <- (w + .symMulmod(R[l, ], Binv[i, l], p)) %% p
    as.integer(w)
  }
  list(anchors = dp, residueFns = lapply(seq_len(k), makeFn))
}


# Decouple the residual directions the free-column gauge entangles: reduce modulo the
# integer scaling span, then direction i is 1 on its own pivot and 0 on the others'.
# Returns (anchors, residue functions), or the free-column gauge if inseparable.
.symCanonGauge <- function(residualFree, scalRows, P, nz, sc) {
  k <- length(residualFree)
  if (k <= 1L) return(list(anchors = residualFree, residueFns = vector("list", k)))
  reduceScal <- .symReduceModRows(scalRows)
  resid <- function(rp, p, zvals = NULL) reduceScal(
    t(vapply(residualFree, function(fc) .symNullResidues(rp, fc, p),
             integer(nz))), p)
  .symGaugeFromRows(residualFree, resid, resid(sc$ref, P), P, k, nz)
}


# .symCanonGauge in log coordinates: a scaling xi_i = c_i * z_i has the low-degree
# log-residue eta_i = xi_i / z_i, so normalising adds no denominator and the fit stays
# sparse. .symLogcoordBacksub maps eta back to xi.
.symLogcoordGauge <- function(residualFree, scalRows, P, nz, sc, zvals0) {
  k <- length(residualFree)
  if (k == 0L) return(list(anchors = residualFree, residueFns = vector("list", k)))

  # integer scaling weights (the log-residues of the scalings) by a centred lift at the
  # base point; reducing modulo them strips scaling admixture from the eta-rows
  W <- if (nrow(scalRows) == 0L) matrix(0L, 0L, nz) else {
    Wm <- matrix(0L, nrow(scalRows), nz)
    for (c in seq_len(nz)) {
      zc <- as.numeric(zvals0[c]) %% P
      if (zc == 0) next
      wc <- .symMulmod(scalRows[, c], .symInvmod(zc, P), P)
      Wm[, c] <- as.integer(ifelse(wc > P / 2, wc - P, wc))
    }
    Wm
  }
  reduceScal <- .symReduceModRows(W)

  etaRows <- function(rp, p, zvals) {
    R <- t(vapply(residualFree, function(fc) .symNullResidues(rp, fc, p),
                  integer(nz)))
    for (c in seq_len(nz)) {
      zc <- as.numeric(zvals[c]) %% p
      if (zc == 0) { if (any(R[, c] %% p != 0)) return(NULL); next }
      R[, c] <- .symMulmod(R[, c], .symInvmod(zc, p), p)
    }
    reduceScal(R, p)
  }
  .symGaugeFromRows(residualFree, etaRows, etaRows(sc$ref, P, zvals0), P, k, nz)
}


# log-coordinate direction back to original coordinates: xi_c = eta_c * z_c
.symLogcoordBacksub <- function(vector, spy) {
  out <- list()
  for (nm in names(vector)) {
    expr <- paste0("(", vector[[nm]], ")*", nm)
    out[[nm]] <- if (is.null(spy)) expr else .symSimplify(expr, spy)
  }
  out
}


# Decouple the residual directions by minimal support, so a scaling with a parameter
# weight (xi_kinh = -nhill * kinh) stays a sparse cocircuit instead of the dense lift of
# .symCanonGauge. Enumerates small column subsets, skipping integer scalings. Same
# return shape as .symCanonGauge.
.symMinsupportGauge <- function(residualFree, scalRows, P, nz, sc, freeCols,
                                  supportCap = 6L, candCap = 20000L, maxSecs = 20) {
  k <- length(residualFree)
  raw <- list(anchors = residualFree, residueFns = vector("list", k))
  if (k == 0L) return(raw)
  deadline <- Sys.time() + maxSecs

  basisRows <- function(rp, p)
    t(vapply(freeCols, function(fc) .symNullResidues(rp, fc, p), integer(nz)))
  rowComb <- function(c, B, p) {
    w <- numeric(nz)
    for (l in seq_along(c)) if (c[l] %% p != 0)
      w <- (w + .symMulmod(B[l, ], c[l], p)) %% p
    as.integer(w)
  }
  # combination of the rows of B vanishing on columns `keep`; NULL at full row rank
  leftNull <- function(B, keep, p) {
    nr <- nrow(B)
    M <- B[, keep, drop = FALSE]
    rr <- .symRrefModp(cbind(matrix(as.numeric(M) %% p, nr), diag(nr)), p)
    z <- which(rowSums(rr$R[, seq_len(ncol(M)), drop = FALSE] %% p != 0) == 0)
    if (!length(z)) return(NULL)
    rr$R[z[1], ncol(M) + seq_len(nr)] %% p
  }
  # the row-space vector supported within columns S (1-based), or NULL
  supported <- function(B, S, p) {
    c <- leftNull(B, setdiff(seq_len(nz), S), p)
    if (is.null(c)) return(NULL)
    rowComb(c, B, p)
  }
  isScaling <- function(v) .symInSpan(t(scalRows), v, nz, P)

  B0 <- basisRows(sc$ref, P)
  # candidate columns: where any residual direction's free-column residue lives
  cols <- sort(unique(unlist(lapply(residualFree,
    function(fc) which(B0[match(fc, freeCols), ] %% P != 0)))))
  if (length(cols) < 2L) return(raw)

  # A wide direction has no small cocircuit, so the subsets are iterated in place under
  # a global budget; once spent, the rest fall through to the free-column fit.
  nCols <- length(cols)
  budget <- as.integer(candCap)
  nextCombo <- function(idx, n, s) {
    i <- s
    while (i >= 1L && idx[i] == n - s + i) i <- i - 1L
    if (i < 1L) return(NULL)
    idx[i] <- idx[i] + 1L
    j <- i + 1L
    while (j <= s) { idx[j] <- idx[j - 1L] + 1L; j <- j + 1L }
    idx
  }
  found <- list(); sel <- matrix(0L, nz, 0L); iter <- 0L
  for (s in 2:min(nCols, supportCap)) {
    if (length(found) >= k || budget <= 0L) break
    idx <- seq_len(s)
    repeat {
      if (length(found) >= k || budget <= 0L) break
      # wall-clock cap as well; the candidate budget alone can take minutes
      iter <- iter + 1L
      if (iter %% 256L == 0L && Sys.time() > deadline) { budget <- 0L; break }
      S <- cols[idx]
      idx <- nextCombo(idx, nCols, s)
      budget <- budget - 1L
      spanned <- any(vapply(found, function(fc) all(fc$supp %in% S), logical(1)))
      if (!spanned) {
        v <- supported(B0, S, P)
        if (!is.null(v) && !isScaling(v)) {
          base <- cbind(if (nrow(scalRows)) t(scalRows) else matrix(0L, nz, 0L), sel)
          if (!.symInSpan(base, v, nz, P)) {
            supp <- which(v %% P != 0)
            found[[length(found) + 1L]] <- list(supp = supp, anchor = supp[1] - 1L, S = S)
            sel <- cbind(sel, v)
          }
        }
      }
      if (is.null(idx)) break
    }
  }
  # possibly fewer than k; the caller sends the rest to the free-column fit
  if (!length(found)) return(raw)

  fns <- lapply(found, function(fc) {
    S <- fc$S; anchor <- fc$anchor
    # log-residue eta_c = xi_c / z_c, 1 on the anchor: a parameter-weighted scaling
    # then has a constant entry (eta = -nhill) that depends on its weight only
    fn <- function(rp, p, zvals = NULL) {
      v <- supported(basisRows(rp, p), S, p)
      if (is.null(v) || v[anchor + 1L] %% p == 0) return(NULL)
      v <- .symMulmod(v, .symInvmod(v[anchor + 1L] %% p, p), p)
      if (is.null(zvals)) return(as.integer(v))
      za <- as.numeric(zvals[anchor + 1L]) %% p
      if (za == 0) return(NULL)
      eta <- integer(nz)
      for (c in S) {
        zc <- as.numeric(zvals[c]) %% p
        if (zc == 0) { if (v[c] %% p != 0) return(NULL) else next }
        eta[c] <- as.integer(.symMulmod(v[c], .symMulmod(za, .symInvmod(zc, p), p), p))
      }
      eta
    }
    # support fixed to S: a NULL residue marks the sample unusable, not the leaf
    # relevant; the result is certified at a fresh point
    attr(fn, "pinnedSupport") <- TRUE
    fn
  })
  list(anchors = vapply(found, function(fc) fc$anchor, integer(1)),
       residueFns = fns, vectors = sel)
}


# An integer given as a decimal string, modulo p.
.symBigMod <- function(s, p) {
  neg <- startsWith(s, "-")
  d <- sub("^-", "", s)
  r <- if (nchar(d) <= 15L) as.numeric(d) %% p else {
    acc <- 0
    for (ch in strsplit(d, "")[[1]]) acc <- (.symMulmod(acc, 10, p) + as.numeric(ch)) %% p
    acc
  }
  if (neg) (p - r) %% p else r
}

# Evaluator of a straight-line tape (op codes 0 const, 1 add, 2 mul, 3 inv) over the
# leaves, modulo p. Returns function(point, p), NULL where an inverse hits zero.
.symTapeFn <- function(tp, nLeaves) {
  op <- as.integer(tp$op); a <- as.integer(tp$a) + 1L; b <- as.integer(tp$b) + 1L
  cnum <- as.character(tp$cnum); cden <- as.character(tp$cden)
  out <- as.integer(tp$out) + 1L
  isC <- which(op == 0L)
  cpos <- integer(length(op)); cpos[isC] <- seq_along(isC)
  cache <- list()
  function(point, p) {
    key <- as.character(p)
    cv <- cache[[key]]
    if (is.null(cv)) {
      cv <- vapply(isC, function(i) {
        d <- .symBigMod(cden[i], p)
        if (d == 0) NA_real_ else .symMulmod(.symBigMod(cnum[i], p), .symInvmod(d, p), p)
      }, numeric(1))
      cache[[key]] <<- cv
    }
    if (anyNA(cv)) return(NULL)
    v <- numeric(nLeaves + length(op))
    v[seq_len(nLeaves)] <- as.numeric(point[seq_len(nLeaves)]) %% p
    for (i in seq_along(op)) {
      s <- nLeaves + i
      v[s] <- switch(op[i] + 1L,
        cv[cpos[i]],
        (v[a[i]] + v[b[i]]) %% p,
        .symMulmod(v[a[i]], v[b[i]], p),
        if (v[a[i]] == 0) return(NULL) else .symInvmod(v[a[i]], p))
    }
    v[out]
  }
}

# Relation rows e_W - W log(r) grad(tau/L) of the exponential leaves at a point.
.symExpRows <- function(rel, relFn, znames, point, p) {
  vals <- relFn(point, p)
  if (is.null(vals)) return(NULL)
  Ws <- unique(as.character(rel$W))
  M <- matrix(0, length(Ws), length(znames))
  M[cbind(seq_along(Ws), match(Ws, znames))] <- 1
  zc <- match(as.character(rel$z), znames)
  ok <- !is.na(zc)
  wr <- match(as.character(rel$W), Ws)
  for (k in which(ok))
    M[wr[k], zc[k]] <- (M[wr[k], zc[k]] - vals[k]) %% p
  M
}

# Components in the exponential leaves back in the model symbols.
.symExpBacksub <- function(vector, back, sd) {
  nm <- as.list(as.character(back$names)); vl <- as.list(as.character(back$values))
  lapply(vector, function(x) {
    out <- tryCatch(sd$expBacksub(as.character(x), nm, vl),
                    error = function(e) as.character(x))
    if (length(out) != 1L || is.na(out)) as.character(x) else out
  })
}


# ---- symbolic reconstruction: monomials, per-prime lift, recast back-substitution ----

# exponent vectors of `nvar` variables with total degree <= `degree`, by degree;
# enumerated directly, not filtered from the (degree+1)^nvar grid
.symMonoTable <- function(nvar, degree) {
  if (nvar == 0L) return(matrix(0L, 1L, 0L))
  gen <- function(n, d) {
    if (n == 1L) return(lapply(0:d, function(e) e))
    unlist(lapply(0:d, function(e) lapply(gen(n - 1L, d - e), function(s) c(e, s))),
           recursive = FALSE)
  }
  m <- do.call(rbind, gen(nvar, degree))
  m <- m[order(rowSums(m)), , drop = FALSE]
  storage.mode(m) <- "integer"
  dimnames(m) <- NULL
  m
}


.symMonoString <- function(expo, vars) {
  terms <- character(0)
  for (k in seq_along(vars)) if (expo[k] != 0)
    terms <- c(terms, if (expo[k] == 1) vars[k] else paste0(vars[k], "^", expo[k]))
  if (!length(terms)) "1" else paste(terms, collapse = "*")
}


.symPolyString <- function(numc, denc, mons, vars) {
  parts <- character(0)
  for (j in seq_len(nrow(mons))) {
    n <- numc[j]; d <- denc[j]
    if (n == "0") next
    coef <- if (d == "1") n else paste0("(", n, "/", d, ")")
    ms <- .symMonoString(mons[j, ], vars)
    parts <- c(parts, if (ms == "1") coef else paste0(coef, "*", ms))
  }
  if (!length(parts)) "0" else paste(parts, collapse = " + ")
}


# One nullspace entry as a rational function of the relevant variables. The per-prime
# fits must share their free column (gauge). Returns coefficient strings, or NULL if
# no closed form of bounded degree fits.
.symReconstructEntry <- function(sampleU, mons, residues, primes) {
  nMon <- nrow(mons)
  coefRes <- matrix(0L, 2L * nMon, length(primes))
  freeCol <- NULL
  for (j in seq_along(primes)) {
    fit <- symFitRational(sampleU, mons, as.integer(residues[, j]), primes[j])
    if (!identical(fit$status, "ok")) return(NULL)
    if (is.null(freeCol)) freeCol <- fit$freeCol
    else if (fit$freeCol != freeCol) return(NULL)
    coefRes[, j] <- fit$coeffs
  }
  rec <- symRatRecon(coefRes, as.integer(primes))
  if (any(rec$den == "0")) return(NULL)
  num <- seq_len(nMon); den <- nMon + seq_len(nMon)
  list(numCoefN = rec$num[num], numCoefD = rec$den[num],
       denCoefN = rec$num[den], denCoefD = rec$den[den])
}


.symSimplify <- function(expr, spy) {
  out <- tryCatch(as.character(spy$cancel(spy$sympify(expr))),
                  error = function(e) expr)
  if (length(out) != 1L || is.na(out)) expr else out
}


# recast coordinates back in a direction, in Python: E -> base^exp, L -> log(base)
.symRecastBacksub <- function(vector, recast, sd) {
  eN <- as.list(vapply(recast, function(r) as.character(r$E), ""))
  lN <- as.list(vapply(recast, function(r) as.character(r$L), ""))
  bN <- as.list(vapply(recast, function(r) as.character(r$base), ""))
  xN <- as.list(vapply(recast, function(r) as.character(r$exp), ""))
  lapply(vector, function(x) {
    out <- tryCatch(sd$recastBacksub(as.character(x), eN, lN, bN, xN),
                    error = function(e) as.character(x))
    if (length(out) != 1L || is.na(out)) as.character(x) else out
  })
}


# Forward-sampling reconstruction of a residual direction that depends on the resting
# state x*(theta), where each prime's steady-state slice differs and the backward CRT
# is inconsistent. Solving f = 0 linearly for a turnover subset of rates holds at every
# prime, so the resting states become shared sample coordinates (symbols in the report).
# Returns a verified closed form or a support-only fallback.
.symPerprimeForward <- function(f, sc, kcall, kcallFwd, znames, zSlots, leafNames, nz,
                                  scaling, stateColNames, paramNames, recast, sd, spy, ctrl,
                                  physCols, models, realStateNames, solveParamNames, solveHeld,
                                  acIn = NULL, listAnchors = FALSE) {
  # three 31-bit primes (~2^93) suffice for small coefficients; each costs a kernel pass
  rp8 <- unique(as.integer(c(.symPrimes, 2147483563, 2147483549, 2147483543, 2147483497)))
  nPrimeFwd <- max(2L, min(3L, length(rp8) - 1L))
  primes <- rp8[seq_len(nPrimeFwd)]; qv <- rp8[length(primes) + 1L]; P1 <- primes[1]
  dirSupport <- NULL   # this direction's support, set once its anchor is known
  fb <- function(reason) {
    supp <- if (!is.null(dirSupport)) dirSupport
            else { v <- .symNullResidues(sc$ref, f, P1); .symSort(znames[v != 0]) }
    list(support = supp, type = "general", closedForm = FALSE, reason = reason) }
  if (is.null(sd) || is.null(kcallFwd) || .symExpired(ctrl)) {
    if (listAnchors) return(integer(0)); return(fb("forward path unavailable")) }
  rdiag <- nzchar(Sys.getenv("DMOD_SYM_ROBUSTDIAG"))
  stripName <- function(nm) sub("\\|c[0-9]+$", "", nm)
  npt <- length(sc$point0); poolN <- sc$poolNext
  draw <- function(n) { u <- sc$pool(poolN + seq_len(n) - 1L); poolN <<- poolN + n; u }
  freeColsOf <- function(rp) setdiff(0:(nz - 1L), as.integer(rp$pivots))
  isStateCol <- logical(nz); isStateCol[which(znames %in% stateColNames)] <- TRUE
  physSet <- if (is.null(physCols)) (seq_len(nz) - 1L) else as.integer(physCols)
  nS <- length(scaling)
  Wn <- lapply(scaling, function(s) names(s$vector)); Wv <- lapply(scaling, function(s) as.character(unlist(s$vector)))
  tangentsAt <- function(zvals, p) {
    env <- as.list(setNames(as.numeric(zvals), znames)); M <- matrix(0L, nS, nz); drp <- logical(nS)
    for (j in seq_len(nS)) { cols <- match(Wn[[j]], znames)
      for (t in seq_along(cols)) { cc <- cols[t]; if (is.na(cc)) next
        w <- suppressWarnings(as.numeric(Wv[[j]][t]))
        if (is.na(w)) { wv <- .symEvalModq(Wv[[j]][t], env, p, sd)
          if (is.null(wv) || is.na(wv)) { drp[j] <- TRUE; break }; w <- as.numeric(wv) }
        zc <- if (isStateCol[cc]) 1 else (as.numeric(zvals[cc]) %% p)
        M[j, cc] <- as.integer((w %% p * zc) %% p) } }
    M[!drp, , drop = FALSE] }
  reduceRows <- function(B, S, p) { if (nrow(S) == 0L) return(B %% p)
    Sr <- .symRrefModp(S, p)
    for (j in seq_along(Sr$piv)) { fac <- B[, Sr$piv[j] + 1L] %% p
      for (i in seq_len(nrow(B))) if (fac[i] != 0) B[i, ] <- (B[i, ] - .symMulmod(Sr$R[j, ], fac[i], p)) %% p }
    B }
  # the residual direction with ac = 1 and every physical pivot 0, a gauge consistent
  # across points: nullspace mod scalings, physical-first RREF. Defined whether or not
  # `ac` is a pivot here (the forward pivots can differ from the backward ones).
  extract <- function(rp, zvals, p, ac) {
    B <- t(vapply(freeColsOf(rp), function(fc) .symNullResidues(rp, fc, p), integer(nz)))
    Bred <- reduceRows(B, tangentsAt(zvals, p), p)
    Bred <- Bred[apply(Bred, 1L, function(x) any(x %% p != 0)), , drop = FALSE]
    if (nrow(Bred) == 0L) return(NULL)
    rr <- .symPhysRref(Bred, physSet, nz, p)
    pr <- match(ac, rr$piv)
    if (!is.na(pr)) return(as.integer(rr$R[pr, ] %% p))          # ac is a pivot: its row
    if ((ac + 1L) > nz) return(NULL)
    v <- integer(nz); v[ac + 1L] <- 1L                          # ac free: 1 on ac, 0 on pivots
    for (i in seq_along(rr$piv)) v[rr$piv[i] + 1L] <- as.integer((p - rr$R[i, ac + 1L]) %% p)
    v }
  kxF <- function(pt, p, ac, sr) { rp <- kcallFwd(pt, p, sc$NtUsed, sr)
    if (!isTRUE(rp$ok) || rp$rank != sc$rank) return(NULL); extract(rp, pt[zSlots + 1L], p, ac) }

  # 1. anchors: physical-first RREF of the residual null vectors (mod scalings) at a
  # backward-valid point gives each direction a distinct physical pivot. `listAnchors`
  # returns them all; otherwise reconstruct the one at `acIn` (or the first).
  ptbB <- NULL
  for (att in seq_len(200L)) { cand <- draw(npt)
    r <- kcall(cand, P1, sc$NtUsed); if (isTRUE(r$ok) && r$rank == sc$rank) { ptbB <- cand; break } }
  if (is.null(ptbB)) { if (listAnchors) return(integer(0)); return(fb("forward: no backward base for support")) }
  rp0 <- kcall(ptbB, P1, sc$NtUsed)
  B0 <- t(vapply(freeColsOf(rp0), function(fc) .symNullResidues(rp0, fc, P1), integer(nz)))
  Bred0 <- reduceRows(B0, tangentsAt(ptbB[zSlots + 1L], P1), P1)
  Bred0 <- Bred0[apply(Bred0, 1L, function(x) any(x %% P1 != 0)), , drop = FALSE]
  if (nrow(Bred0) == 0L) { if (listAnchors) return(integer(0)); return(fb("forward: empty residual")) }
  rr0 <- .symPhysRref(Bred0, physSet, nz, P1)
  anchors <- rr0$piv[rr0$piv %in% physSet]
  if (listAnchors) return(anchors)
  if (!length(anchors)) return(fb("forward: no physical residual anchor"))
  ac <- if (!is.null(acIn)) as.integer(acIn) else anchors[1]
  arow <- match(ac, rr0$piv)
  if (is.na(arow)) return(fb("forward: requested anchor is not a residual pivot"))
  v0 <- as.integer(.symMulmod(rr0$R[arow, ], .symInvmod(rr0$R[arow, ac + 1L] %% P1, P1), P1))
  dirSupport <- .symSort(unique(stripName(znames[which(v0 %% P1 != 0)])))
  suppNames <- unique(stripName(znames[c(ac, setdiff(which(v0 %% P1 != 0) - 1L, ac)) + 1L]))

  # 2. turnover transversal (avoids the direction's support + the recast coords)
  keepFree <- unique(c(suppNames, grep("^_E_|^_L_", leafNames, value = TRUE)))
  tsv <- setNames(as.list(((seq_along(realStateNames) * 7919L + 11L) %% (P1 - 1L)) + 1L), realStateNames)
  tpv <- setNames(as.list(((seq_along(solveParamNames) * 104729L + 3L) %% (P1 - 1L)) + 1L), solveParamNames)
  ts <- tryCatch(sd$solveForwardModular(models[[1]], realStateNames, solveParamNames, tsv, tpv, P1,
                   forcings = if (length(solveHeld)) solveHeld else NULL,
                   keepFree = as.list(keepFree)), error = function(e) NULL)
  if (is.null(ts) || !isTRUE(ts$ok)) return(fb(paste("forward: no turnover transversal;",
                                                     if (is.null(ts)) "solve error" else ts$why)))
  solveRates <- as.character(ts$solveRates)
  if (rdiag) message(sprintf("[fwd] anchor=%s support=%d transversal(%d)={%s}",
                             znames[ac + 1L], length(suppNames), length(solveRates),
                             paste(solveRates, collapse = ",")))

  # 3. a forward base valid at every prime
  base0 <- NULL
  for (att in seq_len(80L)) { cand <- draw(npt)
    if (all(vapply(primes, function(pp) !is.null(kxF(cand, pp, ac, solveRates)), logical(1)))) { base0 <- cand; break } }
  if (is.null(base0)) return(fb("forward: no base valid at all primes"))

  # 4. relevance over every free leaf but the solved rates (parameters, _L_, _E_, states);
  # a missed leaf would make an entry base0-specific (step 6 catches that)
  cand <- setdiff(seq_along(leafNames), match(solveRates, leafNames))
  cand <- cand[!is.na(cand)]
  v0f <- kxF(base0, P1, ac, solveRates)
  physSupp <- setdiff(intersect(which(v0f %% P1 != 0) - 1L, physSet), ac)
  relBy <- replicate(length(physSupp), integer(0), simplify = FALSE)
  for (li in cand) { if (.symExpired(ctrl)) return(fb("forward: timeout in relevance"))
    ch <- rep(FALSE, length(physSupp))
    for (dv in c(3L, 7L, 11L)) { pt <- as.numeric(base0); pt[li] <- pt[li] + dv
      v <- kxF(pt, P1, ac, solveRates); if (is.null(v)) next
      ch <- ch | ((v[physSupp + 1L] - v0f[physSupp + 1L]) %% P1 != 0) }
    for (i in which(ch)) relBy[[i]] <- c(relBy[[i]], li) }
  unionRel <- sort(unique(unlist(relBy)))
  if (!length(unionRel)) return(fb("forward: residual constant in every leaf"))
  maxRel <- max(vapply(relBy, length, integer(1)))
  if (rdiag) message(sprintf("[fwd] physSupp=%d unionRel=%d maxRel=%d", length(physSupp), length(unionRel), maxRel))
  if (nzchar(Sys.getenv("DMOD_SYM_FWDREL"))) {
    for (i in seq_along(physSupp))
      message(sprintf("[fwdrel] %-26s : %s", stripName(znames[physSupp[i] + 1L]),
                      paste(leafNames[relBy[[i]]], collapse = ", ")))
    return(fb("fwdrel diagnostic"))
  }
  # the forward path frees the states and both recast partners, so the direction cap is
  # relaxed; each entry is still bounded by relevanceCapSparse
  fwdCapDir <- max(as.integer(ctrl$relevanceCapDir), 48L)
  if (length(unionRel) > fwdCapDir || maxRel > ctrl$relevanceCapSparse)
    return(fb(sprintf("forward: couples %d leaves (entry up to %d)", length(unionRel), maxRel)))

  # 5. one sample bank for all primes, topped up per fit degree on demand (each point
  # costs a kernel pass, so low-degree directions stay cheap)
  dCap <- max(1L, min(3L, as.integer(ctrl$degreeCap)))
  needOf <- function(k, d) 2L * nrow(.symMonoTable(k, d)) + 20L
  bankU <- matrix(0L, 0L, length(unionRel)); bankV <- lapply(primes, function(.) matrix(0L, 0L, nz))
  fillBankTo <- function(n) {
    tries <- 0L
    while (nrow(bankU) < n && tries < 12L * n) { if (.symExpired(ctrl)) break
      tries <- tries + 1L; pt <- as.numeric(base0); pt[unionRel] <- draw(length(unionRel))
      vs <- lapply(primes, function(pp) kxF(pt, pp, ac, solveRates))
      if (any(vapply(vs, is.null, logical(1)))) next
      bankU <<- rbind(bankU, pt[unionRel])
      for (jp in seq_along(primes)) bankV[[jp]] <<- rbind(bankV[[jp]], vs[[jp]]) }
    if (rdiag && nrow(bankU) >= n) message(sprintf("[fwd] bank at %d points (all %d primes)", nrow(bankU), length(primes)))
    nrow(bankU) >= n }

  recEntry <- function(reli, c9) {
    reliCols <- match(reli, unionRel); vars <- leafNames[reli]
    for (d in 0:dCap) {
      mons <- .symMonoTable(length(reli), d); nMon <- nrow(mons)
      if (!fillBankTo(needOf(length(reli), d))) next
      need <- nrow(bankU)
      refFree <- NULL; coefRes <- matrix(0L, 2L * nMon, length(primes)); ok <- TRUE
      for (jp in seq_along(primes)) { pp <- primes[jp]
        sU <- matrix(as.integer(bankU[seq_len(need), reliCols, drop = FALSE]), need)
        rv <- as.integer(bankV[[jp]][seq_len(need), c9 + 1L])
        fit <- symFitRational(sU, matrix(as.integer(mons), nMon), rv, pp)
        if (!identical(fit$status, "ok")) { ok <- FALSE; break }
        raw <- as.numeric(fit$coeffs); if (is.null(refFree)) refFree <- fit$freeCol
        dn <- raw[refFree + 1L] %% pp; if (dn == 0) { ok <- FALSE; break }
        coefRes[, jp] <- as.integer(.symMulmod(raw, .symInvmod(dn, pp), pp)) }
      if (!ok) next
      rec <- tryCatch(sd$symRatReconBig(coefRes, as.integer(primes)), error = function(e) NULL)
      if (is.null(rec) || any(rec$den == "0")) next
      numI <- seq_len(nMon); denI <- nMon + seq_len(nMon)
      numStr <- .symPolyString(rec$num[numI], rec$den[numI], mons, vars)
      denStr <- .symPolyString(rec$num[denI], rec$den[denI], mons, vars)
      return(list(expr = if (denStr == "1") numStr else paste0("(", numStr, ")/(", denStr, ")"), reli = vars)) }
    NULL }

  entries <- list(); entryReli <- list()
  entries[[znames[ac + 1L]]] <- "1"; entryReli[[znames[ac + 1L]]] <- character(0)
  for (i in seq_along(physSupp)) { if (.symExpired(ctrl)) return(fb("forward: timeout in fitting"))
    c9 <- physSupp[i]; nm <- znames[c9 + 1L]; reli <- relBy[[i]]
    if (!length(reli)) {
      if (!fillBankTo(1L)) return(fb("forward: empty bank"))
      resc <- vapply(seq_along(primes), function(jp) bankV[[jp]][1, c9 + 1L], integer(1))
      rc <- tryCatch(sd$symRatReconBig(matrix(as.integer(resc), 1L), as.integer(primes)), error = function(e) NULL)
      if (is.null(rc) || rc$den[1] == "0") return(fb("forward: a constant entry could not be lifted"))
      entries[[nm]] <- if (rc$den[1] == "1") rc$num[1] else paste0(rc$num[1], "/", rc$den[1])
      entryReli[[nm]] <- character(0); next }
    rr <- recEntry(reli, c9)
    if (rdiag) message(sprintf("[fwd]   entry %-24s nvar=%d -> %s", stripName(nm), length(reli),
                               if (is.null(rr)) sprintf("NO FIT (deg>%d)", dCap) else "ok"))
    if (is.null(rr)) return(fb(sprintf("forward: entry %s not a bounded rational (deg>%d)", nm, dCap)))
    entries[[nm]] <- if (is.null(spy)) rr$expr else .symSimplify(rr$expr, spy); entryReli[[nm]] <- rr$reli }

  # 6. verify at a fresh prime and a point with every leaf redrawn, so the form is
  # certified independent of base0 and of the leaves judged irrelevant
  ptv <- NULL
  for (t in seq_len(400L)) { pt <- draw(npt)
    if (!is.null(kxF(pt, qv, ac, solveRates))) { ptv <- pt; break } }
  if (is.null(ptv)) return(fb("forward: no fresh verification point"))
  vver <- kxF(ptv, qv, ac, solveRates)
  for (nm in names(entries)) { col <- match(nm, znames); rl <- entryReli[[nm]]
    pred <- tryCatch(sd$evalRationalMod(entries[[nm]], as.list(rl),
              as.list(as.integer(ptv[match(rl, leafNames)])), qv), error = function(e) NULL)
    okv <- !is.null(pred) && ((as.integer(pred) - vver[col]) %% qv == 0)
    if (rdiag) message(sprintf("[fwd]   verify %-24s -> %s", stripName(nm), if (okv) "ok" else "FAIL"))
    if (!okv) return(fb("forward: reconstructed direction failed fresh-point verification")) }
  if (rdiag) message("[fwd] VERIFIED (base-independent): closing direction (states are resting-level symbols)")
  list(support = .symSort(names(entries)), vector = entries, type = "general", closedForm = TRUE)
}


# Every residual direction by the forward path: one anchor set, then one
# .symPerprimeForward() per anchor. Returns the directions (closed or support-only), or
# NULL without anchors.
.symPerprimeForwardMulti <- function(residualFree, sc, kcall, kcallFwd, znames, zSlots,
                                        leafNames, nz, scaling, stateColNames, paramNames,
                                        recast, sd, spy, ctrl, physCols, models,
                                        realStateNames, solveParamNames, solveHeld) {
  if (!length(residualFree)) return(list())
  # the caller adopts the set only if every direction closes, replacing both the peel
  # and the per-column results
  fwd1 <- function(acIn, listAnchors)
    .symPerprimeForward(residualFree[1], sc, kcall, kcallFwd, znames, zSlots, leafNames,
                          nz, scaling, stateColNames, paramNames, recast, sd, spy, ctrl,
                          physCols, models, realStateNames, solveParamNames, solveHeld,
                          acIn = acIn, listAnchors = listAnchors)
  anchors <- fwd1(NULL, TRUE)
  if (!length(anchors)) return(NULL)
  lapply(anchors, function(a) fwd1(a, FALSE))
}


# rational expression at integer `env` modulo q, exact in Python; NA when the
# denominator vanishes mod q
.symEvalModq <- function(expr, env, q, sd) {
  v <- tryCatch(sd$evalRationalMod(as.character(expr), as.list(names(env)),
                                   as.list(as.numeric(unlist(env))),
                                   as.integer(q)),
                error = function(err) NULL)
  if (is.null(v)) NA_integer_ else as.integer(v)
}


# ---- verification: nullspace membership, fresh-prime re-check, saturation guard ------

# Check a reconstructed direction (free column `f`) against the null vector at the
# base point modulo `.symVerifyPrime`. FALSE only on a definite mismatch.
.symVerifyDirection <- function(entry, f, znames, leafNames, point0, NtUsed,
                                  kcall, sd, residueFn = NULL) {
  if (is.null(sd) || is.null(entry$vector)) return(TRUE)
  q <- .symVerifyPrime
  rq <- kcall(point0, q, as.integer(NtUsed))
  if (!isTRUE(rq$ok)) return(TRUE)
  aq <- if (is.null(residueFn)) .symNullResidues(rq, f, q) else residueFn(rq, q, NULL)
  if (is.null(aq)) return(TRUE)
  env <- setNames(as.list(as.numeric(point0[seq_along(leafNames)])), leafNames)
  for (nm in names(entry$vector)) {
    col <- match(nm, znames)
    if (is.na(col)) next
    got <- .symEvalModq(entry$vector[[nm]], env, q, sd)
    if (is.na(got)) return(TRUE)
    if (as.integer(got) %% q != as.integer(aq[col]) %% q) return(FALSE)
  }
  TRUE
}


# Gauge-independent check that a direction lies in the nullspace at a fresh point, which
# catches base-point values baked into a representative. Inconclusive is TRUE.
.symVerifyInNullspace <- function(entry, f, znames, leafNames, point0, NtUsed,
                                     kcall, pool, poolNext, nz, sd) {
  if (is.null(sd) || is.null(entry$vector)) return(TRUE)
  q <- .symVerifyPrime
  pt <- as.numeric(point0)
  pt[] <- pool(poolNext + seq_along(pt) - 1L)
  rq <- kcall(pt, q, as.integer(NtUsed))
  if (!isTRUE(rq$ok)) return(TRUE)
  env <- setNames(as.list(pt[seq_along(leafNames)]), leafNames)
  v <- integer(nz); v[f + 1L] <- 1L
  for (nm in names(entry$vector)) {
    col <- match(nm, znames)
    if (is.na(col)) next
    got <- .symEvalModq(entry$vector[[nm]], env, q, sd)
    if (is.na(got)) return(TRUE)
    v[col] <- as.integer(got %% q)
  }
  freeCols <- setdiff(seq_len(nz) - 1L, as.integer(rq$pivots))
  recon <- numeric(nz)
  for (fcol in freeCols) {
    vf <- v[fcol + 1L] %% q
    if (vf != 0) recon <- (recon + .symMulmod(.symNullResidues(rq, fcol, q), vf, q)) %% q
  }
  all((recon - v) %% q == 0)
}


# Strict check for the per-prime (coupled steady-state) path, where a random point rarely
# has a modular steady state and the lenient checks would pass a spurious fit. Retries
# fresh points until one solves with the same pivots; an exhausted budget fails.
.symVerifyPerprime <- function(entry, f, znames, leafNames, point0, NtUsed,
                                 kcall, pivots, pool, poolNext, nz, sd,
                                 relLeaves = NULL, tries = 80L) {
  if (is.null(sd) || is.null(entry$vector)) return(FALSE)
  q <- .symVerifyPrime
  pn <- poolNext
  # perturb only the relevant leaves, as the sampling did; a fully random point almost
  # never has a modular steady state
  pertIdx <- if (is.null(relLeaves)) seq_along(point0) else relLeaves
  for (t in seq_len(tries)) {
    pt <- as.numeric(point0)
    if (length(pertIdx)) pt[pertIdx] <- pool(pn + seq_along(pertIdx) - 1L)
    pn <- pn + length(pertIdx) + 1L
    rq <- kcall(pt, q, as.integer(NtUsed))
    # an all-constant direction (no relevant leaf) is checked once at the base point
    if (!length(pertIdx) && !isTRUE(rq$ok)) return(FALSE)
    if (!isTRUE(rq$ok) || !identical(as.integer(rq$pivots), as.integer(pivots))) next
    env <- setNames(as.list(pt[seq_along(leafNames)]), leafNames)
    v <- integer(nz); v[f + 1L] <- 1L; bad <- FALSE
    for (nm in names(entry$vector)) {
      col <- match(nm, znames); if (is.na(col)) next
      got <- .symEvalModq(entry$vector[[nm]], env, q, sd)
      if (is.na(got)) { bad <- TRUE; break }
      v[col] <- as.integer(got %% q)
    }
    if (bad) next
    freeCols <- setdiff(seq_len(nz) - 1L, as.integer(rq$pivots))
    recon <- numeric(nz)
    for (fcol in freeCols) {
      vf <- v[fcol + 1L] %% q
      if (vf != 0) recon <- (recon + .symMulmod(.symNullResidues(rq, fcol, q), vf, q)) %% q
    }
    return(all((recon - v) %% q == 0))
  }
  FALSE
}


# ---- event gaps: rank along a line in the gap lengths --------------------------------

# Series rows S (rows x nzIn*N) with column c scaled by scale[c] and moved to column
# colMap[c] of a width nzOut*N block.
.symSeriesEmbed <- function(S, N, colMap, nzOut, p, scale = NULL) {
  out <- matrix(0, nrow(S), nzOut * N)
  for (c in seq_along(colMap)) {
    v <- S[, (c - 1L) * N + seq_len(N), drop = FALSE]
    if (!is.null(scale)) v[] <- .symMulmod(v, scale[c], p)
    out[, (colMap[c] - 1L) * N + seq_len(N)] <- v
  }
  out
}

# Constant rows as series rows (every entry at eps^0).
.symSeriesConst <- function(M, N) {
  out <- matrix(0, nrow(M), ncol(M) * N)
  if (nrow(M)) out[, (seq_len(ncol(M)) - 1L) * N + 1L] <- M
  out
}

# Attach the series rank of the stacked series blocks to a reduced kernel result. A
# kernel that moves with the gap lengths and is polynomial in them replaces the rows
# by their value at the actual gaps.
.symSeriesStack <- function(res, sblocks, nz, N, p) {
  S <- do.call(rbind, sblocks)
  storage.mode(S) <- "integer"
  sr <- symSeriesRank(S, as.integer(nz), as.integer(N), p, integer(0),
                      atOneBelow = as.integer(res$rank))
  res$rankS <- sr$rank; res$S <- sr$S; res$N <- N; res$atOne <- !is.null(sr$R)
  if (res$atOne) {
    res$R <- sr$R; res$pivots <- as.integer(sr$pivots); res$rank <- length(sr$pivots)
  }
  res
}

# The rank that decides identifiability: over the gap series where there is one.
.symRankOf <- function(r) if (is.null(r$rankS)) r$rank else r$rankS

# Progress of a kernel result in the saturation: both the rank of the stacked
# coefficient rows (whose nullspace the reconstruction uses) and the series rank.
.symRankScore <- function(r) r$rank + .symRankOf(r)


# Find a generic base point of maximal rank and raise the Lie order (and the gap series
# order Mtot) until the rank saturates. kcall(point, p, Nt, Mtot) returns the kernel list
# (ok, R, pivots, rank, dim); `blockCall` the same per condition, where the Lie order is
# decided; `budget` the codimension of the specialisation (NA: no certificate); `maxM`
# caps the gap order. Returns NULL without a usable point, else the reference
# reduction, the rank and the orders used.
.symSaturateCertify <- function(kcall, nLeaves, nz, maxM = 0L,
                                warm = function(pts, primes) invisible(),
                                probeBlock = 1L, blockCall = NULL,
                                budget = NA_integer_) {
  P <- .symPrimes[1]
  pool <- .symPool()
  point0 <- pool(seq_len(nLeaves))
  poolNext <- nLeaves + 1L

  # A flat step means saturation for one unspecialised condition, so the order is decided
  # per block and the stack built once at the largest. A specialisation can hide at most
  # `budget` growths, so flat steps are counted cumulatively; past the budget the rank is
  # final (vignette("Symmetries")). With `budget` NA the plateau is a heuristic, backed
  # by verify = TRUE. DMOD_SYM_LIEPLATEAU caps the flat steps, DMOD_SYM_LIEPLATEAU_BLOCK
  # forces the per-block value (0: stacked rule), DMOD_SYM_LIEDIAG traces.
  plateauNeed <- max(1L, as.integer(Sys.getenv("DMOD_SYM_LIEPLATEAU", "3")))
  blockOverride <- suppressWarnings(
    as.integer(Sys.getenv("DMOD_SYM_LIEPLATEAU_BLOCK", NA_character_)))
  lieDiag <- nzchar(Sys.getenv("DMOD_SYM_LIEDIAG"))
  perBlock <- length(blockCall) > 1L && (is.na(blockOverride) || blockOverride > 0L)
  budget <- if (is.null(budget) || length(budget) != 1L) NA_integer_ else as.integer(budget)
  # the budget covers one filtration: a block, or the stack of a single condition
  useBudget <- perBlock || length(blockCall) == 1L
  needUsed <- if (!is.na(blockOverride) && blockOverride > 0L) blockOverride
              else if (!useBudget || is.na(budget)) plateauNeed
              else min(plateauNeed, budget + 1L)
  certified <- useBudget && !is.na(budget) && needUsed >= budget + 1L

  # Raise the Lie order until the rank is flat for `need` orders. Returns the last kernel
  # and `grew`, the last order that raised the rank. From `from` > 1 the rank is
  # monotone, so equal ends of [from, from + need] settle the block without the middle.
  scanNt <- function(call1, point, Mtot, from, need, label) {
    if (from > 1L) {
      lo <- call1(point, P, from, Mtot)
      if (!isTRUE(lo$ok)) return(NULL)
      hi <- call1(point, P, from + need, Mtot)
      if (!isTRUE(hi$ok)) return(NULL)
      if (.symRankScore(lo) == .symRankScore(hi)) {
        if (lieDiag) message("[liediag] ", label, " Mtot=", Mtot, " rank ", lo$rank,
                             " flat over Lie order ", from, "-", from + need,
                             " (nz=", nz, ")")
        return(list(res = hi, Nt = from + need, grew = from))
      }
    }
    Nt <- max(1L, from); prev <- -1L; flat <- 0L; grew <- Nt
    res <- NULL; ranks <- integer(0)
    repeat {
      r <- call1(point, P, Nt, Mtot)
      if (!isTRUE(r$ok)) return(NULL)
      res <- r; ranks <- c(ranks, as.integer(.symRankOf(r)))
      # flat steps accumulate and are not reset by a growth
      if (.symRankScore(r) == prev) flat <- flat + 1L else grew <- Nt
      if (.symRankOf(r) >= nz || (flat >= need && Nt >= 2L) || Nt > nz + 1L) break
      prev <- .symRankScore(r); Nt <- Nt + 1L
    }
    if (lieDiag) message("[liediag] ", label, " Mtot=", Mtot, " ranks from Lie order ",
                         max(1L, from), ": ", paste(ranks, collapse = ","),
                         " (nz=", nz, ")")
    list(res = res, Nt = Nt, grew = grew)
  }

  saturateNt <- function(point, Mtot) {
    if (!perBlock) return(scanNt(kcall, point, Mtot, 1L, needUsed, "stacked"))
    orders <- integer(length(blockCall)); Nt <- 1L; driver <- 1L
    for (bi in seq_along(blockCall)) {
      sb <- scanNt(blockCall[[bi]], point, Mtot, Nt, needUsed, paste0("block ", bi))
      if (is.null(sb)) return(NULL)
      orders[bi] <- sb$grew
      if (sb$grew > Nt) { Nt <- sb$grew; driver <- bi }
    }
    r <- kcall(point, P, Nt, Mtot)
    if (!isTRUE(r$ok)) return(NULL)
    if (lieDiag) message("[liediag] block Lie orders: ", paste(orders, collapse = ","),
                         " -> stacked order ", Nt, ", rank ", r$rank)
    list(res = r, Nt = Nt, blockOrders = orders, blockDriver = driver)
  }

  # several generic points may be tried before one admits a steady-state point
  # over GF(p) for every condition (saturate the Lie order at gap order 0 first)
  sat <- NULL
  for (attempt in 1:50) {
    # after the first failure, warm the solves of the next candidates in one parallel
    # batch; they are the pool draws the loop makes anyway, so the result is unchanged
    if (probeBlock > 1L && attempt >= 2L && ((attempt - 2L) %% probeBlock) == 0L)
      warm(lapply(seq_len(probeBlock) - 1L, function(k)
             if (k == 0L) point0
             else pool(poolNext + (k - 1L) * nLeaves + seq_len(nLeaves) - 1L)),
           rep(list(P), probeBlock))
    sat <- saturateNt(point0, 0L)
    if (!is.null(sat)) break
    point0 <- pool(poolNext + seq_len(nLeaves) - 1L); poolNext <- poolNext + nLeaves
  }
  if (is.null(sat)) return(NULL)
  NtUsed <- sat$Nt; MtotUsed <- 0L; saturatedM <- TRUE

  # raise the gap order until the rank stops growing (exact generic-timing rank);
  # if the cap is hit while still growing, the truncated rank is conservative
  if (maxM > 0L) repeat {
    if (.symRankOf(sat$res) >= nz) break
    if (MtotUsed >= maxM) { saturatedM <- FALSE; break }
    satM <- saturateNt(point0, MtotUsed + 1L)
    if (is.null(satM)) break
    MtotUsed <- MtotUsed + 1L
    if (.symRankScore(satM$res) <= .symRankScore(sat$res)) break
    sat <- satM; NtUsed <- sat$Nt
  }

  # a kernel that moves with the gap lengths is read at the actual gaps once it is known
  # as a polynomial in them, which takes further gap orders
  while (.symRankOf(sat$res) < sat$res$rank && MtotUsed < maxM) {
    r <- kcall(point0, P, NtUsed, MtotUsed + 1L)
    if (!isTRUE(r$ok)) break
    MtotUsed <- MtotUsed + 1L; sat$res <- r
  }

  # cross-prime rank check at the same point; the probes are warmed as one batch
  warm(rep(list(point0), length(.symPrimes) - 1L), as.list(.symPrimes[-1]))
  rankMax <- .symRankScore(sat$res)
  for (pj in .symPrimes[-1]) {
    rj <- kcall(point0, pj, NtUsed, MtotUsed)
    if (isTRUE(rj$ok)) rankMax <- max(rankMax, .symRankScore(rj))
  }
  while (.symRankScore(sat$res) < rankMax) {
    point0 <- pool(poolNext + seq_len(nLeaves) - 1L); poolNext <- poolNext + nLeaves
    sat <- saturateNt(point0, MtotUsed)
    if (is.null(sat)) return(NULL)
    NtUsed <- sat$Nt
  }
  list(ref = sat$res, NtUsed = NtUsed, MtotUsed = MtotUsed, saturatedM = saturatedM,
       point0 = point0, pool = pool, poolNext = poolNext,
       blockOrders = sat$blockOrders, blockDriver = sat$blockDriver,
       budget = budget, plateau = needUsed, certified = certified,
       rank = sat$res$rank, rankS = .symRankOf(sat$res), pivots = sat$res$pivots)
}


# Nullspace basis of the reference reduction: one column per free coordinate.
.symNullspaceBasis <- function(ref, freeCols, P) {
  if (!length(freeCols)) return(matrix(0L, ref$dim, 0L))
  vapply(freeCols, function(fc) .symNullResidues(ref, fc, P), integer(ref$dim))
}


# Classify a closed-form direction as "scaling" (removed by fixing one coordinate) or
# "general" (removed by reparametrisation) on its canonical poly-primitive generator, so
# both engines agree and a disguised scaling is recognised. A symbolic (Hill) weight is
# kept. Updates $type, $vector (integer weights for a scaling) and $degree (-1 if not
# polynomial), which selects the certifyPoly candidates.
.symClassifyDirection <- function(d, sd) {
  if (is.null(d$vector)) return(d)
  if (isTRUE(d$type == "scaling")) {
    wint <- suppressWarnings(as.integer(as.character(unlist(d$vector))))
    if (anyNA(wint)) return(d)                       # Hill / symbolic weight: keep
    lit <- setNames(as.list(paste0(as.character(unlist(d$vector)), "*",
                                   names(d$vector))), names(d$vector))
  } else {
    lit <- setNames(as.list(as.character(unlist(d$vector))), names(d$vector))
  }
  cls <- tryCatch(sd$classifyDirection(lit), error = function(e) NULL)
  if (is.null(cls) || is.null(cls$type)) return(d)
  d$type <- as.character(cls$type)
  d$degree <- if (!is.null(cls$degree)) as.integer(cls$degree) else NULL
  if (identical(d$type, "scaling")) {
    d$vector <- lapply(cls$weights, function(w) as.character(w))
    d$support <- .symSort(names(cls$weights))
  } else {
    d$vector <- lapply(cls$components, function(x) as.character(x))
    d$support <- .symSort(names(cls$components))
  }
  d
}


# classify every closed-form direction; support-only ones are left alone
.symRelabelDirections <- function(nonId, sd) {
  if (is.null(sd) || !length(nonId)) return(nonId)
  lapply(nonId, function(d) if (is.null(d$vector)) d
                            else .symClassifyDirection(d, sd))
}


# reconstControl(certifyPoly = TRUE): sets $certified on each general polynomial
# direction of degree 1 to degreeCap in the span of the polynomial engine's generators
# (a Lie point symmetry). FALSE does not make it less non-identifiable.
.symCertifyPoly <- function(nonId, modelLines, obsLines, forcings, fixed,
                              parameters, ctrl, sd) {
  cap <- if (is.null(ctrl$degreeCap)) 4L else as.integer(ctrl$degreeCap)
  isCand <- function(d) !isTRUE(d$type == "scaling") && !is.null(d$degree) &&
    !is.na(d$degree) && d$degree >= 1L && d$degree <= cap
  cand <- Filter(isCand, nonId)
  if (!length(cand) || is.null(sd)) return(nonId)
  deg <- if (!is.null(ctrl$certifyPolyDeg)) as.integer(ctrl$certifyPolyDeg)
         else max(vapply(cand, function(d)
           if (is.null(d$degree)) 1L else as.integer(d$degree), integer(1)))
  gens <- tryCatch({
    reticulate::py_capture_output(
      r <- sd$symmetryDetectiondMod(
        model = modelLines, observation = obsLines,
        ansatz = "multi", pMax = as.integer(max(1L, deg)),
        inputs = if (length(forcings)) forcings else NULL, fixed = fixed,
        allTrafos = TRUE, lieOrder = 0L, exact = TRUE, verify = FALSE,
        backend = "symengine",
        parameters = if (length(parameters)) parameters else NULL, method = "liesym"))
    r
  }, error = function(e) NULL)
  genVecs <- if (is.null(gens)) list()
             else Filter(Negate(is.null), lapply(gens, function(g) g$infinitesimals))
  lapply(nonId, function(d) {
    if (!isCand(d)) return(d)
    d$certified <- length(genVecs) > 0L &&
      isTRUE(tryCatch(sd$certifyInSpan(d$vector, genVecs)$certified,
                      error = function(e) FALSE))
    d
  })
}


# Report on the physical coordinates when the kernel ran over auxiliary ones too
# (`auxNames`: per-condition state columns, recast atoms E and L). Drops directions
# without physical support, recomputes dim/rank/identifiable from the projected
# nullspace and moves auxiliary entries to `auxField`. `flagJoint` marks physical
# entries that still reference an auxiliary coordinate.
.symReportPhysical <- function(result, znames, auxNames, sc, nz, P, auxField,
                                 flagJoint = FALSE) {
  physCoords <- setdiff(znames, auxNames)
  keep <- vapply(result$nonIdentifiable, function(d)
    any(d$support %in% physCoords), logical(1))
  result$nonIdentifiable <- result$nonIdentifiable[keep]
  physCols <- which(znames %in% physCoords)
  N <- .symNullspaceBasis(sc$ref, setdiff(0:(nz - 1L), sc$pivots), P)
  physNull <- if (length(physCols) && ncol(N))
    .symRrefModp(N[physCols, , drop = FALSE], P)$rank else 0L
  # with directions that depend on the event gaps, project the series kernel instead:
  # its part with vanishing physical entries is the kernel of the auxiliary columns
  if (!is.null(sc$rankS) && sc$rankS < sc$rank) {
    auxCols <- which(!(znames %in% physCoords)) - 1L
    rAux <- if (length(auxCols))
      symSeriesRank(sc$ref$S, as.integer(nz), as.integer(sc$ref$N), P, auxCols)$rank
      else 0L
    physNull <- (nz - sc$rankS) - (length(auxCols) - rAux)
  }
  result$dim <- length(physCoords)
  result$coordinates <- physCoords
  result$rank <- as.integer(length(physCoords) - physNull)
  result$identifiable <- (physNull == 0L)
  result$nonIdentifiable <- lapply(result$nonIdentifiable, function(d) {
    d$support <- .symSort(intersect(d$support, physCoords))
    if (!is.null(d$vector)) {
      isAux <- names(d$vector) %in% auxNames
      d[[auxField]] <- d$vector[isAux]
      phys <- d$vector[!isAux]
      if (flagJoint && any(vapply(phys, function(e)
            length(intersect(getSymbols(as.character(e)), auxNames)) > 0L, logical(1))))
        d$jointForm <- TRUE
      d$vector <- phys
    }
    d
  })
  result
}


# Peel the scalings common to all conditions: project each onto the coordinates z, keep
# its tangent w_c * z_c at the base point if in the nullspace and independent of those
# taken. Returns the scalings and their tangent span.
.symPeelScalings <- function(scalRes, znames, nz, zval, P, N, sd = NULL) {
  scaling <- list()
  Bmat <- matrix(0L, nz, 0L)
  inSpan <- function(M, x) .symInSpan(M, x, nz, P)
  # a symbolic weight (xi_kinh = -nhill*kinh) is evaluated at the base point for the
  # tangent and reported verbatim
  env <- as.list(setNames(as.numeric(zval), znames))
  evalW <- function(w) {
    wi <- suppressWarnings(as.numeric(w))
    if (!is.na(wi) && wi == round(wi)) return(as.numeric(wi %% P))
    if (is.null(sd)) return(NA_real_)
    v <- .symEvalModq(w, env, P, sd)
    if (is.null(v) || is.na(v)) NA_real_ else as.numeric(v)
  }
  for (d in scalRes$nonIdentifiable) {
    cols <- match(names(d$vector), znames)
    keep <- !is.na(cols)
    if (!any(keep)) next
    wsym <- as.character(unlist(d$vector))[keep]
    wv <- vapply(wsym, evalW, numeric(1))            # weights evaluated mod P
    if (anyNA(wv)) next
    v <- integer(nz); v[cols[keep]] <- as.integer(wv %% P)
    if (all(v == 0)) next
    tangent <- as.integer((as.numeric(v) * zval) %% P)
    if (all(tangent == 0) || !inSpan(N, tangent) || inSpan(Bmat, tangent)) next
    Bmat <- cbind(Bmat, tangent)
    supp <- names(d$vector)[keep]
    ord <- order(supp, method = "radix")
    scaling[[length(scaling) + 1L]] <- list(
      support = supp[ord],
      vector = setNames(as.list(wsym[ord]), supp[ord]),
      type = "scaling", closedForm = TRUE)
  }
  list(scaling = scaling, Bmat = Bmat)
}


# Expand a scaling onto the per-condition joint coordinates: a state weight goes to each
# of its K columns, a held pivot's weight to its initial-value parameter (heldParamOf).
.symJointExpandScal <- function(scalRes, stateBase, Kc, heldParamOf = character(0)) {
  scalRes$nonIdentifiable <- lapply(scalRes$nonIdentifiable, function(d) {
    vec <- list()
    for (nm in names(d$vector)) {
      w <- d$vector[[nm]]
      if (nm %in% names(heldParamOf))
        vec[[heldParamOf[[nm]]]] <- w
      else if (nm %in% stateBase)
        for (m in seq_len(Kc)) vec[[paste0(nm, "|c", m)]] <- w
      else vec[[nm]] <- w
    }
    d$vector <- vec
    d$support <- names(vec)
    d
  })
  scalRes
}


# ---- parallel steady-state solve pool (PSOCK workers) --------------------------------

# Steady-state solve on a PSOCK worker, through the engine in option
# `dMod.sym.worker_sd`; arguments are plain data. Returns the solve list or NULL.
.symSolveWorker <- function(job, cargs) {
  sd <- getOption("dMod.sym.worker_sd")
  if (is.null(sd)) return(NULL)
  tryCatch(sd$solveSteadyStateModular(
    model = job$model, stateNames = cargs$stateNames,
    paramNames = cargs$paramNames, paramVals = job$paramVals, prime = job$p,
    forcings = cargs$forcings, t0events = job$t0events,
    recast = cargs$recast, lVals = job$lVals, jointMode = TRUE,
    heldStates = job$heldVals),
    error = function(e) NULL)
}


# PSOCK pool of Python engines for the steady-state solves where mclapply cannot fork.
# Workers run Python only. Returns the cluster, or NULL (serial) when n <= 1 or a
# worker cannot start its interpreter.
.symMakeSolveCluster <- function(n) {
  n <- as.integer(n)
  if (is.na(n) || n <= 1L) return(NULL)
  if (!requireNamespace("parallel", quietly = TRUE)) return(NULL)
  cl <- tryCatch(parallel::makeCluster(n), error = function(e) NULL)
  if (is.null(cl)) return(NULL)
  # initialiser per worker, reparented to globalenv so it serialises without this
  # frame; returns TRUE or an error string
  codeDir <- system.file("code", package = "dMod2")
  pyPath  <- Sys.getenv("RETICULATE_PYTHON")
  initFun <- function(codeDir, pyPath) {
    tryCatch({
      if (!requireNamespace("reticulate", quietly = TRUE)) stop("no reticulate")
      if (nzchar(pyPath)) {
        Sys.setenv(RETICULATE_PYTHON = pyPath)
        suppressWarnings(suppressMessages(
          tryCatch(reticulate::use_python(pyPath, required = FALSE),
                   error = function(e) NULL)))
      }
      sysmod <- reticulate::import("sys", convert = TRUE)
      if (!(codeDir %in% sysmod$path)) sysmod$path <- c(codeDir, sysmod$path)
      # an option persists across tasks on the same node
      options(dMod.sym.worker_sd =
                reticulate::import("symmetryDetection", convert = TRUE))
      TRUE
    }, error = function(e) conditionMessage(e))
  }
  environment(initFun) <- globalenv()
  res <- tryCatch(parallel::clusterCall(cl, initFun, codeDir, pyPath),
                  error = function(e) list(conditionMessage(e)))
  ok <- all(vapply(res, isTRUE, logical(1)))
  if (!isTRUE(ok)) {
    if (nzchar(Sys.getenv("DMOD_SYM_TIMING")))
      message("[sym] solve pool init failed: ",
              paste(unique(vapply(res, function(r)
                if (isTRUE(r)) "" else as.character(r)[1], character(1))), collapse = "; "))
    tryCatch(parallel::stopCluster(cl), error = function(e) NULL)
    return(NULL)
  }
  cl
}


# ---- observability engine: build [Obs; df], peel scalings, reconstruct residuals -----

# Multi-condition analytic observability on the per-condition tapes of
# compileObservabilityTapeMulti over shared coordinates. The rows of all conditions
# are stacked, so the verdict and the directions refer to the common nullspace.
.observability_analytic_multi <- function(multi, spy = NULL,
                                          closedForm = FALSE, sd = NULL, cores = 1,
                                          equilZeroStates = character(0),
                                          t0events = list(), nConditions = NULL,
                                          chainOf = NULL, nGaps = 0L,
                                          implicitSteadyState = FALSE,
                                          control = reconstControl(), verify = FALSE,
                                          codimSpec = NA_integer_) {
  ctrl <- control
  jointSS <- isTRUE(multi$jointSteadyState) && isTRUE(implicitSteadyState)
  # ==== parallelism: fork axis vs. kernel threads ===================================
  # `cores` feeds two nested axes: the per-point fork (coresGLp) takes the whole
  # budget, the kernel threads (coresCall) run serial inside it and in full elsewhere
  cores <- as.integer(max(1L, cores))
  coresGLp <- cores
  # per-point kcalls with their steady-state solves forked on unix, serial elsewhere;
  # a crashed child counts as a failed solve
  parMap <- if (coresGLp > 1L && .Platform$OS.type == "unix")
    function(xs, f) lapply(
      parallel::mclapply(xs, f, mc.cores = coresGLp, mc.preschedule = TRUE),
      function(o) if (inherits(o, "try-error")) NULL else o)
    else function(xs, f) lapply(xs, f)
  # warmSolves fills the steady-state solve cache in parallel on every platform (fork
  # on unix, a PSOCK pool elsewhere); the kernel then runs on the cached solves
  usePool <- coresGLp > 1L
  cl <- NULL                      # PSOCK pool, stood up lazily by solveMap below
  poolTried <- FALSE
  warmOff <- FALSE                # set inside the parMap fork: no nested forking
  # one registration covers the lazily created pool: `cl` is read at exit time
  on.exit(if (!is.null(cl)) tryCatch(parallel::stopCluster(cl),
                                     error = function(e) NULL), add = TRUE)
  # replaced in the jointSS block
  warmSolves <- function(pts, primes, conds = NULL) invisible()
  nLeaves <- as.integer(multi$nLeaves)
  nStates <- as.integer(multi$nStates)
  zSlots <- as.integer(multi$zSlots)
  znames <- as.character(multi$znames)
  nz <- length(znames)
  leafNames <- as.character(multi$leafNames)
  ssConstraint <- isTRUE(multi$equilibrate)
  nSegmentTapes <- length(multi$tapes)
  # later events: the kernel propagates across each gap as a power series in the gap
  # length, segments form one chain per condition, maxM caps the series order
  hasGaps <- isTRUE(nGaps > 0L) && !is.null(chainOf)
  chainGroups <- if (hasGaps) split(seq_along(multi$tapes), chainOf) else NULL
  firstOfChain <- if (hasGaps) !duplicated(chainOf) else rep(TRUE, length(multi$tapes))
  maxM <- if (hasGaps) ctrl$gapOrderCap else 0L
  # the kernel threads over conditions/segments, so cap at their count; 1 inside the
  # per-point fork
  nKernelUnits <- if (hasGaps) length(chainGroups) else length(multi$tapes)
  coresCall <- as.integer(max(1L, min(cores, nKernelUnits)))
  # recast coordinates E = base^exp, L = log base extend the leaf space for the fit and
  # are back-substituted, giving closed forms with a free exponent
  recast <- list()
  nAug <- nLeaves
  leafNamesAug <- leafNames
  # auxiliary leaves (per-condition states, recast E/L), exempt from the relevance cap
  auxLeaves <- integer(0)
  # transient recast (free exponent without equilibrate): E and L are free-initial-value
  # leaves tied to (base, exp) by stacked relation rows
  recastTransient <- FALSE
  recastAtomNames <- character(0)

  # ==== marshal the compiled tapes and the recast coordinate metadata ===============
  tapeFields <- function(t) {
    out <- list(
      op = as.integer(t$op), a = as.integer(t$a), b = as.integer(t$b),
      cnum = as.character(t$cnum), cden = as.character(t$cden),
      stateSlots = as.integer(t$stateSlots), fOut = as.integer(t$fOut),
      gOut = as.integer(t$gOut), icLeaf = as.integer(t$icLeaf),
      icNum = as.character(t$icNum), icDen = as.character(t$icDen))
    # IC tape seeding the initial values and their duals (free or carried values, doses,
    # resets; for a joint equilibrate anchor the identity seed filled with the resting state)
    if (!is.null(t$icOp))
      out <- c(out, list(
        icOp = as.integer(t$icOp), icA = as.integer(t$icA),
        icB = as.integer(t$icB), icCnum = as.character(t$icCnum),
        icCden = as.character(t$icCden), icOut = as.integer(t$icOut)))
    # a later segment's left-boundary state-dose event map (replace/add/multiply by
    # a parametric value), applied by the kernel to the propagated state
    if (!is.null(t$evVarIdx))
      out <- c(out, list(
        evOp = as.integer(t$evOp), evA = as.integer(t$evA),
        evB = as.integer(t$evB), evCnum = as.character(t$evCnum),
        evCden = as.character(t$evCden), evOut = as.integer(t$evOut),
        evVarIdx = as.integer(t$evVarIdx), evMethod = as.integer(t$evMethod)))
    # the segment's left boundary time (value and duals)
    if (!is.null(t$tmOp))
      out <- c(out, list(
        tmOp = as.integer(t$tmOp), tmA = as.integer(t$tmA), tmB = as.integer(t$tmB),
        tmCnum = as.character(t$tmCnum), tmCden = as.character(t$tmCden),
        tmOut = as.integer(t$tmOut)))
    out
  }
  tapes <- lapply(multi$tapes, tapeFields)

  # batched per-chunk evaluator for the coupled+gap reconstruction loop (set inside
  # the joint block below when hasGaps); NULL keeps the serial parMap fallback.
  kchunk <- NULL

  # kernel per tape or chain, where .symSaturateCertify() decides the Lie order; stacked
  # constant rows (df tangency, recast relations) stay out
  obsBlockCalls <- function() {
    zs <- zSlots; nl <- nLeaves; ns <- nStates
    if (hasGaps)
      lapply(chainGroups, function(idx) { force(idx)
        function(point, p, Nt, Mtot = 0L)
          symObsNullChain(list(tapes[idx]), nl, ns, zs, as.integer(point), p,
                          as.integer(Nt), as.integer(Mtot), 1L) })
    else
      lapply(seq_along(tapes), function(i) { force(i)
        function(point, p, Nt, Mtot = 0L)
          symObsNullMulti(list(tapes[[i]]), nl, ns, zs, as.integer(point), p,
                          as.integer(Nt), 1L) })
  }
  blockCall <- NULL

  if (ssConstraint) {
    # equilibrate: states stay coordinates seeded at x* per (point, prime); f = 0 enters
    # as stacked df tangency rows
    stateNames <- as.character(multi$stateNames)
    paramNames <- as.character(multi$paramNames)
    forcings <- if (is.null(multi$forcings)) character(0) else as.character(multi$forcings)
    # forcings and forced-zero states are held at zero in the f = 0 solve
    solveHeld <- union(forcings, as.character(equilZeroStates))
    models <- lapply(multi$tapes, function(t) as.character(t$constraintModel))
    # only the equilibrate-seeded segments are solved for a steady state; with
    # gaps only the first segment of each chain anchors (the rest are propagated)
    isEquilTape <- vapply(multi$tapes,
                          function(t) length(t$constraintModel) > 0L, logical(1)) &
                   firstOfChain
    w <- nz + 1L
    # recast atoms: per base a generic coordinate plus L = log base; E is generic and the
    # base solved, or inverted. realStateNames are the names solved for.
    recast <- multi$powerRecast
    if (is.null(recast)) recast <- list()
    realStateNames <- if (length(recast)) as.character(multi$realStateNames) else stateNames
    genName <- function(r) if (isTRUE(r$inverted)) as.character(r$base) else as.character(r$E)
    genNames <- vapply(recast, genName, "")
    lNames <- unique(vapply(recast, function(r) as.character(r$L), ""))
    # held pivots reach the solve as frozen states, never also as parameters
    heldSolveNames <- if (length(multi$heldStateParams))
      as.character(unlist(multi$heldStateParams)) else character(0)
    solveParamNames <- setdiff(c(paramNames, genNames), heldSolveNames)
    # in joint mode E and L are already free-state leaves
    if (length(recast) && !jointSS) {
      leafNamesAug <- c(leafNames, genNames, lNames)
      nAug <- length(leafNamesAug)
    }
    ssWhy <- NULL              # last steady-state failure reason, for diagnostics
    if (jointSS) {
      # joint determining system: observability over (x, theta) stacked with the
      # tangency rows df_rest = [Jx | Jt]; the nullspace of [Obs; df] is the symmetry space
      slotOfName <- function(nm) { i <- match(nm, leafNames)
        if (is.na(i)) NULL else i }
      # held pivots (reduceCQ = FALSE): each moiety pivot is frozen at the residue of its
      # initial-value parameter (heldParamOf), the other states are solved
      heldParamOf <- if (length(multi$heldStateParams))
        unlist(multi$heldStateParams) else character(0)   # named: pivot -> p_e
      heldNames <- setdiff(intersect(names(heldParamOf), realStateNames), solveHeld)
      jointPV <- function(point, p) {
        val <- function(nm) { s <- slotOfName(nm)
          if (is.null(s)) 0 else as.numeric(point[s]) %% p }
        pv <- as.list(setNames(vapply(solveParamNames, val, numeric(1)), solveParamNames))
        lv <- if (length(lNames))
          as.list(setNames(vapply(lNames, val, numeric(1)), lNames)) else NULL
        # freeze value of a pivot = residue of ITS initial-value parameter leaf
        hv <- if (length(heldNames))
          as.list(setNames(vapply(heldNames,
            function(nm) val(heldParamOf[[nm]]), numeric(1)), heldNames)) else NULL
        list(paramVals = pv, lVals = lv, heldVals = hv)
      }
      # ---- per-condition state coordinates of the joint system ---------------------
      # Each condition has its own resting state, so every state coordinate (states, E,
      # L) is per condition; a shared column would miss non-scaling directions. States
      # and E are log-normalised so a scaling weight is equal across columns; L is additive.
      perCondCols <- which(znames %in% as.character(multi$zStateNames))   # states + E + L
      logCols <- which(znames %in% setdiff(as.character(multi$zStateNames), lNames))  # states + E
      jointStateSlot <- zSlots[logCols] + 1L        # 1-based leaf slot to read x*_c
      equilConds <- which(isEquilTape)              # one anchor tape per condition
      Kc <- length(equilConds)
      sharedIdx <- setdiff(seq_len(nz), perCondCols)   # real params only
      nShared <- length(sharedIdx); nSt <- length(perCondCols)
      nzWide <- nShared + Kc * nSt
      # local znames-column (1..nz) -> wide column, for the mi-th equilibrate condition
      wideCols <- lapply(seq_len(Kc), function(mi) {
        wc <- integer(nz)
        wc[sharedIdx] <- seq_len(nShared)
        wc[perCondCols] <- nShared + (mi - 1L) * nSt + seq_len(nSt)
        wc
      })
      # wide labels and leaf slots: state coordinates duplicated per condition as "x|c<m>"
      stateBase <- znames[perCondCols]
      znamesWide <- c(znames[sharedIdx],
                      unlist(lapply(seq_len(Kc), function(mi)
                        paste0(stateBase, "|c", mi))))
      zSlotsWide <- c(zSlots[sharedIdx], rep(zSlots[perCondCols], times = Kc))
      zStateNamesWide <- if (nSt) znamesWide[(nShared + 1L):nzWide] else character(0)
      # per-condition copies; the outer names are widened below and closures would see that
      nzL <- nz; znamesL <- znames; zSlotsL <- zSlots
      # ---- recast relation rows tying E = base^exp and L = log base ----------------
      # without the algebraic ties E and L would be spuriously free. Per condition:
      #   E:  xi_E/E = exp * xi_base/base + log(base) * xi_exp   (log-weight columns)
      #   L:  xi_L   = xi_base/base
      # with exp and log(base) evaluated at the sample point.
      recastRel <- lapply(recast, function(rc) list(
        baseCol = match(as.character(rc$base), znamesL),
        ECol    = match(as.character(rc$E),    znamesL),
        LCol    = match(as.character(rc$L),    znamesL),
        expCol  = match(as.character(rc$exp),  znamesL),
        expSlot = slotOfName(as.character(rc$exp)),
        LSlot   = slotOfName(as.character(rc$L)),
        # a parameter base (e.g. a Km under the exponent) is not log-normalised, so its
        # coefficient carries 1/base
        baseSlot  = slotOfName(as.character(rc$base)),
        baseParam = as.character(rc$base) %in% paramNames))
      # df columns allowed to map to no coordinate (held at 0 or fixed); any other
      # unmapped column is a lost constraint and traced under DMOD_JOINT_DIAG
      dfDroppable <- unique(c(as.character(forcings), as.character(equilZeroStates),
                              setdiff(as.character(leafNames), znames)))
      # df constraint rows for one condition over the nz znames columns
      dfRowsCond <- function(sol, p) {
        Jx <- sol$dfJx; Jt <- sol$dfJt
        scn <- as.character(sol$dfStateCols); pcn <- as.character(sol$dfParamCols)
        rows <- vector("list", length(Jx))
        for (i in seq_along(Jx)) {
          row <- numeric(nzL); xi <- as.numeric(Jx[[i]])
          for (j in seq_along(scn)) { col <- match(scn[j], znamesL)
            if (!is.na(col)) row[col] <- (row[col] + xi[j]) %% p
            else if (!scn[j] %in% dfDroppable && nzchar(Sys.getenv("DMOD_JOINT_DIAG")))
              message("[jointdiag] df drops non-fixed state column: ", scn[j]) }
          for (th in pcn) {
            # a held pivot's column maps to its initial-value parameter
            thz <- if (th %in% names(heldParamOf)) heldParamOf[[th]] else th
            col <- match(thz, znamesL)
            if (!is.na(col)) row[col] <- (row[col] + as.numeric(Jt[[th]])[i]) %% p
            else if (!th %in% dfDroppable && nzchar(Sys.getenv("DMOD_JOINT_DIAG")))
              message("[jointdiag] df drops non-fixed param column: ", th) }
          rows[[i]] <- row
        }
        if (!length(rows)) matrix(0, 0, nzL) else do.call(rbind, rows)
      }
      # ---- the coupled steady-state solve and its cache ----------------------------
      # A solve depends only on the solve parameters, log coordinates and held pivots,
      # so it is memoised on them (failures too); the parallel warm pool fills this cache.
      jointSolveCache <- new.env(parent = emptyenv())
      # keyed on the substituted model, so conditions differing only outside f (a dose, a
      # scale) share a solve; in jointMode the t0 events do not enter the result
      modelKey <- vapply(models, function(m) paste0(m, collapse = "\n"), character(1))
      if (nzchar(Sys.getenv("DMOD_SYM_TIMING")))
        message(sprintf("[sym] %d equilibrate condition(s), %d distinct steady state(s)",
                        Kc, length(unique(modelKey[equilConds]))))
      solveKey <- function(p, pv, ci)
        paste(modelKey[[ci]], p, paste0(unlist(pv$paramVals), collapse = ","),
              paste0(unlist(pv$lVals), collapse = ","),
              paste0(unlist(pv$heldVals), collapse = ","), sep = "|")
      # constant (point-independent) arguments of every solve, shipped once to the pool
      solveConst <- list(stateNames = realStateNames, paramNames = solveParamNames,
                         forcings = if (length(solveHeld)) solveHeld else NULL,
                         recast = if (length(recast)) recast else NULL)
      # serial solve cost, which decides whether a batch repays starting the PSOCK pool
      solveSecs <- 0; solveN <- 0L
      solveRaw <- function(p, pv, ci) {
        evC <- if (ci <= length(t0events) && length(t0events[[ci]]))
          t0events[[ci]] else NULL
        .tSolve <- Sys.time()
        out <- tryCatch(sd$solveSteadyStateModular(
          model = models[[ci]], stateNames = realStateNames,
          paramNames = solveParamNames, paramVals = pv$paramVals, prime = p,
          forcings = if (length(solveHeld)) solveHeld else NULL, t0events = evC,
          recast = if (length(recast)) recast else NULL, lVals = pv$lVals,
          jointMode = TRUE,
          heldStates = if (length(heldNames)) pv$heldVals else NULL),
          error = function(e) NULL)
        solveSecs <<- solveSecs + as.numeric(Sys.time() - .tSolve, units = "secs")
        solveN <<- solveN + 1L
        # per-solve trace, also from forked children
        if (nzchar(Sys.getenv("DMOD_SYM_SOLVEDIAG")))
          message(sprintf("[sym]   solve cond %d, pid %d: %.1fs", ci, Sys.getpid(),
                          as.numeric(Sys.time() - .tSolve, units = "secs")))
        out
      }
      # ---- filling the solve cache in parallel (both platforms) --------------------
      # Parallel map behind warmSolves: mclapply on unix, else a PSOCK pool started on
      # the first batch whose measured cost repays it and then reused. NULL: stay serial.
      # DMOD_SYM_SOLVEPOOL forces the pool on any platform.
      forcePool <- nzchar(Sys.getenv("DMOD_SYM_SOLVEPOOL"))
      mapKind <- "serial"          # which map the last fill used (DMOD_SYM_TIMING)
      solveMap <- function(jobs) {
        mapKind <<- "serial"
        if (!requireNamespace("parallel", quietly = TRUE)) return(NULL)
        n <- as.integer(min(coresGLp, length(jobs)))
        # a persistent pool keeps each worker's compiled models across waves, so it is
        # preferred once a measured batch repays its start-up; until then unix forks
        if (is.null(cl) && (usePool || forcePool) && !poolTried &&
            (forcePool || (solveN > 0L && length(jobs) * (solveSecs / solveN) >= 5))) {
          poolTried <<- TRUE
          cl <<- .symMakeSolveCluster(coresGLp)
        }
        if (!is.null(cl)) {
          # reparented so it serialises without the dMod namespace
          worker <- .symSolveWorker
          environment(worker) <- globalenv()
          mapKind <<- sprintf("pool x%d", coresGLp)
          # sorted by model with chunk.size = 1, so a model returns to the worker
          # holding its compile
          ord <- order(vapply(jobs, function(j) j$mkey, character(1)),
                       method = "radix")
          res <- tryCatch(parallel::parLapply(cl, jobs[ord], worker, cargs = solveConst,
                                              chunk.size = 1L), error = function(e) NULL)
          if (is.null(res) || length(res) != length(jobs)) return(NULL)
          out <- vector("list", length(jobs)); out[ord] <- res
          return(out)
        }
        # fork on unix, dynamically scheduled since solve times vary by orders of magnitude
        if (.Platform$OS.type == "unix" && !forcePool) {
          mapKind <<- sprintf("fork x%d", n)
          return(parallel::mclapply(jobs, function(job) solveRaw(job$p, job$pv, job$ci),
                                    mc.cores = n, mc.preschedule = FALSE))
        }
        NULL
      }
      # Fill the solve cache for a batch of points with the distinct uncached jobs. An
      # unanswered job stays uncached and is solved serially on demand.
      warmSolves <- function(pts, primes, conds = NULL) {
        if (warmOff || coresGLp <= 1L) return(invisible())
        if (is.null(conds)) conds <- equilConds
        jobs <- list(); seen <- new.env(parent = emptyenv())
        for (idx in seq_along(pts)) {
          p <- primes[[idx]]; pt <- pts[[idx]]
          pv <- jointPV(pt, p)                      # independent of the condition
          for (ci in conds) {
            key <- solveKey(p, pv, ci)
            if (!is.null(jointSolveCache[[key]]) || !is.null(seen[[key]])) next
            seen[[key]] <- TRUE
            evC <- if (ci <= length(t0events) && length(t0events[[ci]]))
              t0events[[ci]] else NULL
            jobs[[length(jobs) + 1L]] <- list(key = key, ci = ci, p = p, pv = pv,
              mkey = modelKey[[ci]],           # worker affinity: same model, same node
              model = models[[ci]], paramVals = pv$paramVals, lVals = pv$lVals,
              heldVals = if (length(heldNames)) pv$heldVals else NULL,
              t0events = evC)
          }
        }
        if (length(jobs) < 2L) {
          if (nzchar(Sys.getenv("DMOD_SYM_TIMING")) && length(jobs))
            message("[sym] warm fill: 1 solve, serial (nothing to spread)")
          return(invisible())
        }
        .tFill <- Sys.time()
        res <- solveMap(jobs)
        if (nzchar(Sys.getenv("DMOD_SYM_TIMING")))
          message(sprintf("[sym] warm fill: %d solve(s), %s, %.1fs%s", length(jobs),
                          mapKind, as.numeric(Sys.time() - .tFill, units = "secs"),
                          if (is.null(res)) " -> unavailable, staying serial" else ""))
        if (is.null(res) || length(res) != length(jobs)) return(invisible())
        for (i in seq_along(jobs)) {
          r <- res[[i]]
          if (inherits(r, "try-error")) r <- NULL
          # ok = FALSE is cached as a failure; NULL (worker error) is left for a serial retry
          if (!is.null(r)) jointSolveCache[[jobs[[i]]$key]] <-
            if (isTRUE(r$ok)) r else list(ok = FALSE)
        }
        invisible()
      }
      # solve one equilibrate condition and return its solution plus the point with
      # this condition's on-manifold x* written into the state leaves (cache-backed)
      jointSolveCond <- function(point, p, ci) {
        pv <- jointPV(point, p)
        key <- solveKey(p, pv, ci)
        sol <- jointSolveCache[[key]]
        if (is.null(sol)) {
          sol <- solveRaw(p, pv, ci)
          jointSolveCache[[key]] <- if (is.null(sol) || !isTRUE(sol$ok)) list(ok = FALSE) else sol
          if (is.null(sol) || !isTRUE(sol$ok)) return(NULL)
        } else if (!isTRUE(sol$ok)) return(NULL)
        ptc <- as.numeric(point[seq_len(nLeaves)])
        # seed the pre-event resting value: the IC tape applies the events, and df is
        # linearised at the same point
        vb <- sol$valBy
        for (k in seq_along(realStateNames)) {
          s <- slotOfName(realStateNames[k]); v <- vb[[realStateNames[k]]]
          if (!is.null(s) && !is.null(v)) ptc[s] <- as.numeric(v) %% p
        }
        list(sol = sol, ptc = ptc)
      }
      # ---- the forward solve variant (states chosen, rates solved) -----------------
      # states and free parameters from the point, f = 0 solved linearly for the turnover
      # rates `solveRates`; valid at every prime. Same {sol, ptc} as the backward solve.
      jointSolveCondFwd <- function(point, p, ci, solveRates) {
        rd <- function(nm) { s <- slotOfName(nm); if (is.null(s)) 0 else as.numeric(point[s]) %% p }
        sv <- as.list(setNames(vapply(realStateNames, rd, numeric(1)), realStateNames))
        pv <- as.list(setNames(vapply(solveParamNames, rd, numeric(1)), solveParamNames))
        sol <- tryCatch(sd$solveForwardModular(models[[ci]], realStateNames, solveParamNames,
                          sv, pv, p, forcings = if (length(solveHeld)) solveHeld else NULL,
                          solveRates = solveRates), error = function(e) NULL)
        if (is.null(sol) || !isTRUE(sol$ok)) return(NULL)
        ptc <- as.numeric(point[seq_len(nLeaves)])
        for (r in names(sol$rates)) { s <- slotOfName(r)
          if (!is.null(s)) ptc[s] <- as.numeric(sol$rates[[r]]) %% p }
        for (nm in realStateNames) { s <- slotOfName(nm); v <- sol$valBy[[nm]]
          if (!is.null(s) && !is.null(v)) ptc[s] <- as.numeric(v) %% p }
        list(sol = sol, ptc = ptc)
      }
      # One condition's [Obs; df] blocks at its own resting state, state columns
      # log-normalised by x*_c, embedded in the wide space. Shared by kcall4 and kchunk.
      # NULL: degenerate point.
      oneCondBlocks <- function(mi, obs, sc0, p) {
        if (!isTRUE(obs$ok)) return(NULL)
        oR <- matrix(as.numeric(obs$R), nrow = obs$rank, ncol = nzL)
        dR <- dfRowsCond(sc0$sol, p)
        # a resting value 0 mod p degenerates its column; reject the point
        xvals <- as.numeric(sc0$ptc[jointStateSlot]) %% p
        if (any(xvals == 0)) {
          if (nzchar(Sys.getenv("DMOD_SYM_FWDDIAG")))
            message("[fwddiag] zero log-normal coord at slots ",
                    paste(jointStateSlot[xvals == 0], collapse = ","))
          return(NULL) }
        for (m in seq_along(logCols)) {
          xv <- xvals[m]
          oR[, logCols[m]] <- .symMulmod(oR[, logCols[m]], xv, p)
          if (nrow(dR)) dR[, logCols[m]] <- .symMulmod(dR[, logCols[m]], xv, p)
        }
        # state columns go to block mi, shared columns overlap
        wc <- wideCols[[mi]]
        bl <- list()
        eO <- matrix(0, nrow(oR), nzWide); eO[, wc] <- oR
        bl[[length(bl) + 1L]] <- eO
        if (nrow(dR)) {
          eD <- matrix(0, nrow(dR), nzWide); eD[, wc] <- dR
          bl[[length(bl) + 1L]] <- eD
        }
        # recast relation rows at this condition's point
        for (rc in recastRel) {
          expv <- as.numeric(sc0$ptc[rc$expSlot]) %% p
          Lv   <- as.numeric(sc0$ptc[rc$LSlot]) %% p
          # 1/base factor for a non-log-normalised parameter base (state base: 1)
          bScale <- 1
          if (isTRUE(rc$baseParam)) {
            bv <- as.numeric(sc0$ptc[rc$baseSlot]) %% p
            if (bv == 0) return(NULL)                  # degenerate point, resample
            bScale <- .symInvmod(bv, p)
          }
          r1 <- numeric(nzWide)                       # E = base^exp
          r1[wc[rc$ECol]]    <- 1
          r1[wc[rc$baseCol]] <- (p - .symMulmod(expv, bScale, p)) %% p   # -exp/base
          r1[wc[rc$expCol]]  <- (-Lv) %% p
          r2 <- numeric(nzWide)                       # L = log(base)
          r2[wc[rc$LCol]]    <- 1
          r2[wc[rc$baseCol]] <- (p - bScale) %% p                         # -1/base
          bl[[length(bl) + 1L]] <- rbind(r1, r2)
        }
        # the same blocks over the gap series: the chain's series rows scaled and
        # embedded like oR, every other block constant
        if (!is.null(obs$S)) {
          scl <- rep(1, nzL); scl[logCols] <- xvals
          S <- matrix(as.numeric(obs$S), nrow(obs$S), ncol(obs$S))
          cst <- if (length(bl) > 1L) do.call(rbind, bl[-1]) else matrix(0, 0, nzWide)
          attr(bl, "series") <- rbind(.symSeriesEmbed(S, obs$N, wc, nzWide, p, scl),
                                      .symSeriesConst(cst, obs$N))
        }
        bl
      }
      # one condition's observability kernel at its seeded resting state: a chain of
      # segments across post-t0 event gaps, or the single equilibrate segment.
      condObs <- function(ci, ptc, p, Nt, Mtot) {
        if (hasGaps)
          symObsNullChain(list(tapes[chainGroups[[chainOf[ci]]]]), nLeaves, nStates,
                          zSlotsL, as.integer(ptc %% p), p, as.integer(Nt),
                          as.integer(Mtot), 1L)
        else
          symObsNullMulti(list(tapes[[ci]]), nLeaves, nStates, zSlotsL,
                          as.integer(ptc %% p), p, as.integer(Nt), 1L)
      }
      kcall4 <- function(point, p, Nt, Mtot = 0L, solveFn = jointSolveCond) {
        # Warm this point's Kc solves in parallel (backward solve only), but only after
        # the first condition solved serially: a degenerate point fails there.
        if (missing(solveFn) && Kc > 1L &&
            !is.null(jointSolveCond(point, p, equilConds[1])))
          warmSolves(list(point), list(p))
        blocks <- list(); sblocks <- list()
        for (mi in seq_len(Kc)) {
          ci <- equilConds[mi]
          sc0 <- solveFn(point, p, ci)
          if (is.null(sc0)) { ssWhy <<- "joint solve failed"; return(list(ok = FALSE)) }
          b <- oneCondBlocks(mi, condObs(ci, sc0$ptc, p, Nt, Mtot), sc0, p)
          if (is.null(b)) return(list(ok = FALSE))
          sblocks <- c(sblocks, list(attr(b, "series")))
          blocks <- c(blocks, b)
        }
        # stacked reduction over GF(p); only the pivot rows are needed
        rr <- symRrefMod(do.call(rbind, blocks), p)
        res <- list(ok = TRUE, R = rr$R,
                    pivots = as.integer(rr$piv), rank = as.integer(rr$rank), dim = nzWide)
        if (hasGaps) res <- .symSeriesStack(res, sblocks, nzWide, Mtot + 1L, p)
        res
      }
      # one condition's observability rows for the saturation, without the constant
      # df and recast rows
      blockCall <- lapply(seq_len(Kc), function(mi) { force(mi)
        function(point, p, Nt, Mtot = 0L) {
          sc0 <- jointSolveCond(point, p, equilConds[mi])
          if (is.null(sc0)) { ssWhy <<- "joint solve failed"; return(list(ok = FALSE)) }
          condObs(equilConds[mi], sc0$ptc, p, Nt, Mtot)
        } })
      # ---- the batched twin of the serial per-point loop ---------------------------
      # Coupled + gap path: solve every (point, condition) seed, evaluate all chain
      # kernels in one OpenMP batch, then assemble and reduce per point. Identical to
      # looping kcall4.
      if (hasGaps) {
        chainsList <- lapply(seq_len(Kc), function(mi)
          tapes[chainGroups[[chainOf[equilConds[mi]]]]])
        # no per-point fork here, so the batch takes the full `cores`
        coresChunk <- cores
        kchunk <- function(pointList, primeVec, Nt) {
          warmSolves(pointList, primeVec)
          nP <- length(pointList)
          perCond <- vector("list", nP)   # per point: list of Kc sc0, or NULL if any fails
          seedRows <- list(); evChain <- integer(0); evPrime <- numeric(0)
          for (i in seq_len(nP)) {
            pt <- pointList[[i]]; pp <- primeVec[[i]]
            per <- vector("list", Kc); okAll <- TRUE
            for (mi in seq_len(Kc)) {
              sc0 <- jointSolveCond(pt, pp, equilConds[mi])
              if (is.null(sc0)) { okAll <- FALSE; break }
              per[[mi]] <- sc0
            }
            if (!okAll) next
            perCond[[i]] <- per
            for (mi in seq_len(Kc)) {
              seedRows[[length(seedRows) + 1L]] <- as.integer(per[[mi]]$ptc %% pp)
              evChain <- c(evChain, mi - 1L); evPrime <- c(evPrime, pp)
            }
          }
          kr <- if (length(seedRows))
            symObsNullChainSeedBatch(chainsList, as.integer(evChain),
              do.call(rbind, seedRows), as.numeric(evPrime), nLeaves, nStates,
              zSlotsL, as.integer(Nt), as.integer(MtotUsed),
              as.integer(min(coresChunk, length(evChain))))
            else list()
          out <- vector("list", nP); e <- 0L
          for (i in seq_len(nP)) {
            if (is.null(perCond[[i]])) { out[[i]] <- list(ok = FALSE); next }
            blocks <- list(); sblocks <- list(); bad <- FALSE
            for (mi in seq_len(Kc)) {
              e <- e + 1L
              b <- oneCondBlocks(mi, kr[[e]], perCond[[i]][[mi]], primeVec[[i]])
              if (is.null(b)) bad <- TRUE
              else { blocks <- c(blocks, b); sblocks <- c(sblocks, list(attr(b, "series"))) }
            }
            out[[i]] <- if (bad) list(ok = FALSE) else {
              rr <- symRrefMod(do.call(rbind, blocks), primeVec[[i]])
              .symSeriesStack(list(ok = TRUE, R = rr$R, pivots = as.integer(rr$piv),
                                   rank = as.integer(rr$rank), dim = nzWide),
                              sblocks, nzWide, MtotUsed + 1L, primeVec[[i]]) }
          }
          out
        }
      }
      # from here on the analysis runs over the wide per-condition coordinate space;
      # the kernel/df closures above keep the local names (nzL/znamesL/zSlotsL)
      nz <- nzWide; znames <- znamesWide; zSlots <- zSlotsWide
      multi$zStateNames <- zStateNamesWide
      jointStateBase <- stateBase        # original state names, for the scaling peel
      jointKc <- Kc
      # per-condition state and recast leaves are auxiliary and exempt from the
      # relevance gate; the dense-fit cap grows by the state count to match
      auxLeaves <- which(leafNamesAug %in% setdiff(leafNames, paramNames))
      ctrl$relevanceCap <- ctrl$relevanceCap + nSt
    }
  } else if (length(multi$powerRecast) || length(multi$expAtoms)) {
    # Transient recast (free exponent without equilibrate): E and L are free leaves tied
    # by the linearised relations
    #   E:  xi_E - (exp*E/base) xi_base - (E*log base) xi_exp = 0
    #   L:  xi_L - (1/base) xi_base = 0
    # exact at a generic point since base, log(base) and base^exp are algebraically
    # independent. An exponential leaf W = r^(tau/L) adds xi_W - W log(r) grad(tau/L) . xi.
    recastTransient <- TRUE
    recast <- multi$powerRecast
    if (is.null(recast)) recast <- list()
    for (i in seq_along(recast)) recast[[i]]$inverted <- FALSE
    recastAtomNames <- c(as.character(multi$recastAtomNames),
                         as.character(multi$expAtoms))
    expRel <- multi$expRelation
    expRelFn <- if (length(multi$expAtoms)) .symTapeFn(expRel, nLeaves) else NULL
    auxLeaves <- which(leafNamesAug %in% recastAtomNames)
    # ---- assembling the joint kernel call kcall4 -----------------------------------
    slotOfL <- function(nm) { i <- match(nm, leafNames); if (is.na(i)) NA_integer_ else i }
    colOfN  <- function(nm) { i <- match(nm, znames);    if (is.na(i)) NA_integer_ else i }
    recastRelT <- lapply(recast, function(rc) list(
      Ecol = colOfN(rc$E), baseCol = colOfN(rc$base), expCol = colOfN(rc$exp),
      Lcol = colOfN(rc$L), Eslot = slotOfL(rc$E), baseSlot = slotOfL(rc$base),
      expSlot = slotOfL(rc$exp), Lslot = slotOfL(rc$L)))
    # solveFn unused; the signature matches the joint-branch kcall4
    kcall4 <- function(point, p, Nt, Mtot = 0L, solveFn = jointSolveCond) {
      pt <- as.integer(point)
      o <- if (hasGaps)
        symObsNullChain(lapply(chainGroups, function(idx) tapes[idx]), nLeaves,
                        nStates, zSlots, pt, p, as.integer(Nt), as.integer(Mtot), coresCall)
        else symObsNullMulti(tapes, nLeaves, nStates, zSlots, pt, p, as.integer(Nt), coresCall)
      if (!isTRUE(o$ok)) return(list(ok = FALSE))
      oR <- matrix(as.numeric(o$R), nrow = o$rank, ncol = nz)
      rel <- list(); seenL <- integer(0)
      for (rc in recastRelT) {
        base0 <- as.numeric(point[rc$baseSlot]) %% p
        if (base0 == 0) return(list(ok = FALSE))          # degenerate point, resample
        invb  <- .symInvmod(base0, p)
        E0    <- as.numeric(point[rc$Eslot])  %% p
        exp0  <- as.numeric(point[rc$expSlot]) %% p
        L0    <- as.numeric(point[rc$Lslot])  %% p
        r1 <- numeric(nz)                                 # E = base^exp
        r1[rc$Ecol]    <- 1
        r1[rc$baseCol] <- (p - .symMulmod(exp0, .symMulmod(E0, invb, p), p)) %% p
        r1[rc$expCol]  <- (p - .symMulmod(E0, L0, p)) %% p
        rel[[length(rel) + 1L]] <- r1
        if (!(rc$Lcol %in% seenL)) {                      # L = log(base), once per base
          seenL <- c(seenL, rc$Lcol)
          r2 <- numeric(nz)
          r2[rc$Lcol]    <- 1
          r2[rc$baseCol] <- (p - invb) %% p
          rel[[length(rel) + 1L]] <- r2
        }
      }
      if (!is.null(expRelFn)) {
        eM <- .symExpRows(expRel, expRelFn, znames, pt, p)
        if (is.null(eM)) return(list(ok = FALSE))         # degenerate point, resample
        rel <- c(rel, lapply(seq_len(nrow(eM)), function(i) eM[i, ]))
      }
      relM <- do.call(rbind, rel)
      rr <- symRrefMod(rbind(oR, relM), p)
      res <- list(ok = TRUE, R = rr$R, pivots = as.integer(rr$piv),
                  rank = as.integer(rr$rank), dim = nz)
      if (hasGaps)
        res <- .symSeriesStack(res, list(matrix(as.numeric(o$S), nrow(o$S), ncol(o$S)),
                                         .symSeriesConst(relM, o$N)), nz, o$N, p)
      res
    }
    blockCall <- obsBlockCalls()
  } else {
    # solveFn unused; the signature matches the joint-branch kcall4
    kcall4 <- function(point, p, Nt, Mtot = 0L, solveFn = jointSolveCond) {
      pt <- as.integer(point)
      if (hasGaps)
        symObsNullChain(lapply(chainGroups, function(idx) tapes[idx]), nLeaves,
                        nStates, zSlots, pt, p, as.integer(Nt), as.integer(Mtot), coresCall)
      else
        symObsNullMulti(tapes, nLeaves, nStates, zSlots, pt, p, as.integer(Nt), coresCall)
    }
    blockCall <- obsBlockCalls()
  }

  # the saturation loop rejects a point on its first condition, so only that solve is
  # prefetched; kcall4 fills the rest once the point passes
  warmProbe <- if (jointSS)
    function(pts, primes) warmSolves(pts, primes, conds = equilConds[1])
    else function(pts, primes) invisible()
  # recast rows add to the specialisation; a gap chain gets no budget
  satBudget <- if (hasGaps || is.na(codimSpec)) NA_integer_
               else as.integer(codimSpec) + 2L * length(recast) +
                    as.integer(if (is.null(multi$expCodim)) 0L else multi$expCodim)
  sc <- .symSaturateCertify(kcall4, nAug, nz, maxM, warm = warmProbe,
                            probeBlock = max(1L, min(8L, coresGLp)),
                            blockCall = blockCall, budget = satBudget)
  if (is.null(sc)) {
    if (ssConstraint && !is.null(ssWhy))
      warning("symmetryDetection(): no steady-state point over the finite field ",
              "after the generic-point retries (", ssWhy, "). The equilibrate ",
              "constraint could not be evaluated.", call. = FALSE)
    return(NULL)
  }
  # joint mode: log-normalised state columns have unit value, so their slots are set
  # to 1 for the peel; the solved point is kept for the verify guard
  point0Solved <- sc$point0
  if (jointSS)
    for (col in which(znames %in% as.character(multi$zStateNames)))
      sc$point0[zSlots[col] + 1L] <- 1
  # the reconstruction samples at the saturated gap order, baked into a 3-arg kcall
  MtotUsed <- if (is.null(sc$MtotUsed)) 0L else as.integer(sc$MtotUsed)
  # ==== the kernel drivers: per-point kcall, batched kbatch, forward kcallFwd =======
  kcall <- function(point, p, Nt) kcall4(point, p, Nt, MtotUsed)
  # forward kernel (joint mode): the same [Obs; df] with the forward solve for the
  # turnover rates `solveRates`, so all primes share one slice
  kcallFwd <- if (jointSS) function(point, p, Nt, solveRates)
      kcall4(point, p, Nt, MtotUsed, function(pt, pp, cc) jointSolveCondFwd(pt, pp, cc, solveRates))
    else NULL
  # kernel over many (point, prime) pairs: one OpenMP batch on the plain path; the
  # constraint, gap and recast paths need per-point seeding or extra rows and loop kcall
  canBatch <- !ssConstraint && !hasGaps && !recastTransient
  kbatch <- function(pointList, primeVec, Nt) {
    if (canBatch && length(pointList) > 0L) {
      M <- do.call(rbind, lapply(pointList, function(pp)
        as.integer(pp[seq_len(nLeaves)])))
      symObsNullBatch(tapes, nLeaves, nStates, zSlots, M, as.numeric(primeVec),
                      as.integer(Nt), cores)
    } else if (!is.null(kchunk) && length(pointList) > 0L &&
               !nzchar(Sys.getenv("DMOD_SYM_NOCHUNK"))) {
      # coupled + gap path in one OpenMP batch; DMOD_SYM_NOCHUNK forces the serial loop
      kchunk(pointList, primeVec, Nt)
    } else {
      # warm the solves, then the per-point kernel; inside a unix fork the kernel runs
      # serial, restored on exit
      warmSolves(pointList, primeVec)
      if (coresGLp > 1L && .Platform$OS.type == "unix") {
        ccSaved <- coresCall; coresCall <<- 1L
        # no nested warm fill inside a child
        wsSaved <- warmOff; warmOff <<- TRUE
        on.exit({coresCall <<- ccSaved; warmOff <<- wsSaved}, add = TRUE)
      }
      parMap(seq_along(pointList),
             function(i) kcall(pointList[[i]], primeVec[[i]], Nt))
    }
  }
  if (hasGaps && isFALSE(sc$saturatedM))
    warning("symmetryDetection(): the gap power-series order hit the cap (",
            maxM, ") before the rank stabilised; the reported rank is a sound ",
            "(conservative) lower bound. Raise reconstControl(gapOrderCap=) for ",
            "the exact generic-timing verdict.", call. = FALSE)

  result <- list(method = "observability", engine = "analytic",
                 conditions = if (is.null(nConditions)) nSegmentTapes
                              else as.integer(nConditions),
                 coordinates = znames,
                 segments = nSegmentTapes, gapOrderUsed = MtotUsed,
                 identifiable = (sc$rankS == nz), rank = as.integer(sc$rankS),
                 dim = as.integer(nz), lieOrderUsed = as.integer(sc$NtUsed),
                 lieOrderDriver = sc$blockDriver, lieBudget = sc$budget,
                 liePlateau = sc$plateau, lieCertified = isTRUE(sc$certified),
                 nonIdentifiable = list())
  class(result) <- "symmetrydetection"
  # joint mode reports in parameter space; full wide rank means all identifiable (the
  # other case is reprojected below)
  if (jointSS) {
    physParams0 <- setdiff(znames, as.character(multi$zStateNames))
    result$dim <- length(physParams0)
    if (sc$rankS == nz) result$rank <- length(physParams0)
  }
  # transient recast reports without the atoms E and L
  if (recastTransient) {
    physCoords0 <- setdiff(znames, recastAtomNames)
    result$dim <- length(physCoords0)
    result$coordinates <- physCoords0
    if (sc$rankS == nz) result$rank <- length(physCoords0)
  }
  # directions depending on the time between events (series rank below the stacked
  # rank) have no closed form and are reported by their support
  gapDirs <- list()
  if (sc$rankS < sc$rank) {
    S <- sc$ref$S
    supp <- symSeriesRank(S, as.integer(nz), as.integer(sc$ref$N), .symPrimes[1],
                          integer(0), support = TRUE)$support
    gapDirs <- rep(list(list(support = .symSort(znames[supp + 1L]), type = "general",
                             closedForm = FALSE,
                             reason = "depends on the time between events")),
                   sc$rank - sc$rankS)
  }
  if (sc$rankS == nz)
    return(result)
  # joint mode reports in parameter space, transient recast without its atoms
  reportPhysical <- function(result) {
    if (jointSS)
      result <- .symReportPhysical(result, znames, as.character(multi$zStateNames),
                                     sc, nz, .symPrimes[1], "stateVector", flagJoint = TRUE)
    if (recastTransient)
      result <- .symReportPhysical(result, znames, recastAtomNames,
                                     sc, nz, .symPrimes[1], "recastVector")
    result
  }
  if (sc$rank == nz) {
    result$nonIdentifiable <- gapDirs
    return(reportPhysical(result))
  }

  P <- .symPrimes[1]
  freeCols <- setdiff(0:(nz - 1L), sc$pivots)
  N <- .symNullspaceBasis(sc$ref, freeCols, P)

  # ==== peel the exact scalings common to every condition ===========================
  # closed form from the integer kernel; their span is excluded before reconstruction
  scaling <- list(); Bmat <- matrix(0L, nz, 0L)
  if (!is.null(sd)) {
    # candidates from each segment's regime, validated against the nullspace N
    modelLines <- lapply(multi$tapes, function(t) as.character(t$modelLines))
    obsLines   <- lapply(multi$tapes, function(t) as.character(t$obsLines))
    # duplicate regimes add nothing to the intersection of scaling lattices
    regimeKey <- vapply(seq_along(modelLines), function(i)
      paste(c(modelLines[[i]], "\x1f", obsLines[[i]]), collapse = "\n"),
      character(1))
    keepRegime <- !duplicated(regimeKey)
    modelLines <- modelLines[keepRegime]
    obsLines   <- obsLines[keepRegime]
    inputs <- if (!is.null(multi$forcings)) as.character(multi$forcings) else NULL
    # baked leaves, against the original names in joint mode (znames is wide there)
    fixed <- if (jointSS)
      setdiff(leafNames, c(setdiff(znames, as.character(multi$zStateNames)), jointStateBase))
      else setdiff(leafNames, znames)
    # recast atoms impose c_E = exp*c_base, giving parameter-weighted (Hill) scalings
    # over Q(exp) without the rational fit
    scalRes <- tryCatch(sd$scalingSymmetriesMulti(
      perCondModel = modelLines, perCondObs = obsLines,
      inputs = if (length(inputs)) inputs else NULL,
      fixed = if (length(fixed)) fixed else NULL,
      recast = if (length(recast)) recast else NULL), error = function(e) NULL)
    if (!is.null(scalRes)) {
      # state weights onto every per-condition column
      if (jointSS) scalRes <- .symJointExpandScal(scalRes, jointStateBase, jointKc,
                                                     heldParamOf)
      peel <- .symPeelScalings(scalRes, znames, nz, sc$point0[zSlots + 1L], P, N,
                                 sd = sd)
      scaling <- peel$scaling; Bmat <- peel$Bmat
    }
  }

  scalCols <- ncol(Bmat)          # scaling tangents (fixed before the residual fit)

  # ==== residual directions: what the scalings do not span ==========================
  residualFree <- integer(0)
  for (fc in freeCols) {
    bf <- .symNullResidues(sc$ref, fc, P)
    if (!.symInSpan(Bmat, bf, nz, P)) {
      residualFree <- c(residualFree, fc); Bmat <- cbind(Bmat, bf) }
  }
  scalRows <- if (scalCols > 0L) t(Bmat[, seq_len(scalCols), drop = FALSE])
              else matrix(0L, 0L, nz)

  # ==== closed-form reconstruction of the residual directions =======================
  if (isTRUE(closedForm)) {
    .t0 <- Sys.time()
    to <- if (is.null(ctrl$timeout)) Inf else ctrl$timeout
    ctrl$deadline <- if (is.finite(to)) .t0 + to else NULL
    .tlog <- if (nzchar(Sys.getenv("DMOD_SYM_TIMING")))
      function(msg) message(sprintf("[sym %6.1fs] %s",
                                    as.numeric(Sys.time() - .t0, units = "secs"), msg))
      else function(msg) invisible()
    # the solve fill has run since the saturation loop; report its mode
    if (jointSS && length(residualFree))
      .tlog(sprintf("solve fill: %s",
                    if (coresGLp <= 1L) "serial"
                    else if (.Platform$OS.type == "unix")
                      sprintf("fork, %d cores", coresGLp)
                    else sprintf("PSOCK pool, up to %d workers", coresGLp)))
    .tlog(sprintf("start: %d residual direction(s)", length(residualFree)))
    # Relevance probe shared by all directions: one kernel per one-leaf perturbation,
    # retried until the pivots match (a pivot shift would mark the leaf relevant everywhere).
    probeNext <- sc$poolNext
    relProbe <- vector("list", nAug)
    pending <- seq_len(nAug)
    if (length(residualFree)) for (att in seq_len(ctrl$probeRetries)) {
      if (!length(pending) || .symExpired(ctrl)) break
      perts <- lapply(pending, function(li) {
        pert <- sc$point0
        pert[li] <- sc$pool(probeNext + (li - 1L) * ctrl$probeRetries + (att - 1L))
        pert
      })
      # independent per-leaf kernels in one batch
      cands <- kbatch(perts, rep(P, length(perts)), sc$NtUsed)
      resolved <- logical(length(pending))
      for (j in seq_along(pending)) {
        li <- pending[j]; pert <- perts[[j]]; cand <- cands[[j]]
        relProbe[[li]] <- list(
          rp = cand, zvals = if (length(zSlots)) pert[zSlots + 1L] else NULL,
          pertval = pert[li])
        if (!is.null(cand) && isTRUE(cand$ok) &&
            identical(as.integer(cand$pivots), as.integer(sc$pivots)))
          resolved[j] <- TRUE
      }
      pending <- pending[!resolved]
    }
    poolNext <- probeNext + nAug * ctrl$probeRetries
    .tlog("relevance probe done")

    zvals0 <- if (length(zSlots)) sc$point0[zSlots + 1L] else NULL
    # physical columns (0-based) the per-prime joint fit reconstructs; auxiliary ones skip
    physColsPP <- if (jointSS)
      which(!(znames %in% as.character(multi$zStateNames))) - 1L else NULL
    # one direction in a gauge (free-column when rfn is NULL): interpolate, verify,
    # else downgrade to support-only
    reconstructOne <- function(fc, rfn, logCoords = FALSE, sharedBank = NULL,
                               fastOnly = FALSE, perPrime = FALSE) {
      if (.symExpired(ctrl)) {
        v <- if (!is.null(rfn)) rfn(sc$ref, P, zvals0) else .symNullResidues(sc$ref, fc, P)
        if (is.null(v)) v <- .symNullResidues(sc$ref, fc, P)
        return(list(support = .symSort(znames[v != 0]), type = "general",
                    closedForm = FALSE,
                    reason = "reconstruction time budget (reconstControl(timeout=)) exceeded"))
      }
      # coupled steady state: a perturbation solves at only some primes, so samples are
      # collected per prime and lifted by CRT (free-column gauge only)
      dir <- if (isTRUE(perPrime))
        .symInterpolatePerprime(fc, sc$ref, sc$pivots, znames, zSlots,
                                  leafNamesAug, nAug, sc$point0, sc$pool, poolNext,
                                  sc$NtUsed, kcall, kbatch, spy, relProbe, ctrl,
                                  auxLeaves = auxLeaves, physCols = physColsPP)
        else
        .symInterpolateDirection(fc, sc$ref, sc$pivots, znames, zSlots,
                                   leafNamesAug, nAug, sc$point0, sc$pool,
                                   poolNext, sc$NtUsed, kcall, spy, relProbe,
                                   rfn, ctrl, kbatch, sharedBank, fastOnly,
                                   auxLeaves = auxLeaves)
      poolNext <<- dir$poolNext
      e <- dir$entry
      # relevant leaves for the strict verifier, removed from the reported direction
      ppRel <- e$relevantLeaves; e$relevantLeaves <- NULL
      # log-gauge entries eta = xi / z back to xi before verifying
      if (logCoords && isTRUE(e$closedForm))
        e$vector <- .symLogcoordBacksub(e$vector, spy)
      ok <- isTRUE(e$closedForm) && (
        if (isTRUE(perPrime))
          # coupled path: strict verification at a genuine steady-state point
          .symVerifyPerprime(e, fc, znames, leafNamesAug, sc$point0, sc$NtUsed,
                               kcall, sc$pivots, sc$pool, poolNext, nz, sd,
                               relLeaves = ppRel)
        else if (logCoords)
          .symVerifyInNullspace(e, fc, znames, leafNamesAug, sc$point0,
                                   sc$NtUsed, kcall, sc$pool, poolNext, nz, sd)
        else
          .symVerifyDirection(e, fc, znames, leafNamesAug, sc$point0,
                                sc$NtUsed, kcall, sd, rfn) &&
          (is.null(rfn) ||
           .symVerifyInNullspace(e, fc, znames, leafNamesAug, sc$point0,
                                    sc$NtUsed, kcall, sc$pool, poolNext, nz, sd)))
      if (isTRUE(e$closedForm) && !ok) {
        v <- if (!is.null(rfn)) rfn(sc$ref, P, zvals0) else .symNullResidues(sc$ref, fc, P)
        if (is.null(v)) v <- .symNullResidues(sc$ref, fc, P)
        e <- list(support = .symSort(znames[v != 0]), type = "general",
                  closedForm = FALSE,
                  reason = paste("a closed form was reconstructed but failed",
                                 "verification at a fresh prime"))
      }
      if (isTRUE(e$closedForm) && length(recast) && !is.null(sd))
        e$vector <- .symRecastBacksub(e$vector, recast, sd)
      if (isTRUE(e$closedForm) && length(multi$expBack$names) && !is.null(sd)) {
        e$vector <- .symExpBacksub(e$vector, multi$expBack, sd)
        v <- tryCatch(sd$dropMonomialContent(e$vector), error = function(err) NULL)
        if (!is.null(v)) e$vector <- lapply(v[names(e$vector)], as.character)
      }
      e
    }

    allClosed <- function(rec) all(vapply(rec, function(e) isTRUE(e$closedForm),
                                          logical(1)))

    # peel minimal-support cocircuits that close by the cheap log-coordinate read-off
    # (e.g. a Hill-weighted scaling); only the rest reach the wide free-column fit
    peeled <- list()
    msPeel <- .symMinsupportGauge(residualFree, scalRows, P, nz, sc, freeCols,
                                    candCap = ctrl$minsupportCandCap)
    if (length(msPeel$anchors) && !is.null(msPeel$vectors)) {
      peelVec <- matrix(0L, nz, 0L)
      for (gi in seq_along(msPeel$anchors)) {
        if (is.null(msPeel$residueFns[[gi]])) next
        e <- reconstructOne(msPeel$anchors[gi], msPeel$residueFns[[gi]],
                            logCoords = TRUE, fastOnly = TRUE)
        if (isTRUE(e$closedForm)) {
          peeled[[length(peeled) + 1L]] <- e
          peelVec <- cbind(peelVec, msPeel$vectors[, gi])
        }
      }
      if (length(peeled)) {
        # keep a basis completion of the residual set independent of the peeled span
        keep <- integer(0); span <- peelVec
        for (fc in residualFree) {
          bf <- .symNullResidues(sc$ref, fc, P)
          if (.symInSpan(span, bf, nz, P)) next
          keep <- c(keep, fc); span <- cbind(span, bf)
        }
        residualFree <- keep
      }
    }
    .tlog(sprintf("peel done: %d peeled, %d remaining", length(peeled), length(residualFree)))

    # minimal-support gauge over the remaining residual set, in log coordinates;
    # `fastOnly` reads a monomial off the base point and probe without sampling
    ms <- .symMinsupportGauge(residualFree, scalRows, P, nz, sc, freeCols,
                                candCap = ctrl$minsupportCandCap)
    msTry <- function(fastOnly) {
      if (length(ms$anchors) != length(residualFree) ||
          all(vapply(ms$residueFns, is.null, logical(1)))) return(NULL)
      rec <- lapply(seq_along(ms$anchors), function(gi)
        reconstructOne(ms$anchors[gi], ms$residueFns[[gi]], logCoords = TRUE,
                       fastOnly = fastOnly))
      if (allClosed(rec)) rec else NULL
    }

    # the no-sampling read-off first; it closes parameter-weighted scalings cheaply
    interp <- msTry(TRUE)
    .tlog(if (is.null(interp)) "ms fastOnly: no close, going dense"
          else "ms fastOnly: all closed")
    if (is.null(interp)) {
      # free-column gauge with one dense sample bank over the union of relevant leaves,
      # shared by all directions
      metas <- lapply(residualFree, function(fc)
        .symDirectionRelevance(fc, sc$ref, sc$pivots, zSlots, nAug, sc$point0,
                                 relProbe, NULL))
      needs <- vapply(metas, function(m) .symDenseNeed(m$relByEntry, ctrl)$maxNeed,
                      integer(1))
      denseDir <- which(needs > 0L)
      .tlog(sprintf("relevance/need: %d dense dir(s), needs=[%s]",
                    length(denseDir), paste(needs, collapse = ",")))
      bank <- NULL
      if (length(denseDir)) {
        unionLeaves <- sort(unique(unlist(lapply(metas[denseDir], function(m) m$relevant))))
        bk <- .symBuildSharedBank(unionLeaves, max(needs[denseDir]), sc$point0,
                                     sc$pool, poolNext, kbatch, sc$pivots, sc$NtUsed,
                                     length(.symPrimes), ctrl)
        poolNext <- bk$poolNext
        if (isTRUE(bk$ok)) bank <- bk
        .tlog(sprintf("shared bank built (%d leaves, ok=%s)",
                      length(unionLeaves), isTRUE(bk$ok)))
      }
      # an equilibrate model rarely fills the all-prime bank (a point seldom solves at
      # every prime), so its directions fall back to per-prime reconstruction
      bankMissing <- ssConstraint && length(denseDir) > 0L && is.null(bank)
      interp <- lapply(seq_along(residualFree), function(ii) {
        fc <- residualFree[ii]
        r <- if (bankMissing) reconstructOne(fc, NULL, perPrime = TRUE)
             else reconstructOne(fc, NULL, sharedBank = bank)
        if (!isTRUE(r$closedForm) && ssConstraint && !bankMissing)
          r <- reconstructOne(fc, NULL, perPrime = TRUE)
        .tlog(sprintf("dense dir %d/%d done (closed=%s)", ii, length(residualFree),
                      isTRUE(r$closedForm)))
        r
      })

      # Coupled recast directions the pivot-pinned gauge cannot sample (perturbing the
      # recast leaf shifts the pivots): the forward path, adopted only if it closes every
      # residual direction; otherwise the support-only result stays.
      if (ssConstraint && !is.null(sd) && !is.null(kcallFwd) && length(recast) &&
          !allClosed(interp)) {
        fwd <- tryCatch(
          .symPerprimeForwardMulti(residualFree, sc, kcall, kcallFwd, znames, zSlots,
            leafNamesAug, nz, scaling, as.character(multi$zStateNames),
            as.character(multi$paramNames), recast, sd, spy, ctrl, physColsPP,
            models, realStateNames, solveParamNames, solveHeld),
          error = function(e) NULL)
        fwdClosed <- if (is.null(fwd)) list()
                     else Filter(function(e) isTRUE(e$closedForm), fwd)
        if (length(fwdClosed) == length(residualFree) && length(fwdClosed)) {
          interp <- lapply(fwdClosed, function(e) {
            if (length(recast)) e$vector <- .symRecastBacksub(e$vector, recast, sd)
            e })
          .tlog(sprintf("forward-multi closed %d residual direction(s)", length(fwdClosed)))
        }
      }

      # rescue directions the free-column gauge entangles: retry each open one in the
      # canonical gauge, then in log coordinates; each closed form is verified separately
      rescueEach <- function(gauge, logCoords = FALSE) {
        if (length(gauge$anchors) != length(residualFree)) return(invisible())
        for (ii in which(!vapply(interp, function(e) isTRUE(e$closedForm), logical(1)))) {
          if (is.null(gauge$residueFns[[ii]]) || .symExpired(ctrl)) next
          rec <- reconstructOne(gauge$anchors[ii], gauge$residueFns[[ii]],
                                logCoords = logCoords)
          if (isTRUE(rec$closedForm)) interp[[ii]] <<- rec
        }
        invisible()
      }
      # the rescues need the all-prime bank, so not for a coupled steady-state model
      if (!bankMissing) {
        if (!allClosed(interp))
          rescueEach(.symCanonGauge(residualFree, scalRows, P, nz, sc))
        if (!allClosed(interp))
          rescueEach(.symLogcoordGauge(residualFree, scalRows, P, nz, sc, zvals0), TRUE)
        # final rescue: the minimal-support gauge with kernel sampling, for a pinned
        # entry that is a bounded-degree rational rather than a bare monomial.
        if (!allClosed(interp)) {
          .tlog("rescues exhausted, trying ms sampling")
          msrec <- msTry(FALSE)
          if (!is.null(msrec)) interp <- msrec
        }
      }
    }
    .tlog(sprintf("reconstruction done: %d/%d closed",
                  sum(vapply(interp, function(e) isTRUE(e$closedForm), logical(1))),
                  length(interp)))
    result$nonIdentifiable <- c(scaling, peeled, interp, gapDirs)
  } else {
    support <- lapply(residualFree, function(fc) {
      v <- .symNullResidues(sc$ref, fc, P)
      list(support = .symSort(znames[v != 0]), type = "general", closedForm = FALSE)
    })
    result$nonIdentifiable <- c(scaling, support, gapDirs)
  }
  # ==== report in the physical coordinate space =====================================
  # Joint mode reports in parameter space: the state components, fixed by the parameters
  # on the resting manifold, move to $stateVector. Transient recast drops only E and L
  # (already back-substituted); state initial values stay physical coordinates.
  result <- reportPhysical(result)
  # ==== the saturation guard (verify = TRUE) and the return value ===================
  # where the budget did not certify the Lie order, check that the rank does not grow
  # further up (no second analysis)
  if (isTRUE(verify) && !isTRUE(sc$certified))
    result$verification <- tryCatch(
      .symSzSaturationGuard(kcall, point0Solved, sc$NtUsed, sc$rankS),
      error = function(e) list(ok = NA, method = "saturation guard",
                               reason = conditionMessage(e)))
  result
}


# ---- reconstruction internals: rational-entry interpolation & sampling ---------------

# Laurent candidate exponents: numerator monomials (degree <= dNum) minus single
# denominator monomials (degree <= dDen)
.symLaurentCandidates <- function(nvar, dNum, dDen) {
  num <- .symMonoTable(nvar, dNum)
  den <- .symMonoTable(nvar, dDen)
  cand <- do.call(rbind, lapply(seq_len(nrow(den)),
                                function(d) sweep(num, 2L, den[d, ], "-")))
  cand <- unique(cand)
  storage.mode(cand) <- "integer"
  cand
}


# Lift one Ben-Or-Tiwari polynomial recovered per prime to the rationals: require
# the supports (exponent sets) to agree across primes, then Chinese-remainder the
# coefficients. Returns list(exps, num, den) of coefficient strings, or NULL.
.symBotReconcile <- function(perPrime) {
  np <- length(perPrime)
  key <- function(m) apply(m, 1L, paste, collapse = ",")
  k1 <- key(perPrime[[1]]$exps)
  t <- length(k1)
  coefMat <- matrix(0L, t, np)
  coefMat[, 1L] <- perPrime[[1]]$coeffs
  for (pj in seq_len(np)[-1]) {
    kp <- key(perPrime[[pj]]$exps)
    if (length(kp) != t || !setequal(kp, k1)) return(NULL)
    coefMat[, pj] <- perPrime[[pj]]$coeffs[match(k1, kp)]
  }
  rec <- symRatRecon(coefMat, as.integer(.symPrimes))
  if (any(rec$den == "0")) return(NULL)
  list(exps = perPrime[[1]]$exps, num = rec$num, den = rec$den)
}


# Laurent terms as one entry string over a common monomial denominator
.symLaurentAssemble <- function(perPrime, reli, leafNames, nvar) {
  rc <- .symBotReconcile(perPrime)
  if (is.null(rc)) return(NULL)
  expRows <- rc$exps
  t <- nrow(expRows)
  mu <- apply(expRows, 2L, function(col) max(0L, -min(col)))
  vars <- leafNames[reli]
  numStr <- .symPolyString(rc$num, rc$den,
                             expRows + matrix(mu, t, nvar, byrow = TRUE), vars)
  denStr <- .symMonoString(mu, vars)
  if (denStr == "1") numStr else paste0("(", numStr, ")/(", denStr, ")")
}


# One wide entry by sparse Laurent interpolation: geometric samples until the
# Ben-Or-Tiwari term count stabilises, terms identified at the smallest degrees that
# fit, coefficients lifted across primes. NULL if not a bounded Laurent polynomial.
.symSparseEntry <- function(reli, supportCol, f, point0, leafNames, NtUsed,
                              kcall, pivots, residueFn = NULL,
                              ctrl = reconstControl(), zSlots = NULL) {
  nvar <- length(reli)
  bases <- .symSieve(nvar)
  np <- length(.symPrimes)
  maxLen <- 2L * ctrl$termCap + 2L

  # grow geometric samples until the term count (BM order) stabilises
  seqs <- lapply(seq_len(np), function(.) integer(0))
  cur <- lapply(seq_len(np), function(.) rep(1, nvar))
  have <- 0L
  repeat {
    if (.symExpired(ctrl)) return(NULL)
    target <- min(if (have == 0L) 16L else 2L * have, maxLen)
    for (pj in seq_len(np)) {
      p <- .symPrimes[pj]
      for (k in (have + 1L):target) {
        if (.symExpired(ctrl)) return(NULL)
        pt <- point0; pt[reli] <- cur[[pj]]
        rp <- kcall(pt, p, NtUsed)
        if (!isTRUE(rp$ok) || !identical(as.integer(rp$pivots), as.integer(pivots)))
          return(NULL)
        zv <- if (is.null(zSlots)) NULL else pt[zSlots + 1L]
        nvp <- if (is.null(residueFn)) .symNullResidues(rp, f, p) else residueFn(rp, p, zv)
        if (is.null(nvp)) return(NULL)
        seqs[[pj]][k] <- nvp[supportCol + 1L]
        cur[[pj]] <- (cur[[pj]] * bases) %% p
      }
    }
    have <- target
    if (2L * symBMorder(seqs[[1]], .symPrimes[1]) < have) break
    # no stable term count: not a bounded Laurent polynomial
    if (have >= maxLen) return(NULL)
  }

  # smallest degrees first; dDen = 0 is the polynomial case
  for (dDen in 0:ctrl$laurentDegDen) for (dNum in seq_len(ctrl$laurentDegNum)) {
    if (choose(nvar + dNum, nvar) * choose(nvar + dDen, nvar) > ctrl$laurentCandCap)
      next
    candM <- matrix(as.integer(.symLaurentCandidates(nvar, dNum, dDen)), ncol = nvar)
    perPrime <- vector("list", np)
    okAll <- TRUE
    for (pj in seq_len(np)) {
      monoRes <- symMonoResidues(candM, as.integer(bases), .symPrimes[pj])
      res <- symSparsePoly(seqs[[pj]], candM, as.integer(monoRes), .symPrimes[pj])
      if (!identical(res$status, "ok") || res$nterms == 0L) { okAll <- FALSE; break }
      perPrime[[pj]] <- res
    }
    if (!okAll) next
    out <- .symLaurentAssemble(perPrime, reli, leafNames, nvar)
    if (!is.null(out)) return(out)
  }
  NULL
}


# Reconstruct one wide entry with a general (multi-term) denominator. Along a ray from a
# generic shift s the entry is a univariate rational in t; a Cauchy fit normalised to
# B(0) = 1 evaluates N(u)/D(s) and D(u)/D(s) at any u, both sparse polynomials, recovered
# by Ben-Or-Tiwari on the geometric schedule and lifted across primes. Returns the entry
# string, or NULL when no bounded-degree rational form fits.
.symGeneralRationalEntry <- function(reli, supportCol, f, point0, leafNames,
                                        NtUsed, kcall, pivots, residueFn = NULL,
                                        ctrl = reconstControl(), zSlots = NULL) {
  nvar <- length(reli)
  np <- length(.symPrimes)
  bases <- .symSieve(nvar)
  s <- as.numeric(point0[reli])

  sampleR <- function(relvals, p) {
    pt <- point0; pt[reli] <- relvals %% p
    rp <- kcall(pt, p, NtUsed)
    if (!isTRUE(rp$ok) || !identical(as.integer(rp$pivots), as.integer(pivots)))
      return(NULL)
    zv <- if (is.null(zSlots)) NULL else pt[zSlots + 1L]
    nvp <- if (is.null(residueFn)) .symNullResidues(rp, f, p) else residueFn(rp, p, zv)
    if (is.null(nvp)) return(NULL)
    nvp[supportCol + 1L]
  }
  # value of (N(u)/D(s), D(u)/D(s)) at u over prime p, via the ray Cauchy fit
  evalND <- function(uvec, p, dN, dD) {
    v <- (uvec - s) %% p
    tn <- 0:(dN + dD + 4L)
    rv <- lapply(tn, function(tj) sampleR((s + tj * v) %% p, p))
    if (any(vapply(rv, is.null, logical(1)))) return(NULL)
    res <- symCauchyEval(as.integer(tn), as.integer(unlist(rv)), dN, dD, p)
    if (!identical(res$status, "ok")) return(NULL)
    c(res$N, res$D)
  }

  # numerator/denominator degrees, probed once at a generic point
  P1 <- .symPrimes[1]
  u0 <- vapply(seq_len(nvar), function(i) (bases[i] * bases[i] + 3) %% P1, numeric(1))
  deg <- NULL
  for (tot in 2:(ctrl$generalDegNum + ctrl$generalDegDen)) {
    if (.symExpired(ctrl)) return(NULL)
    for (a in seq.int(max(1L, tot - ctrl$generalDegDen), min(tot - 1L, ctrl$generalDegNum))) {
      bD <- tot - a
      if (bD < 1L || bD > ctrl$generalDegDen) next
      if (!is.null(evalND(u0, P1, a, bD))) { deg <- c(a, bD); break }
    }
    if (!is.null(deg)) break
  }
  if (is.null(deg)) return(NULL)
  dN <- deg[1]; dD <- deg[2]
  if (choose(nvar + dN, nvar) + choose(nvar + dD, nvar) > ctrl$laurentCandCap)
    return(NULL)
  candN <- matrix(as.integer(.symMonoTable(nvar, dN)), ncol = nvar)
  candD <- matrix(as.integer(.symMonoTable(nvar, dD)), ncol = nvar)

  # geometric Ben-Or-Tiwari sampling of both N/D(s) and D/D(s), grown until both
  # term counts stabilise
  maxLen <- 2L * ctrl$termCap + 2L
  Nseq <- lapply(seq_len(np), function(.) integer(0)); Dseq <- Nseq
  cur <- lapply(seq_len(np), function(.) rep(1, nvar))
  have <- 0L
  repeat {
    if (.symExpired(ctrl)) return(NULL)
    target <- min(if (have == 0L) 16L else 2L * have, maxLen)
    for (pj in seq_len(np)) {
      p <- .symPrimes[pj]
      for (k in (have + 1L):target) {
        if (.symExpired(ctrl)) return(NULL)
        e <- evalND(cur[[pj]], p, dN, dD)
        if (is.null(e)) return(NULL)
        Nseq[[pj]][k] <- e[1]; Dseq[[pj]][k] <- e[2]
        cur[[pj]] <- (cur[[pj]] * bases) %% p
      }
    }
    have <- target
    if (2L * symBMorder(Nseq[[1]], P1) < have &&
        2L * symBMorder(Dseq[[1]], P1) < have) break
    if (have >= maxLen) return(NULL)
  }

  fitPoly <- function(seqList, candM) {
    perPrime <- vector("list", np)
    for (pj in seq_len(np)) {
      mr <- symMonoResidues(candM, as.integer(bases), .symPrimes[pj])
      res <- symSparsePoly(seqList[[pj]], candM, as.integer(mr), .symPrimes[pj])
      if (!identical(res$status, "ok") || res$nterms == 0L) return(NULL)
      perPrime[[pj]] <- res
    }
    .symBotReconcile(perPrime)
  }
  Nrc <- fitPoly(Nseq, candN); if (is.null(Nrc)) return(NULL)
  Drc <- fitPoly(Dseq, candD); if (is.null(Drc)) return(NULL)
  vars <- leafNames[reli]
  numS <- .symPolyString(Nrc$num, Nrc$den, Nrc$exps, vars)
  denS <- .symPolyString(Drc$num, Drc$den, Drc$exps, vars)
  if (denS == "0") return(NULL)
  if (denS == "1") numS else paste0("(", numS, ")/(", denS, ")")
}


# Relevance of one direction from the shared probe: its base nullspace vector, the
# support columns to fit, and the leaves moving the direction and each entry, so each
# entry is fit over its own variables only.
.symDirectionRelevance <- function(f, ref, pivots, zSlots, nLeaves, point0,
                                     relProbe, residueFn) {
  P <- .symPrimes[1]
  zvals0 <- if (is.null(zSlots)) NULL else point0[zSlots + 1L]
  nv <- function(rp, p, zvals = zvals0) if (is.null(residueFn)) .symNullResidues(rp, f, p)
                        else residueFn(rp, p, zvals)
  # A pivot-shifted probe is unusable, not proof of relevance. The support-pinned and
  # free-column gauges skip it and rely on verification; canonical/log gauges count it.
  optimistic <- isTRUE(attr(residueFn, "pinnedSupport")) || is.null(residueFn)
  base_nv <- nv(ref, P)
  if (is.null(base_nv)) base_nv <- .symNullResidues(ref, f, P)
  supportCols <- setdiff(which(base_nv != 0) - 1L, f)
  relevant <- integer(0)
  relByEntry <- replicate(length(supportCols), integer(0), simplify = FALSE)
  for (li in seq_len(nLeaves)) {
    rp <- relProbe[[li]]$rp; zvp <- relProbe[[li]]$zvals
    nvp <- if (!is.null(rp) && isTRUE(rp$ok) &&
               identical(as.integer(rp$pivots), as.integer(pivots))) nv(rp, P, zvp)
           else NULL
    if (is.null(nvp)) {
      if (optimistic) next
      relevant <- c(relevant, li)
      relByEntry <- lapply(relByEntry, function(s) c(s, li))
      next
    }
    if (any(nvp != base_nv)) relevant <- c(relevant, li)
    for (i in which(nvp[supportCols + 1L] != base_nv[supportCols + 1L]))
      relByEntry[[i]] <- c(relByEntry[[i]], li)
  }
  list(base_nv = base_nv, supportCols = supportCols, relevant = relevant,
       relByEntry = relByEntry)
}


# Points (over the union of relevant leaves) needed for the dense fit of the widest
# entry, and whether any entry is dense at all.
.symDenseNeed <- function(relByEntry, ctrl) {
  isDense <- function(r) length(r) >= 1L && length(r) <= ctrl$relevanceCap
  denseRel <- max(0L, vapply(relByEntry,
                             function(r) if (isDense(r)) length(r) else 0L, integer(1)))
  anyDense <- any(vapply(relByEntry, isDense, logical(1)))
  maxNeed <- if (anyDense)
    2L * choose(denseRel + ctrl$degreeCap, denseRel) - 1L + ctrl$sampleSlack else 0L
  list(maxNeed = as.integer(maxNeed), anyDense = anyDense)
}


# Dense sample bank shared by the free-column directions: pivot-consistent points over
# the union of relevant leaves, the kernel result stored per point and prime. `ok` is
# FALSE when it cannot be filled.
.symBuildSharedBank <- function(union, needPts, point0, pool, poolNext, kbatch,
                                   pivots, NtUsed, nP, ctrl = NULL) {
  nu <- length(union)
  U <- matrix(0L, 0L, nu); points <- list(); rps <- list()
  tries <- 0L
  .bd <- nzchar(Sys.getenv("DMOD_SYM_BANKDIAG"))
  .nSolveFail <- 0L; .nPivMismatch <- 0L; .nGood <- 0L
  while (length(points) < needPts && tries < 30L * needPts) {
    # bail early when no candidate solves at every prime (coupled equilibrate case);
    # the caller falls back to per-prime reconstruction
    if (length(points) == 0L && tries >= min(2L * needPts, 48L)) {
      if (.bd) message(sprintf("[bankdiag/bail] union=%d needPts=%d tries=%d good=0 solveFail=%d pivMismatch=%d",
                               nu, needPts, tries, .nSolveFail, .nPivMismatch))
      break
    }
    if (!is.null(ctrl) && .symExpired(ctrl)) {
      if (.bd) message(sprintf("[bankdiag/expired] union=%d needPts=%d tries=%d good=%d solveFail=%d pivMismatch=%d",
                               nu, needPts, tries, .nGood, .nSolveFail, .nPivMismatch))
      return(list(ok = FALSE, poolNext = poolNext))
    }
    chunk <- max(1L, min(needPts - length(points) + 5L, 30L * needPts - tries))
    # small chunks until the first point lands, so the bail fires early
    if (length(points) == 0L) chunk <- min(chunk, 24L)
    # bounded chunks under a deadline, so the expiry check fires in time
    if (!is.null(ctrl) && !is.null(ctrl$deadline)) chunk <- min(chunk, 96L)
    tries <- tries + chunk
    uv <- lapply(seq_len(chunk), function(ci) {
      u <- pool(poolNext + seq_len(nu) - 1L); poolNext <<- poolNext + nu; u
    })
    cand <- lapply(uv, function(u) { pt <- point0; pt[union] <- u; pt })
    res <- kbatch(rep(cand, each = nP), rep(.symPrimes, times = chunk), NtUsed)
    for (ci in seq_len(chunk)) {
      if (length(points) >= needPts) break
      good <- TRUE; rpc <- vector("list", nP)
      for (j in seq_len(nP)) {
        rp <- res[[(ci - 1L) * nP + j]]
        if (!isTRUE(rp$ok)) { good <- FALSE; if (.bd) .nSolveFail <- .nSolveFail + 1L; break }
        if (!identical(as.integer(rp$pivots), as.integer(pivots))) {
          good <- FALSE; if (.bd) .nPivMismatch <- .nPivMismatch + 1L; break
        }
        rpc[[j]] <- rp
      }
      if (!good) next
      if (.bd) .nGood <- .nGood + 1L
      U <- rbind(U, uv[[ci]])
      points[[length(points) + 1L]] <- cand[[ci]]
      rps[[length(rps) + 1L]] <- rpc
    }
  }
  if (.bd) message(sprintf("[bankdiag] union=%d needPts=%d tries=%d good=%d solveFail=%d pivMismatch=%d",
                           nu, needPts, tries, .nGood, .nSolveFail, .nPivMismatch))
  if (length(points) < needPts) return(list(ok = FALSE, poolNext = poolNext))
  list(ok = TRUE, U = U, points = points, rps = rps, union = union,
       poolNext = poolNext)
}


# Each entry of a support-pinned direction as a Laurent monomial coeff * prod(leaf^exp):
# exponents from the base-to-probe ratios at one prime, the coefficient lifted across
# primes at the base point. NULL if an entry is not such a monomial.
.symPinnedMonomials <- function(f, pivots, znames, leafNames, point0,
                                  NtUsed, kcall, nv, base_nv, supportCols,
                                  relByEntry, relProbe, spy, maxExp = 6L) {
  P <- .symPrimes[1]
  powmod <- function(b, e, p) {
    r <- 1; b <- b %% p
    while (e > 0) { if (e %% 2 == 1) r <- .symMulmod(r, b, p); e <- e %/% 2
                    if (e > 0) b <- .symMulmod(b, b, p) }
    r
  }
  # smallest a in [-maxExp, maxExp] with base^a == target (mod p), or NULL
  findExp <- function(base, target, p) {
    target <- target %% p
    if (target == 1) return(0L)
    acc <- 1
    for (a in seq_len(maxExp)) { acc <- .symMulmod(acc, base, p)
      if (acc == target) return(a) }
    acc <- 1; bi <- .symInvmod(base, p)
    for (a in seq_len(maxExp)) { acc <- .symMulmod(acc, bi, p)
      if (acc == target) return(-a) }
    NULL
  }

  entries <- list()
  for (i in seq_along(supportCols)) {
    ci <- supportCols[i]; reli <- relByEntry[[i]]
    e0 <- base_nv[ci + 1L] %% P
    if (e0 == 0) return(NULL)
    exps <- integer(length(reli))
    for (j in seq_along(reli)) {
      l <- reli[j]; pr <- relProbe[[l]]
      if (is.null(pr$rp) || !isTRUE(pr$rp$ok) ||
          !identical(as.integer(pr$rp$pivots), as.integer(pivots))) return(NULL)
      ev <- nv(pr$rp, P, pr$zvals)
      if (is.null(ev)) return(NULL)
      bval <- as.numeric(point0[l]) %% P; pval <- as.numeric(pr$pertval) %% P
      if (bval == 0 || pval == 0) return(NULL)
      a <- findExp(.symMulmod(pval, .symInvmod(bval, P), P),
                   .symMulmod(ev[ci + 1L] %% P, .symInvmod(e0, P), P), P)
      if (is.null(a)) return(NULL)
      exps[j] <- a
    }
    # coefficient entry(base) / prod(leaf^exp), lifted over the pivot-consistent primes
    monoDenom <- function(pj) {
      d <- 1
      for (j in seq_along(reli)) {
        bv <- as.numeric(point0[reli[j]]) %% pj
        d <- .symMulmod(d, if (exps[j] >= 0) powmod(bv, exps[j], pj)
                            else powmod(.symInvmod(bv, pj), -exps[j], pj), pj)
      }
      d
    }
    residP <- .symMulmod(e0, .symInvmod(monoDenom(P), P), P)
    usePrimes <- P; resid <- residP
    for (pj in .symPrimes[-1]) {
      rp <- kcall(point0, pj, NtUsed)
      if (is.null(rp) || !isTRUE(rp$ok) ||
          !identical(as.integer(rp$pivots), as.integer(pivots))) next
      nvp <- nv(rp, pj)
      if (is.null(nvp)) next
      usePrimes <- c(usePrimes, pj)
      resid <- c(resid, .symMulmod(nvp[ci + 1L] %% pj, .symInvmod(monoDenom(pj), pj), pj))
    }
    rec <- symRatRecon(matrix(as.integer(resid), 1L), as.integer(usePrimes))
    if (rec$den[1] == "0") return(NULL)
    coef <- if (rec$den[1] == "1") rec$num[1] else paste0(rec$num[1], "/(", rec$den[1], ")")
    # assemble coeff * prod(leaf^exp) as a Laurent monomial string
    num <- character(0); den <- character(0)
    for (j in seq_along(reli)) {
      v <- leafNames[reli[j]]; a <- exps[j]
      if (a == 0) next
      e <- abs(a)
      term <- if (e == 1) v else paste0(v, "^", e)
      if (a > 0) num <- c(num, term) else den <- c(den, term)
    }
    numS <- paste(c(if (coef != "1" || !length(num)) coef, num), collapse = "*")
    expr <- if (!length(den)) numS
            else paste0("(", numS, ")/(", paste(den, collapse = "*"), ")")
    entries[[znames[ci + 1L]]] <- if (is.null(spy)) expr else .symSimplify(expr, spy)
  }
  entries[[znames[f + 1L]]] <- "1"
  list(support = .symSort(names(entries)), vector = entries, type = "general",
       closedForm = TRUE)
}


# Per-prime reconstruction for the coupled steady-state path, where a point rarely solves
# at every prime: each prime fits its own points and the coefficients are lifted by CRT
# (the function is the same, only the points differ). Primes that never solve are
# dropped. Same return as .symInterpolateDirection; anchor column `f` is 1.
.symInterpolatePerprime <- function(f, ref, pivots, znames, zSlots, leafNames,
                                      nLeaves, point0, pool, poolNext, NtUsed,
                                      kcall, kbatch, spy, relProbe, ctrl,
                                      auxLeaves = integer(0), physCols = NULL) {
  nP <- length(.symPrimes)
  zvalsOf <- function(pt) if (is.null(zSlots)) NULL else pt[zSlots + 1L]
  nv <- function(rp, p) .symNullResidues(rp, f, p)
  rel <- .symDirectionRelevance(f, ref, pivots, zSlots, nLeaves, point0,
                                  relProbe, NULL)
  base_nv <- rel$base_nv; supportCols <- rel$supportCols
  relByEntry <- rel$relByEntry
  # joint mode fits only the physical columns; the auxiliary ones are projected out later
  if (!is.null(physCols)) {
    keep <- which(supportCols %in% physCols)
    supportCols <- supportCols[keep]; relByEntry <- relByEntry[keep]
  }
  relevant <- sort(unique(unlist(relByEntry)))
  fallback <- function(reason) list(poolNext = poolNext,
    entry = list(support = .symSort(znames[base_nv != 0]), type = "general",
                 closedForm = FALSE, reason = reason))
  timedOut <- "reconstruction time budget (reconstControl(timeout=)) exceeded"
  if (.symExpired(ctrl)) return(fallback(timedOut))
  ex <- function(v) if (length(auxLeaves)) setdiff(v, auxLeaves) else v
  maxRelPhys <- max(0L, vapply(relByEntry, function(e) length(ex(e)), integer(1)))
  if (nzchar(Sys.getenv("DMOD_SYM_SPARSEDIAG")))
    message(sprintf("perprime f=%d: support=%d relevant=%d maxRelPhys=%d",
                    f, length(supportCols), length(relevant), maxRelPhys))
  if (length(ex(relevant)) > ctrl$relevanceCapDir || maxRelPhys > ctrl$relevanceCapSparse)
    return(fallback(sprintf("direction couples %d parameters (an entry up to %d)",
                            length(ex(relevant)), maxRelPhys)))

  # a constant entry lifts from the base point over the primes where it solves
  constEntry <- function(col) {
    vals <- integer(0); prs <- numeric(0)
    for (pj in .symPrimes) {
      rp <- kcall(point0, pj, NtUsed)
      if (!isTRUE(rp$ok) || !identical(as.integer(rp$pivots), as.integer(pivots))) next
      v <- nv(rp, pj); if (is.null(v)) next
      vals <- c(vals, v[col + 1L]); prs <- c(prs, pj)
    }
    if (!length(prs)) return(NULL)
    rec <- symRatRecon(matrix(as.integer(vals), 1L), as.integer(prs))
    if (rec$den[1] == "0") return(NULL)
    if (rec$den[1] == "1") rec$num[1] else paste0(rec$num[1], "/", rec$den[1])
  }
  if (!length(relevant)) {
    entries <- list()
    for (i in seq_along(supportCols)) {
      e <- constEntry(supportCols[i])
      if (is.null(e)) return(fallback("a constant entry could not be lifted"))
      entries[[znames[supportCols[i] + 1L]]] <- e
    }
    entries[[znames[f + 1L]]] <- "1"
    return(list(poolNext = poolNext, entry = list(support = .symSort(names(entries)),
      vector = entries, type = "general", closedForm = TRUE,
      relevantLeaves = integer(0))))
  }

  # per-prime point banks over the relevant leaves, grown on demand
  nrel <- length(relevant)
  bankU <- lapply(seq_len(nP), function(.) matrix(0L, 0L, nrel))
  bankR <- lapply(seq_len(nP), function(.) matrix(0L, 0L, length(supportCols)))
  triesP <- integer(nP); deadP <- logical(nP)
  ensurePrime <- function(j, need) {
    p <- .symPrimes[j]
    while (nrow(bankU[[j]]) < need && !deadP[j]) {
      if (.symExpired(ctrl)) break
      have <- nrow(bankU[[j]])
      chunk <- max(1L, min(need - have + 4L, 128L))
      uv <- lapply(seq_len(chunk), function(.) {
        u <- pool(poolNext + seq_len(nrel) - 1L); poolNext <<- poolNext + nrel; u })
      cand <- lapply(uv, function(u) { pt <- point0; pt[relevant] <- u; pt })
      res <- kbatch(cand, rep(p, chunk), NtUsed)
      for (ci in seq_len(chunk)) {
        rp <- res[[ci]]
        if (!isTRUE(rp$ok) || !identical(as.integer(rp$pivots), as.integer(pivots))) next
        vv <- nv(rp, p); if (is.null(vv)) next
        bankU[[j]] <<- rbind(bankU[[j]], uv[[ci]])
        bankR[[j]] <<- rbind(bankR[[j]], vv[supportCols + 1L])
      }
      triesP[j] <<- triesP[j] + chunk
      # drop a prime that never yields a usable point
      if (nrow(bankU[[j]]) == 0L && triesP[j] >= 8L * need + 32L) { deadP[j] <<- TRUE; break }
      if (triesP[j] >= 80L * need + 256L) break
    }
    nrow(bankU[[j]]) >= need
  }

  entries <- list()
  for (i in seq_along(supportCols)) {
    reli <- relByEntry[[i]]
    if (!length(reli)) {
      e <- constEntry(supportCols[i])
      if (is.null(e)) return(fallback("a constant entry could not be lifted"))
      entries[[znames[supportCols[i] + 1L]]] <- e
      next
    }
    cols_i <- match(reli, relevant)
    fitted <- NULL
    for (degree in 0:ctrl$degreeCap) {
      if (.symExpired(ctrl)) return(fallback(timedOut))
      mons <- .symMonoTable(length(reli), degree)
      need <- 2L * nrow(mons) - 1L + ctrl$sampleSlack
      if (need > ctrl$perprimeCap) break
      live <- which(vapply(seq_len(nP), function(j) !deadP[j] && ensurePrime(j, need),
                           logical(1)))
      if (length(live) < ctrl$perprimeMinPrimes) next
      # fit per prime, then normalise every prime on the first fit's free column so
      # the CRT combines one representative
      raw <- vector("list", length(live)); refFree <- NULL; okAll <- TRUE
      for (jj in seq_along(live)) { j <- live[jj]
        su <- matrix(as.integer(bankU[[j]][seq_len(need), cols_i, drop = FALSE]), need)
        rv <- as.integer(bankR[[j]][seq_len(need), i])
        fit <- symFitRational(su, matrix(as.integer(mons), nrow(mons)), rv, .symPrimes[j])
        if (!identical(fit$status, "ok")) { okAll <- FALSE; break }
        raw[[jj]] <- as.numeric(fit$coeffs)
        if (is.null(refFree)) refFree <- fit$freeCol
      }
      if (!okAll) next
      coefRes <- matrix(0L, 2L * nrow(mons), length(live))
      for (jj in seq_along(live)) { j <- .symPrimes[live[jj]]
        d <- raw[[jj]][refFree + 1L] %% j
        if (d == 0) { okAll <- FALSE; break }              # refFree not in support here
        coefRes[, jj] <- as.integer(.symMulmod(raw[[jj]], .symInvmod(d, j), j))
      }
      if (!okAll) next
      rec <- symRatRecon(coefRes, as.integer(.symPrimes[live]))
      if (any(rec$den == "0")) next
      nMon <- nrow(mons); numI <- seq_len(nMon); denI <- nMon + seq_len(nMon)
      fitted <- list(numN = rec$num[numI], numD = rec$den[numI],
                     denN = rec$num[denI], denD = rec$den[denI], mons = mons)
      if (nzchar(Sys.getenv("DMOD_SYM_SPARSEDIAG")))
        message(sprintf("  perprime entry %d/%d col=%s reli=%d: fit at degree %d (%d primes, %d pts/prime)",
                        i, length(supportCols), znames[supportCols[i] + 1L],
                        length(reli), degree, length(live), need))
      break
    }
    if (is.null(fitted)) {
      if (nzchar(Sys.getenv("DMOD_SYM_SPARSEDIAG")))
        message(sprintf("  perprime entry %d/%d col=%s reli=%d: FAILED (degreeCap/perprimeCap)",
                        i, length(supportCols), znames[supportCols[i] + 1L], length(reli)))
      return(fallback(paste("a joint entry could not be fit per-prime within",
                            "degreeCap/perprimeCap; raise them or the direction is",
                            "not a bounded-degree rational")))
    }
    vars <- leafNames[reli]
    numStr <- .symPolyString(fitted$numN, fitted$numD, fitted$mons, vars)
    denStr <- .symPolyString(fitted$denN, fitted$denD, fitted$mons, vars)
    expr <- if (denStr == "1") numStr else paste0("(", numStr, ")/(", denStr, ")")
    if (!is.null(spy)) expr <- .symSimplify(expr, spy)
    entries[[znames[supportCols[i] + 1L]]] <- expr
  }
  entries[[znames[f + 1L]]] <- "1"
  list(poolNext = poolNext, entry = list(support = .symSort(names(entries)),
    vector = entries, type = "general", closedForm = TRUE,
    relevantLeaves = relevant))
}


.symInterpolateDirection <- function(f, ref, pivots, znames, zSlots, leafNames,
                                       nLeaves, point0, pool, poolNext, NtUsed,
                                       kcall, spy, relProbe = NULL,
                                       residueFn = NULL, ctrl = reconstControl(),
                                       kbatch = NULL, sharedBank = NULL,
                                       fastOnly = FALSE, auxLeaves = integer(0)) {
  if (is.null(kbatch))
    kbatch <- function(pl, pv, Nt) Map(function(pp, pr) kcall(pp, pr, Nt), pl, pv)
  P <- .symPrimes[1]
  zvalsOf <- function(pt) if (is.null(zSlots)) NULL else pt[zSlots + 1L]
  zvals0 <- zvalsOf(point0)
  # `nv`: the direction's nullspace vector at a kernel result (free-column residue, or
  # residueFn's row, which may need the z-values); `f` is the gauge column fixed to 1
  nv <- function(rp, p, zvals = zvals0) if (is.null(residueFn)) .symNullResidues(rp, f, p)
                        else residueFn(rp, p, zvals)
  rel <- .symDirectionRelevance(f, ref, pivots, zSlots, nLeaves, point0,
                                  relProbe, residueFn)
  base_nv <- rel$base_nv; supportCols <- rel$supportCols
  relevant <- rel$relevant; relByEntry <- rel$relByEntry
  fallback <- function(reason = NULL) {
    list(poolNext = poolNext,
         entry = list(support = .symSort(znames[base_nv != 0]), type = "general",
                      closedForm = FALSE, reason = reason))
  }
  timedOut <- "reconstruction time budget (reconstControl(timeout=)) exceeded"
  if (.symExpired(ctrl)) return(fallback(timedOut))

  # a constant entry comes from the base-point kernel, evaluated once over all primes
  constBase <- NULL
  constEntry <- function(col) {
    if (is.null(constBase)) {
      kr <- kbatch(rep(list(point0), length(.symPrimes)), .symPrimes, NtUsed)
      constBase <<- lapply(seq_along(.symPrimes), function(k) {
        rp <- kr[[k]]
        if (!isTRUE(rp$ok) || !identical(as.integer(rp$pivots), as.integer(pivots)))
          return(NULL)
        nv(rp, .symPrimes[k])
      })
    }
    res <- vapply(seq_along(.symPrimes), function(k) {
      nvp <- constBase[[k]]
      if (is.null(nvp)) NA_real_ else nvp[col + 1L]
    }, numeric(1))
    if (anyNA(res)) return(NULL)
    rec <- symRatRecon(matrix(as.integer(res), 1L), as.integer(.symPrimes))
    if (rec$den[1] == "0") return(NULL)
    if (rec$den[1] == "1") rec$num[1] else paste0(rec$num[1], "/", rec$den[1])
  }

  if (!length(relevant)) {
    entries <- list()
    for (i in seq_along(supportCols)) {
      e <- constEntry(supportCols[i])
      if (is.null(e)) return(fallback("a constant entry could not be lifted from its residues"))
      entries[[znames[supportCols[i] + 1L]]] <- e
    }
    entries[[znames[f + 1L]]] <- "1"
    return(list(poolNext = poolNext,
                entry = list(support = .symSort(names(entries)), vector = entries,
                             type = "general", closedForm = TRUE)))
  }

  # pinned support: read Laurent-monomial entries (eta = -nhill) off the base point and
  # probes without resampling, which rarely keeps the pivots for a recast exponent
  if (isTRUE(attr(residueFn, "pinnedSupport")) && !is.null(relProbe)) {
    mono <- .symPinnedMonomials(f, pivots, znames, leafNames, point0,
                                  NtUsed, kcall, nv, base_nv, supportCols,
                                  relByEntry, relProbe, spy)
    if (!is.null(mono)) return(list(poolNext = poolNext, entry = mono))
    # fastOnly: no kernel sampling
    if (isTRUE(fastOnly))
      return(fallback("pinned entry is not a bounded-degree monomial"))
  }

  maxRel <- max(0L, vapply(relByEntry, length, integer(1)))
  # the relevance caps count parameters only, not the auxiliary leaves
  ex <- function(v) if (length(auxLeaves)) setdiff(v, auxLeaves) else v
  relPhys <- ex(relevant)
  maxRelPhys <- max(0L, vapply(relByEntry, function(e) length(ex(e)), integer(1)))
  if (nzchar(Sys.getenv("DMOD_SYM_SPARSEDIAG")))
    message(sprintf("interp f=%d: relevant=%d (phys %d) maxRel=%d (phys %d) supportCols=%d (capDir=%d capSparse=%d)",
                    f, length(relevant), length(relPhys), maxRel, maxRelPhys,
                    length(supportCols), ctrl$relevanceCapDir, ctrl$relevanceCapSparse))
  if (length(relPhys) > ctrl$relevanceCapDir || maxRelPhys > ctrl$relevanceCapSparse)
    return(fallback(sprintf(
      paste("direction couples %d parameters (an entry up to %d), above",
            "relevanceCapDir=%d / relevanceCapSparse=%d; raise the relevant cap"),
      length(relPhys), maxRelPhys, ctrl$relevanceCapDir, ctrl$relevanceCapSparse)))
  nrel <- length(relevant)
  nP <- length(.symPrimes)
  dn <- .symDenseNeed(relByEntry, ctrl)
  sampleU <- matrix(0L, 0L, nrel)
  rstore <- lapply(seq_along(supportCols), function(i) matrix(0L, 0L, nP))
  if (dn$anyDense) {
    maxNeed <- dn$maxNeed
    if (!is.null(sharedBank)) {
      # project the shared bank onto this direction's leaves and residues
      nAcc <- nrow(sharedBank$U)
      if (nAcc < maxNeed)
        return(fallback("the shared bank holds fewer points than this direction needs"))
      sampleU <- sharedBank$U[, match(relevant, sharedBank$union), drop = FALSE]
      rstore <- lapply(seq_along(supportCols), function(i) matrix(0L, nAcc, nP))
      for (a in seq_len(nAcc)) {
        zvc <- zvalsOf(sharedBank$points[[a]])
        for (j in seq_len(nP)) {
          nvp <- nv(sharedBank$rps[[a]][[j]], .symPrimes[j], zvc)
          if (is.null(nvp)) return(fallback("shared-bank residue extraction failed"))
          for (i in seq_along(supportCols)) rstore[[i]][a, j] <- nvp[supportCols[i] + 1L]
        }
      }
    } else {
      # sample in chunks over all primes at once, keeping pivot-consistent candidates
      tries <- 0L
      while (nrow(sampleU) < maxNeed && tries < 30L * maxNeed) {
        if (.symExpired(ctrl)) return(fallback(timedOut))
        # bail early when perturbations admit no valid steady state (e.g. a held recast
        # E inconsistent with the perturbed parameters)
        if (nrow(sampleU) == 0L && tries >= 2L * maxNeed)
          return(fallback("no valid steady state under perturbation of the relevant leaves"))
        chunk <- max(1L, min(maxNeed - nrow(sampleU) + ctrl$sampleSlack,
                             30L * maxNeed - tries))
        if (!is.null(ctrl$deadline)) chunk <- min(chunk, 96L)
        tries <- tries + chunk
        uv <- lapply(seq_len(chunk), function(ci) {
          u <- pool(poolNext + seq_len(nrel) - 1L); poolNext <<- poolNext + nrel
          u
        })
        cand <- lapply(uv, function(u) { pt <- point0; pt[relevant] <- u; pt })
        res <- kbatch(rep(cand, each = nP), rep(.symPrimes, times = chunk), NtUsed)
        for (ci in seq_len(chunk)) {
          if (nrow(sampleU) >= maxNeed) break
          good <- TRUE
          rrow <- matrix(0L, length(supportCols), nP)
          zvc <- zvalsOf(cand[[ci]])
          for (j in seq_len(nP)) {
            rp <- res[[(ci - 1L) * nP + j]]
            if (!isTRUE(rp$ok) || !identical(as.integer(rp$pivots), as.integer(pivots))) {
              good <- FALSE; break
            }
            nvp <- nv(rp, .symPrimes[j], zvc)
            if (is.null(nvp)) { good <- FALSE; break }
            rrow[, j] <- nvp[supportCols + 1L]
          }
          if (!good) next
          sampleU <- rbind(sampleU, uv[[ci]])
          for (i in seq_along(supportCols)) rstore[[i]] <- rbind(rstore[[i]], rrow[i, ])
        }
      }
      if (nrow(sampleU) < maxNeed)
        return(fallback(paste("dense sampling could not collect enough",
                              "pivot-consistent points; raise sampleSlack or probeRetries")))
    }
  }

  # reconstruct each entry over its own relevant variables: a constant from the
  # base point, a narrow entry by the dense fit, a wide one by sparse Laurent
  entries <- list()
  for (i in seq_along(supportCols)) {
    reli <- relByEntry[[i]]
    if (!length(reli)) {
      e <- constEntry(supportCols[i])
      if (is.null(e)) return(fallback("a constant entry could not be lifted from its residues"))
      entries[[znames[supportCols[i] + 1L]]] <- e
      next
    }
    if (length(reli) > ctrl$relevanceCap) {
      # single-monomial denominator (Laurent) first, then a general denominator
      e <- .symSparseEntry(reli, supportCols[i], f, point0, leafNames, NtUsed,
                             kcall, pivots, residueFn, ctrl, zSlots)
      if (is.null(e))
        e <- .symGeneralRationalEntry(reli, supportCols[i], f, point0,
                                         leafNames, NtUsed, kcall, pivots,
                                         residueFn, ctrl, zSlots)
      if (is.null(e) && .symExpired(ctrl)) return(fallback(timedOut))
      if (is.null(e)) return(fallback(sprintf(
        paste("an entry couples %d variables and the sparse fit hit its caps",
              "(laurentDegNum=%d, generalDegNum=%d, generalDegDen=%d, termCap=%d);",
              "raise these, or raise relevanceCap=%d to use the dense fit"),
        length(reli), ctrl$laurentDegNum, ctrl$generalDegNum, ctrl$generalDegDen,
        ctrl$termCap, ctrl$relevanceCap)))
      if (!is.null(spy)) e <- .symSimplify(e, spy)
      entries[[znames[supportCols[i] + 1L]]] <- e
      next
    }
    cols_i <- match(reli, relevant)
    vars_i <- leafNames[reli]
    fitted <- NULL
    for (degree in 0:ctrl$degreeCap) {
      mons <- .symMonoTable(length(reli), degree)
      need <- 2L * nrow(mons) - 1L + ctrl$sampleSlack
      rec <- .symReconstructEntry(
        matrix(as.integer(sampleU[seq_len(need), cols_i, drop = FALSE]), need),
        matrix(as.integer(mons), nrow(mons)),
        rstore[[i]][seq_len(need), , drop = FALSE], .symPrimes)
      if (!is.null(rec)) { fitted <- list(rec = rec, mons = mons); break }
    }
    if (is.null(fitted))
      return(fallback(sprintf("an entry exceeded degreeCap=%d; raise it",
                              ctrl$degreeCap)))
    numStr <- .symPolyString(fitted$rec$numCoefN, fitted$rec$numCoefD,
                               fitted$mons, vars_i)
    denStr <- .symPolyString(fitted$rec$denCoefN, fitted$rec$denCoefD,
                               fitted$mons, vars_i)
    expr <- if (denStr == "1") numStr else paste0("(", numStr, ")/(", denStr, ")")
    if (!is.null(spy)) expr <- .symSimplify(expr, spy)
    entries[[znames[supportCols[i] + 1L]]] <- expr
  }

  entries[[znames[f + 1L]]] <- "1"
  list(poolNext = poolNext,
       entry = list(support = .symSort(names(entries)), vector = entries,
                    type = "general", closedForm = TRUE))
}


# Condition grid, events and initial values as the per-segment inputs of
# compileObservabilityTapeMulti, one chain per condition split at the event times after
# t0 (the earliest event time, else 0). Returns flat per-segment lists plus
# chainOf/posInChain; later segments get a dummy ic0, their state is propagated.
.symResolveConditions <- function(conditions, events, initial, symbols, states,
                                    constStates = character(0),
                                    forcings = character(0), t0 = NULL,
                                    equilibrate = FALSE,
                                    condSubs = NULL, condInitial = NULL,
                                    nCondObs = 0L) {
  # per-condition trafo or observation lists set the condition count without a grid
  nGrid <- if (is.null(conditions)) 0L else nrow(as.data.frame(conditions))
  Kcond <- max(1L, length(condSubs), length(condInitial), nGrid, nCondObs)
  if (nGrid && length(condSubs) && length(condSubs) != nGrid)
    stop("symmetryDetection(): the per-condition `trafo` list length (",
         length(condSubs), ") must match the condition grid rows (", nGrid, ").",
         call. = FALSE)
  if (is.null(conditions)) conditions <- data.frame(row.names = as.character(seq_len(Kcond)))
  else conditions <- as.data.frame(conditions, stringsAsFactors = FALSE)
  conds <- rownames(conditions)
  if (is.null(conds) || !length(conds))
    conds <- as.character(seq_len(max(1L, nrow(conditions))))
  cols <- colnames(conditions)
  subCols <- intersect(cols, symbols)
  initial <- if (is.null(initial)) NULL else as.eqnvec(initial)
  dynStates <- setdiff(states, constStates)
  # a grid column naming a dynamic state sets its initial value, as elsewhere in dMod;
  # a constant state (rhs = 0) is substituted
  icCols  <- intersect(subCols, dynStates)
  subCols <- setdiff(subCols, icCols)

  ev <- if (is.null(events)) NULL else as.data.frame(events, stringsAsFactors = FALSE)
  tstr <- if (is.null(ev) || !nrow(ev)) character(0) else trimws(as.character(ev$time))
  tnum <- suppressWarnings(as.numeric(tstr))
  if (is.null(t0)) t0 <- if (length(tnum) && any(!is.na(tnum)))
    min(tnum, na.rm = TRUE) else 0
  t0 <- as.numeric(t0)
  evVar    <- if (is.null(ev) || !nrow(ev)) character(0) else as.character(ev$var)
  evMethod <- if (is.null(ev) || !nrow(ev)) character(0) else as.character(ev$method)
  isSwitch <- evVar %in% constStates

  cell <- function(ci, col) as.character(conditions[ci, col])
  resolve <- function(val, ci) {
    val <- as.character(val)
    v <- if (val %in% cols) cell(ci, val) else val
    # this condition's trafo substitutions apply to event values too
    if (length(condSubs) && ci <= length(condSubs) && length(condSubs[[ci]])) {
      sub <- condSubs[[ci]]
      v <- replaceSymbols(names(sub), unlist(sub), v)
    }
    v
  }

  # Segment of each event: 1 for S0 (a numeric time <= t0), else the rank of its time
  # among the distinct later times. A time given in parameters takes its place in the
  # order of the eventlist, which must then list the numeric times ascending.
  isT0 <- !is.na(tnum) & tnum <= t0
  key <- ifelse(is.na(tnum), tstr, as.character(tnum))
  postKeys <- if (!any(is.na(tnum))) as.character(sort(unique(tnum[!isT0])))
              else unique(key[!isT0])
  if (any(is.na(tnum))) {
    kn <- suppressWarnings(as.numeric(postKeys))
    if (is.unsorted(kn[!is.na(kn)], strictly = TRUE))
      stop("symmetryDetection(): with an event time given in parameters, the ",
           "eventlist sets the order of the events and must list the numeric times ",
           "ascending.", call. = FALSE)
  }
  segOf <- ifelse(isT0, 1L, 1L + match(key, postKeys))
  segTimes <- c(as.character(t0), postKeys)   # left endpoints of S0, S1, ...
  nSeg <- length(segTimes)
  # dummy seed for a propagation-seeded later segment: not a free unknown
  dummyIc0 <- setNames(as.list(rep("0", length(dynStates))), dynStates)

  subsList <- list(); ic0List <- list(); ev0List <- list(); segEvList <- list()
  timeList <- list()
  equilList <- logical(0); chainOf <- integer(0); posInChain <- integer(0)
  for (k in seq_along(conds)) {
    subsBase <- list()
    for (col in subCols) subsBase[[col]] <- cell(k, col)
    # a per-condition trafo overrides the grid substitutions
    if (length(condSubs) && !is.null(condSubs[[k]]))
      for (nm in names(condSubs[[k]])) subsBase[[nm]] <- condSubs[[k]][[nm]]
    # initial values from the trafo list, else the shared `initial`
    initK <- if (length(condInitial)) condInitial[[k]] else initial
    initK <- if (is.null(initK)) NULL else as.eqnvec(initK)
    for (j in seq_len(nSeg)) {
      # latest switch value per constant state; a forcing is 0 before its first event
      subs <- subsBase
      for (cs in constStates) {
        si <- which(isSwitch & evVar == cs & evMethod == "replace" & segOf <= j)
        if (length(si)) subs[[cs]] <- resolve(ev$value[si[which.max(segOf[si])]], k)
        else if (cs %in% forcings) subs[[cs]] <- "0"
      }
      if (j == 1L) {
        # S0: t0 events compose the initial value (with equilibrate, a dose on the
        # resting state)
        leIdx <- which(!isSwitch & isT0)
        ev0 <- lapply(leIdx, function(i) list(
          var = evVar[i], method = evMethod[i], value = resolve(ev$value[i], k)))
        ic0 <- list()
        for (X in dynStates) {
          isForcing <- X %in% forcings
          hasInit <- !is.null(initK) && X %in% names(initK)
          hasCell <- X %in% icCols
          ops <- leIdx[evVar[leIdx] == X]
          if (!isForcing && !hasInit && !hasCell && !length(ops)) next
          # a per-condition trafo wins over a grid cell, as it does for parameters
          cur <- if (isForcing) "0" else if (hasInit) as.character(initK[[X]])
                 else if (hasCell) resolve(cell(k, X), k) else X
          for (i in ops) cur <- .symComposeEvent(cur, evMethod[i], resolve(ev$value[i], k))
          ic0[[X]] <- cur
        }
      } else {
        # later segment: propagated state with dummy ic0; the kernel applies the doses
        ev0 <- list(); ic0 <- dummyIc0
        evIdx <- which(!isSwitch & segOf == j)
        segEv <- lapply(evIdx, function(i) list(
          var = evVar[i], method = evMethod[i], value = resolve(ev$value[i], k)))
      }
      subsList[[length(subsList) + 1L]] <- subs
      ic0List[[length(ic0List) + 1L]] <- ic0
      ev0List[[length(ev0List) + 1L]] <- ev0
      segEvList[[length(segEvList) + 1L]] <- if (j == 1L) list() else segEv
      timeList[[length(timeList) + 1L]] <- resolve(segTimes[j], k)
      equilList <- c(equilList, isTRUE(equilibrate))
      chainOf <- c(chainOf, k); posInChain <- c(posInChain, j)
    }
  }
  list(subs = subsList, ic0 = ic0List, segEquil = equilList, events0 = ev0List,
       segEvents = segEvList, chainOf = chainOf, posInChain = posInChain,
       conditions = conds, nConditions = length(conds), nGaps = nSeg - 1L,
       times = timeList)
}

# A log parameter theta is analysed in X = base^theta and reported in theta: X ->
# base^theta, xi_theta = xi_X / (X log(base)). A scaling of X becomes a translation.
.symLogParamBack <- function(res, lp, sd) {
  if (is.null(res) || !length(lp)) return(res)
  X  <- vapply(lp, function(e) as.character(e$X), "")
  th <- vapply(lp, function(e) as.character(e$theta), "")
  b  <- vapply(lp, function(e) as.character(e$base), "")
  ren <- function(v) { i <- match(v, X); v[!is.na(i)] <- th[i[!is.na(i)]]; v }
  res$coordinates <- ren(res$coordinates)
  res$nonIdentifiable <- lapply(res$nonIdentifiable, function(d) {
    d$support <- .symSort(ren(d$support))
    vec <- d$vector
    if (is.null(vec)) return(d)
    if (isTRUE(d$type == "scaling")) {
      if (!any(names(vec) %in% X)) return(d)
      vec <- setNames(lapply(names(vec), function(k) .symScalingComponent(vec[[k]], k)),
                      names(vec))
      d$type <- "general"
    }
    d$vector <- setNames(lapply(names(vec), function(k)
      as.character(sd$logParamBacksub(as.character(vec[[k]]), k, X, th, b))),
      ren(names(vec)))
    d
  })
  res
}


# compose one event onto a current initial-value expression string
.symComposeEvent <- function(cur, method, value) {
  switch(as.character(method),
    replace  = value,
    add      = paste0("(", cur, ") + (", value, ")"),
    multiply = paste0("(", cur, ") * (", value, ")"),
    stop("unsupported event method '", method, "'"))
}


# the operator d/d<var> with U+2202, escaped to keep the source ASCII for R CMD check
.symPartial <- function(var) paste0("\u2202/\u2202", var)

# Label of the component along a coordinate, eta(P): Unicode has no subscripts for
# capitals or '_', and eta_k_d would not show where the name ends.
.symEta <- function(var) paste0("\u03b7(", var, ")")

# ---- result display: the shared print()/summary() renderer ---------------------------

# the generator component xi_i of a scaling from its weight w and coordinate:
# w = 1 -> "var", w = -1 -> "-var", else "w*var"
.symScalingComponent <- function(w, var) {
  w <- as.character(w)
  if (w == "1") var
  else if (w == "-1") paste0("-", var)
  else paste0(w, "*", var)
}

# one signed generator term "<coef> d/d<var>": pull a leading '-' out for the
# join, elide a unit coefficient, and parenthesise a multi-term (sum) coefficient.
.symSignedTerm <- function(xi, var) {
  xi <- trimws(xi)
  neg <- startsWith(xi, "-")
  mag <- if (neg) trimws(sub("^-", "", xi)) else xi
  bare <- gsub("\\*\\*", "^", mag)                    # exponents are not sums
  needParen <- grepl("[+]", bare) || grepl(".[-]", bare)
  coef <- if (mag == "1") ""
          else paste0(if (needParen) paste0("(", mag, ")") else mag, " ")
  list(neg = neg, text = paste0(coef, .symPartial(var)))
}

# join signed terms into "a - b + c"
.symJoinGenerator <- function(signed) {
  if (!length(signed)) return("0")
  out <- paste0(if (signed[[1]]$neg) "-" else "", signed[[1]]$text)
  for (k in seq_along(signed)[-1])
    out <- paste0(out, if (signed[[k]]$neg) " - " else " + ", signed[[k]]$text)
  out
}

# the symmetry-class tag, carrying the total degree of a general direction that has one
.symClassTag <- function(d) {
  if (isTRUE(d$type == "general") && !is.null(d$degree) && isTRUE(d$degree >= 0L))
    sprintf("general (%d)", d$degree)
  else d$type
}

# one line per generator: sum_i xi_i d/dz_i and its class tag; $generator holds the
# components xi_i for every class
.symDirectionLine <- function(d) {
  if (is.null(d$generator))
    return(paste0("[", d$type, ", support only] involves: ",
                  paste(d$support, collapse = ", ")))
  signed <- lapply(names(d$generator),
                   function(k) .symSignedTerm(as.character(d$generator[[k]]), k))
  tag <- .symClassTag(d)
  if (isTRUE(d$certified))  tag <- paste0(tag, ", certified")
  if (isFALSE(d$verified))  tag <- paste0(tag, ", unverified")
  paste0(.symJoinGenerator(signed), "   [", tag, "]")
}

# plural group heading of a class; other labels (the polynomial engine's
# "scaling, translation", ...) are capitalised
.symTypeLabel <- function(t) switch(t,
  scaling = "Scalings", general = "General",
  paste0(toupper(substring(t, 1, 1)), substring(t, 2)))

# a non-negative integer as subscript digits U+2080..U+2089 (escaped; one column each)
.symSubscript <- function(n) {
  sub <- c("\u2080", "\u2081", "\u2082", "\u2083", "\u2084",
           "\u2085", "\u2086", "\u2087", "\u2088", "\u2089")
  paste(sub[as.integer(strsplit(as.character(n), "")[[1]]) + 1L], collapse = "")
}

# left-justify by display columns; sprintf pads by bytes, which multibyte glyphs break
.symLjust <- function(x, w) paste0(x, strrep(" ", pmax(0L, w - nchar(x))))


# ---- one generator as a component table ----------------------------------------

# Split a component at its top-level binary '+'/'-' for wrapping; a sign after '*',
# '(' or '**', or inside "1e-5", is not a cut.
.symSplitTerms <- function(s) {
  if (!nzchar(s)) return("0")
  ch <- strsplit(s, "")[[1]]
  depth <- 0L; start <- 1L; out <- character(0)
  for (i in seq_along(ch)) {
    if (ch[i] == "(") depth <- depth + 1L
    else if (ch[i] == ")") depth <- depth - 1L
    else if (depth == 0L && (ch[i] == "+" || ch[i] == "-") && i > start) {
      j <- i - 1L
      while (j >= 1L && ch[j] == " ") j <- j - 1L
      prev <- if (j >= 1L) ch[j] else ""
      if (grepl("^[A-Za-z0-9_.)]$", prev) &&
          !(prev %in% c("e", "E") && j > 1L && grepl("^[0-9.]$", ch[j - 1L]))) {
        out <- c(out, paste(ch[start:(i - 1L)], collapse = "")); start <- i
      }
    }
  }
  trimws(c(out, paste(ch[start:length(ch)], collapse = "")))
}

# every term as "<sign> <magnitude>", so the signs sit in one column and a wrapped
# continuation line still opens with its operator
.symSignedTerms <- function(s) {
  vapply(.symSplitTerms(s), function(x)
    paste0(if (startsWith(x, "-")) "- " else "+ ", trimws(sub("^[+-]", "", x))),
    character(1), USE.NAMES = FALSE)
}

# One generator as a table, one coordinate per line: eta(i) left, the component right.
# A wide component wraps at its top-level '+'/'-'; a single wide term overflows.
.symFormatGenerator <- function(generator, indent, width = getOption("width")) {
  nm <- names(generator)
  if (!length(nm)) return(paste0(indent, "0"))
  ops <- vapply(nm, .symEta, character(1), USE.NAMES = FALSE)
  ow  <- max(nchar(ops))
  pad <- strrep(" ", nchar(indent) + ow + 4L)
  room <- max(24L, width - nchar(pad))
  unlist(lapply(seq_along(nm), function(i) {
    rows <- character(0); cur <- ""
    for (t in .symSignedTerms(trimws(as.character(generator[[i]])))) {
      if (!nzchar(cur)) cur <- t
      else if (nchar(cur) + 1L + nchar(t) <= room) cur <- paste(cur, t)
      else { rows <- c(rows, cur); cur <- t }
    }
    rows <- c(rows, cur)
    paste0(c(paste0(indent, .symLjust(ops[i], ow), " :  "), rep(pad, length(rows) - 1L)),
           rows)
  }), use.names = FALSE)
}

# directions grouped by type with labels X1, X2, ...; shared by the listing and the
# reduction report
.symOrdered <- function(object) {
  syms <- object$symmetries
  if (!length(syms)) return(list(syms = list(), labels = character(0)))
  types <- vapply(syms, function(d) as.character(d$type), character(1))
  syms <- syms[order(match(types, c("scaling", "general")), seq_along(syms))]
  list(syms = syms,
       labels = vapply(seq_along(syms), function(i) paste0("X", .symSubscript(i)),
                       character(1)))
}

# certificate and verification notes for the label line
.symFlags <- function(d) {
  f <- character(0)
  if (isTRUE(d$certified)) f <- c(f, "certified")
  if (isFALSE(d$verified)) f <- c(f, "unverified")
  paste(f, collapse = ", ")
}

# the display components of a direction: the factored form when finalisation reached
# one, else the canonical components
.symShown <- function(d) if (!is.null(d$display)) d$display else d$generator

# Print the directions grouped by class under a legend, each a label line with notes
# and its component table. Verbose adds the finite transformation and the reason a
# closed form was missed.
.symCatGenerators <- function(object, verbose = FALSE, width = getOption("width")) {
  o <- .symOrdered(object)
  if (!length(o$syms)) return(invisible())
  cat("Generators  X = \u03a3\u1d62 \u03b7(i) \u2202\u1d62\n\n")
  lw <- max(nchar(o$labels)); cur <- ""
  for (i in seq_along(o$syms)) {
    d <- o$syms[[i]]
    if (d$type != cur) { cur <- d$type; cat(.symTypeLabel(cur), ":\n", sep = "") }
    lab <- paste0("  ", .symLjust(o$labels[i], lw), "  ")
    ind <- strrep(" ", nchar(lab))
    notes <- .symFlags(d)
    if (is.null(d$generator)) {
      cat(sub(" +$", "", paste0(lab, notes, if (nzchar(notes)) ", ", "support only")),
          "\n", ind, paste(d$support, collapse = ", "), "\n", sep = "")
      if (isTRUE(verbose) && !is.null(d$reason))
        cat(ind, "reason: ", d$reason, "\n", sep = "")
      cat("\n")
      next
    }
    cat(sub(" +$", "", paste0(lab, notes)), "\n", sep = "")
    cat(.symFormatGenerator(.symShown(d), ind, width), sep = "\n")
    cat("\n")
    if (isTRUE(verbose) && !is.null(d$transformation)) {
      tr <- d$transformation
      cat(ind, "transformation: ",
          paste(vapply(names(tr), function(k) paste0(k, " -> ", tr[[k]]), character(1)),
                collapse = ",  "), "\n", sep = "")
    }
  }
}

# ---- the reduction report ------------------------------------------------------

# Integer weight matrix W of the scalings, one row each. Fixing coordinates S removes
# every scaling iff rank(W[, S]) equals the number of scalings.
.symWeights <- function(syms) {
  isScal <- vapply(syms, function(d) isTRUE(d$type == "scaling") && !is.null(d$weights),
                   logical(1))
  if (!any(isScal)) return(NULL)
  cols <- unique(unlist(lapply(syms[isScal], function(d) names(d$weights))))
  W <- matrix(0, nrow = sum(isScal), ncol = length(cols), dimnames = list(NULL, cols))
  k <- 0L
  for (i in which(isScal)) {
    k <- k + 1L; w <- syms[[i]]$weights
    W[k, names(w)] <- suppressWarnings(as.numeric(unlist(w)))
  }
  if (anyNA(W)) return(NULL)
  list(W = W, rows = which(isScal))
}

.symRank <- function(M) if (!length(M)) 0L else qr(M, tol = 1e-7)$rank

# the coordinates a direction lives on, whether or not it reached a closed form
.symCoords <- function(d) if (!is.null(d$generator)) names(d$generator) else d$support

# How to remove the directions: a scaling by fixing one coordinate at any value, a
# general direction by reparametrisation. With `fixed`, reports the weight rank.
.symCatReduction <- function(object, fixed = NULL, width = getOption("width")) {
  o <- .symOrdered(object)
  if (!length(o$syms)) return(invisible())
  syms <- o$syms; labels <- o$labels
  isScal <- vapply(syms, function(d) isTRUE(d$type == "scaling"), logical(1))
  cat("Reduction:\n")

  if (is.null(fixed)) {
    say <- function(...) cat(strwrap(paste0(...), width = max(40L, width - 2L),
                                     prefix = "  "), sep = "\n")
    if (any(isScal))
      say(paste(labels[isScal], collapse = ", "),
          if (sum(isScal) == 1L)
            " is a scaling: gauge it by holding any one of its coordinates fixed,"
          else " are scalings: gauge them by holding one coordinate of each fixed,",
          " at any value.")
    if (any(!isScal))
      say(paste(labels[!isScal], collapse = ", "),
          if (sum(!isScal) == 1L) " is not a scaling: reparametrise onto its"
          else " are not scalings: reparametrise onto their",
          " invariants with  symmetryReduction(<result>).")
    return(invisible())
  }

  fixed <- unique(as.character(fixed))
  wm <- .symWeights(syms)
  nScal <- if (is.null(wm)) 0L else nrow(wm$W)
  cat("  candidate: ", .symPlural(length(fixed), "coordinate", "coordinates"), "\n", sep = "")
  if (nScal) {
    inW <- intersect(fixed, colnames(wm$W))
    r <- .symRank(wm$W[, inW, drop = FALSE])
    cat(sprintf("  scalings:  %d of %d removed (weight rank %d / %d)%s\n",
                r, nScal, r, nScal,
                if (r == nScal) "" else sprintf(", %d survive", nScal - r)))
    red <- character(0); acc <- character(0)
    for (cn in inW)
      if (.symRank(wm$W[, c(acc, cn), drop = FALSE]) == .symRank(wm$W[, acc, drop = FALSE]))
        red <- c(red, cn) else acc <- c(acc, cn)
    if (length(red))
      cat("             redundant (raise no rank): ", paste(red, collapse = ", "),
          "\n", sep = "")
  }
  if (any(!isScal))
    cat("  general:   ", paste(labels[!isScal], collapse = ", "),
        " (fixing does not settle these; reparametrise instead)\n", sep = "")
  unknown <- setdiff(fixed, unique(unlist(lapply(syms, .symCoords))))
  if (length(unknown))
    cat("  no effect: ", paste(unknown, collapse = ", "),
        " (not a coordinate of any direction)\n", sep = "")
  invisible()
}

# "<n> <singular|plural>"
.symPlural <- function(n, one, many) sprintf("%d %s", n, if (n == 1L) one else many)

# result section of print() and summary(): verdict, generators, optionally the reduction
.symCatResult <- function(object, verbose = FALSE, fixed = NULL,
                          width = getOption("width"), fixing = FALSE) {
  m <- object$method; isObs <- m == "observability"
  n <- length(object$symmetries)
  if (isObs && isTRUE(object$identifiable)) {
    cat(sprintf("Result:  structurally locally identifiable (rank %d / %d)\n",
                object$rank, object$dim))
    return(invisible())
  }
  if (isObs) {
    dirs <- .symPlural(n, "non-identifiable direction", "non-identifiable directions")
    cat(sprintf("Result:  rank %d / %d  |  %s\n\n", object$rank, object$dim, dirs))
  } else if (m == "scaling") {
    cat(sprintf("Result:  %s (exact integer kernel)\n\n",
                .symPlural(n, "scaling symmetry", "scaling symmetries")))
  } else if (n == 0L) {
    cat(sprintf("Result:  no Lie-symmetry generator found up to pMax=%s\n",
                format(object$info$settings$pMax)),
        "         (a non-exhaustive search; not a proof of identifiability)\n", sep = "")
    return(invisible())
  } else {
    cat(sprintf("Result:  %s\n\n",
                .symPlural(n, "polynomial Lie-symmetry generator",
                            "polynomial Lie-symmetry generators")))
  }
  .symCatGenerators(object, verbose, width)
  # the reduction for summary() or an explicit `fixed`
  if (isTRUE(fixing) || !is.null(fixed)) .symCatReduction(object, fixed, width)
}


# ---- finalisation: normalise every engine's raw result into the public object -----
# One raw direction ($vector, or $infinitesimals from the polynomial engine) as a
# public $symmetries element with $generator and, for a scaling, $weights.
.symPublicSymmetry <- function(d) {
  gen <- if (!is.null(d$vector)) d$vector else d$infinitesimals
  weights <- NULL
  # a scaling's integer weights (xi_i = w_i z_i) expand to components and are kept
  if (isTRUE(d$type == "scaling") && !is.null(d$vector) && is.null(d$infinitesimals)) {
    weights <- d$vector
    gen <- setNames(as.list(vapply(names(d$vector),
             function(k) .symScalingComponent(d$vector[[k]], k), character(1))),
             names(d$vector))
  }
  # an eqnvec in R's power syntax, the right-hand side of the flow dz_i/deps = xi_i(z)
  gen <- if (is.null(gen)) NULL else
    as.eqnvec(setNames(gsub("\\*\\*", "^", vapply(gen, as.character, character(1))),
                       names(gen)))
  structure(list(
    type           = sub("^Type: ", "", as.character(d$type)),
    generator      = gen,
    weights        = weights,
    degree         = d$degree,
    support        = if (!is.null(d$support)) d$support else names(gen),
    # explicit iff a full generator exists, whichever engine produced it
    explicit       = !is.null(gen),
    reason         = d$reason,
    certified      = isTRUE(d$certified),
    transformation = d$transformation,
    verified       = if (is.null(d$verified)) NA else isTRUE(d$verified)
  ), class = "symmetrygenerator")
}

# Factored display form of the components, when shorter; $generator stays canonical.
# Skipped without sympy.
.symDisplayForm <- function(syms) {
  if (!length(syms)) return(syms)
  spy <- tryCatch(reticulate::import("sympy", convert = TRUE), error = function(e) NULL)
  if (is.null(spy)) return(syms)
  lapply(syms, function(d) {
    if (is.null(d$generator) || isTRUE(d$type == "scaling")) return(d)
    d$display <- lapply(d$generator, function(x) {
      x <- as.character(x)
      if (nchar(x) > 4000L) return(gsub("\\*\\*", "^", x))
      y <- tryCatch(as.character(spy$factor(spy$sympify(gsub("\\^", "**", x)))),
                    error = function(e) x)
      if (length(y) != 1L || is.na(y) || nchar(y) >= nchar(x)) y <- x
      gsub("\\*\\*", "^", y)
    })
    d
  })
}


# Complete form f * X of a general generator, f = D^2 / (1 + sum_i Q_i^2) with
# Q_i = D^2 xi_i, divided by z_i on a positive coordinate: same orbits and
# invariants, flow defined for all s in R. A scaling is complete already.
.symComplete <- function(syms, positive = TRUE, coordinates = NULL) {
  if (!length(syms)) return(syms)
  if (is.character(positive) && length(coordinates)) {
    unknown <- setdiff(positive, as.character(coordinates))
    if (length(unknown))
      warning("symmetryDetection(): `positive` names no coordinate: ",
              paste(unknown, collapse = ", "), "; ignored.", call. = FALSE)
  }
  spy <- tryCatch(reticulate::import("sympy", convert = TRUE), error = function(e) NULL)
  lapply(syms, function(d) {
    if (is.null(d$generator)) return(d)
    d$completeGenerator <- d$generator
    d$factor            <- "1"
    if (isTRUE(d$type == "scaling")) return(d)
    xi  <- vapply(d$generator, as.character, character(1))
    pos <- if (isTRUE(positive)) names(xi) else if (isFALSE(positive)) character(0)
           else intersect(names(xi), positive)
    D2  <- paste0("(", .symDenominator(xi, spy), ")^2")
    Q   <- ifelse(names(xi) %in% pos,
                  paste0(D2, "*(", xi, ")/", names(xi)), paste0(D2, "*(", xi, ")"))
    f   <- .symTidy(paste0(D2, "/(1 + ", paste0("(", Q, ")^2", collapse = " + "), ")"), spy)
    d$completeGenerator <- as.eqnvec(setNames(paste0("(", f, ")*(", xi, ")"), names(xi)))
    d$factor            <- f
    d
  })
}

# common denominator of the components, "1" when they are polynomial or sympy
# cannot tell
.symDenominator <- function(xi, spy) {
  if (is.null(spy)) return("1")
  tryCatch({
    comp   <- gsub("\\^", "**", xi)
    locals <- .symRedLocals(comp, spy)
    dens   <- lapply(comp, function(x)
      spy$fraction(spy$cancel(spy$together(spy$sympify(x, locals = locals))))[[2]])
    gsub("\\*\\*", "^",
         as.character(Reduce(function(a, b) spy$lcm(a, b), dens, spy$Integer(1L))))
  }, error = function(e) "1")
}

# an expression in R's power syntax, cancelled by sympy when that makes it shorter
.symTidy <- function(x, spy) {
  if (is.null(spy)) return(x)
  y <- tryCatch({
    e <- spy$sympify(gsub("\\^", "**", x), locals = .symRedLocals(x, spy))
    gsub("\\*\\*", "^", as.character(spy$factor(spy$cancel(e))))
  }, error = function(e) x)
  if (length(y) == 1L && !is.na(y) && nchar(y) < nchar(x)) y else x
}

# The public "symmetrydetection" object: verdict at top level, computation in $info.
# `identifiable` is NA for the non-exhaustive scaling and polynomial engines.
.symFinalize <- function(raw, method, settings, call, elapsed = NA_real_,
                          coordinates = NULL) {
  isObs   <- method == "observability"
  rawSyms <- if (isObs || method == "scaling") raw$nonIdentifiable else raw
  if (is.null(rawSyms)) rawSyms <- list()
  syms    <- .symComplete(.symDisplayForm(lapply(rawSyms, .symPublicSymmetry)),
                          settings$positive %||% TRUE, coordinates)

  rank <- if (!is.null(raw$rank)) as.integer(raw$rank) else NA_integer_
  dim  <- if (!is.null(raw$dim))  as.integer(raw$dim)  else NA_integer_
  identifiable <- if (isObs) isTRUE(rank == dim) else NA

  engine <- switch(method,
    scaling    = "integer-kernel",
    polynomial = "lie-ansatz",
    if (identical(raw$engine, "symbolic")) "symbolic" else "modular")

  info <- list(
    engine       = engine,
    lieOrderUsed = raw$lieOrderUsed,
    lieOrderDriver = raw$lieOrderDriver,
    lieBudget = raw$lieBudget,
    liePlateau = raw$liePlateau,
    lieCertified = raw$lieCertified,
    gapOrderUsed = raw$gapOrderUsed,
    conditions   = raw$conditions,
    segments     = raw$segments,
    # all coordinates of the analysis, for symmetryReduction()'s identity trafo
    coordinates  = if (!is.null(raw$coordinates)) as.character(raw$coordinates)
                   else coordinates,
    settings     = settings,
    elapsed      = elapsed,
    verification = raw$verification)

  structure(list(
    method       = method,
    identifiable = identifiable,
    rank         = rank,
    dim          = dim,
    symmetries   = syms,
    info         = info,
    call         = call), class = "symmetrydetection")
}

# ---- the single report renderer shared by print() and summary() -----------------
# a compact "k=v" join of the settings relevant to `method` (more when verbose)
.symSettingsLine <- function(s, method, verbose) {
  keys <- switch(method,
    observability =
      if (verbose) c("reduceCQ", "equilibrate", "reconstruct", "verify",
                     "symEngine", "certifyPoly", "degreeCap")
      else c("reduceCQ", "equilibrate", "reconstruct"),
    scaling    = if (verbose) "reduceCQ" else character(0),
    polynomial = if (verbose) c("ansatz", "pMax", "polyBackend")
                 else c("ansatz", "pMax"))
  keys <- keys[vapply(keys, function(k) !is.null(s[[k]]), logical(1))]
  if (!length(keys)) return(NULL)
  paste(vapply(keys, function(k) paste0(k, "=", format(s[[k]])), character(1)),
        collapse = ", ")
}

# the Schwartz-Zippel saturation-guard verdict, one line (+ detail when verbose)
.symGuardLines <- function(v, verbose) {
  if (!is.list(v)) return(NULL)
  head <- if (isTRUE(v$ok))
      sprintf("saturation guard: PASSED (%s)", v$reason)
    else if (isFALSE(v$ok))
      sprintf("saturation guard: FAILED (%s; directions may be over-reported)", v$reason)
    else
      sprintf("saturation guard: inconclusive (%s)",
              if (!is.null(v$reason)) v$reason else "unavailable")
  out <- head
  if (isTRUE(verbose) && !is.null(v$kernelRank))
    out <- c(out, sprintf("  kernel rank %d (extended %d), checked to Lie order %d%s",
                          v$kernelRank, v$kernelRankExtended, v$ordersChecked,
                          if (!is.na(v$growAt)) sprintf(", grew at %d", v$growAt) else ""))
  out
}

.symReport <- function(object, verbose = FALSE, fixed = NULL,
                       width = getOption("width")) {
  m    <- object$method
  info <- object$info
  s    <- info$settings
  isObs <- m == "observability"
  bar  <- strrep("-", 60)

  engLabel <- switch(as.character(info$engine),
    modular          = "modular (GF(p) + CRT)",
    symbolic         = "symbolic (pure sympy)",
    `integer-kernel` = "integer kernel (exact)",
    `lie-ansatz`     = paste0("Lie ansatz (", if (!is.null(s$polyBackend)) s$polyBackend
                                              else "symengine", ")"),
    as.character(info$engine))

  cat(bar, "\n", sep = "")
  cat(sprintf("symmetryDetection  |  method: %s   engine: %s\n", m, engLabel))
  cat(bar, "\n", sep = "")

  # ---- computation report ----
  comp <- character(0)
  if (isObs && !is.null(info$lieOrderUsed))
    comp <- c(comp, sprintf("Lie order %d (gap order %d)%s",
                            info$lieOrderUsed,
                            if (!is.null(info$gapOrderUsed)) info$gapOrderUsed else 0L,
                            if (!is.null(info$lieOrderDriver))
                              sprintf(", set by condition %d", info$lieOrderDriver) else ""))
  if (isObs && !is.null(info$liePlateau))
    comp <- c(comp, if (isTRUE(info$lieCertified))
        sprintf("saturation: certified (plateau %d > codimension %d)",
                info$liePlateau, info$lieBudget)
      else sprintf("saturation: provisional (plateau %d, %s)", info$liePlateau,
                   if (is.null(info$lieBudget) || is.na(info$lieBudget))
                     "no codimension bound available"
                   else sprintf("codimension %d would certify", info$lieBudget)))
  if (!is.null(info$conditions) && info$conditions > 1L)
    comp <- c(comp, paste0(.symPlural(info$conditions, "condition", "conditions"), ", ",
                           .symPlural(info$segments, "segment", "segments")))
  setline <- .symSettingsLine(s, m, verbose)
  if (!is.null(setline)) comp <- c(comp, paste0("settings: ", setline))
  comp <- c(comp, .symGuardLines(info$verification, verbose))
  if (!is.null(info$elapsed) && is.finite(info$elapsed) && info$elapsed >= 0.05)
    comp <- c(comp, sprintf("elapsed: %.1fs", info$elapsed))
  if (length(comp)) {
    cat("Computation:\n")
    for (l in comp) cat("  ", l, "\n", sep = "")
  }
  cat("\n")
  .symCatResult(object, verbose, fixed, width, fixing = TRUE)
  invisible(object)
}

# print() shows the result section only; summary() adds header and computation
#' @export
print.symmetrydetection <- function(x, fixed = NULL, width = getOption("width"), ...) {
  .symCatResult(x, fixed = fixed, width = width)
  invisible(x)
}

#' @export
summary.symmetrydetection <- function(object, verbose = FALSE, fixed = NULL,
                                      width = getOption("width"), ...)
  .symReport(object, verbose, fixed, width)


# A result of the chart L = log(a + b*v) (logArgChart) reported in v, see
# logArgBackVec; numeric weights make a scaling.
.symLogArgBack <- function(raw, la, sd) {
  if (is.null(la) || is.null(raw) || !is.list(raw)) return(raw)
  ren <- function(x) if (!length(x)) x else replaceSymbols(la$L, la$v, as.character(x))
  back <- function(vec, transform = FALSE) {
    if (is.null(vec) || !length(vec)) return(vec)
    out <- sd$logArgBackVec(lapply(vec, function(x)
      gsub("^", "**", as.character(x), fixed = TRUE)), la$map, transform)
    lapply(out, function(x) gsub("**", "^", x, fixed = TRUE))
  }
  fixDir <- function(d) {
    if (!is.list(d)) return(d)
    abSyms <- getSymbols(unlist(lapply(la$map, function(m) c(m$a, m$b))))
    hit <- any(c(names(d$vector), names(d$infinitesimals)) %in% c(la$L, abSyms)) ||
      any(grepl(paste0("\\b(", paste(la$L, collapse = "|"), ")\\b"),
                unlist(c(d$vector, d$infinitesimals))))
    if (!hit) return(d)
    if (!is.null(d$vector)) {
      vec <- d$vector
      if (isTRUE(d$type == "scaling"))
        vec <- setNames(lapply(names(vec), function(k)
          .symScalingComponent(vec[[k]], k)), names(vec))
      vec <- back(vec)
      w <- tryCatch(sd$scalingWeights(lapply(vec, function(x) gsub("^", "**", x, fixed = TRUE))),
                    error = function(e) NULL)
      if (!is.null(w)) {
        d$vector <- lapply(w, as.character); d$type <- "scaling"; d$degree <- 1L
      } else {
        d$vector <- vec; d$type <- "general"
        if (!is.null(d$degree)) d$degree <- -1L
      }
      d$support <- .symSort(names(d$vector))
    } else if (!is.null(d$support)) d$support <- .symSort(ren(d$support))
    if (!is.null(d$infinitesimals)) d$infinitesimals <- back(d$infinitesimals)
    if (!is.null(d$transformation)) d$transformation <- back(d$transformation, TRUE)
    d
  }
  if (!is.null(raw$coordinates)) raw$coordinates <- ren(raw$coordinates)
  if (!is.null(raw$nonIdentifiable))
    raw$nonIdentifiable <- lapply(raw$nonIdentifiable, fixDir)
  else if (is.null(raw$rank) && is.null(raw$method)) raw <- lapply(raw, fixDir)
  raw
}
