## =====================================================================
##  harness.R: time and score the objective of a PEtab problem.
## =====================================================================

bench_source <- function(dir = "benchmarks/R") {
  for (f in c("petab.R", "tiers.R", "plots.R", "report.R")) {
    p <- file.path(dir, f)
    if (file.exists(p)) source(p)
  }
  invisible(TRUE)
}


## ---------------------------------------------------------------------
##  Timing
## ---------------------------------------------------------------------

##  The minimum over `bursts`, each burst long enough (`target` seconds) that
##  a sub-millisecond call is not measured as clock noise. The machines
##  scatter; the minimum is the number that repeats.
bench_time <- function(fn, bursts = 5L, target = 0.25) {
  once <- system.time(fn())[["elapsed"]]
  n <- max(1L, ceiling(target / max(once, 1e-3)))
  min(vapply(seq_len(bursts), function(i)
    system.time(for (j in seq_len(n)) fn())[["elapsed"]] / n, 0))
}


## ---------------------------------------------------------------------
##  The arms
## ---------------------------------------------------------------------

##  What a reverse call returns says which direction answered it. A wrapper
##  that drops `sweep` would otherwise time the forward mode as the adjoint.
bench_reverse <- function(obj, p) {
  r <- obj(p, deriv = TRUE, sweep = "reverse")
  if (!identical(attr(r, "sweep"), "reverse"))
    stop("sweep = \"reverse\" was answered by ",
         if (is.null(attr(r, "sweep"))) "the forward mode" else attr(r, "sweep"),
         call. = FALSE)
  r
}

##  The reference gradient: the forward sensitivities at a tight tolerance.
bench_reference <- function(pet, tol, ref_atol = 1e-12, ref_rtol = 1e-10) {
  petab_tolerance(pet, ref_atol, ref_rtol)
  on.exit(petab_tolerance(pet, tol, tol))
  pet$obj(pet$bestfit, deriv = TRUE, hessian = FALSE)$gradient
}

##  max |g - g_ref| / max |g_ref|, and the CVODES-style score
##  max |g - g_ref| / (rtol |g_ref| + gradtol), at most 1 within tolerance.
bench_error <- function(g, ref, rtol, gradtol) {
  nm <- intersect(names(ref), names(g))
  d <- abs(g[nm] - ref[nm])
  c(rel = max(d) / max(abs(ref[nm])), tau = max(d / (rtol * abs(ref[nm]) + gradtol)))
}

##  Every arm of one problem. `asa` is a second import on the Sundials backend.
bench_problem <- function(pet, entry, tol, bursts, target, asa = NULL) {
  obj <- pet$obj; p <- pet$bestfit
  ref <- bench_reference(pet, tol)
  gradtol <- tol * 1e-3 * max(abs(ref))
  arms <- list(
    value   = function() obj(p, deriv = FALSE),
    forward = function() obj(p, deriv = TRUE, hessian = FALSE),
    reverse = function() bench_reverse(obj, p))
  out <- list()
  for (a in names(arms)) {
    r <- arms[[a]]()
    out[[a]] <- list(sec = bench_time(arms[[a]], bursts, target),
                     err = if (a == "value") c(rel = NA, tau = NA)
                           else bench_error(r$gradient, ref, tol, gradtol))
  }
  dMod2::controls(pet$x, name = "optionsReverse") <- list(refine = TRUE, gradtol = gradtol)
  r <- bench_reverse(obj, p)
  out$refine <- list(sec = bench_time(function() bench_reverse(obj, p), bursts, target),
                     err = bench_error(r$gradient, ref, tol, gradtol))
  dMod2::controls(pet$x, name = "optionsReverse") <- NULL
  if (!is.null(asa)) {
    oa <- asa$obj
    r <- bench_reverse(oa, p)
    out$asa <- list(sec = bench_time(function() bench_reverse(oa, p), bursts, target),
                    err = bench_error(r$gradient, ref, tol, gradtol))
  }
  split <- bench_split(pet, p)
  data.frame(model = entry$short, problem = entry$name, n_theta = length(p),
             n_states = length(pet$reactions$states), n_conditions = length(pet$dataList),
             compile_s = attr(pet, "compile_s"), value = obj(p, deriv = FALSE)$value,
             arm = names(out), sec = vapply(out, `[[`, 0, "sec"),
             cost = vapply(out, `[[`, 0, "sec") / out$value$sec,
             err_rel = vapply(out, function(o) o$err[["rel"]], 0),
             err_tau = vapply(out, function(o) o$err[["tau"]], 0),
             solve_frac = unname(split[names(out)]),
             gradtol = gradtol, stringsAsFactors = FALSE, row.names = NULL)
}


## ---------------------------------------------------------------------
##  Solve against R
## ---------------------------------------------------------------------

##  The share of a call spent in cppDE, from wrapping its two entry points for
##  the duration of a 2 s loop; a nested call counts once. The rest is dMod2's R.
bench_split <- function(pet, p, seconds = 2) {
  ns <- asNamespace("cppDE")
  acc <- new.env()
  acc$depth <- 0L
  wrap <- function(nm) {
    f0 <- get(nm, envir = ns)
    f1 <- function(...) {
      acc$depth <- acc$depth + 1L
      t0 <- proc.time()[["elapsed"]]
      on.exit({
        acc$depth <- acc$depth - 1L
        if (acc$depth == 0L) acc$t <- acc$t + proc.time()[["elapsed"]] - t0
      })
      f0(...)
    }
    unlockBinding(nm, ns); assign(nm, f1, envir = ns); lockBinding(nm, ns)
    f0
  }
  orig <- lapply(c(".batchRun", "solveODE"), wrap)
  on.exit(for (i in 1:2) {
    nm <- c(".batchRun", "solveODE")[i]
    unlockBinding(nm, ns); assign(nm, orig[[i]], envir = ns); lockBinding(nm, ns)
  })
  frac <- function(fn) {
    fn(); acc$t <- 0; n <- 0L
    t0 <- proc.time()[["elapsed"]]
    while (proc.time()[["elapsed"]] - t0 < seconds) { fn(); n <- n + 1L }
    acc$t / (proc.time()[["elapsed"]] - t0)
  }
  c(value   = frac(function() pet$obj(p, deriv = FALSE)),
    forward = frac(function() pet$obj(p, deriv = TRUE, hessian = FALSE)),
    reverse = frac(function() pet$obj(p, deriv = TRUE, sweep = "reverse")),
    refine = NA, asa = NA)
}
