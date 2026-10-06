## =====================================================================
##  tiers.R: the problem sets, and what each problem exercises.
## =====================================================================

##  A tier is an explicit list, so a run is comparable to the last one. Each
##  larger tier contains the smaller ones.

BENCH_TIERS <- list(
  tiny   = c("Boehm", "Raia", "Lucarelli"),
  medium = c("Boehm", "Raia", "Lucarelli", "Bachmann", "Isensee"),
  full   = c("Boehm", "Raia", "Lucarelli", "Bachmann", "Isensee", "Chen", "Lang")
)

##  Traits rather than size: each names a path through dMod2 that a problem
##  takes, and the tier coverage says which paths a run did not see.
BENCH_TRAITS <- c(
  conditions = "many experimental conditions (>= 10)",
  events     = "events or switching times",
  preeq      = "preequilibration (Pimpl steady state)",
  parameters = "many estimated parameters (>= 50)",
  klu        = "sparse linear solver (KLU)",
  errormodel = "estimated error model parameters"
)

##  The traits of an imported problem, from the problem itself.
bench_traits <- function(pet, entry) {
  meta <- attr(pet, "petab_meta")
  tabs <- tryCatch(dMod2::readPetabTables(entry$yaml), error = function(e) NULL)
  meas <- tabs$measurements
  preeq <- !is.null(meas$preequilibrationConditionId) &&
    any(nzchar(stats::na.omit(as.character(meas$preequilibrationConditionId))))
  ev <- meta$eventsSource
  noise <- unlist(meta$obs_meta$noise)
  noise_syms <- unique(unlist(lapply(noise, function(s)
    all.vars(tryCatch(parse(text = s)[[1L]], error = function(e) NULL)))))
  ode <- pet$odemodel
  klu <- isTRUE(attr(ode$extended, "sparse")) || isTRUE(attr(ode$reversed, "sparse"))
  c(conditions = length(pet$dataList) >= 10L,
    events     = (!is.null(ev) && NROW(ev) > 0L) || length(meta$switches) > 0L,
    preeq      = preeq,
    parameters = length(pet$bestfit) >= 50L,
    klu        = klu,
    errormodel = any(noise_syms %in% names(pet$bestfit)) ||
      any(grepl("^noiseParameter", noise_syms)))
}

##  Problems per trait over a run, with a warning for a trait nobody covers.
bench_coverage <- function(traits) {
  if (!length(traits)) return(invisible(NULL))
  m <- do.call(rbind, traits)
  cov <- colSums(m)
  cat("\nTrait coverage\n")
  for (k in names(BENCH_TRAITS))
    cat(sprintf("  %-11s %2d  %s%s\n", k, cov[[k]], BENCH_TRAITS[[k]],
                if (cov[[k]] == 0L) "  <- none" else ""))
  if (any(cov == 0L))
    warning("no problem in this run covers: ",
            paste(names(cov)[cov == 0L], collapse = ", "), call. = FALSE)
  invisible(cov)
}
