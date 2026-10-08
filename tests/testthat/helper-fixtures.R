# Compiled fixtures shared by the test files. testthat sources helpers once per
# file and process, so the cache lives in globalenv and holds one build per file.

.dmod_fx_cache <- function() {
  if (!exists("..dmod_fx_cache..", envir = globalenv(), inherits = FALSE))
    assign("..dmod_fx_cache..", new.env(parent = emptyenv()),
           envir = globalenv())
  get("..dmod_fx_cache..", envir = globalenv())
}

.dmod_fx_workdir <- function() {
  cache <- .dmod_fx_cache()
  if (is.null(cache$workdir)) {
    cache$workdir <- file.path(tempdir(), "dmod_fx")
    dir.create(cache$workdir, recursive = TRUE, showWarnings = FALSE)
  }
  cache$workdir
}

.dmod_with_fx_workdir <- function(expr) {
  oldwd <- setwd(.dmod_fx_workdir())
  on.exit(setwd(oldwd), add = TRUE)
  force(expr)
}


## ---- Registration ---------------------------------------------------------

# A test file names, before first use, what it links into the fixture's shared
# object: `extra` maps the uncompiled fixture to a named list of fn objects,
# `multicond` adds the four-condition chain. Everything then builds in one go.
fx_register <- function(extra = NULL, multicond = FALSE) {
  cache <- .dmod_fx_cache()
  if (!is.null(cache$decay))
    stop("fx_register() has to run before the fixture is built")
  if (!is.null(extra)) cache$extraFns <- c(cache$extraFns, list(extra))
  if (multicond) cache$wantMulticond <- TRUE
  invisible(NULL)
}

# The registered extra objects, compiled.
fx_extra <- function() {
  fx_decay_compiled()
  .dmod_fx_cache()$extra
}


## ---- Linear decay ----------------------------------------------------------

# One-state decay dA/dt = -k*A observed as y = A, with an identity and a log
# trafo for condition C1. Elements: m, xfn, gfn, pfn_id, pfn_log, prd_id,
# prd_log, outerpars_id, outerpars_log.
fx_decay_compiled <- function() {
  cache <- .dmod_fx_cache()
  if (!is.null(cache$decay)) return(cache$decay)

  .dmod_with_fx_workdir({
    reactions <- addReaction(eqnlist(), from = "A", to = "", rate = "k*A",
                             description = "linear decay")
    m <- odemodel(reactions, modelname = "fx_decay", compile = FALSE)
    xfn <- Xs(m)
    gfn <- Y(c(y = "A"), f = xfn, condition = NULL, attachInput = FALSE,
             modelname = "fx_decay_obs", compile = FALSE)
    pfn_id <- P(eqnvec(A = "A", k = "k"), condition = "C1",
                modelname = "fx_decay_p_id", compile = FALSE)
    pfn_log <- P(eqnvec(A = "exp(A_log)", k = "exp(k_log)"), condition = "C1",
                 modelname = "fx_decay_p_log", compile = FALSE)
    base <- list(m = m, xfn = xfn, gfn = gfn, pfn_id = pfn_id, pfn_log = pfn_log)

    mc <- if (isTRUE(cache$wantMulticond)) .fx_multicond_parts(base)
    extra <- do.call(c, lapply(cache$extraFns, function(f) f(base)))
    objs <- c(list(xfn, gfn, pfn_id, pfn_log), mc[c("gfn", "pfn")], unname(extra))
    do.call(compile, c(Filter(Negate(is.null), objs),
                       list(output = "fx_decay_all", cores = test_cores())))

    cache$extra <- extra
    if (!is.null(mc)) cache$decay_mc <- .fx_multicond_finish(base, mc)
    cache$decay <- c(base, list(
      prd_id        = gfn * xfn * pfn_id,
      prd_log       = gfn * xfn * pfn_log,
      outerpars_id  = c(A = 1.0, k = 0.5),
      outerpars_log = c(A_log = 0, k_log = log(0.5))))
  })

  cache$decay
}


## ---- Four conditions -------------------------------------------------------

# The decay chain branched over C1..C4 with a scale s_<C>_log of its own per
# condition, so every condition has a different derivative basis. Two or more
# conditions put every leaf on its batched entry. Elements: conditions, m, xfn,
# gfn (y = s*A), pfn, prd, outerpars.
.fx_multicond_parts <- function(base) {
  conds <- paste0("C", 1:4)
  gfn <- Y(c(y = "s*A"), f = base$xfn, condition = NULL, attachInput = FALSE,
           modelname = "fx_mc_obs", compile = FALSE)
  tree <- data.frame(s_log = paste0("s_", conds, "_log"),
                     row.names = conds, stringsAsFactors = FALSE)
  pfn <- P(branch(eqnvec(A = "exp(A_log)", k = "exp(k_log)", s = "exp(s_log)"),
                  table = tree, apply = "insert"),
           method = "explicit", modelname = "fx_mc_p", compile = FALSE)
  list(conditions = conds, gfn = gfn, pfn = pfn)
}

.fx_multicond_finish <- function(base, mc) {
  list(conditions = mc$conditions, m = base$m, xfn = base$xfn, gfn = mc$gfn,
       pfn = mc$pfn, prd = mc$gfn * base$xfn * mc$pfn,
       outerpars = c(A_log = 0, k_log = log(0.5),
                     setNames(seq(0, 0.3, length.out = 4),
                              paste0("s_", mc$conditions, "_log"))))
}

fx_decay_multicond_compiled <- function() {
  cache <- .dmod_fx_cache()
  if (!is.null(cache$decay_mc)) return(cache$decay_mc)
  base <- fx_decay_compiled()
  if (!is.null(cache$decay_mc)) return(cache$decay_mc)
  .dmod_with_fx_workdir({
    mc <- .fx_multicond_parts(base)
    compile(mc$gfn, mc$pfn, output = "fx_mc_all", cores = test_cores())
    cache$decay_mc <- .fx_multicond_finish(base, mc)
  })
  cache$decay_mc
}


## ---- Data ------------------------------------------------------------------

# Noisy decay data for one condition, closed-form A(t) plus Gaussian noise.
fx_decay_data <- function(pars  = c(A = 1.0, k = 0.5),
                          times = seq(0, 10, by = 1),
                          sigma = 0.05,
                          condition = "C1",
                          seed  = 1L) {
  df <- make_noisy_data(
    truth_fn  = function(t, p) p["A"] * exp(-p["k"] * t),
    pars = pars, times = times, name = "y", sigma = sigma,
    condition = condition, seed = seed)
  as.datalist(df, splitBy = "condition")
}

# Noisy decay data for several conditions, one parameter set each.
fx_decay_data_multi <- function(parslist = list(C1 = c(A = 1.0, k = 0.5),
                                                C2 = c(A = 1.0, k = 1.0)),
                                times = seq(0, 10, by = 1),
                                sigma = 0.05,
                                seed = 1L) {
  out <- do.call(rbind, lapply(seq_along(parslist), function(i) {
    set.seed(seed + i - 1L)
    p <- parslist[[i]]
    data.frame(name = "y", time = times,
               value = p["A"] * exp(-p["k"] * times) +
                 rnorm(length(times), 0, sigma),
               sigma = sigma, condition = names(parslist)[i],
               stringsAsFactors = FALSE)
  }))
  as.datalist(out, splitBy = "condition")
}

# Decay data with a lower limit of quantification that censors the late rows.
fx_decay_data_bloq <- function(pars  = c(A = 1.0, k = 0.5),
                               times = seq(0, 10, by = 0.5),
                               sigma = 0.05,
                               lloq  = 0.1,
                               condition = "C1",
                               seed  = 1L) {
  df <- make_noisy_data(
    truth_fn  = function(t, p) p["A"] * exp(-p["k"] * t),
    pars = pars, times = times, name = "y", sigma = sigma,
    condition = condition, lloq = lloq, seed = seed)
  as.datalist(df, splitBy = "condition")
}
