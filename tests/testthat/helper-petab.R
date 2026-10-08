# Shared fixtures of the PEtab tests (test-petab*.R).

# PEtabTests/ is .Rbuildignore'd and found through DMOD_PETABTESTS, which
# setup.R sets by walking up from the package root; CI may export it.
.petab_repo_dir <- function() {
  p <- Sys.getenv("DMOD_PETABTESTS", unset = "")
  if (nzchar(p) && dir.exists(p)) normalizePath(p, winslash = "/") else ""
}

# Skip integration tests that need importSbml() when the libsbml virtualenv
# is missing. The probe ships inside tests/, unlike PEtabTests/, so R CMD check
# reaches it. Absolute, since tests switch to tempdir().
.libsbml_probe_file <- normalizePath(file.path("fixtures", "petab_probe_model.xml"),
                                     winslash = "/", mustWork = FALSE)

.libsbml_works <- function() {
  if (!file.exists(.libsbml_probe_file)) return(FALSE)
  isTRUE(tryCatch({
    res <- suppressWarnings(importSbml(.libsbml_probe_file))
    !is.null(res$reactions)
  }, error = function(e) FALSE))
}

.petab_boehm_yaml <- function()
  file.path(system.file("extdata/petab_boehm", package = "dMod2"), "Boehm.yaml")

# Links imported problems into one shared object. Their sources must agree on
# preprocessor macros, see compile().
.petab_compile <- function(problems, output) {
  objs <- do.call(c, lapply(problems, function(pp) list(pp$prd, pp$e)))
  do.call(compile, c(Filter(Negate(is.null), objs),
                     list(output = output, cores = test_cores())))
}

# Every problem the file imports as published, imported once and linked into
# one shared object. A failed import is kept as its error for the test to raise.
.petab_published <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    d <- file.path(tempdir(), "dmod_petab_published")
    dir.create(d, showWarnings = FALSE)
    owd <- setwd(d); on.exit(setwd(owd))

    spec <- function(yaml, backend, quiet = FALSE)
      list(yaml = yaml, backend = backend, quiet = quiet)
    specs <- list(boehm = spec(.petab_boehm_yaml(), "deSolve"))
    petab_dir <- .petab_repo_dir()
    if (nzchar(petab_dir)) {
      for (id in sprintf("%04d", 1:16)) {
        if (!file.exists(file.path(petab_dir, id, paste0("_", id, "_solution.yaml")))) next
        specs[[paste0("v1_", id)]] <- spec(file.path(petab_dir, id, paste0("_", id, ".yaml")),
                                           "deSolve", quiet = as.integer(id) >= 7L)
      }
      for (case in c("0001", "0002", "0009", "0016", "0024", "0030"))
        specs[[paste0("v2_", case)]] <- spec(
          file.path(petab_dir, "v2", case, paste0("_", case, ".yaml")),
          if (case %in% c("0001", "0002", "0009")) "deSolve" else "cppDE")
    }
    specs <- Filter(function(s) file.exists(s$yaml), specs)

    problems <- Map(function(key, s) tryCatch({
      imp <- function() importPEtab(s$yaml, backend = s$backend, modelname = key,
                                    compile = FALSE)
      if (s$quiet) suppressWarnings(imp()) else imp()
    }, error = identity), names(specs), specs)
    .petab_compile(Filter(function(pp) !inherits(pp, "error"), problems),
                   "petab_published")
    cache <<- problems
    cache
  }
})

.petab_case <- function(key) {
  pp <- .petab_published()[[key]]
  if (inherits(pp, "error")) stop(pp)
  pp
}

# Boehm under cppDE, exported by two tests. Its KLU macros keep it out of the
# shared object above.
.petab_boehm_cppDE <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    d <- file.path(tempdir(), "dmod_petab_boehm")
    dir.create(d, showWarnings = FALSE)
    owd <- setwd(d); on.exit(setwd(owd))
    cache <<- importPEtab(.petab_boehm_yaml(), backend = "cppDE",
                          modelname = "v1rtA", cores = test_cores())
    cache
  }
})

# The models the hand-built and native-export tests evaluate, built once and
# linked into one shared object. None of them needs libsbml.
.petab_native <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    d <- file.path(tempdir(), "dmod_petab_native")
    dir.create(d, showWarnings = FALSE)
    owd <- setwd(d); on.exit(setwd(owd))

    ab <- eqnlist() %>%
      addReaction("A", "B", rate = "k1*A", description = "fwd") %>%
      addReaction("B", "A", rate = "k2*B", description = "rev")
    a <- eqnlist() %>%
      addReaction("A", "B", rate = "k*A", description = "fwd")
    x_ab <- Xs(odemodel(ab, modelname = "nat_ab", backend = "deSolve", compile = FALSE))
    x_a  <- Xs(odemodel(a,  modelname = "nat_a",  backend = "deSolve", compile = FALSE))

    obs_ab <- eqnvec(obs_a = "A", obs_b = "B")
    obs_a  <- eqnvec(obs_a = "A")
    g_0  <- Y(g = c(obs_a = "A"), f = ab, attachInput = FALSE,
              modelname = "nat_obs_0")
    g_ab <- Y(obs_ab, f = x_ab, condition = NULL, attachInput = FALSE,
              modelname = "nat_obs_ab")
    g_a  <- Y(obs_a, f = x_a, condition = NULL, attachInput = FALSE,
              modelname = "nat_obs_a")

    innerpars <- getParameters(x_ab)
    tr_0 <- structure(innerpars, names = innerpars)
    tr_0["A"] <- "a0"
    tr_0["B"] <- "b0"
    p_0 <- P(tr_0, condition = "c0", modelname = "nat_par_0")
    p_rt1 <- P(as.eqnvec(c(A = "10^(A)", B = "10^(B)",
                           k1 = "10^(K1)", k2 = "10^(K2)")),
               condition = "c1", modelname = "rt1_par")
    p_rt2 <- P(as.eqnvec(c(A = "10^(A)", B = "10^(B)", k = "10^(K)")),
               condition = "closed", modelname = "rt2_par_c") +
             P(as.eqnvec(c(A = "10^(A)", B = "10^(B)", k = "10^(K_OPEN)")),
               condition = "open", modelname = "rt2_par_o")
    p_sig <- P(as.eqnvec(c(A = "10^(A)", B = "10^(B)", k = "10^(K)")),
               condition = "c1", modelname = "rt_sig_par")
    p_rt3 <- P(as.eqnvec(c(A = "10^(A)", B = "10^(B)", k = "10^(K + 5)")),
               condition = "c1", modelname = "rt3_par")

    compile(x_ab, x_a, g_0, g_ab, g_a, p_0, p_rt1, p_rt2, p_sig, p_rt3,
            output = "petab_native", cores = test_cores())
    cache <<- list(ab = ab, a = a, x_ab = x_ab, x_a = x_a,
                   obs_ab = obs_ab, obs_a = obs_a,
                   g_0 = g_0, g_ab = g_ab, g_a = g_a,
                   p_0 = p_0, p_rt1 = p_rt1, p_rt2 = p_rt2,
                   p_sig = p_sig, p_rt3 = p_rt3)
    cache
  }
})


# A module whose preequilibration leaves some states unmoved.
.petab_ac_module <- function()
  eqnlist() |>
    addReaction("AC",   "ACF",  "kf*Fsk*AC",        compartment = "cell") |>
    addReaction("ACF",  "AC",   "kr*ACF",           compartment = "cell") |>
    addReaction("AC",   "pAC",  "kp*AC",            compartment = "cell") |>
    addReaction("pAC",  "AC",   "kdp*pAC",          compartment = "cell") |>
    addReaction("",     "cAMP", "ks*(AC + xi*ACF)", compartment = "cell") |>
    addReaction("cAMP", "",     "kd*cAMP*PDE",      compartment = "cell") |>
    addReaction("PDE",  "pPDE", "kpp*cAMP*PDE",     compartment = "cell") |>
    addReaction("pPDE", "PDE",  "kpd*pPDE",         compartment = "cell")

# exportPEtab writes an outer parameter named like a state as init_<state>.
.init_ids <- function(p) {
  hit <- names(p) %in% c("A", "B")
  names(p)[hit] <- paste0("init_", names(p)[hit])
  p
}
