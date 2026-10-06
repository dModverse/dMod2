## PEtab importer / exporter: the published test cases and Boehm.

## --- libsbml-dependent integration tests ----------------------------------

test_that("PEtab test cases 0001-0006 import and produce solution-matching llh", {

  withr::local_dir(tempdir())
  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ not found -- set DMOD_PETABTESTS to the repo directory")
  if (!.libsbml_works())   skip("libsbml virtualenv not available")

  for (id in sprintf("%04d", 1:6)) {
    sol_path  <- file.path(petab_dir, id, paste0("_", id, "_solution.yaml"))
    if (!file.exists(sol_path)) next

    petab <- .petab_case(paste0("v1_", id))
    sol <- yaml::read_yaml(sol_path)

    out <- petab$obj(petab$bestfit, deriv = FALSE)
    # dMod normL2 returns -2*log L (chi2 + log normaliser); compare against
    # -2 * sol$llh so all 6 cases share the same metric.
    expect_lt(abs(out$value - (-2 * sol$llh)),
              max(0.01, abs(2 * sol$tol_llh)),
              label = paste0("case ", id, " -2*llh"))
  }
})



test_that("PEtab Stage-2 test cases 0007-0016 produce solution-matching llh", {

  withr::local_dir(tempdir())
  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ not found -- set DMOD_PETABTESTS to the repo directory")
  if (!.libsbml_works())   skip("libsbml virtualenv not available")

  # 0007 log10, 0008 replicates, 0009/0010 preequilibration, 0011-0013 init
  # overrides, 0014/0015 noise overrides, 0016 log.
  for (id in sprintf("%04d", 7:16)) {

    sol_path  <- file.path(petab_dir, id, paste0("_", id, "_solution.yaml"))
    if (!file.exists(sol_path)) next

    petab <- .petab_case(paste0("v1_", id))
    sol <- yaml::read_yaml(sol_path)
    out <- petab$obj(petab$bestfit, deriv = FALSE)
    expect_lt(abs(out$value - (-2 * sol$llh)),
              max(0.01, abs(2 * sol$tol_llh)),
              label = paste0("case ", id, " -2*llh"))
  }
})



test_that("two-condition roundtrip preserves objective value", {

  withr::local_dir(tempdir())
  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ not found -- set DMOD_PETABTESTS to the repo directory")
  if (!.libsbml_works())   skip("libsbml virtualenv not available")

  # Two conditions and InitialAssignments, which the roundtrip has to keep.
  petab1 <- .petab_case("v1_0002")
  v1 <- petab1$obj(petab1$bestfit, deriv = FALSE)$value

  out_dir <- file.path(tempdir(), "petab_roundtrip")
  # Exported as v1: a v2 export of a v1 import loses state initials that live
  # only on the trafo.
  yaml2 <- exportPEtabObject(petab1, out_dir, modelID = "rt_out",
                             formatVersion = "1", overwrite = TRUE)

  petab2 <- importPEtab(yaml2, backend = "deSolve",
                        modelname = "rt_back", cores = test_cores())
  v2 <- petab2$obj(petab2$bestfit, deriv = FALSE)$value

  expect_true(is.finite(v1))
  expect_true(is.finite(v2))
  expect_lt(abs(v1 - v2), 1e-6)

  unlink("rt_*"); unlink("*.c"); unlink("*.cpp")
  unlink("*.o"); unlink("*.so")
})



## --- real-world benchmark: Boehm_JProteomeRes2014 -------------------------
##
## Published JAK/STAT5 benchmark: log10 scales, <power/> MathML, an assignment
## rule input and per-observable noise symbols. At the published optimum
## -log L reproduces the benchmark value 138.22 (Hass et al. 2019).

test_that("the bundled Boehm problem imports and matches the published optimum", {

  withr::local_dir(tempdir())
  if (!.libsbml_works())  skip("libsbml virtualenv not available")

  # The package ships this problem, so the check does not depend on a
  # third-party fixture tree being present.
  petab <- .petab_case("boehm")

  # Imported problem shape:
  expect_equal(length(petab$bestfit), 9L)
  expect_setequal(names(attr(petab, "petab_meta")$obs_meta$obs),
                  c("pSTAT5A_rel", "pSTAT5B_rel", "rSTAT5A_rel"))
  # All estimated parameters are on log10 scale per parameters.tsv:
  scales <- attr(petab$bestfit, "petab_scales")
  expect_true(all(scales == "log10"))
  # AssignmentRule for BaF3_Epo must have been inlined → not in `fixed`:
  expect_false("BaF3_Epo" %in% names(attr(petab, "petab_meta")$fixed))

  out <- petab$obj(petab$bestfit, deriv = FALSE)

  # Published optimum: -log L = 138.22 (Hass et al. 2019, "Benchmark
  # problems for dynamic modeling of intracellular processes"). dMod's
  # normL2 returns -2*log L, so we compare against ~276.44.
  expect_lt(abs(out$value - 2 * 138.22), 0.5,
            label = "Boehm -2*logL at published optimum")
})



test_that("exportPEtabObject v2 writes nominalValue verbatim (no parameterScale linearisation)", {
  if (!.libsbml_works()) skip("libsbml virtualenv not available")
  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ not found -- set DMOD_PETABTESTS to the repo directory")

  withr::local_dir(tempdir())
  pp <- .petab_case("v1_0001")
  td <- tempfile("petab_v2_lin_"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  # A v2 export of a log10-scaled v1 problem is silent: the trafo already
  # carries the scale as 10^(...), which v2 keeps.
  expect_silent(
    exportPEtabObject(pp, dir = td, formatVersion = "2.0.0",
                      overwrite = TRUE))

  par_path <- list.files(td, pattern = "^parameters_.*\\.tsv$",
                         full.names = TRUE)
  par_df <- read.delim(par_path, stringsAsFactors = FALSE, na.strings = "")
  expect_false("parameterScale" %in% colnames(par_df))
  # nominalValue is the internal pouter.
  est <- par_df[par_df$estimate == "true", , drop = FALSE]
  expect_equal(est$nominalValue,
               unname(pp$bestfit[est$parameterId]))
})



test_that("exportPEtabObject v2 writes long-format conditions and experiments", {
  if (!.libsbml_works()) skip("libsbml virtualenv not available")
  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ not found -- set DMOD_PETABTESTS to the repo directory")

  withr::local_dir(tempdir())
  # 0001 has a single condition; trivial v2 export should produce one
  # experimentId row.
  petab <- .petab_case("v1_0001")
  td <- tempfile("petab_v2_out_"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  yamlPath <- exportPEtabObject(petab, dir = td, formatVersion = "2.0.0",
                                  overwrite = TRUE)

  cond_path <- list.files(td, pattern = "^conditions_.*\\.tsv$",
                          full.names = TRUE)
  expect_length(cond_path, 1L)
  cond <- read.delim(cond_path, stringsAsFactors = FALSE, na.strings = "")
  expect_setequal(colnames(cond), c("conditionId", "targetId", "targetValue"))

  expect_length(list.files(td, pattern = "^experiments_.*\\.tsv$"), 1L)
  m <- yaml::read_yaml(yamlPath)
  expect_identical(m$format_version, "2.0.0")
  expect_true("model_files" %in% names(m))
  expect_equal(m$model_files[[1]]$language, "sbml")
})



test_that("v2 export → v2 import roundtrips the objective on case 0001", {
  if (!.libsbml_works()) skip("libsbml virtualenv not available")
  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ not found -- set DMOD_PETABTESTS to the repo directory")

  withr::local_dir(tempdir())
  pp1 <- .petab_case("v1_0001")
  td <- tempfile("v2_rt_"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  yamlPath <- exportPEtabObject(pp1, dir = td, formatVersion = "2.0.0",
                                 overwrite = TRUE)

  setwd(td)
  pp2 <- importPEtab(yamlPath, backend = "deSolve", compile = TRUE,
                     modelname = "v2rt_0001", cores = test_cores())

  v1 <- pp1$obj(pp1$bestfit)$value
  v2 <- pp2$obj(pp2$bestfit)$value
  expect_equal(v1, v2, tolerance = 1e-6)
})



test_that("v2 PEtab test cases 0001/0002/0009 import and match published llh", {
  withr::local_dir(tempdir())
  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ not found -- set DMOD_PETABTESTS to the repo directory")
  v2_dir <- file.path(petab_dir, "v2")
  if (!dir.exists(v2_dir)) skip("PEtabTests/v2/ not present")
  if (!.libsbml_works()) skip("libsbml virtualenv not available")

  for (case in c("0001", "0002", "0009")) {
    yamlPath <- file.path(v2_dir, case, paste0("_", case, ".yaml"))
    if (!file.exists(yamlPath)) next
    sol_path  <- file.path(v2_dir, case, paste0("_", case, "_solution.yaml"))
    sol <- yaml::read_yaml(sol_path)
    res <- tryCatch({
      pp <- .petab_case(paste0("v2_", case))
      pp$obj(pp$bestfit)$value
    }, error = function(e) {
      message("v2 case ", case, " import error: ", conditionMessage(e))
      NA_real_
    })
    expect_equal(res, -2 * as.numeric(sol$llh), tolerance = 1e-3,
                 info = sprintf("v2 case %s", case))
  }
})



test_that("v2 experiment periods and promoted event targets match the published llh", {
  withr::local_dir(tempdir())
  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ not found -- set DMOD_PETABTESTS to the repo directory")
  v2_dir <- file.path(petab_dir, "v2")
  if (!dir.exists(v2_dir)) skip("PEtabTests/v2/ not present")
  if (!.libsbml_works()) skip("libsbml virtualenv not available")

  # 0016 and 0030 switch condition mid-run and resize a compartment from an
  # SBML event. 0023 is left out: it fires an event during preequilibration,
  # which is the steady state of the autonomous system here.
  for (case in c("0016", "0030")) {
    yamlPath <- file.path(v2_dir, case, paste0("_", case, ".yaml"))
    if (!file.exists(yamlPath)) next
    sol <- yaml::read_yaml(file.path(v2_dir, case, paste0("_", case, "_solution.yaml")))
    res <- tryCatch({
      pp <- .petab_case(paste0("v2_", case))
      pp$obj(pp$bestfit, deriv = FALSE)$value
    }, error = function(e) {
      message("v2 case ", case, " import error: ", conditionMessage(e))
      NA_real_
    })
    expect_equal(res, -2 * as.numeric(sol$llh),
                 tolerance = max(1e-4, abs(as.numeric(sol$tol_llh) / sol$llh)),
                 info = sprintf("v2 case %s", case))
  }
})



test_that("v2 priors add the truncated log density to the objective", {
  withr::local_dir(tempdir())
  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ not found -- set DMOD_PETABTESTS to the repo directory")
  yamlPath <- file.path(petab_dir, "v2", "0024", "_0024.yaml")
  if (!file.exists(yamlPath)) skip("v2 case 0024 not present")
  if (!.libsbml_works()) skip("libsbml virtualenv not available")
  sol <- yaml::read_yaml(file.path(petab_dir, "v2", "0024", "_0024_solution.yaml"))

  pp <- .petab_case("v2_0024")
  res <- pp$obj(pp$bestfit, deriv = FALSE)

  # The objective is -2 log posterior, so subtracting the prior part leaves
  # the likelihood. Only a declared distribution contributes: `p1` has
  # bounds alone, which dMod hands to the fit rather than to the objective.
  declared <- setdiff(names(sol$log_prior), "p1")
  expect_equal(unname(attr(res, "prior")),
               -2 * sum(unlist(sol$log_prior[declared])), tolerance = 1e-6)
  expect_equal(unname(res$value - attr(res, "prior")),
               -2 * as.numeric(sol$llh), tolerance = 1e-3)
})



test_that("v2 export round-trips a mid-run condition switch on a compartment", {
  withr::local_dir(tempdir())
  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ not found -- set DMOD_PETABTESTS to the repo directory")
  yamlPath <- file.path(petab_dir, "v2", "0030", "_0030.yaml")
  if (!file.exists(yamlPath)) skip("v2 case 0030 not present")
  if (!.libsbml_works()) skip("libsbml virtualenv not available")

  wd <- tempfile("v2_rt_0030_"); dir.create(wd)
  setwd(wd)
  first <- .petab_case("v2_0030")
  td <- file.path(wd, "export"); dir.create(td)
  exported <- exportPEtabObject(first, dir = td, formatVersion = "2.0.0",
                                overwrite = TRUE)

  # The compartment is an event target, so it imports as a state and has to
  # go back out with its size and its non-constant flag intact.
  sbml <- paste(readLines(list.files(td, pattern = "\\.xml$", full.names = TRUE)),
                collapse = "")
  expect_match(sbml, "<compartment id=\"C\"[^>]*size=\"4\"")
  expect_match(sbml, "<compartment id=\"C\"[^>]*constant=\"false\"")

  setwd(td)
  second <- importPEtab(exported, backend = "cppDE", compile = TRUE,
                        modelname = "v2rt0030B", cores = test_cores())
  withr::local_dir(tempdir())
  expect_equal(second$obj(second$bestfit, deriv = FALSE)$value,
               first$obj(first$bestfit, deriv = FALSE)$value,
               tolerance = 1e-6)
})



test_that("a v1 export keeps the preequilibration condition", {
  withr::local_dir(tempdir())
  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ not found -- set DMOD_PETABTESTS to the repo directory")
  if (!.libsbml_works()) skip("libsbml virtualenv not available")

  wd <- tempfile("v1_rt_0009_"); dir.create(wd)
  setwd(wd)
  first <- .petab_case("v1_0009")
  td <- file.path(wd, "export"); dir.create(td)
  exported <- exportPEtabObject(first, dir = td, formatVersion = "1",
                                overwrite = TRUE, modelID = "rt0009")

  # k1 is set per condition, so it belongs to conditions.tsv only.
  cond <- read.delim(file.path(td, "conditions_rt0009.tsv"),
                     stringsAsFactors = FALSE)
  expect_equal(setNames(cond$k1, cond$conditionId),
               c(c0 = 0.8, preeq_c0 = 0.3)[cond$conditionId])
  pars <- read.delim(file.path(td, "parameters_rt0009.tsv"),
                     stringsAsFactors = FALSE)
  expect_false("k1" %in% pars$parameterId)
  # A single sigma per observable stays a plain noise formula.
  obs <- read.delim(file.path(td, "observables_rt0009.tsv"),
                    stringsAsFactors = FALSE)
  expect_equal(as.character(obs$noiseFormula), "0.5")

  setwd(td)
  second <- importPEtab(exported, backend = "cppDE", modelname = "v1rt0009B",
                        cores = test_cores())
  withr::local_dir(tempdir())
  expect_equal(second$obj(second$bestfit, deriv = FALSE)$value,
               first$obj(first$bestfit, deriv = FALSE)$value,
               tolerance = 1e-6)
})



test_that("a v1 export refuses a condition switch during the simulation", {
  withr::local_dir(tempdir())
  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ not found -- set DMOD_PETABTESTS to the repo directory")
  if (!file.exists(file.path(petab_dir, "v2", "0030", "_0030.yaml")))
    skip("v2 case 0030 not present")
  if (!.libsbml_works()) skip("libsbml virtualenv not available")

  expect_error(
    exportPEtabObject(.petab_case("v2_0030"), dir = tempfile("v1_sw_"),
                      formatVersion = "1"),
    "cannot express condition changes")
})
