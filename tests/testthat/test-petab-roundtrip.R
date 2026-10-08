## PEtab importer / exporter: round trips of imported problems.

test_that("states a preequilibration cannot move keep their initial values", {

  withr::local_dir(tempdir())
  if (!.libsbml_works()) skip("libsbml virtualenv not available")

  # Isensee_JCB2018 pattern: the preequilibration keeps AC + pAC and
  # PDE + pPDE split as they start, while the steady-state equations alone
  # leave the split open.
  d <- file.path(tempdir(), "petab_invariant")
  dir.create(d, showWarnings = FALSE)
  est <- c(ks = 0.6, kd = 2, kf = 1.5, kr = 0.5, xi = 3)
  exportSbml(.petab_ac_module(),
             parameters = c(est, Fsk = 0, kp = 0, kdp = 0, kpp = 0, kpd = 0, cell = 1),
             inits = c(AC = 1, pAC = 0, ACF = 0, cAMP = 0.3, PDE = 1, pPDE = 0),
             filepath = file.path(d, "model.xml"), modelID = "invariant")
  tsv <- function(df, f)
    utils::write.table(df, file.path(d, f), sep = "\t", quote = FALSE, row.names = FALSE)
  tsv(data.frame(parameterId = names(est), parameterScale = "log10", lowerBound = 1e-3,
                 upperBound = 1e3, nominalValue = est, estimate = 1), "parameters.tsv")
  tsv(data.frame(conditionId = c("ctrl", "stim"), Fsk = c(0, 2)), "conditions.tsv")
  tsv(data.frame(observableId = c("obs_cAMP", "obs_AC"),
                 observableFormula = c("cAMP", "AC + ACF"), noiseFormula = 0.1),
      "observables.tsv")
  times <- c(0, 1, 5)
  tsv(data.frame(observableId = rep(c("obs_cAMP", "obs_AC"), each = 3),
                 preequilibrationConditionId = "ctrl", simulationConditionId = "stim",
                 time = times, measurement = 1), "measurements.tsv")
  writeLines(c("format_version: 1", "parameter_file: parameters.tsv", "problems:",
               "- condition_files:", "  - conditions.tsv",
               "  measurement_files:", "  - measurements.tsv",
               "  observable_files:", "  - observables.tsv",
               "  sbml_files:", "  - model.xml"), file.path(d, "problem.yaml"))

  pp <- importPEtab(file.path(d, "problem.yaml"), backend = "deSolve",
                    modelname = "petab_invariant",
                    options = list(atol = 1e-12, rtol = 1e-10))
  pred <- pp$prd(times, pp$bestfit, fixed = attr(pp, "petab_meta")$fixed,
                 deriv = FALSE)[[1]]

  ks <- est[["ks"]]; kd <- est[["kd"]]; xi <- est[["xi"]]
  a  <- 2 * est[["kf"]] / (2 * est[["kf"]] + est[["kr"]]); b <- 2 * est[["kf"]] + est[["kr"]]
  c_inf <- ks * (1 + (xi - 1) * a) / kd
  c_b   <- -ks * (xi - 1) * a / (kd - b)
  camp  <- c_inf + c_b * exp(-b * times) + (ks / kd - c_inf - c_b) * exp(-kd * times)
  expect_equal(unname(pred[, "obs_AC"]), rep(1, 3), tolerance = 1e-8)
  expect_equal(unname(pred[, "obs_cAMP"]), camp, tolerance = 1e-8)
})



test_that("exportPEtab keeps log parametrisations free of self-references", {

  withr::local_dir(tempdir())
  if (!.libsbml_works()) skip("libsbml virtualenv not available")

  reactions <- eqnlist() %>%
    addReaction("A", "B", "k1*A", "A to B") %>%
    addReaction("B", "A", "k2*B", "B to A")
  x <- Xs(odemodel(reactions, modelname = "selfref_ode", backend = "deSolve"),
          condition = "C1")
  obs <- eqnvec(y = "log(B)")
  g <- Y(obs, f = x, attachInput = TRUE, modelname = "selfref_obs",
         compile = TRUE)
  truth <- log(c(A = 1, B = 0.5, k1 = 0.25, k2 = 0.25))
  times <- c(0.5, 1, 2, 4)

  for (wrap in c("exp", "exp10")) {
    p <- eqnvec() %>%
      define("x~x", x = getParameters(g, x)) %>%
      insert(sprintf("x~%s(x)", wrap), x = .currentSymbols) %>%
      P(modelname = paste0("selfref_p_", wrap), condition = "C1",
        compile = TRUE)
    pars <- if (wrap == "exp") truth else truth / log(10)
    prd <- g * x * p
    pr <- prd(times, pars, deriv = FALSE)[["C1"]]
    data <- as.datalist(data.frame(
      name = "y", time = times, value = pr[match(times, pr[, "time"]), "y"],
      sigma = 0.1, condition = "C1"))
    obj_native <- normL2(data, prd)

    for (fv in c("1", "2.0.0")) {
      out_dir <- file.path(tempdir(), paste0("petab_selfref_", wrap, fv))
      yaml_out <- suppressWarnings(suppressMessages(exportPEtab(
        data, reactions, obs, p, pars, formatVersion = fv,
        dir = out_dir, overwrite = TRUE)))

      sbml <- readLines(file.path(out_dir, "dMod_export.xml"))
      expect_false(any(grepl("initialAssignment", sbml)))
      pt <- read.delim(file.path(out_dir, "parameters_dMod_export.tsv"),
                       stringsAsFactors = FALSE)
      nominal <- setNames(pt$nominalValue, pt$parameterId)
      expect_equal(nominal[c("init_A", "init_B", "k1", "k2")],
                   c(init_A = 1, init_B = 0.5, k1 = 0.25, k2 = 0.25))
      if (fv == "1")
        expect_true(all(pt$parameterScale ==
                        if (wrap == "exp") "log" else "log10"))

      cond <- read.delim(file.path(out_dir, "conditions_dMod_export.tsv"),
                         stringsAsFactors = FALSE)
      cond_map <- if (fv == "1") unlist(cond[1, c("A", "B")])
                  else setNames(cond$targetValue, cond$targetId)
      expect_equal(cond_map[c("A", "B")], c(A = "init_A", B = "init_B"))

      petab <- suppressWarnings(importPEtab(
        yaml_out, backend = "deSolve",
        modelname = paste0("selfref_imp_", wrap, substr(fv, 1, 1))))
      bf <- petab$bestfit
      v_petab <- petab$obj(bf, deriv = FALSE)$value
      expect_equal(v_petab, obj_native(pars, deriv = FALSE)$value,
                   tolerance = 1e-4)
    }
  }

  unlink("selfref_*"); unlink("*.c"); unlink("*.cpp")
  unlink("*.o"); unlink("*.so")
})



test_that("a v1 problem survives an export round trip in both formats", {
  if (!.libsbml_works()) skip("libsbml virtualenv not available")
  withr::local_dir(tempdir())

  wd <- tempfile("v1_rt_"); dir.create(wd)
  setwd(wd)
  first  <- .petab_boehm_cppDE()
  before <- first$obj(first$bestfit)$value

  # v1 keeps `parameterScale`, v2 has no such column and takes linear values,
  # so the exporter has to invert the scale for v2 and only for v2.
  versions <- c("1", "2.0.0")
  second <- lapply(versions, function(version) {
    td <- file.path(wd, paste0("export_", sub("\\.", "", version)))
    dir.create(td)
    exported <- suppressWarnings(
      exportPEtabObject(first, dir = td, formatVersion = version,
                        overwrite = TRUE))
    importPEtab(exported, backend = "cppDE", compile = FALSE,
                modelname = paste0("v1rtB", sub("\\.", "", version)))
  })
  # Both are the same model and compile under the same KLU macros.
  .petab_compile(second, "v1rtB")
  for (i in seq_along(versions))
    expect_equal(second[[i]]$obj(second[[i]]$bestfit)$value, before,
                 tolerance = 1e-4, info = paste("formatVersion", versions[i]))
  withr::local_dir(tempdir())
})



test_that("a v1 conditionName column is metadata, not a condition target", {
  if (!.libsbml_works()) skip("libsbml virtualenv not available")
  withr::local_dir(tempdir())

  wd <- tempfile("v1_condname_"); dir.create(wd)
  setwd(wd)
  petab <- .petab_boehm_cppDE()
  td <- file.path(wd, "export"); dir.create(td)
  suppressWarnings(exportPEtabObject(petab, dir = td, formatVersion = "2.0.0",
                                     overwrite = TRUE))
  cond <- utils::read.delim(
    list.files(td, pattern = "^conditions_.*\\.tsv$", full.names = TRUE),
    stringsAsFactors = FALSE)
  withr::local_dir(tempdir())

  expect_false("conditionName" %in% cond$targetId)
})



test_that("a v1 export names conditions by id in every table", {
  if (!.libsbml_works()) skip("libsbml virtualenv not available")
  withr::local_dir(tempdir())

  # Boehm keys its conditions by name; the tables on disk store the id.
  wd <- tempfile("v1_condid_"); dir.create(wd)
  setwd(wd)
  petab <- .petab_boehm_cppDE()
  td <- file.path(wd, "export"); dir.create(td)
  exportPEtabObject(petab, dir = td, formatVersion = "1", overwrite = TRUE,
                    modelID = "condid")
  cond <- utils::read.delim(file.path(td, "conditions_condid.tsv"),
                            stringsAsFactors = FALSE)
  meas <- utils::read.delim(file.path(td, "measurements_condid.tsv"),
                            stringsAsFactors = FALSE)
  withr::local_dir(tempdir())

  expect_true(all(meas$simulationConditionId %in% cond$conditionId))
})



test_that("importPEtab builds a reverse sweep through every piece", {
  withr::local_dir(tempdir())
  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ not found -- set DMOD_PETABTESTS to the repo directory")
  if (!.libsbml_works()) skip("libsbml virtualenv not available")

  wd <- tempfile("petab_rev_"); dir.create(wd)
  setwd(wd)
  pp <- importPEtab(file.path(petab_dir, "0001", "_0001.yaml"), backend = "cppDE",
                    derivMode = c("forward", "reverse"), modelname = "petab_rev",
                    cores = test_cores(), options = list(atol = 1e-12, rtol = 1e-10),
                    optionsSens = list(atol = 1e-10, rtol = 1e-8))
  fwd <- pp$obj(pp$bestfit)
  rev <- pp$obj(pp$bestfit, sweep = "reverse")
  withr::local_dir(tempdir())
  expect_equal(rev$value, fwd$value, tolerance = 1e-6)
  expect_equal(rev$gradient[names(fwd$gradient)], fwd$gradient, tolerance = 1e-4)
})

test_that("importPEtab(sparse =) pins the linear solver of the model", {
  if (!.libsbml_works()) skip("libsbml virtualenv not available")
  withr::local_dir(tempdir())
  auto <- .petab_boehm_cppDE()
  d <- file.path(tempdir(), "dmod_petab_boehm_dense")
  dir.create(d, showWarnings = FALSE)
  dense <- importPEtab(.petab_boehm_yaml(), backend = "cppDE", modelname = "v1dense",
                       sparse = FALSE, cores = test_cores(), outdir = d)
  expect_false(isTRUE(attr(dense$odemodel$extended, "sparse")))
  a <- auto$obj(auto$bestfit)
  b <- dense$obj(dense$bestfit)
  expect_equal(b$value, a$value, tolerance = 1e-6)
  expect_equal(b$gradient, a$gradient, tolerance = 1e-4)
})
