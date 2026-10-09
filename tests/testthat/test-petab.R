## PEtab importer / exporter: parsers, normalisers and native exports.

## --- interpreter resolution -----------------------------------------------

test_that(".dmod_libsbml_python resolves a usable interpreter", {

  # Deliberately not behind .libsbml_works(): that helper turns every failure
  # into FALSE, so a resolver returning an empty path made the libsbml tests
  # skip instead of fail. Nothing starts Python before the call either: the bug
  # lives in the uninitialised state, so a guard that initialises first would
  # let the old code pass here.
  skip_if_not_installed("reticulate")
  withr::local_envvar(c(DMOD_LIBSBML_PYTHON = NA, DMOD_LIBSBML_OK = NA))

  py <- dMod2:::.dmod_libsbml_python()
  expect_true(nzchar(py))
  expect_true(file.exists(py))
})



## --- pure parser unit tests (no SBML) -------------------------------------

test_that(".petab_parse_parameters splits estimated / fixed and tracks scales", {

  # PEtab v1: nominalValue / lowerBound / upperBound are written on the
  # linear scale regardless of parameterScale. The parser pre-transforms
  # estimated parameters and bounds to the parameter scale (dMod's pouter
  # convention); fixed parameters stay on the linear scale because the
  # trafo's scale chain rule only wraps estimated outer parameters.
  df <- data.frame(
    parameterId    = c("a", "b", "c"),
    parameterScale = c("lin", "log10", "log"),
    lowerBound     = c(0, 1e-3, 1e-5),
    upperBound     = c(10, 1e3, 1e5),
    nominalValue   = c(1.0, 100, exp(2)),
    estimate       = c(1L, 1L, 0L),
    stringsAsFactors = FALSE
  )
  pm <- dMod2:::.petab_parse_parameters(df)

  expect_equal(names(pm$pouter), c("a", "b"))
  # a (lin)   = 1.0
  # b (log10) = log10(100) = 2  -- pouter on parameter scale
  expect_equal(unname(pm$pouter), c(1.0, 2.0))
  expect_equal(names(pm$fixed),  c("c"))
  # c is fixed → stays on linear scale (no scale chain rule wraps it).
  expect_equal(unname(pm$fixed["c"]), exp(2))
  expect_equal(pm$scales[["a"]], "lin")
  expect_equal(pm$scales[["b"]], "log10")
  expect_equal(pm$scales[["c"]], "log")
  # lower["b"] = log10(1e-3) = -3
  expect_equal(unname(pm$lower["b"]), -3)
  expect_equal(unname(pm$upper["b"]), 3)
})



test_that(".petab_parse_observables defaults to lin/normal and parses noise", {

  df <- data.frame(
    observableId      = c("o1", "o2"),
    observableFormula = c("A", "B + offset"),
    noiseFormula      = c("0.5", "1"),
    stringsAsFactors  = FALSE
  )
  om <- dMod2:::.petab_parse_observables(df)
  expect_equal(unname(om$obs_trafo),  c("lin", "lin"))
  expect_equal(unname(om$noise_dist), c("normal", "normal"))
  expect_equal(unname(om$noise),      c("0.5", "1"))
})



test_that(".petab_parse_conditions classifies columns as init / parameter", {

  # Case 0002 shape: a0 is in conditions and is a parameter symbol that also
  # parameterises species A's initial. We expect "parameter" classification.
  df <- data.frame(conditionId = c("c0", "c1"),
                   a0          = c(0.8, 0.9),
                   stringsAsFactors = FALSE)
  ci <- dMod2:::.petab_parse_conditions(df,
          sbml_states       = c("A", "B"),
          sbml_compartments = "compartment",
          sbml_pars         = c("a0", "b0", "k1", "k2", "compartment"))
  expect_equal(ci$col_kind[["a0"]], "parameter")
  expect_equal(ci$override_cols, "a0")

  # init kind: column name is a state itself
  df2 <- data.frame(conditionId = "c0", A = 0.5,
                    stringsAsFactors = FALSE)
  ci2 <- dMod2:::.petab_parse_conditions(df2,
           sbml_states = c("A", "B"),
           sbml_compartments = character(),
           sbml_pars = character())
  expect_equal(ci2$col_kind[["A"]], "init")
})



test_that(".petab_parse_measurements unfolds per-row observableParameters", {

  obs_meta <- list(
    obs   = c(obs_a = "observableParameter1_obs_a * A"),
    noise = c(obs_a = "1"),
    obs_trafo  = c(obs_a = "lin"),
    noise_dist = c(obs_a = "normal")
  )

  # Case 0006 shape: same simulation condition, two different obs param
  # values across two rows.
  df <- data.frame(
    observableId          = c("obs_a", "obs_a"),
    simulationConditionId = c("c0", "c0"),
    time                  = c(0, 10),
    measurement           = c(0.7, 0.1),
    observableParameters  = c("10", "15"),
    stringsAsFactors = FALSE
  )
  mi <- dMod2:::.petab_parse_measurements(df, obs_meta)
  expect_equal(nrow(mi$sub_cond_map), 2L)
  expect_true(all(grepl("^c0__", mi$sub_cond_map$sub_condition)))
  # data is partitioned across sub-conditions:
  expect_equal(sort(unique(mi$data$condition)),
               sort(mi$sub_cond_map$sub_condition))

  # Case 0001 shape: no sub-condition splitting, single condition.
  df2 <- data.frame(
    observableId          = c("obs_a", "obs_a"),
    simulationConditionId = c("c0", "c0"),
    time                  = c(0, 10),
    measurement           = c(0.7, 0.1),
    stringsAsFactors = FALSE
  )
  mi2 <- dMod2:::.petab_parse_measurements(df2, obs_meta)
  expect_equal(nrow(mi2$sub_cond_map), 1L)
  expect_equal(mi2$sub_cond_map$sub_condition, "c0")
})



test_that("readPEtabYaml resolves manifest paths correctly", {

  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ not found -- set DMOD_PETABTESTS to the repo directory")

  y <- readPEtabYaml(file.path(petab_dir, "0001", "_0001.yaml"))
  expect_equal(y$formatVersion, 1L)
  expect_true(file.exists(y$problems[[1]]$sbmlFile))
  expect_true(file.exists(y$problems[[1]]$measurementFile))
})



test_that("readPEtabTables returns the expected slots for v1", {
  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ not found -- set DMOD_PETABTESTS to the repo directory")

  tabs <- readPEtabTables(file.path(petab_dir, "0001", "_0001.yaml"))
  expect_named(tabs, c("parameters", "conditions", "measurements",
                       "observables", "experiments", "mapping",
                       "sbmlPath", "sbmlPaths", "formatVersion"))
  expect_s3_class(tabs$parameters,   "data.frame")
  expect_s3_class(tabs$conditions,   "data.frame")
  expect_s3_class(tabs$measurements, "data.frame")
  expect_s3_class(tabs$observables,  "data.frame")
  expect_null(tabs$experiments)   # v1 has no experiments table
  expect_null(tabs$mapping)
  expect_identical(tabs$formatVersion, 1L)
})

test_that("readPEtabYaml and readPEtabTables read the bundled Boehm problem", {
  skip_if_not_installed("yaml")
  yaml <- system.file("extdata/petab_boehm/Boehm.yaml", package = "dMod2")

  y <- readPEtabYaml(yaml)
  expect_identical(y$formatVersion, 1L)
  expect_identical(normalizePath(y$baseDir), normalizePath(dirname(yaml)))
  pr <- y$problems[[1]]
  files <- unlist(pr[c("sbmlFile", "conditionFile", "measurementFile", "observableFile")])
  expect_true(all(file.exists(c(y$parameterFile, files))))
  expect_null(pr$experimentFile)

  tabs <- readPEtabTables(yaml)
  expect_identical(tabs$sbmlPath, pr$sbmlFile)
  expect_equal(nrow(tabs$parameters), 11L)
  expect_equal(sum(tabs$parameters$estimate), 9)
  expect_equal(nrow(tabs$conditions), 1L)
  expect_equal(nrow(tabs$measurements), 48L)
  expect_identical(tabs$observables$observableId,
                   c("pSTAT5A_rel", "pSTAT5B_rel", "rSTAT5A_rel"))
  expect_true(all(tabs$measurements$observableId %in% tabs$observables$observableId))
})



## --- end-to-end fixture test (no SBML import required) -------------------
##
## We hand-build the eqnlist that matches PEtab test case 0001's SBML model
## and verify the trafo+objective machinery against the published solution.
## This avoids a libsbml dependency on every test run.

test_that("hand-built case-0001 fixture produces solution-matching llh", {

  # Reaction network identical to PEtabTests/0001/_model.xml after libsbml
  # would have inlined the kinetic law's compartment factor, i.e. with a unit
  # compartment.
  nat <- .petab_native()
  x <- nat$x_ab
  g <- nat$g_0
  p <- nat$p_0

  # Data exactly matching _measurements.tsv.
  data <- as.datalist(data.frame(
    name      = c("obs_a", "obs_a"),
    time      = c(0, 10),
    value     = c(0.7, 0.1),
    sigma     = c(0.5, 0.5),
    condition = c("c0", "c0"),
    stringsAsFactors = FALSE
  ))

  prd <- g * x * p
  obj <- normL2(data, prd)

  # Nominal pouter from _parameters.tsv
  pouter <- c(a0 = 1.0, b0 = 0.0, k1 = 0.8, k2 = 0.6)

  out <- obj(pouter, deriv = FALSE)

  # PEtab _0001_solution.yaml gives llh = -0.8475016971318833 and
  # chi2 = 0.7918379836848569. dMod's normL2 returns the *full* Gaussian
  # negative log-likelihood multiplied by 2 (i.e. -2*log L), which equals
  # chi2 + sum(log(2*pi*sigma^2)) per data point. We compare against -2*llh.
  expect_lt(abs(out$value - (-2 * -0.8475016971318833)), 0.001)
})



test_that("exportSbml emits InitialAssignment for symbolic state initials", {

  withr::local_dir(tempdir())
  if (!.libsbml_works()) skip("libsbml virtualenv not available")

  reactions <- eqnlist()
  reactions <- addReaction(reactions, "A", "B", "k1*A", "fwd",
                           compartment = "compartment")
  reactions <- addReaction(reactions, "B", "A", "k2*B", "rev",
                           compartment = "compartment")

  # Mixed inits: A is symbolic (→ InitialAssignment), B is numeric.
  inits <- c(A = "a0", B = "0")
  pars  <- c(a0 = 0.8, k1 = 0.8, k2 = 0.6, compartment = 1.0)

  out_xml <- file.path(tempdir(), "ia_export.xml")
  exportSbml(reactions, parameters = pars, inits = inits,
              filepath = out_xml, modelID = "ia_export")

  xml_text <- readLines(out_xml, warn = FALSE)
  expect_true(any(grepl("<initialAssignment", xml_text, fixed = TRUE)),
              info = "no <initialAssignment> emitted for symbolic init")
  expect_true(any(grepl("symbol=\"A\"", xml_text)),
              info = "InitialAssignment for A missing")

  unlink(out_xml)
})



## --- trafo-aware exportPEtab: pure-R helper unit tests --------------------
##
## The strip + classify decomposer should be unit-testable without libsbml
## because it operates only on character RHSes and named eqnvecs.

test_that(".petab_invariant_states finds the states an equilibration cannot move", {
  el   <- .petab_ac_module()
  init <- c(AC = "1", pAC = "0", ACF = "0.0", cAMP = "0.3", PDE = "1", pPDE = "0")
  inv  <- dMod2:::.petab_invariant_states(el, init, c("Fsk", "kp", "kdp", "kpp", "kpd"))
  expect_setequal(inv$zero, c("pAC", "ACF", "pPDE"))
  expect_setequal(inv$frozen, c("AC", "PDE"))
  expect_equal(inv$idle, c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE))

  red <- dMod2:::.petab_reduce_network(el, c(inv$zero, inv$frozen), inv$idle,
                                       c("Fsk", "kp", "kdp", "kpp", "kpd", inv$zero))
  expect_equal(red$states, "cAMP")
  expect_setequal(setdiff(getSymbols(red$rates), "cAMP"), c("ks", "AC", "xi", "kd", "PDE"))

  # With the input on, ACF is fed and AC moves; PDE still cannot.
  inv <- dMod2:::.petab_invariant_states(el, init, c("kp", "kdp", "kpp", "kpd"))
  expect_setequal(inv$zero, c("pAC", "pPDE"))
  expect_setequal(inv$frozen, "PDE")
  # A state that starts above 0 is never pinned at 0.
  inv <- dMod2:::.petab_invariant_states(el, replace(init, "pAC", "0.2"),
                                         c("Fsk", "kp", "kdp", "kpp", "kpd"))
  expect_setequal(inv$zero, c("ACF", "pPDE"))
  expect_setequal(inv$frozen, c("AC", "pAC", "PDE"))
})


test_that(".petab_strip_param_scale compensates the chain rule per-occurrence", {
  # Clean wrap stays clean (importer chain rule re-wraps it)
  expect_equal(
    dMod2:::.petab_strip_param_scale("10^(K_REFLUX)", c(K_REFLUX = "log10")),
    "K_REFLUX")
  # Steady-state-like product of clean wraps
  expect_equal(
    dMod2:::.petab_strip_param_scale(
      "10^(TCA_CELL) * 10^(K_EXPORT_CANA) / 10^(K_REFLUX)",
      c(TCA_CELL = "log10", K_EXPORT_CANA = "log10", K_REFLUX = "log10")),
    "TCA_CELL * K_EXPORT_CANA/K_REFLUX")
  # Compound expression: bare KM inside is compensated with log10(KM)
  expect_equal(
    dMod2:::.petab_strip_param_scale("10^(KM + 5)", c(KM = "log10")),
    "10^(log10(KM) + 5)")
  # Mixed wrap: clean wrap stripped, bare occurrence compensated
  expect_equal(
    dMod2:::.petab_strip_param_scale("10^(K) + K", c(K = "log10")),
    "K + log10(K)")
  # exp(.) for log scale, with compensation
  expect_equal(
    dMod2:::.petab_strip_param_scale("exp(K)", c(K = "log")),
    "K")
  expect_equal(
    dMod2:::.petab_strip_param_scale("exp(KM + 5)", c(KM = "log")),
    "exp(log(KM) + 5)")
  # lin parameters unchanged
  expect_equal(
    dMod2:::.petab_strip_param_scale("K1 * K2 + offset",
      c(K1 = "log10", K2 = "log10", offset = "lin")),
    "log10(K1) * log10(K2) + offset")
  # Pure numeric literal -- passes through
  expect_equal(
    dMod2:::.petab_strip_param_scale("0", c()), "0")
})



test_that(".petab_classify_lhs categorizes per-condition RHSes", {
  conds <- c("c1", "c2")
  stripped <- list(
    c1 = c(s = "1", k = "K", state = "0", iden = "iden", v = "K"),
    c2 = c(s = "1", k = "K", state = "0", iden = "iden", v = "K2"))
  # all-numeric constant
  expect_equal(dMod2:::.petab_classify_lhs(stripped, "s", conds),
               list(kind = "const_numeric", value = 1))
  # all-symbolic constant
  expect_equal(dMod2:::.petab_classify_lhs(stripped, "k", conds),
               list(kind = "const_symbolic", formula = "K"))
  # numeric-zero (state init)
  expect_equal(dMod2:::.petab_classify_lhs(stripped, "state", conds),
               list(kind = "const_numeric", value = 0))
  # identity (RHS == LHS)
  expect_equal(dMod2:::.petab_classify_lhs(stripped, "iden", conds),
               list(kind = "identity"))
  # varying
  res <- dMod2:::.petab_classify_lhs(stripped, "v", conds)
  expect_equal(res$kind, "varying")
  expect_equal(res$per_cond, c(c1 = "K", c2 = "K2"))
  # missing
  expect_equal(dMod2:::.petab_classify_lhs(stripped, "absent", conds),
               list(kind = "missing"))
})



test_that(".petab_classify_lhs collapses 10^0 -> 1 via eval_constant", {
  conds <- c("c1", "c2")
  stripped <- list(c1 = c(s = "10^0"), c2 = c(s = "10^0"))
  res <- dMod2:::.petab_classify_lhs(stripped, "s", conds)
  expect_equal(res, list(kind = "const_numeric", value = 1))
})



## --- trafo-aware exportPEtab: native roundtrip ----------------------------


test_that("native exportPEtab roundtrips outer pouter on log10 scale (1-cond)", {

  withr::local_dir(tempdir())
  if (!.libsbml_works()) skip("libsbml virtualenv not available")

  # Tiny 2-state model, single condition.
  nat <- .petab_native()
  reactions <- nat$ab
  obs       <- nat$obs_ab
  x_native  <- nat$x_ab
  g_native  <- nat$g_ab
  p_native  <- nat$p_rt1

  data <- as.datalist(data.frame(
    name = c("obs_a", "obs_b"), time = c(1, 1),
    value = c(0.5, 0.3), sigma = c(1, 1),
    condition = c("c1", "c1"), stringsAsFactors = FALSE))

  pouter <- c(A = -1, B = -1, K1 = -1, K2 = -1)
  obj_native <- normL2(data, g_native * x_native * p_native)

  out_dir <- file.path(tempdir(), "petab_rt1")
  yaml_out <- exportPEtab(
    data = data, reactions = reactions, observables = obs,
    p = p_native, pouter = pouter,
    parameterScale = "log10", modelID = "rt1_export",
    formatVersion = "1", dir = out_dir, overwrite = TRUE)

  petab <- importPEtab(yaml_out, backend = "deSolve",
                       modelname = "rt1_imp", cores = test_cores())

  expect_setequal(names(petab$bestfit), names(.init_ids(pouter)))
  expect_true(all(attr(petab$bestfit, "petab_scales") == "log10"))

  v_native <- obj_native(pouter, deriv = FALSE)$value
  v_petab  <- petab$obj(.init_ids(pouter)[names(petab$bestfit)], deriv = FALSE)$value
  expect_lt(abs(v_native - v_petab), 1e-3)

  unlink("rt1_*"); unlink("*.c"); unlink("*.cpp")
  unlink("*.o"); unlink("*.so")
})



test_that("native exportPEtab roundtrips per-condition k override (2-cond)", {

  withr::local_dir(tempdir())
  if (!.libsbml_works()) skip("libsbml virtualenv not available")

  nat <- .petab_native()
  reactions <- nat$a
  obs       <- nat$obs_a
  x_native  <- nat$x_a
  g_native  <- nat$g_a
  p_native  <- nat$p_rt2

  data <- as.datalist(data.frame(
    name = c("obs_a", "obs_a"), time = c(1, 1),
    value = c(0.5, 0.5), sigma = c(1, 1),
    condition = c("closed", "open"), stringsAsFactors = FALSE))

  pouter <- c(A = -1, B = -1, K = -1, K_OPEN = -0.5)
  obj_native <- normL2(data, g_native * x_native * p_native)

  out_dir <- file.path(tempdir(), "petab_rt2")
  yaml_out <- exportPEtab(
    data = data, reactions = reactions, observables = obs,
    p = p_native, pouter = pouter,
    parameterScale = "log10", modelID = "rt2_export",
    formatVersion = "1", dir = out_dir, overwrite = TRUE)

  # conditions.tsv must have a `k` column distinguishing closed from open.
  cond_df <- read.delim(file.path(out_dir, "conditions_rt2_export.tsv"),
                        stringsAsFactors = FALSE)
  expect_true("k" %in% colnames(cond_df))
  expect_setequal(cond_df$k, c("K", "K_OPEN"))

  petab <- importPEtab(yaml_out, backend = "deSolve",
                       modelname = "rt2_imp", cores = test_cores())
  expect_setequal(names(petab$bestfit), c("init_A", "init_B", "K", "K_OPEN"))

  v_native <- obj_native(pouter, deriv = FALSE)$value
  v_petab  <- petab$obj(.init_ids(pouter)[names(petab$bestfit)], deriv = FALSE)$value
  expect_lt(abs(v_native - v_petab), 1e-3)

  unlink("rt2_*"); unlink("*.c"); unlink("*.cpp")
  unlink("*.o"); unlink("*.so")
})



test_that("exportPEtab errors on undeclared free symbol after strip", {

  withr::local_dir(tempdir())
  if (!.libsbml_works()) skip("libsbml virtualenv not available")

  reactions <- eqnlist() %>%
    addReaction("A", "B", rate = "k*A", description = "fwd")
  m <- odemodel(reactions, modelname = "err1_ode", compile = FALSE,
                backend = "deSolve")
  x <- Xs(m)
  obs <- eqnvec(obs_a = "A")
  g <- Y(obs, f = x, compile = FALSE, modelname = "err1_obs",
         attachInput = FALSE)
  trafo <- as.eqnvec(c(A = "10^(A)", B = "10^(B)",
                       k = "10^(K) + UNDECLARED"))
  p <- P(trafo, condition = "c1", compile = FALSE, modelname = "err1_par")

  data <- as.datalist(data.frame(
    name = "obs_a", time = 1, value = 0.5, sigma = 1, condition = "c1",
    stringsAsFactors = FALSE))

  # exportPEtab emits an informational warning when parameterScale is supplied
  # to a v2 export (it is ignored on disk). The test only cares about the error.
  suppressWarnings(expect_error(
    exportPEtab(data = data, reactions = reactions, observables = obs,
                p = p, pouter = c(A = 0, B = 0, K = -1),
                parameterScale = "log10",
                dir = tempfile("err1_"), overwrite = TRUE),
    "undeclared symbol"))

  unlink("err1_*"); unlink("*.c"); unlink("*.cpp")
  unlink("*.o"); unlink("*.so")
})



test_that("native exportPEtab roundtrips per-row sigma via noiseParameters column", {

  withr::local_dir(tempdir())
  if (!.libsbml_works()) skip("libsbml virtualenv not available")

  nat <- .petab_native()
  reactions <- nat$a
  obs <- nat$obs_a
  x   <- nat$x_a
  g   <- nat$g_a
  p   <- nat$p_sig

  # Three measurements with three different sigmas -- exercise the
  # noiseParameter1_<obsId> placeholder + per-row noiseParameters path.
  data <- as.datalist(data.frame(
    name = "obs_a", time = c(1, 2, 3),
    value = c(0.5, 0.3, 0.2), sigma = c(0.5, 1.0, 2.0),
    condition = "c1", stringsAsFactors = FALSE))

  pouter <- c(A = 0, B = 0, K = -1)
  obj_native <- normL2(data, g * x * p)

  out_dir <- file.path(tempdir(), "petab_rt_sig")
  yaml_out <- exportPEtab(
    data = data, reactions = reactions, observables = obs,
    p = p, pouter = pouter,
    parameterScale = "log10", modelID = "rt_sig_export",
    formatVersion = "1", dir = out_dir, overwrite = TRUE)

  # observables.tsv must declare the placeholder noiseFormula.
  obs_tsv <- read.delim(file.path(out_dir, "observables_rt_sig_export.tsv"),
                        stringsAsFactors = FALSE)
  expect_equal(obs_tsv$noiseFormula, "noiseParameter1_obs_a")

  # measurements.tsv must contain per-row noiseParameters values.
  meas_tsv <- read.delim(file.path(out_dir, "measurements_rt_sig_export.tsv"),
                         stringsAsFactors = FALSE)
  expect_setequal(as.numeric(meas_tsv$noiseParameters), c(0.5, 1.0, 2.0))

  petab <- importPEtab(yaml_out, backend = "deSolve",
                       modelname = "rt_sig_imp", cores = test_cores())
  v_native <- obj_native(pouter, deriv = FALSE)$value
  v_petab  <- petab$obj(.init_ids(pouter)[names(petab$bestfit)], deriv = FALSE)$value
  expect_lt(abs(v_native - v_petab), 1e-3)

  unlink("rt_sig_*"); unlink("*.c"); unlink("*.cpp")
  unlink("*.o"); unlink("*.so")
})



test_that("native exportPEtab roundtrips compound trafos like 10^(KM + 5)", {

  withr::local_dir(tempdir())
  if (!.libsbml_works()) skip("libsbml virtualenv not available")

  # 1-state, 1-reaction with a non-trivial compound mapping for the rate:
  #   k = 10^(K + 5) -- chain-rule "compensation" path, not strippable.
  nat <- .petab_native()
  reactions <- nat$a
  obs <- nat$obs_a
  x   <- nat$x_a
  g   <- nat$g_a
  p   <- nat$p_rt3

  data <- as.datalist(data.frame(
    name = "obs_a", time = 1, value = 0.5, sigma = 1, condition = "c1",
    stringsAsFactors = FALSE))

  pouter <- c(A = 0, B = 0, K = -1)
  obj_native <- normL2(data, g * x * p)

  out_dir <- file.path(tempdir(), "petab_rt3")
  # v1 has no expressions in conditions.tsv; dMod reads them back.
  expect_warning(yaml_out <- exportPEtab(
    data = data, reactions = reactions, observables = obs,
    p = p, pouter = pouter,
    parameterScale = "log10", modelID = "rt3_export",
    formatVersion = "1", dir = out_dir, overwrite = TRUE),
    "only numbers and parameter ids")

  # The conditions.tsv cell for k must contain the compensated form
  # `10^(log10(K) + 5)` so the importer's chain rule reproduces 10^(K+5).
  cond_df <- read.delim(file.path(out_dir, "conditions_rt3_export.tsv"),
                        stringsAsFactors = FALSE)
  expect_match(as.character(cond_df$k[[1L]]), "log10\\(K\\)")

  petab <- importPEtab(yaml_out, backend = "deSolve",
                       modelname = "rt3_imp", cores = test_cores())
  v_native <- obj_native(pouter, deriv = FALSE)$value
  v_petab  <- petab$obj(.init_ids(pouter)[names(petab$bestfit)], deriv = FALSE)$value
  expect_lt(abs(v_native - v_petab), 1e-3)

  unlink("rt3_*"); unlink("*.c"); unlink("*.cpp")
  unlink("*.o"); unlink("*.so")
})



test_that("exportPEtab rejects self-referencing assignments", {
  expect_error(
    dMod2:::.petab_check_self_refs(list(A = "exp(A)"),
                                   data.frame(conditionId = "C1")),
    "initialAssignment A = exp\\(A\\)")
  expect_error(
    dMod2:::.petab_check_self_refs(list(),
                                   data.frame(conditionId = "C1",
                                              k1 = "exp(k1)")),
    "condition k1 = exp\\(k1\\)")
  expect_true(dMod2:::.petab_check_self_refs(
    list(A = "init_A"), data.frame(conditionId = "C1", k1 = "exp(logk1)")))
})



## --- PEtab v2 (no-SBML pure-parser tests) ---------------------------------

test_that(".petab_major_version recognises v1 and v2 strings", {
  expect_identical(dMod2:::.petab_major_version(1L),       1L)
  expect_identical(dMod2:::.petab_major_version("1"),      1L)
  expect_identical(dMod2:::.petab_major_version("1.0.0"),  1L)
  expect_identical(dMod2:::.petab_major_version("2.0.0"),  2L)
  expect_identical(dMod2:::.petab_major_version("2.1.3"),  2L)
  expect_identical(dMod2:::.petab_major_version(NULL),     1L)  # legacy default
  expect_error(dMod2:::.petab_major_version("v2"),
               regexp = "Unrecognised PEtab format_version")
})



test_that(".petab_v2_normalize_tables converts a single-condition v2 problem", {
  tables <- list(
    parameters = data.frame(
      parameterId  = c("k1", "k2", "init_a"),
      lowerBound   = c(1e-5, 1e-5, 0),
      upperBound   = c(1e3,  1e3,  10),
      nominalValue = c(0.1, 0.5, 1.0),
      estimate     = c("true", "true", "false"),
      stringsAsFactors = FALSE),
    observables = data.frame(
      observableId           = c("o1"),
      observableFormula      = c("A * scale + offset"),
      observablePlaceholders = c("scale;offset"),
      noiseFormula           = c("sigma"),
      noiseDistribution      = c("log-normal"),
      noisePlaceholders      = c("sigma"),
      stringsAsFactors = FALSE),
    conditions = data.frame(
      conditionId = c("c1", "c1"),
      targetId    = c("a0", "k_in"),
      targetValue = c("init_a", "0.4"),
      stringsAsFactors = FALSE),
    measurements = data.frame(
      observableId = c("o1", "o1"),
      experimentId = c("exp1", "exp1"),
      time         = c(0, 10),
      measurement  = c(1.0, 0.6),
      observableParameters = c("1.5;0", "1.5;0"),
      noiseParameters      = c("0.1", "0.1"),
      stringsAsFactors = FALSE),
    experiments = data.frame(
      experimentId = c("exp1"),
      time         = c("0"),
      conditionId  = c("c1"),
      stringsAsFactors = FALSE),
    mapping     = NULL,
    sbmlPath   = "ignored.xml",
    formatVersion = 2L)

  out <- dMod2:::.petab_v2_normalize_tables(tables)

  # parameters: parameterScale synthesised; estimate coerced to 1/0.
  expect_true("parameterScale" %in% colnames(out$parameters))
  expect_equal(unique(out$parameters$parameterScale), "lin")
  expect_equal(out$parameters$estimate, c(1L, 1L, 0L))

  # observables: log-normal split into log + normal; placeholders rewritten.
  expect_equal(out$observables$observableTransformation, "log")
  expect_equal(out$observables$noiseDistribution, "normal")
  expect_match(out$observables$observableFormula,
               "observableParameter1_o1.*observableParameter2_o1")
  expect_match(out$observables$noiseFormula,
               "^noiseParameter1_o1$")

  # conditions: long → wide.
  expect_equal(sort(setdiff(colnames(out$conditions), "conditionId")),
               c("a0", "k_in"))
  r <- which(out$conditions$conditionId == "c1")
  expect_equal(out$conditions$a0[r],   "init_a")
  expect_equal(out$conditions$k_in[r], "0.4")

  # measurements: experimentId rewritten.
  expect_equal(out$measurements$simulationConditionId,
               c("c1", "c1"))
  expect_equal(out$measurements$preequilibrationConditionId,
               c("", ""))
  expect_false("experimentId" %in% colnames(out$measurements))
})



test_that(".petab_v2_normalize_tables handles preequilibration via 2-period experiments", {
  tables <- list(
    parameters = data.frame(parameterId = "k", lowerBound = 0, upperBound = 1,
                            nominalValue = 0.5, estimate = "true",
                            stringsAsFactors = FALSE),
    observables = data.frame(observableId = "o1", observableFormula = "A",
                             noiseFormula = "1",
                             noiseDistribution = "normal",
                             stringsAsFactors = FALSE),
    conditions = data.frame(
      conditionId = c("c_pre", "c_sim"),
      targetId    = c("a0", "a0"),
      targetValue = c("5",  "1"),
      stringsAsFactors = FALSE),
    measurements = data.frame(
      observableId = "o1", experimentId = "exp_with_pre",
      time = 5, measurement = 0.7,
      stringsAsFactors = FALSE),
    experiments = data.frame(
      experimentId = c("exp_with_pre", "exp_with_pre"),
      time         = c("-inf",         "0"),
      conditionId  = c("c_pre",        "c_sim"),
      stringsAsFactors = FALSE),
    mapping     = NULL,
    sbmlPath   = "ignored.xml",
    formatVersion = 2L)

  out <- dMod2:::.petab_v2_normalize_tables(tables)
  expect_equal(out$measurements$simulationConditionId,        "c_sim")
  expect_equal(out$measurements$preequilibrationConditionId,  "c_pre")
})



test_that(".petab_v2_normalize_tables turns later periods into switches", {
  tables <- list(
    parameters = data.frame(parameterId = "k", lowerBound = 0, upperBound = 1,
                            nominalValue = 0.5, estimate = "true",
                            stringsAsFactors = FALSE),
    observables = data.frame(observableId = "o1", observableFormula = "A",
                             noiseFormula = "1",
                             noiseDistribution = "normal",
                             stringsAsFactors = FALSE),
    conditions = data.frame(conditionId = c("c1", "c2", "c3"),
                            targetId = c("a0", "a0", "a0"),
                            targetValue = c("1", "2", "3"),
                            stringsAsFactors = FALSE),
    measurements = data.frame(observableId = "o1", experimentId = "e",
                              time = 0, measurement = 1,
                              stringsAsFactors = FALSE),
    experiments = data.frame(experimentId = rep("e", 3),
                             time = c("-inf", "0", "5"),
                             conditionId = c("c1", "c2", "c3"),
                             stringsAsFactors = FALSE),
    mapping = NULL, sbmlPath = "x", formatVersion = 2L)
  out <- dMod2:::.petab_v2_normalize_tables(tables)
  expect_equal(out$measurements$preequilibrationConditionId, "c1")
  expect_equal(out$measurements$simulationConditionId,       "c2")
  expect_equal(unname(out$startTimes["c2"]), 0)
  expect_equal(out$switches[["c2"]],
               data.frame(time = 5, conditionId = "c3",
                          stringsAsFactors = FALSE))
})



test_that(".petab_v2_normalize_tables applies mapping table substitutions", {
  tables <- list(
    parameters = data.frame(parameterId = c("species_a_init", "k"),
                            lowerBound = c(0, 0), upperBound = c(10, 10),
                            nominalValue = c(1, 0.5),
                            estimate = c("false", "true"),
                            stringsAsFactors = FALSE),
    observables = data.frame(observableId = "o1",
                             observableFormula = "species_a",
                             noiseFormula = "1",
                             noiseDistribution = "normal",
                             stringsAsFactors = FALSE),
    conditions = data.frame(conditionId = "c1",
                            targetId = "species_a",
                            targetValue = "species_a_init",
                            stringsAsFactors = FALSE),
    measurements = data.frame(observableId = "o1", experimentId = "exp",
                              time = 0, measurement = 1,
                              stringsAsFactors = FALSE),
    experiments = data.frame(experimentId = "exp", time = "0",
                             conditionId = "c1",
                             stringsAsFactors = FALSE),
    mapping = data.frame(petabEntityId = "species_a",
                         modelEntityId = "A_internal",
                         stringsAsFactors = FALSE),
    sbmlPath = "x", formatVersion = 2L)
  out <- dMod2:::.petab_v2_normalize_tables(tables)
  expect_equal(out$observables$observableFormula, "A_internal")
  # condition target column renamed to model entity name
  expect_true("A_internal" %in% colnames(out$conditions))
})



test_that("readPEtabYaml dispatches v1 vs v2 schema", {
  td <- tempfile("petab_v2_"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  writeLines("conditionId\ttargetId\ttargetValue\nc1\ta0\t1\n",
             file.path(td, "conditions.tsv"))
  writeLines("experimentId\ttime\tconditionId\nexp\t0\tc1\n",
             file.path(td, "experiments.tsv"))
  writeLines("observableId\tobservableFormula\tnoiseFormula\tnoiseDistribution\no1\tA\t1\tnormal\n",
             file.path(td, "observables.tsv"))
  writeLines("observableId\texperimentId\ttime\tmeasurement\no1\texp\t0\t1\n",
             file.path(td, "measurements.tsv"))
  writeLines("parameterId\tlowerBound\tupperBound\tnominalValue\testimate\nk\t0\t1\t0.5\ttrue\n",
             file.path(td, "parameters.tsv"))
  writeLines("<sbml/>", file.path(td, "model.xml"))

  yaml::write_yaml(list(
    format_version    = "2.0.0",
    parameter_files   = list("parameters.tsv"),
    model_files       = list(my_model = list(location = "model.xml",
                                             language = "sbml")),
    observable_files  = list("observables.tsv"),
    measurement_files = list("measurements.tsv"),
    condition_files   = list("conditions.tsv"),
    experiment_files  = list("experiments.tsv")
  ), file.path(td, "problem.yaml"))

  m <- readPEtabYaml(file.path(td, "problem.yaml"))
  expect_identical(m$formatVersion, 2L)
  expect_equal(m$problems[[1]]$modelID, "my_model")
  expect_match(m$problems[[1]]$sbmlFile,        "model\\.xml$")
  expect_match(m$problems[[1]]$experimentFile,  "experiments\\.tsv$")
  expect_null(m$problems[[1]]$mappingFile)
})



test_that("readPEtabYaml errors on non-SBML model language", {
  td <- tempfile("petab_v2_"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  writeLines("dummy", file.path(td, "model.bngl"))
  writeLines("dummy", file.path(td, "p.tsv"))
  writeLines("dummy", file.path(td, "o.tsv"))
  writeLines("dummy", file.path(td, "m.tsv"))
  yaml::write_yaml(list(
    format_version = "2.0.0",
    parameter_files = list("p.tsv"),
    model_files = list(m = list(location = "model.bngl", language = "bngl")),
    observable_files = list("o.tsv"),
    measurement_files = list("m.tsv")
  ), file.path(td, "problem.yaml"))

  expect_error(readPEtabYaml(file.path(td, "problem.yaml")),
               regexp = "SBML")
})



test_that("exportSbml refuses a compartment sized by its own symbol without a value", {
  if (!.libsbml_works()) skip("libsbml virtualenv not available")
  reactions <- eqnlist() %>%
    addReaction("A", "B", "k*A", "A to B") %>%
    assignCompartment(A = "cell", B = "cell", volume = "cell")
  expect_error(exportSbml(reactions, parameters = c(k = 1), inits = c(A = 1, B = 0),
                          filepath = tempfile(fileext = ".xml")),
               "no size for compartment `cell`")
})



test_that("v1 camel-case prior names parse like their v2 spellings", {
  df <- data.frame(parameterId = c("a", "b"), estimate = 1L,
                   objectivePriorType = c("logNormal", "logLaplace"),
                   objectivePriorParameters = c("0;1", "0;2"))
  specs <- dMod2:::.petab_parse_priors(df, c(a = "lin", b = "lin"))
  expect_equal(vapply(specs, `[[`, character(1), "dist"),
               c(a = "log-normal", b = "log-laplace"))
  expect_equal(vapply(specs, `[[`, character(1), "declared"),
               c(a = "logNormal", b = "logLaplace"))
})



test_that("a prior on the parameter scale is truncated on that scale", {
  df <- data.frame(parameterId = "k", estimate = 1L, lowerBound = 1e-5,
                   upperBound = 1e3, objectivePriorType = "parameterScaleNormal",
                   objectivePriorParameters = "-4;3")
  sp <- dMod2:::.petab_parse_priors(df, c(k = "log10"))$k
  expect_equal(c(sp$lower, sp$upper), c(-5, 3))
})



test_that("a prior term honours hessian = FALSE", {
  # It did not, and the argument fell into `...` and was ignored. The cost was
  # not the wasted work: an objective summed with a prior handed back a zero
  # Hessian to a caller that asked for none, so a reverse-swept objective,
  # which cannot produce one, looked as though it had. That is exactly the
  # invariant a caller uses to check the direction actually arrived.
  specs <- list(list(id = "a", dist = "normal", pars = c(0, 1),
                     lower = -Inf, upper = Inf))
  pf <- dMod2:::.petab_prior_objective(specs)
  pars <- c(a = 0.3, b = 1.1)

  full <- pf(pars, deriv = TRUE)
  expect_false(is.null(full$hessian))
  expect_equal(dim(full$hessian), c(2L, 2L))

  none <- pf(pars, deriv = TRUE, hessian = FALSE)
  expect_null(none$hessian)
  expect_equal(none$gradient, full$gradient)
  expect_equal(none$value, full$value)

  # And through the sum, where it mattered: a term that declares `sweep` is
  # asked for a Hessian, one that does not is asked for none, and a zero matrix
  # from the second used to make the total look Hessian-bearing.
  base <- constraintL2(c(a = 0, b = 0), sigma = 1)
  total <- base + pf
  expect_null(total(pars, deriv = TRUE, hessian = FALSE)$hessian)
})



test_that("importSbml takes the initial assignment of a constant parameter", {
  if (!.libsbml_works()) skip("libsbml virtualenv not available")

  # COPASI writes `pcopy := p` as an initial assignment; a condition that sets
  # `p` has to reach the rate through it (Laske_PLOSComputBiol2019).
  f <- tempfile(fileext = ".xml")
  writeLines(c(
    '<?xml version="1.0" encoding="UTF-8"?>',
    '<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">',
    '<model id="ia"><listOfCompartments>',
    '<compartment id="c" size="1" constant="true" spatialDimensions="3"/></listOfCompartments>',
    '<listOfSpecies><species id="A" compartment="c" initialConcentration="2" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/></listOfSpecies>',
    '<listOfParameters><parameter id="p" value="3" constant="true"/>',
    '<parameter id="pcopy" value="99" constant="true"/>',
    '<parameter id="a0" value="99" constant="true"/></listOfParameters>',
    '<listOfInitialAssignments>',
    '<initialAssignment symbol="pcopy"><math xmlns="http://www.w3.org/1998/Math/MathML"><ci> p </ci></math></initialAssignment>',
    '<initialAssignment symbol="a0"><math xmlns="http://www.w3.org/1998/Math/MathML"><ci> A </ci></math></initialAssignment>',
    '</listOfInitialAssignments>',
    '<listOfReactions><reaction id="r1" reversible="false"><listOfReactants>',
    '<speciesReference species="A" stoichiometry="1" constant="true"/></listOfReactants>',
    '<kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><apply><times/><ci> pcopy </ci><ci> a0 </ci><ci> A </ci></apply></math></kineticLaw>',
    '</reaction></listOfReactions></model></sbml>'), f)
  m <- importSbml(f)
  syms <- getSymbols(m$reactions$rates)
  expect_true("p" %in% syms)
  expect_false(any(c("pcopy", "a0") %in% syms))
  expect_false(any(c("pcopy", "a0") %in% names(m$pars)))
})



test_that("exportSbml declares rate species as modifiers", {
  withr::local_dir(tempdir())
  if (!.libsbml_works()) skip("libsbml virtualenv not available")

  # A rate rule imported as a reaction: B enters the rate of A's production
  # without being consumed or produced.
  reactions <- eqnlist() %>%
    addReaction("", "A", "k1*B - k2*A", "A from B") %>%
    addReaction("B", "", "k3*B", "B decay")
  f <- tempfile(fileext = ".xml")
  exportSbml(reactions, parameters = c(k1 = 1, k2 = 1, k3 = 1),
             inits = c(A = 0, B = 1), filepath = f)
  sbml <- paste(readLines(f), collapse = "")
  expect_match(sbml, "<modifierSpeciesReference species=\"B\"")
  expect_no_match(sbml, "<modifierSpeciesReference species=\"A\"")
})



test_that("v2 → v1 → v2 textual normaliser roundtrips a minimal problem", {
  # Round-trip purely at the table level: build a v2 input, normalise to v1
  # shape, write back as v2, normalise again, and compare key invariants.
  v2 <- list(
    parameters = data.frame(
      parameterId = c("k1"), lowerBound = 1e-3, upperBound = 1e3,
      nominalValue = 0.5, estimate = "true", stringsAsFactors = FALSE),
    observables = data.frame(
      observableId = "o1", observableFormula = "A",
      noiseFormula = "1", noiseDistribution = "normal",
      stringsAsFactors = FALSE),
    conditions = data.frame(
      conditionId = "c1", targetId = "a0", targetValue = "3",
      stringsAsFactors = FALSE),
    measurements = data.frame(
      observableId = "o1", experimentId = "e1",
      time = 0, measurement = 1.0, stringsAsFactors = FALSE),
    experiments = data.frame(
      experimentId = "e1", time = "0", conditionId = "c1",
      stringsAsFactors = FALSE),
    mapping = NULL, sbmlPath = "x", formatVersion = 2L)
  out <- dMod2:::.petab_v2_normalize_tables(v2)
  expect_equal(out$measurements$simulationConditionId, "c1")
  expect_equal(out$conditions$a0[out$conditions$conditionId == "c1"], "3")
})



test_that("SBML roundtrip preserves symbolic volumes, reaction frames and amounts", {

  if (!.libsbml_works()) skip("libsbml virtualenv not available")
  withr::local_dir(tempdir())

  f <- eqnlist(
    smatrix = matrix(c(-1, -1, 1, NA, NA, -1), nrow = 2, byrow = TRUE,
                     dimnames = list(NULL, c("L", "R", "C"))),
    states = c("L", "R", "C"), rates = c("k_on*L*R", "k_off*C"),
    description = c("bind", "unbind"),
    compartments = list(ext = "V_ext", cyt = "V_cyt"),
    compartmentOf = c(L = "ext", R = "cyt", C = "cyt"),
    reactionCompartment = c("ext", NA),
    amountStates = "C")

  path <- file.path(tempdir(), "roundtrip.xml")
  exportSbml(f, parameters = c(k_on = 1, k_off = 2, V_ext = 5, V_cyt = 3),
             inits = c(L = 1, R = 2, C = 0), filepath = path)
  xml <- readLines(path)
  expect_true(any(grepl('id="C".*hasOnlySubstanceUnits="true"', xml)))

  g <- importSbml(path)$reactions

  expect_equal(unname(g$compartmentOf[c("L", "R", "C")]), c("ext", "cyt", "cyt"))
  expect_equal(g$compartments$ext$volume, "V_ext")
  expect_equal(g$reactionCompartment[1], "ext")
  expect_equal(g$amountStates, "C")

  # The rate strings are not simplified by the roundtrip, so compare numerically.
  pars <- c(k_on = 1, k_off = 2, V_ext = 5, V_cyt = 3, L = 0.7, R = 0.3, C = 0.5)
  value <- function(el) vapply(as.eqnvec(el)[c("L", "R", "C")],
                               function(e) eval(parse(text = e), as.list(pars)), numeric(1))
  expect_equal(value(g), value(f))
})
