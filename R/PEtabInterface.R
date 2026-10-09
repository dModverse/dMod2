## R/PEtabInterface.R: PEtab v1/v2 import and export on top of the SBML interface.
## The v2 reader translates its tables into the v1 shapes the `.petab_parse_*`
## helpers consume; the exporter writes either version (default "2.0.0").

# Internal: classify a YAML format_version string/number as the major version
# integer dMod cares about. v1 accepts 1 / "1" / "1.0.0"; v2 accepts 2 /
# "2.0.0" / "2.x.y"; everything else is an error.
.petab_major_version <- function(fv) {
  if (length(fv) == 0L) return(1L)  # v1 default for unversioned YAMLs
  s <- as.character(fv)[1L]
  m <- regmatches(s, regexec("^([1-9][0-9]*)", s))[[1L]]
  if (length(m) < 2L)
    stop("Unrecognised PEtab format_version: ", s)
  as.integer(m[2L])
}


#' Deprecated PEtab Readers
#'
#' `readPetabYaml()` and `readPetabTables()` are deprecated, use
#' [readPEtabYaml()] and [readPEtabTables()].
#'
#' @param yamlPath Path to the PEtab YAML file.
#' @return As [readPEtabYaml()] and [readPEtabTables()].
#' @keywords internal
#' @name readPetab-deprecated
NULL

#' @rdname readPetab-deprecated
#' @export
readPetabYaml <- function(yamlPath) {
  warning("'readPetabYaml' is deprecated, use 'readPEtabYaml'.", call. = FALSE)
  readPEtabYaml(yamlPath)
}

#' @rdname readPetab-deprecated
#' @export
readPetabTables <- function(yamlPath) {
  warning("'readPetabTables' is deprecated, use 'readPEtabTables'.",
          call. = FALSE)
  readPEtabTables(yamlPath)
}


## --- low-level YAML / TSV readers ------------------------------------------

#' Read a PEtab YAML File
#'
#' Reads a PEtab v1 or v2 YAML file, as selected by its `format_version`, and
#' resolves the paths of the files it names. Only one problem per file is
#' supported.
#'
#' @param yamlPath Path to the PEtab YAML file.
#' @return A list with `baseDir`, `formatVersion` (integer major version),
#'   `parameterFile`, and a one-element list `problems` whose entry holds the
#'   absolute paths `sbmlFile`, `conditionFile`, `measurementFile` and
#'   `observableFile`. For v2 the entry also holds `experimentFile` and
#'   `mappingFile` (possibly `NULL`), `modelID`, the first model id, and
#'   `models`, the SBML path of every model id.
#' @seealso [readPEtabTables()], [importPEtab()]
#' @examplesIf requireNamespace("yaml", quietly = TRUE)
#' yaml <- system.file("extdata/petab_boehm/Boehm.yaml", package = "dMod2")
#' str(readPEtabYaml(yaml))
#' @export
readPEtabYaml <- function(yamlPath) {

  .require_ns("yaml", "PEtab import")
  yamlPath <- normalizePath(yamlPath, mustWork = TRUE)
  baseDir  <- dirname(yamlPath)
  m <- yaml::read_yaml(yamlPath)

  major <- .petab_major_version(m$format_version)
  if (!major %in% c(1L, 2L))
    stop("PEtab format_version ", m$format_version,
         " not supported (only v1 and v2).")

  resolve <- function(p) normalizePath(file.path(baseDir, p), mustWork = TRUE)

  if (major == 1L) {
    if (length(m$problems) != 1L)
      stop("Only single-problem PEtab YAML is supported (got ",
           length(m$problems), ").")
    prob <- m$problems[[1]]
    pick_one <- function(x, slot) {
      if (length(x) == 0L) stop("YAML problem missing `", slot, "`.")
      if (length(x) > 1L) stop("Only one ", slot, " per problem is supported.")
      x[[1]]
    }
    sbml_path <- resolve(pick_one(prob$sbml_files, "sbml_files"))
    return(list(
      baseDir         = baseDir,
      formatVersion   = 1L,
      parameterFile   = resolve(m$parameter_file),
      problems = list(list(
        sbmlFile        = sbml_path,
        conditionFile   = resolve(pick_one(prob$condition_files, "condition_files")),
        measurementFile = resolve(pick_one(prob$measurement_files, "measurement_files")),
        observableFile  = resolve(pick_one(prob$observable_files, "observable_files")),
        experimentFile  = NULL,
        mappingFile     = NULL,
        modelID         = NA_character_,
        # `model` is a synthetic key for the v1 single-model case; v1 has no
        # modelId column on measurements, so this only ever appears as the
        # default fall-through in `importPEtab()`.
        models          = setNames(sbml_path, "model")
      ))
    ))
  }

  ## --- v2 -----------------------------------------------------------------
  pick_one_top <- function(x, slot, optional = FALSE) {
    if (length(x) == 0L) {
      if (optional) return(NULL)
      stop("YAML missing required `", slot, "`.")
    }
    if (length(x) > 1L)
      stop("Only one ", slot, " per problem is supported (got ",
           length(x), ").")
    if (is.list(x)) x[[1L]] else as.character(x)[1L]
  }

  param_file <- pick_one_top(m$parameter_files, "parameter_files")

  models <- m$model_files
  if (length(models) == 0L)
    stop("YAML missing required `model_files`.")
  if (is.null(names(models)) || any(!nzchar(names(models))))
    stop("`model_files` entries must be a named mapping `<modelId>: { location: ... }`.")

  resolved_models <- character(length(models))
  names(resolved_models) <- names(models)
  for (i in seq_along(models)) {
    entry <- models[[i]]
    mid   <- names(models)[i]
    if (is.null(entry$location))
      stop("model_files entry `", mid, "` missing `location`.")
    language <- tolower(as.character(entry$language %||% "sbml"))
    if (!identical(language, "sbml"))
      stop("Only SBML models are supported by dMod's PEtab importer (got `",
           language, "` for model id `", mid, "`).")
    resolved_models[i] <- resolve(entry$location)
  }

  obs_file  <- pick_one_top(m$observable_files,   "observable_files")
  meas_file <- pick_one_top(m$measurement_files,  "measurement_files")
  cond_file <- pick_one_top(m$condition_files,    "condition_files",  optional = TRUE)
  exp_file  <- pick_one_top(m$experiment_files,   "experiment_files", optional = TRUE)
  map_file  <- pick_one_top(m$mapping_files,      "mapping_files",    optional = TRUE)

  list(
    baseDir         = baseDir,
    formatVersion   = 2L,
    parameterFile   = resolve(param_file),
    problems = list(list(
      sbmlFile        = unname(resolved_models[1L]),
      conditionFile   = if (!is.null(cond_file)) resolve(cond_file) else NULL,
      measurementFile = resolve(meas_file),
      observableFile  = resolve(obs_file),
      experimentFile  = if (!is.null(exp_file)) resolve(exp_file) else NULL,
      mappingFile     = if (!is.null(map_file)) resolve(map_file) else NULL,
      modelID         = names(resolved_models)[1L],
      models          = resolved_models
    ))
  )
}


#' Read the Tables of a PEtab Problem
#'
#' Reads the TSV tables named in a PEtab YAML file into data frames. The SBML
#' model is not read.
#'
#' @param yamlPath Path to the PEtab YAML file.
#' @return A list with the data frames `parameters`, `conditions`,
#'   `measurements`, `observables`, `experiments` and `mapping` (each `NULL`
#'   when the problem has no such table), `sbmlPath`, the path of the first
#'   model, `sbmlPaths`, the SBML path of every model id, and `formatVersion`
#'   (integer).
#' @seealso [readPEtabYaml()], [importPEtab()]
#' @examplesIf requireNamespace("yaml", quietly = TRUE)
#' yaml <- system.file("extdata/petab_boehm/Boehm.yaml", package = "dMod2")
#' tables <- readPEtabTables(yaml)
#' head(tables$measurements)
#' tables$parameters
#' @export
readPEtabTables <- function(yamlPath) {

  m  <- readPEtabYaml(yamlPath)
  pr <- m$problems[[1]]

  list(
    parameters    = .petab_read_tsv(m$parameterFile),
    conditions    = if (!is.null(pr$conditionFile))
                      .petab_read_tsv(pr$conditionFile) else NULL,
    measurements  = .petab_read_tsv(pr$measurementFile),
    observables   = .petab_read_tsv(pr$observableFile),
    experiments   = if (!is.null(pr$experimentFile))
                      .petab_read_tsv(pr$experimentFile) else NULL,
    mapping       = if (!is.null(pr$mappingFile))
                      .petab_read_tsv(pr$mappingFile) else NULL,
    sbmlPath      = pr$sbmlFile,
    sbmlPaths     = pr$models %||% setNames(pr$sbmlFile, pr$modelID),
    formatVersion = m$formatVersion
  )
}


.petab_read_tsv <- function(path) {
  utils::read.delim(path, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE,
                    check.names = FALSE,
                    na.strings = c("", "NA"),
                    strip.white = TRUE)
}


## --- v2 -> v1 normalizer ---------------------------------------------------

# Internal: key conditions by `conditionName` where it is present and unique,
# else by `conditionId`. Returns the tables and the name -> id map the exporter
# uses to write the ids back.
.petab_use_condition_names <- function(tables) {
  df <- tables$conditions
  if (is.null(df) || !"conditionName" %in% colnames(df))
    return(c(tables, list(conditionIds = NULL)))

  ids   <- as.character(df$conditionId)
  names_ <- as.character(df$conditionName)
  usable <- !is.na(names_) & nzchar(trimws(names_)) &
            !names_ %in% ids[!ids %in% names_] &
            !duplicated(names_) & !duplicated(names_, fromLast = TRUE)
  if (!any(usable)) return(c(tables, list(conditionIds = NULL)))

  key <- ifelse(usable, names_, ids)
  map <- setNames(ids, key)

  df$conditionId <- key
  tables$conditions <- df
  for (col in c("simulationConditionId", "preequilibrationConditionId")) {
    if (!col %in% colnames(tables$measurements)) next
    v <- as.character(tables$measurements[[col]])
    hit <- match(v, ids)
    v[!is.na(hit)] <- key[hit[!is.na(hit)]]
    tables$measurements[[col]] <- v
  }
  c(tables, list(conditionIds = map))
}


## Internal: rewrite v2 tables into v1 shapes. The mapping table, if any, is
## applied textually to every string column, so ids match the imported SBML.
.petab_v2_normalize_tables <- function(tables) {

  ## --- 1. mapping table: build petab->model substitution -----------------
  mapping <- tables$mapping
  apply_mapping <- function(s) s
  if (!is.null(mapping) && nrow(mapping) > 0L) {
    if (!all(c("petabEntityId", "modelEntityId") %in% colnames(mapping)))
      stop("mapping.tsv must have columns `petabEntityId` and `modelEntityId`.")
    aliases <- mapping[!is.na(mapping$modelEntityId) &
                       nzchar(mapping$modelEntityId), , drop = FALSE]
    if (nrow(aliases) > 0L) {
      pat <- paste0("\\b", aliases$petabEntityId, "\\b")
      repl <- aliases$modelEntityId
      apply_mapping <- function(s) {
        if (length(s) == 0L) return(s)
        out <- s
        for (i in seq_along(pat))
          out <- gsub(pat[i], repl[i], out, perl = TRUE)
        out
      }
    }
  }

  ## --- 2. parameters: add parameterScale, coerce estimate ----------------
  par_df <- tables$parameters
  if (!"parameterId" %in% colnames(par_df))
    stop("parameters.tsv missing required column `parameterId`.")
  if (!"parameterScale" %in% colnames(par_df))
    par_df$parameterScale <- "lin"
  if ("estimate" %in% colnames(par_df)) {
    e <- par_df$estimate
    if (is.logical(e)) {
      par_df$estimate <- as.integer(e)
    } else {
      etxt <- tolower(trimws(as.character(e)))
      par_df$estimate <- ifelse(etxt %in% c("true", "1"), 1L,
                         ifelse(etxt %in% c("false", "0"), 0L, NA_integer_))
      if (any(is.na(par_df$estimate)))
        stop("parameters.tsv `estimate` must be true/false (or 1/0); got: ",
             paste(unique(e[is.na(par_df$estimate)]), collapse = ", "))
    }
  }

  ## --- 3. observables: split noiseDistribution, rewrite placeholders -----
  obs_df <- tables$observables
  if (!"observableId" %in% colnames(obs_df))
    stop("observables.tsv missing required column `observableId`.")
  if (!"observableFormula" %in% colnames(obs_df))
    stop("observables.tsv missing required column `observableFormula`.")
  # read.delim infers numeric type when a column parses entirely as numbers;
  # strings are substituted into these formulas, so coerce to character first.
  for (col in c("observableFormula", "noiseFormula",
                "observablePlaceholders", "noisePlaceholders")) {
    if (col %in% colnames(obs_df))
      obs_df[[col]] <- as.character(obs_df[[col]])
  }

  v2_dist <- if ("noiseDistribution" %in% colnames(obs_df))
               obs_df$noiseDistribution else rep("normal", nrow(obs_df))
  v2_dist[is.na(v2_dist) | !nzchar(v2_dist)] <- "normal"
  bad <- setdiff(unique(v2_dist),
                 c("normal", "log-normal", "laplace", "log-laplace"))
  if (length(bad))
    stop("Unsupported v2 noiseDistribution(s): ",
         paste(bad, collapse = ", "),
         ". dMod uses L2 likelihoods only; supported v2 values are ",
         "normal, log-normal, laplace, log-laplace.")
  obs_df$observableTransformation <- ifelse(
    v2_dist %in% c("log-normal", "log-laplace"), "log", "lin")
  obs_df$noiseDistribution <- ifelse(
    v2_dist %in% c("laplace", "log-laplace"), "laplace", "normal")

  rewrite_placeholders <- function(formula, ph_str, prefix, obs_id) {
    if (is.na(formula) || !nzchar(formula)) return(formula)
    if (is.na(ph_str) || !nzchar(ph_str))   return(formula)
    parts <- trimws(strsplit(ph_str, ";", fixed = TRUE)[[1L]])
    parts <- parts[nzchar(parts)]
    for (k in seq_along(parts)) {
      pat <- sprintf("\\b%s\\b", parts[k])
      formula <- gsub(pat, sprintf("%s%d_%s", prefix, k, obs_id),
                      formula, perl = TRUE)
    }
    formula
  }
  if ("observablePlaceholders" %in% colnames(obs_df)) {
    obs_df$observableFormula <- vapply(seq_len(nrow(obs_df)), function(i)
      rewrite_placeholders(obs_df$observableFormula[i],
                           obs_df$observablePlaceholders[i],
                           "observableParameter",
                           obs_df$observableId[i]), character(1))
  }
  if ("noiseFormula" %in% colnames(obs_df) &&
      "noisePlaceholders" %in% colnames(obs_df)) {
    obs_df$noiseFormula <- vapply(seq_len(nrow(obs_df)), function(i)
      rewrite_placeholders(obs_df$noiseFormula[i],
                           obs_df$noisePlaceholders[i],
                           "noiseParameter",
                           obs_df$observableId[i]), character(1))
  }
  obs_df$observableFormula <- apply_mapping(obs_df$observableFormula)
  if ("noiseFormula" %in% colnames(obs_df))
    obs_df$noiseFormula <- apply_mapping(obs_df$noiseFormula)

  ## --- 4. conditions: long -> wide --------------------------------------
  cond_v2 <- tables$conditions
  if (is.null(cond_v2) || nrow(cond_v2) == 0L) {
    cond_wide <- data.frame(conditionId = character(0),
                            stringsAsFactors = FALSE)
  } else {
    needed <- c("conditionId", "targetId", "targetValue")
    miss <- setdiff(needed, colnames(cond_v2))
    if (length(miss))
      stop("v2 conditions.tsv missing column(s): ",
           paste(miss, collapse = ", "))
    cond_v2$targetId    <- apply_mapping(as.character(cond_v2$targetId))
    cond_v2$targetValue <- apply_mapping(as.character(cond_v2$targetValue))
    uniq_cids <- unique(cond_v2$conditionId)
    uniq_tids <- unique(cond_v2$targetId)
    cond_wide <- data.frame(conditionId = uniq_cids,
                            stringsAsFactors = FALSE)
    for (tid in uniq_tids)
      cond_wide[[tid]] <- NA_character_
    for (i in seq_len(nrow(cond_v2))) {
      r <- match(cond_v2$conditionId[i], cond_wide$conditionId)
      tid <- cond_v2$targetId[i]
      if (!is.na(cond_wide[r, tid]) &&
          !identical(cond_wide[r, tid], cond_v2$targetValue[i]))
        stop("Duplicate (conditionId,targetId) in v2 conditions.tsv: ",
             cond_v2$conditionId[i], " / ", tid)
      cond_wide[r, tid] <- cond_v2$targetValue[i]
    }
  }

  ## --- 5. experiments -> (simCondId, preeqCondId) per measurement -------
  meas <- tables$measurements
  required_meas <- c("observableId", "time", "measurement")
  miss <- setdiff(required_meas, colnames(meas))
  if (length(miss))
    stop("measurements.tsv missing required column(s): ",
         paste(miss, collapse = ", "))

  # `modelId` is preserved for the multi-model dispatch in importPEtab().
  # Empty/NA cells are filled later (the importer knows the canonical default).

  exp_df <- tables$experiments
  # A measurement without an experiment applies no overrides. It still needs a
  # condition to key on, so one is synthesised: `cond_1`, or the next free
  # number when the problem already uses that id.
  dflt_cond <- local({
    taken <- as.character(cond_wide$conditionId)
    i <- 1L
    repeat {
      cand <- paste0("cond_", i)
      if (!cand %in% taken) return(cand)
      i <- i + 1L
    }
  })

  exp_map <- list()  # experimentId -> list(sim=..., preeq=...)

  exp_ids_in_meas <- if ("experimentId" %in% colnames(meas))
                       unique(meas$experimentId[!is.na(meas$experimentId) &
                                                nzchar(meas$experimentId)])
                     else character(0)

  if (!is.null(exp_df) && nrow(exp_df) > 0L) {
    needed <- c("experimentId", "time", "conditionId")
    miss <- setdiff(needed, colnames(exp_df))
    if (length(miss))
      stop("experiments.tsv missing column(s): ",
           paste(miss, collapse = ", "))
    for (eid in unique(exp_df$experimentId)) {
      sub <- exp_df[exp_df$experimentId == eid, , drop = FALSE]
      times <- suppressWarnings(as.numeric(sub$time))
      # Accept "-inf"/"-Inf" parsing as -Inf
      txt <- tolower(trimws(as.character(sub$time)))
      times[txt %in% c("-inf", "-infinity")] <- -Inf
      times[txt %in% c("inf", "+inf", "infinity")] <- Inf
      ord <- order(times)
      times <- times[ord]
      # An empty conditionId cell reads as NA and means "no overrides".
      cids  <- as.character(sub$conditionId)[ord]
      cids[is.na(cids)] <- ""
      # A leading -inf period is the preequilibration and the first finite
      # period starts the simulation; integration begins at the first
      # measurement time, so a nonzero start needs no further handling.
      preeq <- ""
      has_preeq <- is.infinite(times[1L]) && times[1L] < 0
      if (has_preeq) {
        preeq <- cids[1L]
        times <- times[-1L]
        cids  <- cids[-1L]
      }
      # A period without a conditionId applies no overrides, which is what the
      # synthesised condition stands for.
      dflt <- function(cid) if (nzchar(cid)) cid else dflt_cond
      if (!length(times)) {
        # Preequilibration only: the simulation holds the steady state.
        exp_map[[eid]] <- list(sim = dflt(preeq), preeq = dflt(preeq), t0 = 0)
      } else {
        # Periods after the first switch conditions mid-run; they become
        # events on the model.
        exp_map[[eid]] <- list(
          sim      = dflt(cids[1L]),
          preeq    = if (has_preeq) dflt(preeq) else "",
          t0       = times[1L],
          switches = if (length(times) > 1L)
                       data.frame(time = times[-1L], conditionId = cids[-1L],
                                  stringsAsFactors = FALSE)
                     else NULL)
      }
    }
  }

  if ("experimentId" %in% colnames(meas)) {
    eid <- ifelse(is.na(meas$experimentId), "", meas$experimentId)
    sim   <- character(nrow(meas))
    preeq <- character(nrow(meas))
    start <- numeric(nrow(meas))
    for (i in seq_along(eid)) {
      e <- eid[i]
      if (!nzchar(e)) {
        # v2: an empty experimentId means the model as is; a sentinel condition
        # gives the trafo stage something to key on.
        sim[i]   <- dflt_cond
        preeq[i] <- ""
        start[i] <- 0
        next
      }
      if (is.null(exp_map[[e]]))
        stop("measurements.tsv references unknown experimentId `", e, "`.")
      sim[i]   <- exp_map[[e]]$sim
      preeq[i] <- exp_map[[e]]$preeq
      start[i] <- exp_map[[e]]$t0
    }
    meas$simulationConditionId        <- sim
    meas$preequilibrationConditionId  <- preeq
    meas$experimentId <- NULL

    # With a synthesised default condition, a no-override row is appended to
    # cond_wide so .petab_parse_conditions sees it.
    if (any(c(sim, preeq) == dflt_cond) &&
        !dflt_cond %in% cond_wide$conditionId) {
      new_row <- cond_wide[NA_integer_, , drop = FALSE][1L, , drop = FALSE]
      new_row$conditionId <- dflt_cond
      cond_wide <- rbind(cond_wide, new_row)
    }
  } else {
    # No experimentId column at all; spec says empty means model-as-is.
    start <- numeric(nrow(meas))
    meas$simulationConditionId <- dflt_cond
    meas$preequilibrationConditionId <- ""
    if (!dflt_cond %in% cond_wide$conditionId) {
      new_row <- cond_wide[NA_integer_, , drop = FALSE][1L, , drop = FALSE]
      new_row$conditionId <- dflt_cond
      cond_wide <- rbind(cond_wide, new_row)
    }
  }

  # observableParameters / noiseParameters columns survive unchanged; their
  # placeholder substitution machinery already keys on the v1-style sentinels
  # written into the observable/noise formulas above.
  for (col in c("observableParameters", "noiseParameters")) {
    if (col %in% colnames(meas))
      meas[[col]] <- ifelse(is.na(meas[[col]]), "", as.character(meas[[col]]))
  }

  # One simulation start per condition. Distinct conditions may start at
  # different times, so the objective groups by this below.
  start_times <- tapply(start, meas$simulationConditionId, function(v) v[1L])
  start_times <- setNames(as.numeric(start_times), names(start_times))

  # Mid-run condition switches, keyed by the condition the experiment starts in.
  switch_map <- list()
  for (eid in names(exp_map)) {
    sw <- exp_map[[eid]]$switches
    if (is.null(sw)) next
    key <- exp_map[[eid]]$sim
    if (!is.null(switch_map[[key]]) && !identical(switch_map[[key]], sw))
      stop("Experiments starting in condition `", key,
           "` declare different period schedules; that needs one model per ",
           "schedule, which is not implemented.")
    switch_map[[key]] <- sw
  }

  list(parameters   = par_df,
       observables  = obs_df,
       startTimes   = start_times,
       switches     = switch_map,
       defaultCondition = dflt_cond,
       conditions   = cond_wide,
       measurements = meas,
       sbmlPath      = tables$sbmlPath,
       sbmlPaths     = tables$sbmlPaths,
       formatVersion = 1L)
}


## --- per-table parsers (unit-testable, no SBML side effects) ---------------

# Internal: parameters.tsv -> list(pouter, lower, upper, fixed, scales, priors).
# Estimated values and bounds are on the parameter scale, `fixed` stays linear;
# `priors` is NULL or one record per parameter from .petab_parse_priors().
.petab_parse_parameters <- function(df) {

  required <- c("parameterId", "parameterScale", "lowerBound", "upperBound",
                "nominalValue", "estimate")
  miss <- setdiff(required, colnames(df))
  if (length(miss))
    stop("parameters.tsv missing required column(s): ",
         paste(miss, collapse = ", "))

  scales <- setNames(df$parameterScale, df$parameterId)
  if (any(!scales %in% c("lin", "log", "log10")))
    stop("Unknown parameterScale(s): ",
         paste(unique(scales[!scales %in% c("lin", "log", "log10")]),
               collapse = ", "))

  # PEtab writes nominal values and bounds linear whatever the parameterScale;
  # dMod's outer parameters live on that scale, so they are transformed here.
  to_num <- function(col) suppressWarnings(as.numeric(col))
  apply_fwd_scale <- function(values, ids) {
    sc <- scales[ids]
    out <- values
    log_idx   <- which(sc == "log")
    log10_idx <- which(sc == "log10")
    if (length(log_idx))   out[log_idx]   <- log(values[log_idx])
    if (length(log10_idx)) out[log10_idx] <- log10(values[log10_idx])
    out
  }

  pouter_idx <- which(df$estimate == 1)
  fixed_idx  <- which(df$estimate == 0)

  pouter_ids <- df$parameterId[pouter_idx]
  fixed_ids  <- df$parameterId[fixed_idx]

  pouter <- setNames(apply_fwd_scale(to_num(df$nominalValue[pouter_idx]),
                                     pouter_ids),
                     pouter_ids)
  # `fixed` parameters are passed straight through to the trafo as numeric
  # constants on the *inner* (linear) scale, no scale wrapping in the
  # trafo, so linear values are kept here regardless of parameterScale.
  fixed  <- setNames(to_num(df$nominalValue[fixed_idx]), fixed_ids)
  lower  <- setNames(apply_fwd_scale(to_num(df$lowerBound[pouter_idx]),
                                     pouter_ids),
                     pouter_ids)
  upper  <- setNames(apply_fwd_scale(to_num(df$upperBound[pouter_idx]),
                                     pouter_ids),
                     pouter_ids)

  priors <- .petab_parse_priors(df, scales)

  list(pouter = pouter, lower = lower, upper = upper,
       fixed = fixed, scales = scales, priors = priors)
}


# The time a steady-state ("inf") measurement is evaluated at.
.petab_ss_time <- 1e10


# Internal: one record list(id, dist, pars, lower, upper, ...) per prior, or
# NULL; v2 `priorDistribution` columns win over v1's. A `parameterScale*` prior
# is on the optimised parameter; a plain `normal` on a non-lin scale is rejected.
.petab_parse_priors <- function(df, scales) {

  pick <- function(a, b)
    if (a %in% colnames(df)) a else if (b %in% colnames(df)) b else NA_character_

  dist_col <- pick("priorDistribution",  "objectivePriorType")
  pars_col <- pick("priorParameters",    "objectivePriorParameters")
  if (is.na(dist_col) || is.na(pars_col)) return(NULL)

  pd  <- as.character(df[[dist_col]])
  pp  <- as.character(df[[pars_col]])
  est <- if ("estimate" %in% colnames(df)) df$estimate == 1L else rep(TRUE, nrow(df))
  active <- !is.na(pd) & nzchar(pd) & est
  if (!any(active)) return(NULL)

  n_expected <- c(uniform = 2L, normal = 2L, `log-normal` = 2L, cauchy = 2L,
                  chisquare = 1L, exponential = 1L, gamma = 2L, laplace = 2L,
                  `log-laplace` = 2L, `log-uniform` = 2L, rayleigh = 1L)

  specs <- list()
  for (i in which(active)) {
    pid <- df$parameterId[i]
    d   <- pd[i]
    sc  <- scales[[pid]] %||% "lin"

    scaled <- grepl("^parameterScale", d)
    base   <- if (scaled) tolower(sub("^parameterScale", "", d)) else d
    # v1 spells the log families in camel case.
    v1_log <- c(logNormal = "log-normal", logLaplace = "log-laplace")
    if (base %in% names(v1_log)) base <- v1_log[[base]]
    if (!base %in% names(n_expected))
      stop("Unsupported priorDistribution `", d, "` for parameter `", pid,
           "`. Known: ", paste(names(n_expected), collapse = ", "),
           " (optionally `parameterScale`-prefixed).")
    if (!scaled && base %in% c("normal", "laplace") && !identical(sc, "lin"))
      stop("priorDistribution `", d, "` for parameter `", pid,
           "` (parameterScale = `", sc, "`) is ambiguous: dMod's prior acts ",
           "on the optimizer's parameter (the parameterScale-transformed ",
           "value). Use `parameterScale", sub("^(.)", "\\U\\1", base,
           perl = TRUE), "` to declare a prior on that view, or set ",
           "parameterScale = `lin`.")

    parts <- trimws(strsplit(pp[i], ";", fixed = TRUE)[[1L]])
    vals  <- suppressWarnings(as.numeric(parts))
    if (length(vals) != n_expected[[base]] || anyNA(vals))
      stop("priorParameters for `", pid, "` must be ", n_expected[[base]],
           " `;`-separated number(s) for `", d, "`; got `", pp[i], "`.")

    lo <- if ("lowerBound" %in% colnames(df)) as.numeric(df$lowerBound[i]) else -Inf
    hi <- if ("upperBound" %in% colnames(df)) as.numeric(df$upperBound[i]) else  Inf
    # A prior on the parameter scale is truncated where that scale puts the
    # bounds; the table holds them linear.
    if (scaled && sc %in% c("log", "log10")) {
      to_scale <- if (sc == "log10") log10 else log
      lo <- to_scale(lo); hi <- to_scale(hi)
    }
    specs[[pid]] <- list(id = pid, dist = base, pars = vals,
                         lower = if (is.na(lo)) -Inf else lo,
                         upper = if (is.na(hi))  Inf else hi,
                         declared = d, declaredPars = pp[i])
  }
  if (!length(specs)) return(NULL)
  specs
}


# Internal: log density with first and second derivative for the PEtab priors
# without a dMod constructor. The plain `log-*` spellings are densities in the
# parameter itself and include the change-of-variables term.
.petab_prior_logdens <- function(dist, pars, x) {
  a <- pars[1L]; b <- if (length(pars) > 1L) pars[2L] else NA_real_
  pos <- x > 0
  switch(dist,
    uniform = list(ld = if (x >= a && x <= b) -log(b - a) else -Inf,
                   d1 = 0, d2 = 0),
    `log-normal` = if (!pos) list(ld = -Inf, d1 = 0, d2 = 0) else {
      u <- log(x)
      list(ld = stats::dlnorm(x, a, b, log = TRUE),
           d1 = -1 / x - (u - a) / (b^2 * x),
           d2 = 1 / x^2 - 1 / (b^2 * x^2) + (u - a) / (b^2 * x^2))
    },
    `log-laplace` = if (!pos) list(ld = -Inf, d1 = 0, d2 = 0) else {
      u <- log(x)
      list(ld = -u - log(2 * b) - abs(u - a) / b,
           d1 = -1 / x - sign(u - a) / (b * x),
           d2 = 1 / x^2 + sign(u - a) / (b * x^2))
    },
    `log-uniform` = if (!pos) list(ld = -Inf, d1 = 0, d2 = 0) else
      list(ld = if (x >= a && x <= b) -log(x) - log(log(b) - log(a)) else -Inf,
           d1 = -1 / x, d2 = 1 / x^2),
    stop("Unsupported prior distribution `", dist, "`.")
  )
}


# Internal: probability mass a PEtab prior puts on [lower, upper]. The
# truncation constant does not depend on the parameter, so it shifts the
# objective without touching gradient or Hessian.
.petab_prior_mass <- function(dist, pars, lower, upper) {
  a <- pars[1L]; b <- if (length(pars) > 1L) pars[2L] else NA_real_
  plap <- function(q, m, s)
    ifelse(q < m, 0.5 * exp((q - m) / s), 1 - 0.5 * exp(-(q - m) / s))
  cdf <- switch(dist,
    uniform       = function(q) stats::punif(q, a, b),
    normal        = function(q) stats::pnorm(q, a, b),
    `log-normal`  = function(q) ifelse(q <= 0, 0, stats::plnorm(q, a, b)),
    cauchy        = function(q) stats::pcauchy(q, a, b),
    chisquare     = function(q) ifelse(q <= 0, 0, stats::pchisq(q, a)),
    exponential   = function(q) ifelse(q <= 0, 0, stats::pexp(q, rate = 1 / a)),
    gamma         = function(q) ifelse(q <= 0, 0, stats::pgamma(q, shape = a, scale = b)),
    laplace       = function(q) plap(q, a, b),
    `log-laplace` = function(q) ifelse(q <= 0, 0, plap(log(q), a, b)),
    `log-uniform` = function(q) ifelse(q <= 0, 0, stats::punif(log(q), log(a), log(b))),
    rayleigh      = function(q) ifelse(q <= 0, 0, 1 - exp(-q^2 / (2 * a^2))),
    stop("Unsupported prior distribution `", dist, "`."))
  m <- cdf(upper) - cdf(lower)
  if (!is.finite(m) || m <= 0) 1 else m
}


# Internal: the prior term of a PEtab objective on the -2 log scale, reported in
# `attr(, "prior")`. Families with a dMod constructor go through it; the others
# and the truncation constant, which leaves gradient and Hessian alone, are here.
.petab_prior_objective <- function(specs, attrName = "prior",
                                   condition = NULL) {

  one <- function(sp) {
    id <- sp$id; a <- sp$pars[1L]
    b  <- if (length(sp$pars) > 1L) sp$pars[2L] else NA_real_
    switch(sp$dist,
      normal      = constraintL2(setNames(a, id), sigma = b,
                                 attrName = attrName, condition = condition),
      laplace     = constraintL1(setNames(a, id), sigma = b,
                                 attrName = attrName, condition = condition),
      cauchy      = constraintCauchy(setNames(a, id), sigma = b,
                                     attrName = attrName, condition = condition),
      gamma       = constraintGamma(setNames(a, id), scale = b,
                                    attrName = attrName, condition = condition),
      exponential = constraintExponential(setNames(a, id),
                                          attrName = attrName, condition = condition),
      chisquare   = constraintChisq(setNames(a, id),
                                    attrName = attrName, condition = condition),
      rayleigh    = constraintRayleigh(setNames(a, id),
                                       attrName = attrName, condition = condition),
      NULL)
  }

  parts <- Filter(Negate(is.null), lapply(specs, one))
  rest  <- specs[vapply(specs, function(sp) is.null(one(sp)), logical(1))]

  # `constraintL2` is the penalty form and drops the normalisation of the
  # Gaussian; PEtab reports a normalised density, so the constant comes back
  # here next to the truncation.
  norm_const <- function(sp)
    if (identical(sp$dist, "normal")) 0.5 * (log(2 * pi) + 2 * log(sp$pars[2L]))
    else 0

  ids   <- vapply(specs, `[[`, character(1), "id")
  const <- unname(vapply(specs, function(sp)
    log(.petab_prior_mass(sp$dist, sp$pars, sp$lower, sp$upper)) + norm_const(sp),
    numeric(1)))

  # `hessian` is declared: left to `...` it would be ignored, and a zero Hessian
  # returned to a caller that asked for none makes a reverse sweep treat the
  # whole objective as having one.
  myfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE,
                   hessian = NULL, conditions = condition, env = NULL,
                   cores = getOption("dMod.cores", 1L)) {

    p   <- list(...)[[match.fnargs(list(...), "pars")]]
    all <- c(p, fixed)
    nms <- names(p)
    cv <- .resolveCurvature(deriv, deriv2, hessian, "forward")
    build_hessian <- cv$hessian

    value <- 2 * sum(const[ids %in% names(all)])
    grad  <- setNames(rep(0, length(nms)), nms)
    hess  <- matrix(0, length(nms), length(nms), dimnames = list(nms, nms))

    for (sp in rest) {
      if (!sp$id %in% names(all)) next
      td <- .petab_prior_logdens(sp$dist, sp$pars, as.numeric(all[[sp$id]]))
      value <- value - 2 * td$ld
      if (deriv && sp$id %in% nms) {
        grad[sp$id] <- grad[sp$id] - 2 * td$d1
        if (build_hessian) hess[sp$id, sp$id] <- hess[sp$id, sp$id] - 2 * td$d2
      }
    }

    out <- objlist(value = unname(value),
                   gradient = if (deriv) grad else NULL,
                   hessian  = if (build_hessian) hess else NULL)
    attr(out, attrName) <- out$value
    attr(out, "env") <- env
    out
  }

  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- condition
  attr(myfn, "parameters") <- ids

  Reduce(`+`, parts, myfn)
}


# Internal: one observable parse from several per-model parses. Each model
# inlines its own SBML assignment rules, so the entries differ per model and
# the first one to define an observable wins.
.petab_merge_obs_meta <- function(parts) {
  flds <- c("obs", "obs_trafo", "noise", "noise_dist")
  setNames(lapply(flds, function(f) {
    v <- do.call(c, lapply(parts, `[[`, f))
    v[!duplicated(names(v))]
  }), flds)
}

# Internal: observables.tsv -> list(obs, noise, obs_trafo, noise_dist) keyed by
# observableId. Placeholders stay unsubstituted and a numeric noise stays numeric
# for the normL2 fast path. Defaults are "lin" and "normal".
.petab_parse_observables <- function(df) {

  if (!"observableId" %in% colnames(df))
    stop("observables.tsv missing required column `observableId`.")
  if (!"observableFormula" %in% colnames(df))
    stop("observables.tsv missing required column `observableFormula`.")

  ids <- df$observableId
  if (anyDuplicated(ids))
    stop("Duplicate observableId(s) in observables.tsv: ",
         paste(unique(ids[duplicated(ids)]), collapse = ", "))

  obs <- setNames(df$observableFormula, ids)

  noise_raw <- if ("noiseFormula" %in% colnames(df)) df$noiseFormula
               else rep("1", length(ids))
  noise_raw[is.na(noise_raw)] <- "1"
  names(noise_raw) <- ids

  trafo <- if ("observableTransformation" %in% colnames(df))
             df$observableTransformation else rep("lin", length(ids))
  trafo[is.na(trafo)] <- "lin"
  names(trafo) <- ids

  dist <- if ("noiseDistribution" %in% colnames(df))
            df$noiseDistribution else rep("normal", length(ids))
  dist[is.na(dist)] <- "normal"
  names(dist) <- ids

  if (any(!trafo %in% c("lin", "log", "log10")))
    stop("Unknown observableTransformation(s).")
  if (any(!dist %in% c("normal", "laplace", "log-normal")))
    stop("Unknown noiseDistribution(s).")

  list(obs = obs, noise = noise_raw, obs_trafo = trafo, noise_dist = dist)
}


# Internal: classify each non-id column of conditions.tsv as "init" (a state),
# "compartment" (a volume) or "parameter" (everything else), by matching it
# against the SBML states, compartments and parameters.
.petab_parse_conditions <- function(df, sbml_states = character(),
                                    sbml_compartments = character(),
                                    sbml_pars = character(),
                                    obs_inner = character()) {

  if (!"conditionId" %in% colnames(df))
    stop("conditions.tsv missing required column `conditionId`.")

  cond_ids <- df$conditionId
  if (anyDuplicated(cond_ids))
    stop("Duplicate conditionId(s) in conditions.tsv: ",
         paste(unique(cond_ids[duplicated(cond_ids)]), collapse = ", "))

  override_cols <- setdiff(colnames(df), c("conditionId", "conditionName"))

  col_kind <- vapply(override_cols, function(cn) {
    if (cn %in% sbml_states)              "init"
    else if (cn %in% sbml_compartments)   "compartment"
    else if (cn %in% sbml_pars)           "parameter"
    else if (cn %in% obs_inner)           "parameter"  # observable inner par
    else { warning("Condition column `", cn,
                   "` does not match any SBML or observable symbol; treating as parameter.")
           "parameter" }
  }, character(1))

  grid <- df
  rownames(grid) <- as.character(cond_ids)
  list(grid = grid, col_kind = col_kind, override_cols = override_cols)
}


# Internal: a readable, name-safe label for sub-conditions that differ in their
# observable/noise parameter strings. Punctuation collapses to "_", over-long
# labels get a short md5 tail.
.petab_subcond_label <- function(s) {
  if (length(s) == 0L || all(is.na(s) | s == "")) return(rep("", length(s)))
  vapply(s, function(x) {
    if (is.na(x) || identical(x, "")) return("")
    lab <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
    lab <- gsub("^_+|_+$", "", lab)
    if (!nzchar(lab))
      lab <- substr(digest::digest(x, algo = "md5"), 1L, 8L)
    if (nchar(lab) > 32L) {
      tag <- substr(digest::digest(x, algo = "md5"), 1L, 6L)
      lab <- paste0(substr(lab, 1L, 24L), "_", tag)
    }
    lab
  }, character(1))
}


# Internal: measurements.tsv -> list(data, sub_cond_map, peq_map). Tuples of
# (simulationConditionId, observableParameters, noiseParameters) that vary
# within one condition are split into sub-conditions.
.petab_parse_measurements <- function(df, obs_meta) {

  needed <- c("observableId", "simulationConditionId", "time", "measurement")
  miss <- setdiff(needed, colnames(df))
  if (length(miss))
    stop("measurements.tsv missing required column(s): ",
         paste(miss, collapse = ", "))

  m <- df
  m$preequilibrationConditionId <-
    if ("preequilibrationConditionId" %in% colnames(m))
      ifelse(is.na(m$preequilibrationConditionId), "",
             m$preequilibrationConditionId)
    else rep("", nrow(m))

  # read.delim infers numeric type when a column is all numeric, but PEtab
  # observable/noise parameters can be either numeric or symbol strings.
  # Cast to character so the substitution step has a uniform input type.
  m$observableParameters <-
    if ("observableParameters" %in% colnames(m))
      ifelse(is.na(m$observableParameters), "",
             as.character(m$observableParameters))
    else rep("", nrow(m))

  m$noiseParameters <-
    if ("noiseParameters" %in% colnames(m))
      ifelse(is.na(m$noiseParameters), "",
             as.character(m$noiseParameters))
    else rep("", nrow(m))

  # A noise formula that becomes a finite constant after per-row substitution is
  # the data sigma (normL2 fast path); one with symbols left keeps sigma = NA
  # and goes through the error model.
  sigma <- vapply(seq_len(nrow(m)), function(i) {
    obsId <- m$observableId[i]
    f <- obs_meta$noise[obsId]
    if (is.na(f) || !nzchar(f)) return(NA_real_)
    f <- .petab_substitute_param_string(f, m$observableParameters[i],
                                        prefix = "observableParameter")
    f <- .petab_substitute_param_string(f, m$noiseParameters[i],
                                        prefix = "noiseParameter")
    .petab_eval_constant(f)
  }, numeric(1))

  # Rows whose numeric sigma varies within one (condition, observable) would
  # each compile their own trafo, so they collapse to a uniform literal of the
  # same arity. Agreeing rows keep their strings, which the exporter writes back.
  numeric_noise <- !is.na(sigma) & nzchar(m$noiseParameters)
  if (any(numeric_noise)) {
    key <- paste(m$simulationConditionId, m$preequilibrationConditionId,
                 m$observableId, sep = "\r")
    varying <- names(which(vapply(
      split(m$noiseParameters[numeric_noise], key[numeric_noise]),
      function(v) length(unique(v)) > 1L, logical(1))))
    collapse <- numeric_noise & key %in% varying
    if (any(collapse)) {
      n_parts <- lengths(strsplit(m$noiseParameters[collapse],
                                  ";", fixed = TRUE))
      m$noiseParameters[collapse] <- vapply(
        pmax(n_parts, 1L),
        function(n) paste(rep("1", n), collapse = ";"),
        character(1))
    }
  }

  # Placeholders are observable-specific, so rows of different observables can
  # share one sub-condition. A (simCondId, preeq) group merges when each
  # observable has a single (obsPar, noisePar) tuple, else it splits per tuple.
  obs_lab <- .petab_subcond_label(m$observableParameters)
  noi_lab <- .petab_subcond_label(m$noiseParameters)
  peq_lab <- .petab_subcond_label(m$preequilibrationConditionId)

  m$sub_condition  <- m$simulationConditionId  # default; overwritten below
  obs_subs_by_sub  <- list()
  noi_subs_by_sub  <- list()

  # Group by (simCondId, preeq); decide merge-vs-split per group.
  group_key <- paste(m$simulationConditionId, m$preequilibrationConditionId,
                     sep = "")
  for (gk in unique(group_key)) {
    ix <- which(group_key == gk)
    sc   <- m$simulationConditionId[ix[1L]]
    peq  <- m$preequilibrationConditionId[ix[1L]]
    plab <- peq_lab[ix[1L]]

    # Per observableId in this group: distinct (obsPar, noisePar) tuples.
    by_obs <- split(ix, m$observableId[ix])
    consistent <- all(vapply(by_obs, function(rs) {
      length(unique(paste(m$observableParameters[rs],
                          m$noiseParameters[rs],
                          sep = ""))) <= 1L
    }, logical(1)))

    if (consistent) {
      # ----- merge: one sub-condition for the whole (sc, peq) group -------
      sub_name <- if (nzchar(plab)) paste0(sc, "__", plab) else sc
      m$sub_condition[ix] <- sub_name
      obs_map <- vapply(by_obs, function(rs) m$observableParameters[rs[1L]],
                        character(1))
      noi_map <- vapply(by_obs, function(rs) m$noiseParameters[rs[1L]],
                        character(1))
      obs_map <- obs_map[nzchar(obs_map)]
      noi_map <- noi_map[nzchar(noi_map)]
      if (length(obs_map)) obs_subs_by_sub[[sub_name]] <- obs_map
      if (length(noi_map)) noi_subs_by_sub[[sub_name]] <- noi_map
      next
    }

    # ----- split: one sub-condition per distinct (obsPar, noisePar) tuple -----
    tuple_key <- paste(peq_lab[ix], obs_lab[ix], noi_lab[ix], sep = "")
    keys_here <- unique(tuple_key)
    suffixes <- vapply(keys_here, function(k) {
      jx <- ix[tuple_key == k][1L]
      parts <- c(if (nzchar(peq_lab[jx])) peq_lab[jx],
                 if (nzchar(obs_lab[jx])) obs_lab[jx],
                 if (nzchar(noi_lab[jx])) noi_lab[jx])
      paste(parts, collapse = "_")
    }, character(1))
    if (anyDuplicated(suffixes)) {
      tags <- vapply(keys_here, function(k)
        substr(digest::digest(k, algo = "md5"), 1L, 6L), character(1))
      suffixes <- paste(suffixes, tags, sep = "_")
    }
    for (i in seq_along(keys_here)) {
      jx <- ix[tuple_key == keys_here[i]]
      sub_name <- paste0(sc, "__", suffixes[i])
      m$sub_condition[jx] <- sub_name
      o <- m$observableParameters[jx[1L]]
      n <- m$noiseParameters[jx[1L]]
      if (nzchar(o)) obs_subs_by_sub[[sub_name]] <- c("*" = o)
      if (nzchar(n)) noi_subs_by_sub[[sub_name]] <- c("*" = n)
    }
  }

  # One row per sub-condition; the substitution maps live in the attributes
  # "obs_subs" / "noi_subs" so the data.frame stays serialisable.
  sub_cond_map <- unique(m[, c("simulationConditionId",
                               "preequilibrationConditionId",
                               "sub_condition")])
  rownames(sub_cond_map) <- NULL
  if (anyDuplicated(sub_cond_map$sub_condition))
    stop("Internal: sub-condition map has duplicate sub_condition rows.")
  attr(sub_cond_map, "obs_subs") <- obs_subs_by_sub
  attr(sub_cond_map, "noi_subs") <- noi_subs_by_sub

  # PEtab measurements are linear and the observable transformation belongs to
  # the likelihood. Transforming the data here and the observable in `g` puts
  # the residual on the chosen scale and keeps normL2's fast path.
  trafo_per_row <- obs_meta$obs_trafo[m$observableId]
  val <- as.numeric(m$measurement)
  log_idx   <- which(trafo_per_row == "log")
  log10_idx <- which(trafo_per_row == "log10")
  if (length(log_idx))   val[log_idx]   <- log(val[log_idx])
  if (length(log10_idx)) val[log10_idx] <- log10(val[log10_idx])

  # PEtab's `time = inf` asks for the steady state. The row is placed at
  # .petab_ss_time and mapped back to `inf` on export.
  meas_time <- as.numeric(m$time)
  meas_time[is.infinite(meas_time) & meas_time > 0] <- .petab_ss_time

  data_df <- data.frame(
    name      = m$observableId,
    time      = meas_time,
    value     = val,
    sigma     = as.numeric(sigma),
    condition = m$sub_condition,
    stringsAsFactors = FALSE
  )

  list(data = data_df,
       sub_cond_map = sub_cond_map,
       has_preeq = any(sub_cond_map$preequilibrationConditionId != ""))
}


# Internal: substitute the placeholders `<prefix>K_<id>` in `formula` with the
# K-th entry of the ";"-separated `repls_str`.
.petab_substitute_param_string <- function(formula, repls_str, prefix) {
  if (length(formula) != 1L) {
    return(vapply(formula, .petab_substitute_param_string, character(1),
                  repls_str = repls_str, prefix = prefix))
  }
  if (is.na(formula) || !nzchar(formula)) return(formula)
  if (is.na(repls_str) || !nzchar(repls_str)) return(formula)
  parts <- trimws(strsplit(repls_str, ";", fixed = TRUE)[[1]])
  for (k in seq_along(parts)) {
    pat <- sprintf("\\b%s%d_[A-Za-z][A-Za-z0-9_]*\\b", prefix, k)
    formula <- gsub(pat, parts[k], formula, perl = TRUE)
  }
  formula
}


# Internal: evaluate a formula string to a finite constant in baseenv(), or
# NA_real_ when it is not numeric, not finite or has a free symbol.
.petab_eval_constant <- function(formula) {
  if (is.na(formula) || !nzchar(formula)) return(NA_real_)
  num <- suppressWarnings(as.numeric(formula))
  if (!is.na(num)) return(num)
  v <- tryCatch(eval(parse(text = formula), envir = baseenv()),
                error = function(e) NULL)
  if (is.numeric(v) && length(v) == 1L && is.finite(v)) v else NA_real_
}


# Internal: TRUE if expression `e` is 0 for all values of its other symbols once
# those in `zero` are 0, judged at two fixed positive points. An expression that
# does not evaluate counts as nonzero.
.petab_vanishes <- function(e, zero = character(0)) {
  if (is.null(e)) return(FALSE)
  syms <- all.vars(e)
  env  <- new.env(parent = baseenv())
  env$piecewise <- function(...) {
    a <- list(...); n <- length(a)
    for (i in seq_len(n %/% 2L)) if (isTRUE(a[[2L * i]])) return(a[[2L * i - 1L]])
    if (n %% 2L) a[[n]] else 0
  }
  k <- seq_along(syms)
  for (probe in list(exp(sin(k)), exp(cos(k)))) {
    vals <- setNames(as.list(probe), syms)
    vals[intersect(syms, zero)] <- 0
    v <- tryCatch(eval(e, vals, env), error = function(err) NA_real_)
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v != 0) return(FALSE)
  }
  TRUE
}


# Internal: states an equilibration from `init` cannot move. A reaction is idle
# if its rate vanishes with `zero_pars` and `zero` at 0; `zero` and `frozen` take
# part in idle reactions only. Returns list(zero, frozen, idle).
.petab_invariant_states <- function(reactions, init, zero_pars = character(0)) {
  st <- reactions$states
  S  <- reactions$smatrix
  if (is.null(S) || !length(st))
    return(list(zero = character(0), frozen = character(0), idle = logical(0)))
  S  <- suppressWarnings(matrix(as.numeric(S), nrow(S), ncol(S), dimnames = list(NULL, st)))
  involved <- !is.na(S) & S != 0
  parse1 <- function(s) tryCatch(str2lang(s), error = function(err) NULL)
  rates  <- lapply(as.character(reactions$rates), parse1)
  starts0 <- vapply(st, function(s)
    !is.na(init[s]) && .petab_vanishes(parse1(init[[s]]), zero_pars), logical(1))
  idleWith <- function(z) vapply(rates, .petab_vanishes, logical(1), zero = c(zero_pars, z))
  onlyIdle <- function(idle) colSums(involved[!idle, , drop = FALSE]) == 0
  zero <- st[starts0]
  repeat {
    idle <- idleWith(zero)
    keep <- zero[onlyIdle(idle)[zero]]
    if (length(keep) == length(zero)) break
    zero <- keep
  }
  list(zero = zero, frozen = setdiff(st[onlyIdle(idle)], zero), idle = idle)
}


# Internal: `reactions` without the reactions flagged in `drop_rx` and the
# states `drop_st`, with the symbols `zero` set to 0 in the remaining rates.
.petab_reduce_network <- function(reactions, drop_st, drop_rx, zero = character(0)) {
  if (!length(drop_st) && !any(drop_rx)) return(reactions)
  keep_st <- setdiff(reactions$states, drop_st)
  keep_rx <- !drop_rx
  rates <- reactions$rates[keep_rx]
  z <- intersect(zero, unique(unlist(lapply(rates, getSymbols))))
  if (length(z)) rates <- replaceSymbols(z, rep("0", length(z)), rates)
  vol <- reactions$volumes
  cof <- reactions$compartmentOf
  rc  <- reactions$reactionCompartment
  eqnlist(smatrix = reactions$smatrix[keep_rx, keep_st, drop = FALSE],
          states = keep_st, rates = rates,
          volumes = if (!is.null(vol)) vol[intersect(names(vol), keep_st)],
          description = reactions$description[keep_rx],
          compartments = reactions$compartments,
          compartmentOf = if (!is.null(cof)) cof[intersect(names(cof), keep_st)],
          reactionCompartment = if (!is.null(rc)) rc[keep_rx],
          amountStates = intersect(reactions$amountStates, keep_st))
}


## --- core trafo / observation / objective builders --------------------------

# Internal: the per-condition parfn from outer (estimated and fixed) parameters
# to the inner set: state initial values, `inner_pars` and the observable/noise
# placeholders `obs_inner`. Arguments come from the `.petab_parse_*` helpers.
.petab_build_trafo <- function(sub_cond_map, conditions, col_kind, override_cols,
                               inits, sbml_pars, states, inner_pars,
                               pouter_names, fixed, scales,
                               obs_inner = character(),
                               reactions = NULL,
                               events = NULL,
                               compile = TRUE,
                               modelname = "petab_trafo",
                               cores = 1L, outdir = .dmodOutdir(),
                               derivMode = "forward") {

  # Inner side that the trafo must produce per condition.
  inner_targets <- unique(c(states, inner_pars, obs_inner))

  # Baseline: every inner target maps to itself, states to their SBML initial
  # expression. Condition overrides and placeholder substitutions follow, then
  # the parameter-scale chain rule.
  build_default <- function() {
    base <- setNames(inner_targets, inner_targets)
    for (st in intersect(states, names(inits))) {
      v <- inits[[st]]
      base[st] <- as.character(v)
    }
    base
  }

  apply_scale_chain_rule <- function(tr) {
    log_pars   <- names(scales)[scales == "log"   & names(scales) %in% pouter_names]
    log10_pars <- names(scales)[scales == "log10" & names(scales) %in% pouter_names]
    # The replacement has to be parenthesised as a whole: `^` is right
    # associative, so substituting into `K^2` would give `10^(K^2)` instead of
    # `(10^K)^2`.
    if (length(log_pars))   tr <- repar(tr, "x ~ (exp(x))", x = log_pars)
    if (length(log10_pars)) tr <- repar(tr, "x ~ (10^(x))", x = log10_pars)
    tr
  }

  apply_row_overrides <- function(tr, cond_id, scope) {
    # scope "all" applies every override column, "state-only" only init columns,
    # for preequilibration stages whose parameter overrides are in the rates.
    if (!cond_id %in% rownames(conditions)) return(tr)
    for (cn in override_cols) {
      v <- conditions[cond_id, cn]
      if (is.na(v) || is.null(v)) next
      v <- trimws(as.character(v))
      if (v == "") next
      kind <- col_kind[[cn]]
      if (kind == "init") {
        tr[cn] <- v
      } else if (identical(scope, "all")) {
        tr <- repar(tr, paste0(cn, " ~ ", v))
      }
    }
    tr
  }

  apply_petab_param_subs <- function(tr, obs_subs, noi_subs) {
    tr <- .petab_apply_subs_map(tr, obs_subs, prefix = "observableParameter")
    tr <- .petab_apply_subs_map(tr, noi_subs, prefix = "noiseParameter")
    # A placeholder no measurement row fills cannot reach a residual; it is
    # pinned so it does not become an undeclared outer parameter.
    unfilled <- grepl("^(observable|noise)Parameter[0-9]+_", names(tr)) &
                names(tr) == unname(tr)
    tr[unfilled] <- "1"
    tr
  }

  obs_subs_attr <- attr(sub_cond_map, "obs_subs") %||% list()
  noi_subs_attr <- attr(sub_cond_map, "noi_subs") %||% list()

  # One `P()` per group of conditions, so the generator and `cores` run over the
  # whole list. Conditions with one preequilibration share its equilibration
  # model, which is sound because this builder runs once per SBML model.
  plain_tr <- list(); pre_tr <- list(); post_tr <- list()
  eq_of    <- character(0); eq_model <- list(); eq_net <- list()

  for (sub in sub_cond_map$sub_condition) {

    row_idx <- which(sub_cond_map$sub_condition == sub)
    sim     <- sub_cond_map$simulationConditionId[row_idx]
    peq     <- sub_cond_map$preequilibrationConditionId[row_idx]
    obs_subs <- obs_subs_attr[[sub]]
    noi_subs <- noi_subs_attr[[sub]]

    has_peq <- length(peq) && nzchar(peq)

    if (!has_peq) {
      # ----- single-stage Pexpl path (Stage-1 path) ---------------------
      tr <- build_default()
      tr <- apply_row_overrides(tr, sim, "all")
      tr <- apply_petab_param_subs(tr, obs_subs, noi_subs)
      tr <- apply_scale_chain_rule(tr)
      plain_tr[[sub]] <- tr
      next
    }

    # ----- preequilibration ---------------------------------------------
    # 1. peq parameter overrides go straight into the rates. The scale chain
    #    rule is applied first, since overrides are linear-scale literals.
    peq_reactions <- reactions
    peq_param_subs <- list()
    peq_state_subs <- list()
    if (peq %in% rownames(conditions)) {
      for (cn in override_cols) {
        v <- conditions[peq, cn]
        if (is.na(v) || is.null(v)) next
        v <- trimws(as.character(v))
        if (v == "") next
        kind <- col_kind[[cn]]
        if (kind == "init") peq_state_subs[[cn]] <- v
        else                 peq_param_subs[[cn]] <- v
      }
    }
    if (length(peq_param_subs)) {
      peq_reactions$rates <- replaceSymbols(
        names(peq_param_subs), unlist(peq_param_subs), peq_reactions$rates)
    }

    # A steady state needs an autonomous system, so a time-dependent input is
    # frozen at its value at the equilibration start; otherwise `time` reaches
    # the integrator as an unset parameter.
    peq_start <- 0
    if (any(grepl("\\btime\\b", peq_reactions$rates)))
      peq_reactions$rates <- replaceSymbols("time", as.character(peq_start),
                                            peq_reactions$rates)

    # 2. p_pre: outer -> c(state inits with peq state overrides applied,
    #    inner pars identity, observable/noise placeholders pre-substituted).
    tr_pre <- build_default()
    for (st in names(peq_state_subs)) tr_pre[st] <- peq_state_subs[[st]]

    # States the equilibration cannot move are undetermined by the steady-state
    # equations; they leave the network and pass to p_post with their initial
    # values. SBML constants outside parameters.tsv at 0 idle their reactions.
    if (is.null(eq_net[[peq]])) {
      zero_pars <- setdiff(names(fixed)[!is.na(fixed) & fixed == 0], names(scales))
      inv <- .petab_invariant_states(peq_reactions, tr_pre[states], zero_pars)
      eq_net[[peq]] <- .petab_reduce_network(peq_reactions, c(inv$zero, inv$frozen),
                                             inv$idle, c(zero_pars, inv$zero))
    }
    peq_reactions <- eq_net[[peq]]

    tr_pre <- apply_petab_param_subs(tr_pre, obs_subs, noi_subs)
    tr_pre <- apply_scale_chain_rule(tr_pre)

    # Totals of the conserved quantities from the species initials, over the
    # network without its structurally zero states, as Pimpl derives them.
    peq_totals <- if (length(peq_reactions$states))
      getTotals(.zeroStatesFromSmatrix(peq_reactions)$eqnlist)
    for (tn in names(peq_totals)) {
      sp <- getSymbols(peq_totals[[tn]])
      if (!all(sp %in% names(tr_pre)))
        stop("PEtab: the preequilibration total ", tn, " needs initial values for ",
             paste(setdiff(sp, names(tr_pre)), collapse = ", "), ".", call. = FALSE)
      tr_pre[tn] <- replaceSymbols(sp, paste0("(", tr_pre[sp], ")"), peq_totals[[tn]])
    }

    pre_tr[[sub]] <- tr_pre

    # 3. p_eq: stable steady state of peq_reactions, started from the initial
    #    state. Parameters and placeholders pass through. Shared by all
    #    sub-conditions with this preequilibration condition.
    eq_of[sub] <- peq
    if (is.null(eq_model[[peq]]) && length(peq_reactions$states))
      eq_model[[peq]] <- Pimpl(peq_reactions, condition = NULL,
                               compile = compile, outdir = outdir,
                               modelname = paste(modelname,
                                                 sanitizeConditions(peq),
                                                 "eq", sep = "_"))

    # 4. p_post: states from p_eq, then the simulation row's overrides, which
    #    may change parameters and re-override states.
    tr_post <- setNames(inner_targets, inner_targets)
    post_tr[[sub]] <- apply_row_overrides(tr_post, sim, "all")
  }

  parfns <- list()
  if (length(plain_tr))
    parfns[[length(parfns) + 1L]] <-
      P(plain_tr, compile = compile, outdir = outdir, cores = cores,
        modelname = modelname, derivMode = derivMode)
  for (k in unique(eq_of)) {
    subs   <- names(eq_of)[eq_of == k]
    tag    <- sanitizeConditions(k)
    p_pre  <- P(pre_tr[subs],  compile = compile, outdir = outdir, cores = cores,
                modelname = paste(modelname, tag, "pre",  sep = "_"),
                derivMode = derivMode)
    p_post <- P(post_tr[subs], compile = compile, outdir = outdir, cores = cores,
                modelname = paste(modelname, tag, "post", sep = "_"),
                derivMode = derivMode)
    parfns[[length(parfns) + 1L]] <- if (is.null(eq_model[[k]])) p_post * p_pre
      else p_post * .fnWithConditions(eq_model[[k]], subs) * p_pre
  }

  Reduce(`+`, parfns)
}


# Internal: override every inner target `<prefix>K_*` of a trafo with the K-th
# entry of the ";"-separated `repls`. A non-NULL `obs_id` restricts this to the
# placeholders of that observable.
.petab_apply_param_substitution <- function(tr, repls, prefix, obs_id = NULL) {
  parts <- trimws(strsplit(repls, ";", fixed = TRUE)[[1]])
  inner <- names(tr)
  pat <- if (is.null(obs_id))
           paste0("^", prefix, "([0-9]+)_(.+)$")
         else
           paste0("^", prefix, "([0-9]+)_",
                  gsub("([][\\\\.|()^$*+?{}])", "\\\\\\1", obs_id), "$")
  for (sym in inner) {
    m <- regmatches(sym, regexec(pat, sym))[[1]]
    if (length(m) >= 2L) {
      k <- as.integer(m[2])
      if (k >= 1L && k <= length(parts)) tr[sym] <- parts[k]
    }
  }
  tr
}


# Internal: apply a per-observable substitution map to a trafo. `subs` maps an
# observableId, or "*" for every placeholder, to a ";"-separated replacement.
.petab_apply_subs_map <- function(tr, subs, prefix) {
  if (length(subs) == 0L) return(tr)
  for (key in names(subs)) {
    repl <- subs[[key]]
    if (!nzchar(repl)) next
    tr <- .petab_apply_param_substitution(
      tr, repl, prefix,
      obs_id = if (identical(key, "*")) NULL else key)
  }
  tr
}


# `g` is evaluated on the whole time grid, also where a species is still empty.
# The floor keeps `log` and `log10` finite there, below any measured level.
.PETAB_LOG_FLOOR <- "1e-15"


# States a noise formula names. Relative noise is often written against the
# species (`sigma * Cer`) rather than against the observable, and then the
# observation function has to pass that species through to the error model.
.petab_noise_states <- function(obs_meta, reactions) {
  syms <- unique(unlist(lapply(obs_meta$noise, function(f)
    if (is.na(suppressWarnings(as.numeric(f)))) getSymbols(f) else character(0))))
  intersect(syms, reactions$states)
}


# Build the observation function `g`. A `log`/`log10` observable is wrapped at
# construction, so with the transformed data the residual is on that scale and
# normL2's fast path holds.
.petab_build_observation_fn <- function(obs, obs_trafo, reactions,
                                        compile = TRUE,
                                        modelname = "petab_obs",
                                        outdir = .dmodOutdir(),
                                        keepStates = character(),
                                        derivMode = "forward") {
  obs_eqn <- mapply(function(formula, trafo) {
    if (identical(trafo, "log"))   sprintf("log(%s + %s)",   formula, .PETAB_LOG_FLOOR)
    else if (identical(trafo, "log10")) sprintf("log10(%s + %s)", formula, .PETAB_LOG_FLOOR)
    else formula
  }, obs, obs_trafo[names(obs)], SIMPLIFY = TRUE, USE.NAMES = TRUE)
  # A PEtab model may declare no species at all; the observables are then
  # pure functions of parameters, so Y gets its symbols directly.
  stateless <- length(reactions$states) == 0L
  Y(g          = as.eqnvec(obs_eqn),
    f          = if (stateless) NULL else as.eqnvec(reactions),
    parameters = if (stateless)
                   setdiff(unique(unlist(lapply(obs_eqn, getSymbols))), "time")
                 else NULL,
    attachInput = length(keepStates) > 0L,
    compile    = compile,
    modelname  = modelname,
    outdir     = outdir,
    derivMode  = derivMode)
}


# Build the error model, or NULL when every noise formula is constant after
# substitution. One condition-free Y over the placeholder-form noise formulas
# suffices, since the trafo binds the placeholders per sub-condition.
.petab_build_error_fn <- function(obs_meta, sub_cond_map, reactions,
                                  compile = TRUE, modelname = "petab_err",
                                  outdir = .dmodOutdir(), derivMode = "forward") {

  obs_subs <- attr(sub_cond_map, "obs_subs") %||% list()
  noi_subs <- attr(sub_cond_map, "noi_subs") %||% list()

  # A noise formula still symbolic after per-row substitution needs the error
  # model; a constant one comes from the data column.
  pick_str <- function(map, sub, obsId) {
    m <- map[[sub]]
    if (is.null(m)) return("")
    if (obsId %in% names(m)) return(unname(m[obsId]))
    if ("*" %in% names(m))   return(unname(m["*"]))
    ""
  }
  any_symbolic <- any(vapply(sub_cond_map$sub_condition, function(sub) {
    any(vapply(names(obs_meta$noise), function(obsId) {
      f <- obs_meta$noise[[obsId]]
      f <- .petab_substitute_param_string(
        f, pick_str(obs_subs, sub, obsId), "observableParameter")
      f <- .petab_substitute_param_string(
        f, pick_str(noi_subs, sub, obsId), "noiseParameter")
      is.na(.petab_eval_constant(f))
    }, logical(1)))
  }, logical(1)))
  if (!any_symbolic) return(NULL)

  # States and untransformed observable formulas, so a noise formula may name an
  # observable; PEtab sigma already lives on the transformed scale.
  obs_eqnvec <- as.eqnvec(setNames(unname(unlist(obs_meta$obs)),
                                   names(obs_meta$obs)))
  reactions_eqnvec <- if (length(reactions$states)) as.eqnvec(reactions) else NULL

  Y(g            = as.eqnvec(obs_meta$noise),
    f            = c(reactions_eqnvec, obs_eqnvec),
    states       = c(names(obs_eqnvec), .petab_noise_states(obs_meta, reactions)),
    attachInput  = FALSE,
    compile      = compile,
    modelname    = modelname,
    outdir       = outdir,
    derivMode    = derivMode)
}


# Build the odemodel and Xs() for the chosen backend. `compile` is forwarded so
# linking can wait for one batched compile(); `events` come from importSbml(),
# `options` / `optionsSens` reach Xs() untouched.
.petab_build_odemodel <- function(reactions, backend,
                                  modelname = "petab_model",
                                  compile = TRUE,
                                  events = NULL,
                                  options = NULL, optionsSens = NULL,
                                  sparse = NULL,
                                  deriv = TRUE, derivMode = "forward",
                                  outdir = .dmodOutdir()) {
  # No species means no dynamics: Xt() supplies the time axis and the
  # observables are evaluated from parameters alone.
  if (length(reactions$states) == 0L)
    return(list(odemodel = NULL, x = Xt()))

  # A PEtab experiment may start after t = 0, and its initial values belong at
  # that start. The backends otherwise force 0 into the integration grid,
  # which would apply the initial values there instead.
  args <- list(reactions, modelname = modelname, backend = backend,
               events = events, compile = compile, includeTimeZero = FALSE,
               deriv = deriv, derivMode = derivMode, outdir = outdir)
  if (!is.null(sparse) && backend != "deSolve") args$sparse <- sparse
  m <- do.call(odemodel, args)
  opts <- list(m)
  if (!is.null(options))     opts$options     <- options
  if (!is.null(optionsSens)) opts$optionsSens <- optionsSens
  list(odemodel = m, x = do.call(Xs, opts))
}


# Internal: give a condition-free obsfn explicit `conds`, all routing to one
# kernel, so per-model error functions can be combined with `+.fn`.
.obsfn_with_conditions <- function(fn, conds) {
  if (is.null(fn) || !length(conds)) return(NULL)
  if (!inherits(fn, "obsfn")) stop("not an obsfn")
  .fnWithConditions(fn, conds)
}


# Internal: promote event targets that are not states to states with zero rate,
# since the solvers apply events to states only. The original value becomes the
# initial condition, so every reference resolves to the same symbol.
.petab_promote_event_targets <- function(sbml, events) {
  if (is.null(events)) return(sbml)
  r <- sbml$reactions
  extra <- setdiff(unique(as.character(events$var)), r$states)
  if (!length(extra)) return(sbml)

  # The SBML default becomes the state's initial value. The parameter stays in
  # place: other initial expressions (an amount species is conc * volume) still
  # reference it as an outer symbol.
  init <- vapply(extra, function(nm) {
    if (nm %in% names(sbml$pars))
      return(format(unname(sbml$pars[nm]), digits = 17))
    vol <- r$compartments[[nm]]$volume
    if (!is.null(vol)) return(as.character(vol))
    "0"
  }, character(1))

  states <- c(r$states, extra)
  smat <- r$smatrix
  if (!is.null(smat)) {
    smat <- cbind(smat, matrix(NA_real_, nrow(smat), length(extra)))
    colnames(smat) <- states
  }

  sbml$reactions <- eqnlist(smatrix = smat, states = states, rates = r$rates,
                            description = r$description,
                            compartments = r$compartments,
                            compartmentOf = r$compartmentOf,
                            reactionCompartment = r$reactionCompartment,
                            amountStates = r$amountStates)
  sbml$inits <- c(sbml$inits, setNames(unname(init), extra))
  sbml$promoted <- extra
  sbml
}


# Internal: mid-run condition switches as a dMod eventlist. A period boundary
# assigns every target its condition overrides at the switch time; values may
# reference states, which is what PEtab's `A + 5` style overrides need.
.petab_switch_events <- function(switches, conditions_df) {
  if (!length(switches) || is.null(conditions_df)) return(NULL)
  ev <- NULL
  for (sw in switches) {
    for (i in seq_len(nrow(sw))) {
      row <- conditions_df[conditions_df$conditionId == sw$conditionId[i], ,
                           drop = FALSE]
      if (!nrow(row)) next
      for (tid in setdiff(colnames(row), "conditionId")) {
        v <- row[[tid]][1L]
        if (is.na(v) || !nzchar(trimws(as.character(v)))) next
        ev <- addEvent(ev, var = tid, time = sw$time[i],
                       value = as.character(v), method = "replace")
      }
    }
  }
  ev
}


# Internal: the per-model dMod pieces. `meas_m` is the model's slice of the
# measurements; `obs_meta_full` is reduced to the observables it uses.
.petab_build_model_pieces <- function(sbml, meas_m, conditions_df,
                                      obs_meta_full, param_meta,
                                      modelname, backend, compile,
                                      sub_cond_prefix = "",
                                      start_times = NULL,
                                      switches = NULL,
                                      options = NULL, optionsSens = NULL,
                                      sparse = NULL,
                                      deriv = TRUE, derivMode = "forward", cores = 1L,
                                      outdir = .dmodOutdir()) {

  # `importSbml` renames ids that R cannot parse or C++ reserves; the PEtab
  # tables are renamed with the same map so they keep matching the model.
  if (length(sbml$renamed)) {
    ren <- function(x) {
      if (!length(x)) return(x)
      v <- .renameIds(unlist(x, use.names = FALSE), sbml$renamed)
      if (is.list(x)) setNames(as.list(v), names(x)) else setNames(v, names(x))
    }
    obs_meta_full$obs   <- ren(obs_meta_full$obs)
    obs_meta_full$noise <- ren(obs_meta_full$noise)
    hit <- match(colnames(conditions_df), names(sbml$renamed))
    colnames(conditions_df)[!is.na(hit)] <- unname(sbml$renamed[hit[!is.na(hit)]])
    for (cn in intersect(c("observableParameters", "noiseParameters"),
                         colnames(meas_m)))
      meas_m[[cn]] <- ren(as.character(meas_m[[cn]]))
  }

  obs_used <- unique(meas_m$observableId)
  obs_meta <- list(
    obs        = obs_meta_full$obs[obs_used],
    obs_trafo  = obs_meta_full$obs_trafo[obs_used],
    noise      = obs_meta_full$noise[obs_used],
    noise_dist = obs_meta_full$noise_dist[obs_used]
  )

  # An observable may read a symbol that an SBML assignment rule defines. The
  # rule was inlined into the rates and its left hand side is no longer a
  # parameter, so the observable and noise formulas need it inlined too.
  rules <- sbml$assignmentRules
  if (length(rules)) {
    lhs <- names(rules)
    rhs <- paste0("(", unlist(rules, use.names = FALSE), ")")
    inline <- function(f) {
      if (!length(f)) return(f)
      v <- unlist(f, use.names = FALSE)
      for (it in seq_len(length(lhs) + 1L)) {
        new_v <- replaceSymbols(lhs, rhs, v)
        if (identical(new_v, v)) break
        v <- new_v
      }
      # keep the container the caller handed in; `as.eqnvec()` later on has no
      # method for a list
      if (is.list(f)) setNames(as.list(v), names(f)) else setNames(v, names(f))
    }
    obs_meta$obs   <- inline(obs_meta$obs)
    obs_meta$noise <- inline(obs_meta$noise)
  }

  # SBML events and PEtab period switches share the eventlist, with non-state
  # targets already promoted. The switch comes first so model events see the
  # switched state, as PEtab prescribes.
  sw_events <- .petab_switch_events(switches, conditions_df)
  all_events <- if (is.null(sw_events)) sbml$events
                else if (is.null(sbml$events)) sw_events
                else rbind(sw_events, as.eventlist(sbml$events))
  sbml <- .petab_promote_event_targets(sbml, all_events)

  states    <- sbml$reactions$states
  rate_syms <- unique(unlist(lapply(sbml$reactions$rates, function(r)
                getSymbols(r, exclude = c(states, "time")))))
  # An event may name a parameter no rate mentions, typically its own trigger
  # time. The model still needs it, so it belongs to the inner set.
  ev_syms <- if (!is.null(all_events) && nrow(all_events))
    unique(unlist(lapply(unlist(all_events[intersect(c("time", "value", "root"),
                                                     colnames(all_events))]),
                         function(z) getSymbols(as.character(z),
                                                exclude = c(states, "time")))))
    else character(0)
  inner_pars <- setdiff(unique(c(rate_syms, ev_syms)), states)

  obs_syms   <- unique(unlist(lapply(obs_meta$obs, function(f)
                  getSymbols(f, exclude = c(states, "time")))))
  noise_syms <- unique(unlist(lapply(obs_meta$noise, function(f) {
                  if (is.na(suppressWarnings(as.numeric(f)))) getSymbols(f)
                  else character(0)
                })))
  obs_inner <- unique(c(obs_syms, noise_syms))
  # A noise formula may name an observable (relative noise, `obs*scale`).
  # The error model resolves it from the observable equations, so it is not
  # a parameter the problem has to declare.
  obs_inner <- setdiff(obs_inner,
                       c(states, inner_pars, names(obs_meta$obs), "time"))

  cond_info <- .petab_parse_conditions(
                 conditions_df,
                 sbml_states       = states,
                 sbml_compartments = names(sbml$reactions$compartments %||% list()),
                 sbml_pars         = names(sbml$pars),
                 obs_inner         = obs_inner)
  meas_info <- .petab_parse_measurements(meas_m, obs_meta)

  # Sub-condition keys are not disjoint across models, so they are prefixed with
  # the modelId before per-model trafos and datalists are summed.
  if (nzchar(sub_cond_prefix)) {
    rn <- function(x) paste0(sub_cond_prefix, x)
    meas_info$sub_cond_map$sub_condition <- rn(meas_info$sub_cond_map$sub_condition)
    meas_info$data$condition             <- rn(meas_info$data$condition)
    rename_attr_keys <- function(scm, key) {
      a <- attr(scm, key)
      if (length(a)) names(a) <- rn(names(a))
      attr(scm, key) <- a
      scm
    }
    meas_info$sub_cond_map <- rename_attr_keys(meas_info$sub_cond_map, "obs_subs")
    meas_info$sub_cond_map <- rename_attr_keys(meas_info$sub_cond_map, "noi_subs")
  }

  pouter_names <- names(param_meta$pouter)
  fixed        <- param_meta$fixed

  # Symbols the model needs that parameters.tsv lacks are fixed at their SBML default.
  init_syms <- unique(unlist(lapply(sbml$inits, function(e) {
    if (is.character(e)) getSymbols(e) else character(0)
  })))
  required <- unique(c(inner_pars, obs_inner, init_syms))
  # Excluded: promoted event targets keep their SBML default since other initial
  # expressions read them, `pi` is a constant, and an observable name in a noise
  # formula is resolved by the error model.
  required <- setdiff(required, c(setdiff(states, sbml$promoted), "time", "pi",
                                  names(obs_meta$obs),
                                  pouter_names, names(fixed)))
  is_petab_placeholder <- grepl("^(observable|noise)Parameter[0-9]+_", required)
  required <- required[!is_petab_placeholder]
  if (length(required)) {
    miss <- setdiff(required, names(sbml$pars))
    if (length(miss))
      stop("Symbol(s) referenced by the model but not in parameters.tsv ",
           "and not in SBML defaults: ",
           paste(miss, collapse = ", "))
    extra_fixed <- as.numeric(sbml$pars[required])
    names(extra_fixed) <- required
    fixed <- c(fixed, extra_fixed)
  }

  ode_pair <- .petab_build_odemodel(sbml$reactions, backend = backend,
                                    modelname = paste0(modelname, "_ode"),
                                    compile = compile,
                                    events = all_events,
                                    options = options,
                                    optionsSens = optionsSens,
                                    sparse = sparse,
                                    deriv = deriv, derivMode = derivMode,
                                    outdir = outdir)
  g <- .petab_build_observation_fn(obs_meta$obs, obs_meta$obs_trafo,
                                   sbml$reactions,
                                   compile = compile,
                                   modelname = paste0(modelname, "_obs"),
                                   outdir = outdir,
                                   keepStates = .petab_noise_states(
                                     obs_meta, sbml$reactions),
                                   derivMode = derivMode)
  e <- .petab_build_error_fn(obs_meta, meas_info$sub_cond_map,
                             sbml$reactions, compile = compile,
                             modelname = paste0(modelname, "_err"),
                             outdir = outdir, derivMode = derivMode)
  p <- .petab_build_trafo(
          sub_cond_map  = meas_info$sub_cond_map,
          conditions    = cond_info$grid,
          col_kind      = cond_info$col_kind,
          override_cols = cond_info$override_cols,
          inits         = sbml$inits,
          sbml_pars     = sbml$pars,
          states        = states,
          inner_pars    = inner_pars,
          pouter_names  = pouter_names,
          fixed         = fixed,
          scales        = param_meta$scales,
          obs_inner     = obs_inner,
          reactions = sbml$reactions,
          events = all_events,
          compile = compile, modelname = paste0(modelname, "_trafo"),
          cores = cores, outdir = outdir, derivMode = derivMode)

  dataList <- as.datalist(meas_info$data, splitBy = "condition")

  # Each sub-condition inherits the simulation start of its condition, so the
  # objective can anchor the time grid where the experiment actually begins.
  t0_sub <- setNames(rep(0, length(dataList)), names(dataList))
  if (length(start_times)) {
    cid <- meas_info$sub_cond_map$simulationConditionId[
             match(names(dataList), meas_info$sub_cond_map$sub_condition)]
    known <- !is.na(cid) & cid %in% names(start_times)
    t0_sub[known] <- unname(start_times[cid[known]])
  }
  attr(dataList, "t0") <- t0_sub

  list(
    odemodel     = ode_pair$odemodel,
    x            = ode_pair$x,
    g            = g, e = e, p = p,
    dataList     = dataList,
    fixed        = fixed,
    obs_meta     = obs_meta,
    sub_cond_map = meas_info$sub_cond_map,
    cond_grid    = cond_info$grid,
    col_kind     = cond_info$col_kind,
    sbml         = sbml,
    sub_conds    = unique(meas_info$sub_cond_map$sub_condition)
  )
}


# Build the objective for `{lin, log, log10} * normal` noise: log scales enter
# through `g` and the transformed data, symbolic sigmas through the error model.
# Laplace and log-normal noise are not implemented.
.petab_build_objective <- function(data, prd, errmodel, obs_meta) {

  if (!all(obs_meta$noise_dist == "normal"))
    stop("noiseDistribution(s) ",
         paste(unique(obs_meta$noise_dist[obs_meta$noise_dist != "normal"]),
               collapse = ", "),
         " are not yet supported. Only `normal` is implemented in Stage 2.")

  if (!all(obs_meta$obs_trafo %in% c("lin", "log", "log10")))
    stop("observableTransformation(s) ",
         paste(unique(obs_meta$obs_trafo[!obs_meta$obs_trafo %in%
                                          c("lin", "log", "log10")]),
               collapse = ", "),
         " are not recognised.")

  # Every condition starts at the start time of its experiment: one grid for
  # all when they share it, else a grid per condition beginning there.
  t0 <- attr(data, "t0")
  if (is.null(t0)) t0 <- setNames(rep(0, length(data)), names(data))
  t0 <- t0[names(data)]
  base_obj <- normL2(data = data, x = prd, errmodel = errmodel,
                     times = if (length(unique(t0)) == 1L) unname(t0[1L]) else as.list(t0))

  # PEtab's likelihood is on linear data, so a log/log10 observable adds the
  # data-only Jacobian term 2 * sum log(y) (log(y ln 10) for log10) to -2 log L.
  # Gradient and Hessian are unchanged.
  jac_offset <- .petab_likelihood_offset(data, obs_meta)
  if (jac_offset == 0) return(base_obj)

  # `sweep` is declared rather than left to `...`: callers read formals() to
  # choose a direction, and an undeclared one silently falls back to forward.
  myfn <- function(..., fixed = NULL, deriv = TRUE, env = NULL,
                   sweep = "forward") {
    out <- base_obj(..., fixed = fixed, deriv = deriv, env = env, sweep = sweep)
    out$value <- out$value + jac_offset
    attr_nm <- "data"
    if (!is.null(attr(out, attr_nm)))
      attr(out, attr_nm) <- attr(out, attr_nm) + jac_offset
    out
  }
  for (a in setdiff(names(attributes(base_obj)), "class"))
    attr(myfn, a) <- attr(base_obj, a)
  class(myfn) <- class(base_obj)
  myfn
}


# Internal: data-side log-likelihood offset for {log, log10} observables. The
# transformed data values are mapped back to linear units first.
.petab_likelihood_offset <- function(data, obs_meta) {
  off <- 0
  for (cn in names(data)) {
    df <- data[[cn]]
    for (obsId in unique(df$name)) {
      trafo <- obs_meta$obs_trafo[[obsId]]
      if (is.null(trafo) || identical(trafo, "lin")) next
      val <- df$value[df$name == obsId]
      y_obs <- if (identical(trafo, "log"))   exp(val)
               else if (identical(trafo, "log10")) 10^val
               else val
      if (identical(trafo, "log"))
        off <- off + 2 * sum(log(y_obs))
      else if (identical(trafo, "log10"))
        off <- off + 2 * sum(log(y_obs * log(10)))
    }
  }
  off
}


## --- public top-level entry point ------------------------------------------

#' Import a PEtab Problem
#'
#' Reads a PEtab v1 or v2 problem, a YAML file with an SBML model and TSV
#' tables, and builds the dMod problem: prediction, observation and error
#' functions, parameter transformation, data and objective. The SBML model is
#' read by [importSbml()].
#'
#' Supported: one problem per YAML file with any number of v2 model files;
#' preequilibration (at most two periods per v2 experiment); observable
#' transformations `lin`, `log` and `log10`; the noise distribution `normal`,
#' and `log-normal` in v2; a v2 mapping table. Priors of the distributions
#' `uniform`, `normal`, `log-normal`, `cauchy`, `chisquare`, `exponential`,
#' `gamma`, `laplace`, `log-laplace`, `log-uniform` and `rayleigh`, also with
#' the `parameterScale` prefix, are added to `obj`, truncated to the parameter
#' bounds. Not supported: model languages other than SBML, other noise
#' distributions, and a plain `normal` or `laplace` prior on a parameter whose
#' scale is not `lin`.
#'
#' @param yamlPath Path to the PEtab YAML manifest.
#' @param backend Required, one of `"deSolve"`, `"cppDE"` or `"Sundials"`,
#'   passed to [odemodel()].
#' @param compile Logical. If `TRUE` (default), the generated sources are
#'   compiled. `FALSE` leaves them uncompiled for inspection.
#' @param cores Number of parallel compilation jobs, passed to [compile()].
#'   Default `1`.
#' @param modelname Base name of the generated files. Default `NULL`, the base
#'   name of the YAML file.
#' @param deriv Logical. If `FALSE`, the model is built without sensitivities;
#'   the objective is then evaluated with `deriv = FALSE`. Default `TRUE`.
#' @param derivMode Derivative directions to build, passed to [odemodel()]:
#'   `"forward"` (default), `"reverse"` or `c("forward", "reverse")`.
#'   `"reverse"` needs `backend = "cppDE"` or `"Sundials"`. First order only.
#' @param options,optionsSens Lists passed to [Xs()]: the solver options of
#'   every solve and their overrides for the solves with sensitivities.
#'   Default `NULL`, the backend's defaults.
#' @param sparse `NULL` (default) lets the cppDE or Sundials backend choose a
#'   sparse (KLU) or dense linear solver; `TRUE` or `FALSE` sets it.
#' @param outdir Directory of the generated sources and shared objects, one
#'   source per condition. Default `getOption("dMod.outdir")`, else the
#'   working directory.
#' @param optionsOde Deprecated name of `options`, accepted with a warning.
#' @return A list of class `"petabproblem"` with components `dataList`,
#'   `reactions`, `odemodel`, `g`, `x`, `p`, `e`, `prd` (the composite
#'   `g * x * p`), `obj`, `bestfit`, `parlower` and `parupper`. `obj` has the
#'   PEtab fixed parameters set; `bestfit` is the PEtab `nominalValue` of each
#'   estimated parameter, with the PEtab scales in its attribute
#'   `"petab_scales"`. The attribute `"petab_meta"` of the list holds what
#'   [exportPEtabObject()] needs, among it the fixed parameters as `fixed`.
#' @seealso [exportPEtabObject()], [readPEtabTables()], [importSbml()],
#'   \code{vignette("PEtab", package = "dMod2")}
#' @examplesIf requireNamespace("reticulate", quietly = TRUE) && requireNamespace("rjson", quietly = TRUE) && requireNamespace("yaml", quietly = TRUE) && reticulate::py_module_available("libsbml")
#' \donttest{
#' yaml <- system.file("extdata/petab_boehm/Boehm.yaml", package = "dMod2")
#' petab <- importPEtab(yaml, backend = "cppDE", cores = 1, outdir = tempdir())
#' petab
#'
#' # The objective at the nominal values
#' petab$obj(petab$bestfit)$value
#'
#' # The prediction needs the fixed parameters
#' fixed <- attr(petab, "petab_meta")$fixed
#' times <- seq(0, 240, length.out = 61)
#' plot(petab$prd(times, c(petab$bestfit, fixed)), petab$dataList)
#'
#' # Write the problem back to PEtab
#' yamlOut <- exportPEtabObject(petab, file.path(tempdir(), "boehm"),
#'                              formatVersion = "1")
#' readPEtabTables(yamlOut)$parameters
#' }
#' @export
importPEtab <- function(yamlPath, backend,
                        compile = TRUE, cores = 1L, modelname = NULL,
                        deriv = TRUE, derivMode = "forward",
                        options = NULL, optionsSens = NULL,
                        sparse = NULL, outdir = .dmodOutdir(), optionsOde = NULL) {

  if (!is.null(optionsOde))
    options <- .optionsOdeAlias(options, list(optionsOde = optionsOde),
                                "importPEtab")

  cores <- as.integer(cores)
  if (length(cores) != 1L || is.na(cores) || cores < 1L)
    stop("`cores` must be a single positive integer.")

  if (missing(backend))
    stop("Argument `backend` is required (one of \"deSolve\", \"cppDE\", ",
         "\"Sundials\").")
  backend <- match.arg(backend, c("deSolve", "cppDE", "Sundials"))
  derivMode <- .matchDerivMode(derivMode, c("forward", "reverse"))
  if ("reverse" %in% derivMode && backend == "deSolve")
    stop("derivMode = \"reverse\" needs backend = 'cppDE' or 'Sundials'; ",
         "the deSolve backend goes forward only.", call. = FALSE)


  yamlPath <- normalizePath(yamlPath, mustWork = TRUE)
  derived <- is.null(modelname)
  if (derived)
    modelname <- sub("\\.ya?ml$", "", basename(yamlPath), ignore.case = TRUE)
  modelname <- gsub("[^A-Za-z0-9_]", "_", modelname)
  # Entry points resolve by name, so two problems sharing a model name serve
  # each other's kernels. A name the caller chose stays untouched; the one
  # derived from the file name steps aside for what is already loaded.
  if (derived) {
    base <- modelname
    i <- 1L
    while (modelname %in% names(getLoadedDLLs())) {
      i <- i + 1L
      modelname <- paste0(base, "_", i)
    }
  }

  tables <- readPEtabTables(yamlPath)
  if (identical(tables$formatVersion, 2L))
    tables <- .petab_v2_normalize_tables(tables)
  else
    tables <- .petab_use_condition_names(tables)

  sbml_paths <- tables$sbmlPaths
  if (is.null(sbml_paths) || !length(sbml_paths))
    sbml_paths <- setNames(tables$sbmlPath, "model")
  default_mid <- names(sbml_paths)[1L]

  meas <- tables$measurements
  if (!"modelId" %in% colnames(meas)) {
    if (length(sbml_paths) > 1L)
      stop("measurements.tsv has no `modelId` column but the YAML declares ",
           length(sbml_paths), " models. Add a `modelId` column referencing ",
           "the appropriate model_files key per row.")
    meas$modelId <- default_mid
  } else {
    nas <- is.na(meas$modelId) | !nzchar(meas$modelId)
    if (any(nas)) {
      if (length(sbml_paths) > 1L)
        stop("measurements.tsv has missing `modelId` values; in a multi-model ",
             "problem every measurement row must specify its modelId.")
      meas$modelId[nas] <- default_mid
    }
    bad <- !meas$modelId %in% names(sbml_paths)
    if (any(bad))
      stop("measurements.tsv references unknown modelId(s) not declared in ",
           "the YAML's model_files: ",
           paste(unique(meas$modelId[bad]), collapse = ", "))
  }
  tables$measurements <- meas

  param_meta <- .petab_parse_parameters(tables$parameters)
  obs_meta   <- .petab_parse_observables(tables$observables)

  # Compilation of the generated source files is deferred to a single batched
  # compile() at the end so `cores` controls native-build concurrency and one
  # g++ invocation per sub-condition is avoided.
  per_model <- list()
  for (mid in names(sbml_paths)) {
    meas_m <- meas[meas$modelId == mid, , drop = FALSE]
    if (nrow(meas_m) == 0L) {
      warning("Model `", mid, "` has no measurements; skipping.",
              call. = FALSE)
      next
    }
    meas_m$modelId <- NULL
    # PEtab assigns its condition targets and table parameters itself, ahead
    # of an SBML initial assignment.
    sbml_m <- importSbml(unname(sbml_paths[mid]),
                         keep = c(setdiff(colnames(tables$conditions),
                                          c("conditionId", "conditionName")),
                                  tables$parameters$parameterId))
    suffix <- if (length(sbml_paths) == 1L) ""
              else paste0("__", gsub("[^A-Za-z0-9_]", "_", mid))
    sc_prefix <- if (length(sbml_paths) == 1L) ""
                 else paste0(gsub("[^A-Za-z0-9_]", "_", mid), "__")
    pieces <- .petab_build_model_pieces(
      sbml             = sbml_m,
      meas_m           = meas_m,
      conditions_df    = tables$conditions,
      obs_meta_full    = obs_meta,
      param_meta       = param_meta,
      modelname        = paste0(modelname, suffix),
      backend          = backend,
      compile          = FALSE,
      sub_cond_prefix  = sc_prefix,
      start_times      = tables$startTimes,
      switches         = tables$switches,
      options          = options,
      optionsSens      = optionsSens,
      sparse           = sparse,
      deriv            = deriv,
      derivMode        = derivMode,
      cores            = cores,
      outdir           = outdir)
    pieces$modelID <- mid
    per_model[[mid]] <- pieces
  }
  if (length(per_model) == 0L)
    stop("No measurements matched any declared model.")

  # For multi-model problems each `prd_M = g_M * x_M * p_M` inherits `p_M`'s
  # sub-conds as conditions, so `Reduce(\`+\`, prd_list)` is disjoint.
  multi_model <- length(per_model) > 1L

  if (!multi_model) {
    pm  <- per_model[[1L]]
    odeobj   <- pm$odemodel
    x        <- pm$x
    g        <- pm$g
    e        <- pm$e
    p        <- pm$p
    dataList <- pm$dataList
    fixed    <- pm$fixed
    sbml     <- pm$sbml
    sub_cond_map <- pm$sub_cond_map
    cond_grid    <- pm$cond_grid
    col_kind     <- pm$col_kind
    obs_meta_used <- pm$obs_meta
    prd <- g * x * p
  } else {
    # Compose per-model prd's, then sum via +.fn (mappings disjoint by
    # construction). Each prd_M's conditions come from p_M (multi-cond),
    # so the result has all sub-conds across all models.
    prd_list <- lapply(per_model, function(pm) pm$g * pm$x * pm$p)
    prd <- Reduce(`+`, prd_list)

    # Error model: each per-model `e_M` is a NULL-condition Y (one Y per
    # model). To combine via `+.fn` each gets explicit conditions = M's
    # sub-conds; +.fn then dispatches each sub-cond to the right kernel.
    e_pieces <- Filter(Negate(is.null),
      lapply(per_model, function(pm)
        if (is.null(pm$e)) NULL
        else .obsfn_with_conditions(pm$e, pm$sub_conds)))
    e <- if (length(e_pieces)) Reduce(`+`, e_pieces) else NULL

    # Combined dataList (sub-cond names are disjoint by construction).
    t0_all   <- unlist(lapply(per_model, function(pm) attr(pm$dataList, "t0")))
    dataList <- do.call(c, lapply(per_model, `[[`, "dataList"))
    class(dataList) <- "datalist"
    attr(dataList, "t0") <- t0_all

    # Combined fixed: union across models. Conflicting defaults (the same
    # symbol declared in multiple SBMLs with different values) are an error;
    # add them under different names if they really mean different things.
    fixed <- Reduce(function(a, b) {
      ov <- intersect(names(a), names(b))
      if (length(ov)) {
        diff <- ov[abs(a[ov] - b[ov]) > 1e-12]
        if (length(diff))
          stop("Conflicting SBML default(s) for symbol(s) `",
               paste(diff, collapse = "`, `"), "` across models.")
      }
      c(a, b[setdiff(names(b), names(a))])
    }, lapply(per_model, `[[`, "fixed"))

    # The returned problem shows one representative model; the per-model pieces
    # sit in `attr(., "petab_meta")$models`.
    odeobj <- per_model[[1L]]$odemodel
    sbml   <- per_model[[1L]]$sbml
    x <- prd  # x is conceptually the underlying ODE predictor; `prd`
              # already absorbs the per-model dispatch, so expose that.
    g <- NULL
    p <- NULL
    sub_cond_map <- do.call(rbind, lapply(per_model, `[[`, "sub_cond_map"))
    cond_grid    <- per_model[[1L]]$cond_grid
    col_kind     <- per_model[[1L]]$col_kind
    # The per-model parses have the SBML assignment rules inlined into their
    # observable and noise formulas, which the global parse does not.
    obs_meta_used <- .petab_merge_obs_meta(lapply(per_model, `[[`, "obs_meta"))
  }

  if (isTRUE(compile)) {
    # One shared object, not one per source: a problem with many conditions
    # writes hundreds of files and would run into R's DLL limit.
    if (is.null(e)) dMod2::compile(prd, output = modelname, cores = cores)
    else            dMod2::compile(prd, e, output = modelname, cores = cores)
  }

  raw_obj <- .petab_build_objective(data = dataList, prd = prd, errmodel = e,
                                    obs_meta = obs_meta_used)

  # PEtab priors become a `prior` term added to the objective. It reads
  # parameter values from c(pars, fixed), so this is independent of which
  # parameters end up in `baked_fixed` below.
  if (!is.null(param_meta$priors))
    raw_obj <- raw_obj + .petab_prior_objective(param_meta$priors)

  # SBML defaults that turn out to be inner targets of `p` are set by the trafo,
  # not fixed parameters. They move to `sbml_only_pars`, so a re-export declares
  # them in the SBML but not in parameters.tsv.
  outer_p <- if (multi_model)
               unique(unlist(lapply(per_model, function(pm) getParameters(pm$p))))
             else
               getParameters(p)
  sbml_only_pars <- fixed[setdiff(names(fixed), outer_p)]
  fixed          <- fixed[intersect(names(fixed), outer_p)]

  # Bake `fixed` into the objective so users only pass the estimated bestfit
  # vector. `fixed = NULL` (the default) lets the closure inject the PEtab
  # fixed values; an explicit `fixed = ...` overrides per call.
  baked_fixed <- fixed
  obj <- function(pars, fixed = NULL, ..., sweep = "forward") {
    if (is.null(fixed)) fixed <- baked_fixed
    raw_obj(pars, fixed = fixed, ..., sweep = sweep)
  }
  # Copy the objfn class and attributes over: mstrust() and profile() dispatch
  # on them to reload the shared object inside a worker.
  for (a in setdiff(names(attributes(raw_obj)), "class"))
    attr(obj, a) <- attr(raw_obj, a)
  class(obj) <- class(raw_obj)

  bestfit <- param_meta$pouter
  attr(bestfit, "petab_scales") <- param_meta$scales[names(bestfit)]

  reactions_top <- if (multi_model)
                     lapply(per_model, function(pm) pm$sbml$reactions)
                   else sbml$reactions

  out <- list(
    dataList  = dataList,
    reactions = reactions_top,
    odemodel  = if (multi_model) lapply(per_model, `[[`, "odemodel") else odeobj,
    g         = if (multi_model) lapply(per_model, `[[`, "g") else g,
    x         = if (multi_model) lapply(per_model, `[[`, "x") else x,
    p         = if (multi_model) lapply(per_model, `[[`, "p") else p,
    e         = e,
    prd       = prd,
    obj       = obj,
    bestfit   = bestfit,
    parlower  = param_meta$lower,
    parupper  = param_meta$upper
  )
  attr(out, "petab_meta") <- list(
    fixed          = fixed,
    startTimes     = tables$startTimes,
    eventsSource   = if (multi_model)
                       do.call(rbind, lapply(per_model, function(pm) pm$sbml$eventsSource))
                     else sbml$eventsSource,
    switches       = tables$switches,
    defaultCondition = tables$defaultCondition,
    conditionIds     = tables$conditionIds,
    promoted       = if (multi_model)
                       unique(unlist(lapply(per_model, function(pm) pm$sbml$promoted)))
                     else sbml$promoted,
    sbml_only_pars = sbml_only_pars,
    sbml_pars      = if (multi_model) NULL else sbml$pars,
    inits          = if (multi_model) lapply(per_model, function(pm) pm$sbml$inits)
                     else sbml$inits,
    modelID        = if (multi_model) names(per_model) else modelname,
    sourceYaml     = yamlPath,
    sub_cond_map   = sub_cond_map,
    obs_meta       = obs_meta_used,
    param_meta     = param_meta,
    # The wide condition grid lets `exportPEtabObject` rebuild conditions.tsv;
    # the overrides otherwise exist only inside the per-condition trafo.
    cond_grid      = cond_grid,
    col_kind       = col_kind,
    # Per-model pieces keyed by modelId, NULL for a single-model problem.
    models         = if (multi_model) per_model else NULL
  )
  class(out) <- "petabproblem"
  out
}


#' Print a PEtab Problem
#' @param x A `petabproblem` returned by [importPEtab()].
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.petabproblem <- function(x, ...) {
  meta  <- attr(x, "petab_meta") %||% list()
  scm   <- meta$sub_cond_map
  obs   <- meta$obs_meta$obs
  fixed <- meta$fixed
  models <- meta$models
  mid_label <- if (length(meta$modelID) > 1L)
                 paste(meta$modelID, collapse = ", ")
               else (meta$modelID %||% "")
  cat("<PEtab problem ", mid_label, ">\n", sep = "")
  if (!is.null(meta$sourceYaml))
    cat("  source:        ", meta$sourceYaml, "\n", sep = "")
  if (!is.null(models))
    cat("  models:        ", length(models), " (",
        paste(names(models), collapse = ", "), ")\n", sep = "")
  if (!is.null(scm))
    cat("  conditions:    ", nrow(scm),
        " (", length(unique(scm$simulationConditionId)),
        " sim, ", length(unique(scm$preequilibrationConditionId[
          scm$preequilibrationConditionId != ""])),
        " preeq)\n", sep = "")
  if (!is.null(obs))
    cat("  observables:   ", paste(names(obs), collapse = ", "), "\n",
        sep = "")
  cat("  measurements:  ", sum(vapply(x$dataList, nrow, 0L)), "\n", sep = "")
  cat("  bestfit (n=", length(x$bestfit), "): ",
      paste(names(x$bestfit), collapse = ", "), "\n", sep = "")
  if (length(fixed) > 0L)
    cat("  fixed   (n=", length(fixed), "): ",
        paste(names(fixed), collapse = ", "), "\n", sep = "")
  if (!is.null(meta$param_meta$priors))
    cat("  priors  (n=", length(meta$param_meta$priors), "): ",
        paste(names(meta$param_meta$priors), collapse = ", "), "\n",
        sep = "")
  invisible(x)
}


## --- exporter --------------------------------------------------------------
## The exporter inverts `.petab_build_trafo` on `getEquations(p)`: scale wraps
## are stripped, then each LHS is classified (see .petab_classify_lhs).

# Internal: TRUE if `e` is the symbol `op_sym`, possibly in redundant
# parentheses as in `10^(K)`.
.petab_is_bare <- function(e, op_sym) {
  if (is.symbol(e) && identical(e, op_sym)) return(TRUE)
  if (is.call(e) && length(e) == 2L &&
      identical(e[[1L]], as.symbol("(")) &&
      .petab_is_bare(e[[2L]], op_sym)) return(TRUE)
  FALSE
}

# Internal: scale encoded by `e` if it is a clean wrap of `op_sym`:
# `exp(op)` -> "log", `10^(op)` or `exp10(op)` -> "log10", otherwise NA.
.petab_wrap_scale <- function(e, op_sym) {
  if (!is.call(e)) return(NA_character_)
  head <- e[[1L]]
  if (length(e) == 2L && .petab_is_bare(e[[2L]], op_sym)) {
    if (identical(head, as.symbol("exp")))   return("log")
    if (identical(head, as.symbol("exp10"))) return("log10")
  }
  if (length(e) == 3L && identical(head, as.symbol("^")) &&
      is.numeric(e[[2L]]) && length(e[[2L]]) == 1L && e[[2L]] == 10 &&
      .petab_is_bare(e[[3L]], op_sym)) return("log10")
  NA_character_
}

# Internal: parameter scale of each outer id read off the trafo: "log" if every
# occurrence is a clean `exp(.)` wrap, "log10" for `10^(.)` / `exp10(.)`, else "lin".
.petab_detect_scales <- function(eqs, ids) {
  seen <- setNames(vector("list", length(ids)), ids)
  walk <- function(e) {
    if (is.call(e)) {
      for (op in ids) {
        sc <- .petab_wrap_scale(e, as.symbol(op))
        if (!is.na(sc)) { seen[[op]] <<- c(seen[[op]], sc); return(invisible()) }
      }
      for (i in seq_along(e)[-1L]) walk(e[[i]])
    } else if (is.symbol(e) && as.character(e) %in% ids) {
      seen[[as.character(e)]] <<- c(seen[[as.character(e)]], "lin")
    }
    invisible()
  }
  for (eqv in eqs) for (rhs in as.character(unclass(eqv))) {
    if (!is.na(suppressWarnings(as.numeric(rhs)))) next
    walk(parse(text = rhs, keep.source = FALSE)[[1L]])
  }
  vapply(ids, function(op) {
    u <- unique(seen[[op]])
    if (length(u) == 1L) u else "lin"
  }, character(1))
}

# Internal: SBML and PEtab math know no `exp10`; rewrite it as `10^(.)`.
.petab_expand_exp10 <- function(e) {
  if (!is.call(e)) return(e)
  for (i in seq_along(e)[-1L]) e[[i]] <- .petab_expand_exp10(e[[i]])
  if (identical(e[[1L]], as.symbol("exp10")) && length(e) == 2L)
    return(call("^", 10, call("(", e[[2L]])))
  e
}

# Internal: stop when an SBML initial assignment or a condition entry refers to
# its own target. `conditions = FALSE` skips the table, where a v2 targetValue
# may read the target's current value.
.petab_check_self_refs <- function(inits, cond_grid, conditions = TRUE) {
  bad <- character(0)
  for (st in names(inits)) {
    v <- inits[[st]]
    if (is.character(v) && is.na(suppressWarnings(as.numeric(v))) &&
        st %in% getSymbols(v))
      bad <- c(bad, sprintf("initialAssignment %s = %s", st, v))
  }
  cols <- if (conditions) setdiff(colnames(cond_grid),
                                  c("conditionId", "conditionName", "condition"))
          else character(0)
  for (tid in cols) {
    v <- as.character(cond_grid[[tid]])
    v <- v[!is.na(v) & nzchar(trimws(v)) &
           is.na(suppressWarnings(as.numeric(v)))]
    for (x in unique(v))
      if (tid %in% getSymbols(x))
        bad <- c(bad, sprintf("condition %s = %s", tid, x))
  }
  if (length(bad))
    stop("PEtab export would write self-referencing assignment(s): ",
         paste(bad, collapse = "; "),
         ". Give the outer parameter an id distinct from the model entity.",
         call. = FALSE)
  invisible(TRUE)
}

# Compensate a log/log10-scaled outer parameter `op` in `expr` for the importer's
# chain rule: a direct `10^(op)` / `exp(op)` exponent stays bare, every other
# occurrence is wrapped in `log10(.)` / `log(.)`, so any expression round-trips.
.petab_compensate_chain_rule <- function(expr, op, scale) {
  op_sym  <- as.symbol(op)
  if (scale == "log10") {
    inv_head     <- as.symbol("log10")
  } else {  # "log"
    inv_head     <- as.symbol("log")
  }
  is_clean_wrap <- function(e) identical(.petab_wrap_scale(e, op_sym), scale)
  walk <- function(e) {
    if (is_clean_wrap(e)) return(op_sym)
    if (is.call(e)) {
      for (i in seq_along(e)[-1L]) e[[i]] <- walk(e[[i]])
      return(e)
    }
    if (is.symbol(e) && identical(e, op_sym))
      return(call(as.character(inv_head), op_sym))
    e
  }
  walk(expr)
}

# Compensate a single RHS string for the importer's chain rule, for every
# log/log10-scaled outer parameter at once. Pure numeric literals pass
# through unchanged.
.petab_strip_param_scale <- function(rhs_str, scales) {
  s <- as.character(rhs_str)
  if (!is.na(suppressWarnings(as.numeric(s)))) return(s)
  expr <- tryCatch(parse(text = s, keep.source = FALSE)[[1L]],
                   error = function(e)
                     stop("Cannot parse trafo RHS `", s, "`: ",
                          conditionMessage(e), call. = FALSE))
  log10_pars <- names(scales)[scales == "log10"]
  log_pars   <- names(scales)[scales == "log"]
  for (op in log10_pars)
    expr <- .petab_compensate_chain_rule(expr, op, "log10")
  for (op in log_pars)
    expr <- .petab_compensate_chain_rule(expr, op, "log")
  expr <- .petab_expand_exp10(expr)
  paste(deparse(expr, width.cutoff = 500L), collapse = "")
}

# Apply the strip to every per-condition eqnvec.
.petab_strip_trafo <- function(eqs, scales) {
  lapply(eqs, function(eqv) {
    nm  <- names(eqv)
    raw <- as.character(unclass(eqv))
    out <- vapply(raw,
                  function(rhs) .petab_strip_param_scale(rhs, scales),
                  character(1))
    names(out) <- nm
    out
  })
}

# Classify one LHS across conditions as "missing", "identity" (RHS == LHS in
# all), "const_numeric", "const_symbolic" (same RHS in all) or "varying".
.petab_classify_lhs <- function(stripped_eqs, lhs, conds) {
  rhs <- vapply(conds, function(c) {
    e <- stripped_eqs[[c]]
    v <- if (lhs %in% names(e)) e[[lhs]] else NA_character_
    if (length(v) == 0L || is.na(v)) NA_character_ else as.character(v)
  }, character(1))
  if (any(is.na(rhs))) return(list(kind = "missing"))
  if (all(rhs == lhs)) return(list(kind = "identity"))
  # Try to collapse closed-form numeric expressions (e.g. "10^0" -> 1) so they
  # land as SBML defaults / initialConcentrations rather than conditions.tsv
  # columns referencing literal arithmetic.
  vals <- vapply(rhs, .petab_eval_constant, numeric(1))
  if (length(unique(rhs)) == 1L) {
    if (!is.na(vals[[1L]]))
      return(list(kind = "const_numeric", value = unname(vals[[1L]])))
    return(list(kind = "const_symbolic", formula = unname(rhs[[1L]])))
  }
  # Per-condition values: if every cell evaluates to a constant, store the
  # numeric values; otherwise pass the raw formulas through.
  if (all(!is.na(vals)))
    return(list(kind = "varying",
                per_cond = setNames(as.character(unname(vals)), conds)))
  list(kind = "varying", per_cond = rhs)
}

# Decompose the per-condition trafo into list(conditions_df, inits,
# sbml_extra_pars), the latter being parameters the SBML has to declare so
# condition columns and collapsed inner pars resolve.
.petab_decompose_trafo <- function(eqs, states, inner_pars, obs_inner,
                                   pouter_names, fixed, scales) {

  conds        <- names(eqs)
  # Strip 10^(...) / exp10(...) / exp(...) wraps: parameters.tsv holds the
  # linear value (v1 adds `parameterScale = log10/log`), so the stripped
  # formula is exact in both versions.
  stripped     <- .petab_strip_trafo(eqs, scales)
  inner_targets <- unique(c(states, inner_pars, obs_inner))

  inits          <- list()
  cond_overrides <- list()
  sbml_extra_pars <- numeric(0)

  declare_extra <- function(nm, val) {
    if (!nm %in% c(pouter_names, names(fixed)))
      sbml_extra_pars[[nm]] <<- val
  }

  for (lhs in inner_targets) {
    cls    <- .petab_classify_lhs(stripped, lhs, conds)
    state  <- lhs %in% states

    switch(cls$kind,
      missing = stop(sprintf(
        "Trafo does not cover %s `%s` for at least one condition. Every state and inner parameter must appear on the trafo's LHS.",
        if (state) "state" else "inner parameter", lhs), call. = FALSE),

      identity = {
        if (state) {
          # A state mapped to a symbol of its own name is set by an outer
          # parameter; it goes out as an initialAssignment, which the importer's
          # chain rule turns back into the scaled form.
          if (lhs %in% pouter_names || lhs %in% names(fixed))
            inits[[lhs]] <- lhs
          # Otherwise the post-loop default 0 applies; an undeclared symbol
          # would already have been raised.
        } else {
          declare_extra(lhs, 1)
        }
      },

      const_numeric = {
        if (state) inits[[lhs]] <- cls$value
        else       declare_extra(lhs, cls$value)
      },

      const_symbolic = {
        cond_overrides[[lhs]] <- setNames(rep(cls$formula, length(conds)),
                                          conds)
        if (state) {
          # PEtab sets the initial value from the condition table; the SBML
          # species keeps a numeric placeholder.
          inits[[lhs]] <- 0
        } else {
          # Constant-across-conditions inner_par mapping. Emit as a
          # conditions.tsv column with the same value in every row, and
          # register the inner par as an SBML SId (placeholder default).
          declare_extra(lhs, 1)
        }
      },

      varying = {
        cond_overrides[[lhs]] <- cls$per_cond[conds]
        if (state) {
          # SBML default for a per-condition state init: a placeholder; the
          # importer reads the per-row override anyway.
          v0 <- suppressWarnings(as.numeric(cls$per_cond[[1L]]))
          inits[[lhs]] <- if (!is.na(v0)) v0 else 0
        } else {
          declare_extra(lhs, 1)
        }
      }
    )

    # Free-symbol sanity check on the post-strip RHS.
    if (cls$kind %in% c("const_symbolic", "varying")) {
      rhs_strs <- if (cls$kind == "varying") cls$per_cond else cls$formula
      free <- unique(unlist(lapply(rhs_strs, function(s) {
        if (!is.na(suppressWarnings(as.numeric(s)))) character(0)
        else getSymbols(s, exclude = c(states, "time"))
      })))
      declared <- c(pouter_names, names(fixed), names(sbml_extra_pars))
      miss <- setdiff(free, declared)
      if (length(miss))
        stop(sprintf(
          "Trafo RHS for `%s` references undeclared symbol(s): %s. Add to `pouter` or `fixed`.",
          lhs, paste(miss, collapse = ", ")), call. = FALSE)
    }
  }

  # Every state must have an init recorded (default 0).
  for (st in setdiff(states, names(inits))) inits[[st]] <- 0

  conditions_df <- data.frame(conditionId      = conds,
                              row.names        = conds,
                              stringsAsFactors = FALSE)
  for (cn in names(cond_overrides))
    conditions_df[[cn]] <- unname(cond_overrides[[cn]])

  list(conditions_df   = conditions_df,
       inits           = inits,
       sbml_extra_pars = sbml_extra_pars)
}


#' Export a dMod Problem to PEtab
#'
#' Writes a PEtab v1 or v2 problem to `dir`: the parameter, observable,
#' condition and measurement tables, for v2 also the experiment table, an SBML
#' model and a YAML file. Parameter scales, fixed values and per-condition
#' values are read off the parameter transformation `p`.
#'
#' @section Limitations:
#' Pre-equilibration cannot be expressed through the trafo `p` alone, use
#' [exportPEtabObject()] for those problems. Fixed parameters are always
#' written with `parameterScale = "lin"`.
#'
#' @param data A [datalist] (or a list of data.frames keyed by condition,
#'   each with `name`, `time`, `value` columns; or a long-format data.frame
#'   that [as.datalist()] accepts). Only used for the measurements;
#'   `attr(data, "condition.grid")` is ignored.
#' @param reactions An [eqnlist] describing the ODE network.
#' @param observables Observable formulas keyed by observableId. Accepts a
#'   named character vector, an [eqnvec], or an observation function produced
#'   by [Y()] (formulas read from `attr(observables, "equations")`).
#' @param p A `parfn` produced by [P()]. Required: the symbolic trafo is
#'   the source of truth for parameters.tsv, conditions.tsv, and
#'   `<initialAssignment>` formulas.
#' @param pouter Named numeric vector of estimated outer parameters, on the
#'   chosen `parameterScale`. Names become `parameterId`s in parameters.tsv.
#' @param errors Noise formulas keyed by observableId. Accepts a named
#'   character vector, [eqnvec], or a Y-built error function. If `NULL`,
#'   defaults to `"1"` per observable, or to the column `sigma` of the data
#'   where it is given.
#' @param lower,upper Named numeric vectors of bounds, on the same scale as
#'   `pouter`. Default `NULL`: five decades around `pouter` for log and log10
#'   parameters, `-Inf` and `Inf` for linear ones.
#' @param fixed Named numeric vector of parameters that are not estimated, on
#'   the linear scale. Default `NULL`.
#' @param parameterScale `NULL` (default) reads the scale off `p`: a
#'   parameter entering only as `exp(X)` is `"log"`, one entering only as
#'   `10^(X)` or `exp10(X)` is `"log10"`, anything else `"lin"`. Otherwise a
#'   scalar `"lin"`/`"log"`/`"log10"` (broadcast to all names in `pouter`) or
#'   a named character vector keyed by parameterId.
#' @param observableTransformation `"lin"` (default), `"log"` or `"log10"`,
#'   a single value or one per observableId.
#' @param noiseDistribution `"normal"` (default), `"laplace"` or
#'   `"log-normal"`, a single value or one per observableId. For v2,
#'   `"laplace"` with a log transformation is written as `"log-laplace"`.
#' @param modelID SBML model identifier. Default `"dMod_export"`.
#' @param dir Output directory, created if missing.
#' @param formatVersion `"2.0.0"` (default) or `"1"`.
#' @param overwrite Logical, whether existing files in `dir` are overwritten.
#'   Default `FALSE`.
#' @return Path to the written YAML manifest, invisibly.
#' @details An outer parameter named like a state (e.g. `A ~ exp(A)`) cannot
#'   share the species id in SBML. It is written as `init_<state>` and the
#'   condition table sets the species to it. An outer parameter named
#'   like an inner parameter keeps its name when the trafo is a pure scale
#'   wrap; for any other mapping it is written as `<name>_outer`.
#' @seealso [exportPEtabObject()] for a problem returned by [importPEtab()],
#'   [importPEtab()].
#' @examplesIf requireNamespace("reticulate", quietly = TRUE) && requireNamespace("rjson", quietly = TRUE) && requireNamespace("yaml", quietly = TRUE) && reticulate::py_module_available("libsbml")
#' f <- addReaction(eqnlist(), from = "A", to = "", rate = "k*A")
#' p <- P(eqnvec(A = "exp(logA)", k = "exp(logk)"), condition = "C1",
#'        modelname = "petab_export_p", compile = FALSE, outdir = tempdir())
#' data <- datalist(C1 = data.frame(name = "y", time = 0:4, sigma = 0.1,
#'                                  value = c(2.1, 1.2, 0.7, 0.5, 0.3)))
#' yaml <- exportPEtab(data, f, observables = c(y = "A"), p = p,
#'                     pouter = c(logA = log(2), logk = log(0.5)),
#'                     dir = file.path(tempdir(), "petab_decay"),
#'                     formatVersion = "1")
#' readPEtabTables(yaml)$parameters
#' @export
exportPEtab <- function(data, reactions, observables, p, pouter,
                        errors = NULL,
                        lower = NULL, upper = NULL, fixed = NULL,
                        parameterScale = NULL,
                        observableTransformation = "lin",
                        noiseDistribution = "normal",
                        modelID = "dMod_export",
                        formatVersion = "2.0.0",
                        dir, overwrite = FALSE) {

  ## --- 1. datalist normalisation ------------------------------------------
  if (is.data.frame(data)) data <- as.datalist(data)
  if (is.list(data) && !inherits(data, "datalist")) data <- as.datalist(data)
  if (!inherits(data, "datalist"))
    stop("`data` must be a datalist, list of data.frames, or long-format data.frame.")

  ## --- 2. extract per-condition trafo from p ------------------------------
  if (!is.function(p) || is.null(attr(p, "mappings")))
    stop("`p` must be a parfn produced by P(), needed to decompose the parameter trafo.")
  eqs <- getEquations(p)
  if (!is.list(eqs)) eqs <- list(eqs)
  conds <- names(eqs)
  if (is.null(conds) || any(!nzchar(conds)))
    stop("`getEquations(p)` returned an unnamed list: every condition must have a name.")

  ## --- 3. reactions / observables / errors --------------------------------
  if (!inherits(reactions, "eqnlist"))
    stop("`reactions` must be an eqnlist (the ODE network with stoichiometry).")
  states <- reactions$states

  pull_eqns <- function(obj, what) {
    if (is.function(obj)) {
      eqns <- attr(obj, "equations")
      if (is.null(eqns))
        stop(sprintf(
          "`%s` is a function but has no `equations` attribute.", what))
      if (is.list(eqns)) eqns <- eqns[[1L]]
      return(setNames(as.character(eqns), names(eqns)))
    }
    if (inherits(obj, "eqnvec"))
      return(setNames(as.character(unclass(obj)), names(obj)))
    if (is.character(obj) && !is.null(names(obj))) return(obj)
    stop(sprintf("`%s` must be a named character vector, eqnvec, or Y-built fn.", what))
  }
  obs_eqns <- pull_eqns(observables, "observables")
  if (is.null(names(obs_eqns)) || any(!nzchar(names(obs_eqns))))
    stop("`observables` entries must be named (observableId -> formula).")
  obs_ids <- names(obs_eqns)

  # Per-row data sigmas go out through the `noiseParameter1_<obsId>` placeholder,
  # so they survive a round trip; without a sigma column the noise is "1".
  has_per_row_sigma <- any(vapply(data, function(d) {
    "sigma" %in% colnames(d) && any(!is.na(d$sigma))
  }, logical(1)))
  if (is.null(errors)) {
    err_eqns <- if (has_per_row_sigma)
      setNames(paste0("noiseParameter1_", obs_ids), obs_ids)
    else
      setNames(rep("1", length(obs_ids)), obs_ids)
  } else {
    err_eqns <- pull_eqns(errors, "errors")
  }
  miss <- setdiff(obs_ids, names(err_eqns))
  if (length(miss))
    stop("`errors` is missing entries for: ", paste(miss, collapse = ", "))
  err_eqns <- err_eqns[obs_ids]
  use_per_row_noise <- is.null(errors) && has_per_row_sigma

  ## --- 4. inner_pars / obs_inner (mirror importer logic) ------------------
  rate_syms <- unique(unlist(lapply(reactions$rates, function(r)
                getSymbols(r, exclude = c(states, "time")))))
  inner_pars <- setdiff(rate_syms, states)
  obs_syms <- unique(unlist(lapply(obs_eqns, function(f)
                getSymbols(f, exclude = c(states, "time")))))
  noise_syms <- unique(unlist(lapply(err_eqns, function(f) {
                if (is.na(suppressWarnings(as.numeric(f)))) getSymbols(f)
                else character(0)
              })))
  obs_inner <- setdiff(unique(c(obs_syms, noise_syms)),
                       c(states, inner_pars, "time"))
  # PEtab placeholders are bound per row in measurements.tsv, not inner
  # parameters, so the decomposer must not expect a mapping for them.
  obs_inner <- obs_inner[!grepl(
    "^(observable|noise)Parameter[0-9]+_", obs_inner)]

  ## --- 5. validate pouter / fixed / scales --------------------------------
  pouter_ids <- names(pouter)
  if (is.null(pouter_ids) || any(!nzchar(pouter_ids)))
    stop("`pouter` must be a named numeric vector.")
  if (is.null(fixed)) fixed <- numeric(0)
  if (length(fixed) > 0L && (is.null(names(fixed)) || any(!nzchar(names(fixed)))))
    stop("`fixed` must be a named numeric vector.")
  if (anyDuplicated(c(pouter_ids, names(fixed))))
    stop("Parameter ids overlap between `pouter` and `fixed`.")

  broadcast_scale <- function(s, ids) {
    if (length(s) == 1L && (is.null(names(s)) || !nzchar(names(s))))
      return(setNames(rep(unname(s), length(ids)), ids))
    miss <- setdiff(ids, names(s))
    if (length(miss))
      stop("`parameterScale` is missing entries for: ", paste(miss, collapse = ", "))
    s[ids]
  }
  scales_pouter <- if (is.null(parameterScale))
                     .petab_detect_scales(eqs, pouter_ids)
                   else broadcast_scale(parameterScale, pouter_ids)
  bad <- setdiff(unique(scales_pouter), c("lin", "log", "log10"))
  if (length(bad))
    stop("Unknown parameterScale(s): ", paste(bad, collapse = ", "))

  # Fixed parameters are linear: exportPEtabObject writes them with
  # parameterScale "lin".
  scales_fixed <- setNames(rep("lin", length(fixed)), names(fixed))
  scales_all   <- c(scales_pouter, scales_fixed)

  ## --- 6. bounds ----------------------------------------------------------
  # Default: five decades either side of pouter on log scales, open on lin.
  width <- c(lin = Inf, log = 5 * log(10), log10 = 5)[scales_pouter]
  if (is.null(lower)) lower <- setNames(unname(pouter - width), pouter_ids)
  if (is.null(upper)) upper <- setNames(unname(pouter + width), pouter_ids)
  if (!setequal(names(lower), pouter_ids))
    stop("`lower` must be named like `pouter`.")
  if (!setequal(names(upper), pouter_ids))
    stop("`upper` must be named like `pouter`.")
  lower <- lower[pouter_ids]; upper <- upper[pouter_ids]

  ## --- 6b. outer ids clashing with model entities --------------------------
  # An outer id equal to a state moves to `init_<state>`. One equal to an inner
  # parameter stays shared for a pure scale wrap and moves to `<k>_outer` else.
  stripped <- .petab_strip_trafo(eqs, scales_all)
  model_ids <- c(states, inner_pars, obs_inner)
  outer_ids <- c(pouter_ids, names(fixed))
  renames <- character(0)
  for (op in intersect(outer_ids, model_ids)) {
    if (!op %in% states &&
        all(vapply(stripped, function(e) !op %in% names(e) || e[[op]] == op,
                   logical(1)))) next
    new <- paste0(if (op %in% states) "init_" else "", op,
                  if (op %in% states) "" else "_outer")
    while (new %in% c(outer_ids, model_ids, renames)) new <- paste0(new, "_")
    renames[[op]] <- new
  }
  if (length(renames)) {
    message("exportPEtab: outer parameter(s) renamed to keep SBML ids unique: ",
            paste(names(renames), renames, sep = " -> ", collapse = ", "))
    eqs <- lapply(eqs, function(eqv) {
      out <- replaceSymbols(names(renames), unname(renames),
                            as.character(unclass(eqv)))
      setNames(out, names(eqv))
    })
    rn <- function(x) {
      hit <- names(x) %in% names(renames)
      names(x)[hit] <- renames[names(x)[hit]]
      x
    }
    pouter <- rn(pouter); lower <- rn(lower); upper <- rn(upper)
    scales_pouter <- rn(scales_pouter); fixed <- rn(fixed)
    scales_all <- rn(scales_all)
    pouter_ids <- names(pouter)
  }

  ## --- 7. decompose trafo -------------------------------------------------
  decomp <- .petab_decompose_trafo(
    eqs           = eqs,
    states        = states,
    inner_pars    = inner_pars,
    obs_inner     = obs_inner,
    pouter_names  = pouter_ids,
    fixed         = fixed,
    scales        = scales_all)
  .petab_check_self_refs(decomp$inits, decomp$conditions_df)

  # Decomposer SBML defaults split in two: override targets need an SBML
  # <parameter> but no parameters.tsv row; the rest are fixed outer parameters
  # of the imported p and go to parameters.tsv.
  override_targets <- setdiff(colnames(decomp$conditions_df), "conditionId")
  bound_idx        <- names(decomp$sbml_extra_pars) %in% override_targets
  sbml_only_pars   <- decomp$sbml_extra_pars[bound_idx]
  sbml_unbound     <- decomp$sbml_extra_pars[!bound_idx]
  fixed_full <- c(fixed, sbml_unbound)
  if (anyDuplicated(names(fixed_full)))
    stop("Internal error: duplicate fixed parameter id after decomposition: ",
         paste(names(fixed_full)[duplicated(names(fixed_full))], collapse = ", "))

  ## --- 8. observable metadata --------------------------------------------
  broadcast_obs <- function(x, ids, what, allowed) {
    if (length(x) == 1L && (is.null(names(x)) || !nzchar(names(x))))
      x <- setNames(rep(unname(x), length(ids)), ids)
    else {
      miss <- setdiff(ids, names(x))
      if (length(miss))
        stop(sprintf("`%s` is missing entries for: %s", what,
                     paste(miss, collapse = ", ")))
      x <- x[ids]
    }
    bad <- setdiff(unique(x), allowed)
    if (length(bad))
      stop(sprintf("Unknown %s value(s): %s", what, paste(bad, collapse = ", ")))
    x
  }
  obs_trafo  <- broadcast_obs(observableTransformation, obs_ids,
                              "observableTransformation",
                              c("lin", "log", "log10"))
  noise_dist <- broadcast_obs(noiseDistribution, obs_ids,
                              "noiseDistribution",
                              c("normal", "laplace", "log-normal"))
  obs_meta <- list(obs        = obs_eqns,
                   noise      = err_eqns,
                   obs_trafo  = obs_trafo,
                   noise_dist = noise_dist)

  ## --- 9. drive condition set from p (not from data); split sub-conds -----
  # Per-row sigmas that vary within a condition split it into sub-conditions
  # named `<sim_cond>__<noi_hash>`, as the importer does.
  data_conds <- names(data)
  miss_in_data <- setdiff(conds, data_conds)
  miss_in_p    <- setdiff(data_conds, conds)
  if (length(miss_in_p))
    warning("data has condition(s) not covered by p (dropped): ",
            paste(miss_in_p, collapse = ", "))

  data_filtered    <- list()
  scm_rows         <- list()
  obs_subs_by_sub  <- list()
  noi_subs_by_sub  <- list()

  for (c in intersect(conds, data_conds)) {
    d <- data[[c]]
    if (use_per_row_noise && "sigma" %in% colnames(d) && nrow(d) > 0L) {
      sig_str <- ifelse(is.na(d$sigma), "1",
                        trimws(formatC(d$sigma, digits = 15, format = "g")))
      noi_h   <- .petab_subcond_label(sig_str)
      uniq_h  <- unique(noi_h)
      single  <- length(uniq_h) <= 1L
      for (h in uniq_h) {
        sel  <- which(noi_h == h)
        sub_name <- if (single || h == "") c else paste0(c, "__", h)
        d_sub <- d[sel, , drop = FALSE]
        d_sub$sigma <- NA_real_  # importer rebuilds from noiseParameters
        data_filtered[[sub_name]] <- d_sub
        scm_rows[[length(scm_rows) + 1L]] <- data.frame(
          simulationConditionId       = c,
          preequilibrationConditionId = "",
          sub_condition               = sub_name,
          stringsAsFactors            = FALSE)
        noi_subs_by_sub[[sub_name]] <- c("*" = sig_str[sel[[1L]]])
      }
    } else {
      data_filtered[[c]] <- d
      scm_rows[[length(scm_rows) + 1L]] <- data.frame(
        simulationConditionId       = c,
        preequilibrationConditionId = "",
        sub_condition               = c,
        stringsAsFactors            = FALSE)
    }
  }
  for (c in miss_in_data) {
    data_filtered[[c]] <- data.frame(
      name = character(0), time = numeric(0),
      value = numeric(0), sigma = numeric(0),
      stringsAsFactors = FALSE)
    scm_rows[[length(scm_rows) + 1L]] <- data.frame(
      simulationConditionId       = c,
      preequilibrationConditionId = "",
      sub_condition               = c,
      stringsAsFactors            = FALSE)
  }
  attr(data_filtered, "class") <- attr(data, "class")  # preserve datalist
  sub_cond_map <- do.call(rbind, scm_rows)
  rownames(sub_cond_map) <- NULL
  attr(sub_cond_map, "obs_subs") <- obs_subs_by_sub
  attr(sub_cond_map, "noi_subs") <- noi_subs_by_sub

  ## --- 10. tag pouter with scales attr (consumed by exportPEtabObject) ---
  attr(pouter, "petab_scales") <- scales_pouter

  petab <- list(
    pouter         = pouter,
    lower          = lower,
    upper          = upper,
    fixed          = fixed_full,
    sbml_only_pars = sbml_only_pars,
    obs_meta       = obs_meta,
    sub_cond_map   = sub_cond_map,
    condition.grid = decomp$conditions_df,
    data           = data_filtered,
    reactions      = reactions,
    inits          = decomp$inits,
    modelID        = modelID)

  exportPEtabObject(petab, dir, modelID = modelID,
                    formatVersion = formatVersion,
                    overwrite = overwrite)
}


#' Export an Imported PEtab Problem
#'
#' Writes a problem returned by [importPEtab()] back to PEtab: the parameter,
#' observable, condition and measurement tables, for v2 also the experiment
#' table, the SBML model and a YAML file. Conditions the importer split are
#' merged again. For a problem built in dMod, use [exportPEtab()].
#'
#' Limits:
#' \itemize{
#'   \item Parameter scales are read from `attr(petab$bestfit,
#'     "petab_scales")`; without it they are `"lin"`.
#'   \item v2 has no parameter scales: values and bounds are written on the
#'     linear scale, with a warning.
#'   \item v2 has no `log10` observable transformation: it is written as
#'     `"log-normal"` noise, with a warning.
#'   \item Initial values set per condition only in the parameter
#'     transformation are not written to a v2 condition table; use
#'     `formatVersion = "1"` to keep them.
#'   \item SBML `<algebraicRule>` elements are not imported and therefore not
#'     written.
#' }
#'
#' @param petab A `petabproblem` returned by [importPEtab()], or a list with
#'   the same components.
#' @param dir Output directory, created if missing.
#' @param modelID SBML model identifier. Default `NULL`: the model id of the
#'   imported problem, else `"dMod_export"`.
#' @param formatVersion `"2.0.0"` (default) or `"1"`.
#' @param overwrite Logical. If `FALSE` (default), existing files in `dir` are
#'   an error.
#' @return Path to the written YAML file, invisibly.
#' @seealso [importPEtab()], [exportPEtab()]
#' @inherit importPEtab examples
#' @export
exportPEtabObject <- function(petab, dir, modelID = NULL,
                              formatVersion = "2.0.0",
                              overwrite = FALSE) {

  .require_ns("yaml", "PEtab export")
  stopifnot(inherits(petab, "petabproblem") || is.list(petab))
  major <- .petab_major_version(formatVersion)
  if (!major %in% c(1L, 2L))
    stop("`formatVersion` must be '1' or '2.0.0' (got ", formatVersion, ").")

  meta <- attr(petab, "petab_meta") %||% list()
  fixed        <- meta$fixed        %||% petab$fixed        %||% numeric(0)
  # SBML-only parameters: declared in the SBML model so conditions.tsv
  # targetIds resolve, but excluded from parameters.tsv (their values are
  # determined by per-condition overrides, not the nominal default).
  sbml_only    <- meta$sbml_only_pars %||% petab$sbml_only_pars %||% numeric(0)
  start_times  <- meta$startTimes %||% list()
  switch_map   <- meta$switches   %||% list()
  inits_meta   <- meta$inits        %||% petab$inits
  obs_meta     <- meta$obs_meta     %||% petab$obs_meta
  sub_cond_map <- meta$sub_cond_map %||% petab$sub_cond_map
  param_meta   <- meta$param_meta   %||% petab$param_meta
  # The saved conditions table (from importPEtab or exportPEtab) wins over the
  # data's condition.grid, which lacks the override columns.
  cond_grid <- meta$cond_grid
  if (is.null(cond_grid))
    cond_grid <- petab$condition.grid
  if (is.null(cond_grid))
    cond_grid <- attr(petab$dataList %||% petab$data, "condition.grid")
  bestfit      <- petab$bestfit %||% petab$pouter
  parlower     <- petab$parlower %||% petab$lower
  parupper     <- petab$parupper %||% petab$upper
  data_list    <- petab$dataList %||% petab$data
  reactions    <- petab$reactions
  if (is.null(modelID))
    modelID <- meta$modelID %||% petab$modelID %||% "dMod_export"

  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  paths <- list(
    parameters   = file.path(dir, paste0("parameters_",   modelID, ".tsv")),
    observables  = file.path(dir, paste0("observables_",  modelID, ".tsv")),
    conditions   = file.path(dir, paste0("conditions_",   modelID, ".tsv")),
    experiments  = file.path(dir, paste0("experiments_",  modelID, ".tsv")),
    measurements = file.path(dir, paste0("measurements_", modelID, ".tsv")),
    sbml         = file.path(dir, paste0(modelID, ".xml")),
    yaml         = file.path(dir, paste0(modelID, ".yaml"))
  )
  must_exist <- c("parameters", "observables", "conditions",
                  "measurements", "sbml", "yaml")
  if (major == 2L) must_exist <- c(must_exist, "experiments")
  if (!overwrite) {
    existing <- vapply(paths[must_exist], file.exists, logical(1))
    if (any(existing))
      stop("Output file(s) already exist; pass overwrite = TRUE to replace: ",
           paste(unlist(paths[must_exist])[existing], collapse = ", "))
  }

  # v1 simulates one condition per experiment from t = 0; a period switch or a
  # later start would be lost.
  if (major == 1L && length(switch_map))
    stop("PEtab v1 cannot express condition changes during a simulation ",
         "(experiment periods in ", paste(names(switch_map), collapse = ", "),
         "). Export with formatVersion = \"2.0.0\".", call. = FALSE)
  late <- names(start_times)[vapply(start_times, function(t) t != 0, logical(1))]
  if (major == 1L && length(late))
    stop("PEtab v1 simulates from t = 0, but ", paste(late, collapse = ", "),
         " start later. Export with formatVersion = \"2.0.0\".", call. = FALSE)

  scales <- attr(bestfit, "petab_scales")
  if (is.null(scales))
    scales <- setNames(rep("lin", length(bestfit)), names(bestfit))

  est_ids <- names(bestfit)

  # `pouter` lives on the PEtab parameter scale, both tables on disk hold
  # linear values, so the importer's forward transform is inverted here.
  apply_inv_scale <- function(values, ids) {
    sc <- scales[ids]
    out <- values
    log_idx   <- which(sc == "log")
    log10_idx <- which(sc == "log10")
    if (length(log_idx))   out[log_idx]   <- exp(values[log_idx])
    if (length(log10_idx)) out[log10_idx] <- 10 ^ values[log10_idx]
    out
  }

  # Bounds are round numbers; drop the exp/log round-off.
  bound <- function(b) signif(unname(apply_inv_scale(b[est_ids], est_ids)), 12)

  # Symbols the condition table, the observables and the per-row placeholder
  # strings name. A fixed parameter only they use appears in no SBML
  # element, so it has to stay in parameters.tsv.
  cond_syms <- unique(unlist(lapply(
    setdiff(colnames(cond_grid), c("conditionId", "conditionName", "condition")),
    function(tid) {
      v <- as.character(cond_grid[[tid]])
      v <- v[!is.na(v) & nzchar(trimws(v))]
      unlist(lapply(v, getSymbols))
    })))
  subs_syms <- unlist(lapply(
    c(attr(sub_cond_map, "obs_subs") %||% list(),
      attr(sub_cond_map, "noi_subs") %||% list()),
    function(m) unlist(strsplit(as.character(m), ";", fixed = TRUE))))
  obs_syms <- unique(c(
    unlist(lapply(c(unlist(obs_meta$obs), unlist(obs_meta$noise)), getSymbols)),
    trimws(subs_syms)))

  if (major == 1L) {
    est_df <- data.frame(
      parameterId    = est_ids,
      parameterScale = unname(scales[est_ids]),
      lowerBound     = bound(parlower),
      upperBound     = bound(parupper),
      nominalValue   = unname(apply_inv_scale(bestfit, est_ids)),
      estimate       = 1L,
      stringsAsFactors = FALSE
    )
    # Compartment sizes and species values live in the SBML model, condition
    # targets in conditions.tsv; v1 allows neither in parameters.tsv.
    cond_targets <- setdiff(colnames(cond_grid),
                            c("conditionId", "conditionName", "condition"))
    fixed_v1 <- fixed[!names(fixed) %in% c(names(reactions$compartments),
                                           reactions$states, cond_targets)]
    # An imported problem keeps its SBML constants in the SBML, as the source
    # did: the importer reads a constant outside parameters.tsv as structure
    # (a rate constant at 0 idles its reaction in a preequilibration).
    if (!is.null(param_meta))
      fixed_v1 <- fixed_v1[names(fixed_v1) %in%
                             c(names(param_meta$fixed), cond_syms, obs_syms)]
    fixed_df <- if (length(fixed_v1)) data.frame(
      parameterId    = names(fixed_v1),
      parameterScale = "lin",
      lowerBound     = -Inf,
      upperBound     = Inf,
      nominalValue   = unname(fixed_v1),
      estimate       = 0L,
      stringsAsFactors = FALSE
    ) else NULL
  } else {
    # v2 has no parameterScale column, so non-lin parameters are linearised with
    # one warning; a re-import optimises the linear parameter.
    nonlin <- est_ids[scales[est_ids] != "lin"]
    if (length(nonlin))
      warning("exportPEtabObject: PEtab v2 has no `parameterScale`. ",
              "Linearising ", paste(nonlin, collapse = ", "),
              ". Export to formatVersion = \"1\" to keep the scale.",
              call. = FALSE)

    # v2 keeps non-estimated parameters the SBML declares out of parameters.tsv,
    # except one a conditions.tsv `targetValue` names, which needs declaring.
    est_df <- data.frame(
      parameterId  = est_ids,
      lowerBound   = bound(parlower),
      upperBound   = bound(parupper),
      nominalValue = unname(apply_inv_scale(bestfit, est_ids)),
      estimate     = "true",
      stringsAsFactors = FALSE
    )
    keep_fixed <- intersect(names(fixed), c(cond_syms, obs_syms))
    fixed_df <- if (length(keep_fixed)) data.frame(
      parameterId  = keep_fixed,
      lowerBound   = -Inf,
      upperBound   = Inf,
      nominalValue = unname(fixed[keep_fixed]),
      estimate     = "false",
      stringsAsFactors = FALSE
    ) else NULL
  }

  # Round-trip priors verbatim: the parsed record keeps the declared
  # spelling and parameter string, so a re-export reproduces the source
  # tables rather than collapsing every distribution onto one.
  priors <- param_meta$priors
  if (!is.null(priors)) {
    pd <- rep(NA_character_, nrow(est_df))
    pp <- rep(NA_character_, nrow(est_df))
    hits <- match(names(priors), est_df$parameterId)
    keep <- !is.na(hits)
    pd[hits[keep]] <- vapply(priors[keep], `[[`, character(1), "declared")
    pp[hits[keep]] <- vapply(priors[keep], `[[`, character(1), "declaredPars")
    if (major == 1L) {
      v1_name <- c("log-normal" = "logNormal", "log-laplace" = "logLaplace")
      renamed <- !is.na(pd) & pd %in% names(v1_name)
      pd[renamed] <- v1_name[pd[renamed]]
      bad <- setdiff(stats::na.omit(pd),
                     c("uniform", "normal", "laplace", "logNormal", "logLaplace",
                       "parameterScaleUniform", "parameterScaleNormal",
                       "parameterScaleLaplace"))
      if (length(bad))
        warning("PEtab v1 has no prior distribution ",
                paste(bad, collapse = ", "), ". Export with formatVersion = ",
                "\"2.0.0\" to keep it.", call. = FALSE)
      est_df$objectivePriorType       <- pd
      est_df$objectivePriorParameters <- pp
      if (!is.null(fixed_df)) {
        fixed_df$objectivePriorType       <- NA_character_
        fixed_df$objectivePriorParameters <- NA_character_
      }
    } else {
      # v2 has no parameter scale: a prior on the scaled parameter becomes the
      # matching linear one, log10-normal(mu, sd) as log-normal(mu ln 10, sd ln 10).
      v2_name <- c(logNormal = "log-normal", logLaplace = "log-laplace")
      renamed <- !is.na(pd) & pd %in% names(v2_name)
      pd[renamed] <- v2_name[pd[renamed]]
      on_scale <- which(!is.na(pd) & grepl("^parameterScale", pd))
      for (i in on_scale) {
        base <- tolower(sub("^parameterScale", "", pd[i]))
        sc   <- unname(scales[est_df$parameterId[i]])
        v    <- as.numeric(strsplit(pp[i], ";", fixed = TRUE)[[1L]])
        f    <- if (identical(sc, "log10")) log(10) else 1
        if (identical(sc, "lin")) {
          pd[i] <- base
        } else if (base == "uniform") {
          pd[i] <- "log-uniform"
          v <- exp(v * f)
        } else {
          pd[i] <- paste0("log-", base)
          v <- v * f
        }
        pp[i] <- paste(vapply(v, format, character(1), digits = 15),
                       collapse = ";")
      }
      if (length(on_scale))
        warning("PEtab v2 has no priors on the parameter scale; ",
                paste(est_df$parameterId[on_scale], collapse = ", "),
                " take the matching prior of the linear parameter.",
                call. = FALSE)
      est_df$priorDistribution <- pd
      est_df$priorParameters   <- pp
      if (!is.null(fixed_df)) {
        fixed_df$priorDistribution <- NA_character_
        fixed_df$priorParameters   <- NA_character_
      }
    }
  }

  utils::write.table(rbind(est_df, fixed_df), paths$parameters,
                     sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  ## --- observables.tsv ---------------------------------------------------
  om <- obs_meta

  # Numeric noise lives in the data's `sigma` column and may vary by row; such an
  # observable exports one noise placeholder and each row writes its sigma.
  .sigma_of <- function(o) unlist(lapply(data_list, function(d)
    d$sigma[as.character(d$name) == o]), use.names = FALSE)
  numeric_sigma <- vapply(names(om$obs), function(o) {
    v <- .sigma_of(o)
    length(v) > 0L && all(is.finite(v))
  }, logical(1))

  # v2 has no observable transformation: a log10 observable goes out as
  # `log-normal`, whose residual is on the natural log. The noise is scaled by
  # ln(10) so the same standard deviation describes the same distribution.
  ln10 <- 2.302585092994046
  .rescale_v2 <- major != 1L & om$obs_trafo[names(om$obs)] == "log10"

  # An observable with a single sigma throughout writes it as its noise
  # formula; the per-row placeholder is kept for varying sigmas.
  one_sigma <- vapply(names(om$obs), function(o) {
    v <- unique(.sigma_of(o))
    if (numeric_sigma[[o]] && length(v) == 1L) v else NA_real_
  }, numeric(1))
  has_one <- !is.na(one_sigma)
  numeric_sigma <- numeric_sigma & !has_one
  if (any(has_one))
    om$noise[names(om$obs)[has_one]] <- vapply(
      one_sigma[has_one] * ifelse(.rescale_v2[has_one], ln10, 1),
      format, character(1), digits = 15)

  if (major == 1L) {
    noi_v1 <- unname(om$noise[names(om$obs)])
    # v2 lets a noise formula name an observable; v1 needs its formula.
    noi_v1 <- vapply(noi_v1, function(f) {
      if (is.na(f)) return(f)
      hit <- intersect(getSymbols(f), names(om$obs))
      if (!length(hit)) return(f)
      replaceSymbols(hit, paste0("(", om$obs[hit], ")"), f)
    }, character(1), USE.NAMES = FALSE)
    noi_v1[numeric_sigma] <- paste0("noiseParameter1_", names(om$obs)[numeric_sigma])
    obs_df <- data.frame(
      observableId             = names(om$obs),
      observableFormula        = unname(om$obs),
      observableTransformation = unname(om$obs_trafo[names(om$obs)]),
      noiseFormula             = noi_v1,
      noiseDistribution        = unname(om$noise_dist[names(om$obs)]),
      stringsAsFactors = FALSE
    )
  } else {
    obs_ids <- names(om$obs)
    obs_trafo  <- om$obs_trafo[obs_ids]
    noise_dist <- om$noise_dist[obs_ids]
    if (any(obs_trafo == "log10"))
      warning("PEtab v2 dropped the `log10` observable transformation; ",
              "emitting `log-normal` with the noise scaled by ln(10) for: ",
              paste(obs_ids[obs_trafo == "log10"], collapse = ", "),
              ". Residuals are preserved; the likelihood value shifts by the ",
              "log10/log Jacobian.", call. = FALSE)
    v2_dist <- ifelse(obs_trafo %in% c("log", "log10") &
                      noise_dist == "normal", "log-normal",
              ifelse(obs_trafo == "lin" & noise_dist == "normal", "normal",
              ifelse(obs_trafo %in% c("log", "log10") &
                     noise_dist == "laplace", "log-laplace",
              ifelse(obs_trafo == "lin" & noise_dist == "laplace", "laplace",
                     NA_character_))))
    if (any(is.na(v2_dist)))
      stop("Cannot encode (obs_trafo, noise_dist) into v2 noiseDistribution: ",
           paste(sprintf("%s=(%s,%s)", obs_ids[is.na(v2_dist)],
                         obs_trafo[is.na(v2_dist)],
                         noise_dist[is.na(v2_dist)]), collapse = "; "))
    obs_form <- unname(om$obs[obs_ids])
    noi_form <- unname(om$noise[obs_ids])
    # A symbolic noise keeps its parameters and takes the factor in the formula;
    # a numeric one is written per row and takes it there instead.
    scale_here <- unname(.rescale_v2[obs_ids]) & !unname(numeric_sigma[obs_ids]) &
                  !unname(has_one[obs_ids])
    if (any(scale_here))
      noi_form[scale_here] <- sprintf("(%s) * %.15g", noi_form[scale_here], ln10)
    noi_form[unname(numeric_sigma[obs_ids])] <-
      paste0("noiseParameter1_", obs_ids[unname(numeric_sigma[obs_ids])])
    scan_placeholders <- function(formula, prefix, obs_id) {
      if (is.na(formula) || !nzchar(formula)) return("")
      pat <- sprintf("%s[0-9]+_%s", prefix, obs_id)
      m <- regmatches(formula, gregexpr(pat, formula, perl = TRUE))[[1L]]
      if (!length(m)) return("")
      idx <- as.integer(sub(sprintf("^%s([0-9]+)_.*", prefix), "\\1", m))
      uniq <- m[!duplicated(idx)]
      uniq <- uniq[order(idx[!duplicated(idx)])]
      paste(uniq, collapse = ";")
    }
    obs_ph   <- vapply(seq_along(obs_ids), function(i)
      scan_placeholders(obs_form[i], "observableParameter", obs_ids[i]),
      character(1))
    noise_ph <- vapply(seq_along(obs_ids), function(i)
      scan_placeholders(noi_form[i], "noiseParameter", obs_ids[i]),
      character(1))
    obs_df <- data.frame(
      observableId           = obs_ids,
      observableFormula      = obs_form,
      observablePlaceholders = obs_ph,
      noiseFormula           = noi_form,
      noiseDistribution      = v2_dist,
      noisePlaceholders      = noise_ph,
      stringsAsFactors = FALSE
    )
  }
  utils::write.table(obs_df, paths$observables, sep = "\t", quote = FALSE,
                     row.names = FALSE, na = "")

  dflt_cond <- meta$defaultCondition %||% petab$defaultCondition %||% ""
  # Conditions may be keyed by their human-readable name; the tables on disk
  # use the id, so the map recorded at import is applied on the way out.
  cond_ids <- meta$conditionIds %||% petab$conditionIds
  as_id <- function(cid) {
    if (is.null(cond_ids) || !cid %in% names(cond_ids)) return(cid)
    unname(cond_ids[[cid]])
  }

  ## --- conditions.tsv ---------------------------------------------------
  cond_inits <- list()
  scm <- sub_cond_map
  uniq_sims  <- unique(scm$simulationConditionId)
  uniq_preeq <- setdiff(unique(scm$preequilibrationConditionId), "")
  # Conditions a period switches into are not simulation conditions, but the
  # exported table still has to define them.
  sw_conds <- unlist(lapply(switch_map, function(sw) as.character(sw$conditionId)))
  uniq_conds <- unique(c(uniq_sims, uniq_preeq, sw_conds))

  if (major == 1L) {
    cond_rows <- intersect(uniq_conds, rownames(cond_grid))
    cond_df <- cond_grid[cond_rows, , drop = FALSE]
    cond_df$conditionId <- vapply(cond_rows, as_id, character(1))
    cond_df <- cond_df[, c("conditionId",
                           setdiff(colnames(cond_df), "conditionId")),
                       drop = FALSE]
    # v1 holds only numbers and parameter ids in conditions.tsv. A species set
    # to one expression in every condition takes it as an SBML initial
    # assignment; a per-condition expression has no v1 form.
    is_expr <- function(v) {
      v <- trimws(as.character(v))
      !is.na(v) & nzchar(v) & is.na(suppressWarnings(as.numeric(v))) &
        !grepl("^[A-Za-z_][A-Za-z0-9_]*$", v)
    }
    expr_cols <- character(0)
    for (tid in setdiff(colnames(cond_df),
                        c("conditionId", "conditionName", "condition"))) {
      v <- as.character(cond_df[[tid]])
      if (!any(is_expr(v))) next
      if (tid %in% reactions$states && length(unique(v)) == 1L) {
        cond_inits[[tid]] <- v[[1L]]
        cond_df[[tid]] <- NULL
      } else {
        expr_cols <- c(expr_cols, tid)
      }
    }
    if (length(expr_cols))
      warning("PEtab v1 allows only numbers and parameter ids in the ",
              "condition table, but ", paste(expr_cols, collapse = ", "),
              " take expressions. dMod reads them back; export with ",
              "formatVersion = \"2.0.0\" for other tools.", call. = FALSE)
  } else {
    # long-format: one row per (conditionId, targetId) override.
    # Skip the dMod-internal `condition` column whose values just echo
    # conditionId, it's a `as.datalist` artifact, not a real override.
    rows <- list()
    cg_cols <- setdiff(colnames(cond_grid),
                       c("conditionId", "conditionName", "condition"))
    for (cid in uniq_conds) {
      if (!cid %in% rownames(cond_grid)) next
      for (tid in cg_cols) {
        v <- cond_grid[cid, tid]
        if (is.na(v) || !nzchar(trimws(as.character(v)))) next
        # Self-mapping (value == conditionId) is also a dMod artifact.
        if (identical(as.character(v), cid)) next
        rows[[length(rows) + 1L]] <- data.frame(
          conditionId  = as_id(cid),
          targetId     = tid,
          targetValue  = as.character(v),
          stringsAsFactors = FALSE)
      }
    }
    cond_df <- if (length(rows)) do.call(rbind, rows) else
      data.frame(conditionId = character(0), targetId = character(0),
                 targetValue = character(0), stringsAsFactors = FALSE)
  }
  utils::write.table(cond_df, paths$conditions, sep = "\t", quote = FALSE,
                     row.names = FALSE, na = "")

  ## --- experiments.tsv (v2 only) ----------------------------------------
  # One synthesised experimentId per (sim, preeq) sub-condition tuple.
  uniq_pairs <- unique(scm[, c("simulationConditionId",
                               "preequilibrationConditionId"),
                           drop = FALSE])
  exp_id_for <- function(sim, preeq) {
    if (nzchar(dflt_cond) && identical(sim, dflt_cond)) return("")
    if (!nzchar(preeq)) as_id(sim) else paste0(as_id(preeq), "__", as_id(sim))
  }
  uniq_pairs$experimentId <- mapply(
    exp_id_for,
    uniq_pairs$simulationConditionId,
    uniq_pairs$preequilibrationConditionId)

  if (major == 2L) {
    # The synthesised condition, and any condition without a row in the
    # condition table, is written as an empty conditionId, which v2 reads as
    # "no overrides".
    named <- function(cid) {
      if (nzchar(dflt_cond) && identical(cid, dflt_cond)) return("")
      if (!as_id(cid) %in% cond_df$conditionId) return("")
      as_id(cid)
    }
    exp_rows <- list()
    for (i in seq_len(nrow(uniq_pairs))) {
      sim <- uniq_pairs$simulationConditionId[i]
      pre <- uniq_pairs$preequilibrationConditionId[i]
      t0  <- if (is.null(start_times[[sim]])) 0 else unname(start_times[[sim]])
      sw  <- switch_map[[sim]]
      eid <- uniq_pairs$experimentId[i]
      if (!nzchar(eid)) {
        # The default condition needs an experiment row only when the
        # schedule says more than "start at 0": a preequilibration, a later
        # start, or a period switch.
        if (!nzchar(pre) && t0 == 0 && is.null(sw)) next
        eid <- "experiment_1"
        uniq_pairs$experimentId[i] <- eid
      }
      if (nzchar(pre))
        exp_rows[[length(exp_rows) + 1L]] <- data.frame(
          experimentId = eid, time = "-inf", conditionId = named(pre),
          stringsAsFactors = FALSE)
      exp_rows[[length(exp_rows) + 1L]] <- data.frame(
        experimentId = eid, time = format(t0, digits = 17),
        conditionId = named(sim), stringsAsFactors = FALSE)
      for (k in seq_len(NROW(sw)))
        exp_rows[[length(exp_rows) + 1L]] <- data.frame(
          experimentId = eid, time = format(sw$time[k], digits = 17),
          conditionId = named(sw$conditionId[k]), stringsAsFactors = FALSE)
    }
    exp_df <- if (length(exp_rows)) do.call(rbind, exp_rows) else
      data.frame(experimentId = character(0), time = character(0),
                 conditionId = character(0), stringsAsFactors = FALSE)
    utils::write.table(exp_df, paths$experiments, sep = "\t", quote = FALSE,
                       row.names = FALSE, na = "")
  }

  ## --- measurements.tsv -------------------------------------------------
  scm_obs_subs <- attr(scm, "obs_subs") %||% list()
  scm_noi_subs <- attr(scm, "noi_subs") %||% list()
  pick_subs <- function(sub, obsId, map) {
    m <- map[[sub]]
    if (is.null(m)) return("")
    if (obsId %in% names(m)) return(unname(m[obsId]))
    if ("*" %in% names(m))   return(unname(m["*"]))
    ""
  }
  # PEtab keeps `measurement` on the linear scale and applies the observable
  # transformation inside the likelihood; the datalist holds the transformed
  # values, so invert before writing.
  unscale_obs <- function(value, trafo) {
    out <- value
    lg  <- which(trafo == "log")
    l10 <- which(trafo == "log10")
    if (length(lg))  out[lg]  <- exp(value[lg])
    if (length(l10)) out[l10] <- 10 ^ value[l10]
    out
  }
  meas_rows <- do.call(rbind, lapply(names(data_list), function(sub) {
    df <- data_list[[sub]]
    ix <- match(sub, scm$sub_condition)
    sim <- scm$simulationConditionId[ix]
    pre <- scm$preequilibrationConditionId[ix]
    obs_par <- vapply(df$name, function(o) pick_subs(sub, o, scm_obs_subs),
                      character(1))
    noi_par <- vapply(df$name, function(o) pick_subs(sub, o, scm_noi_subs),
                      character(1))
    # Rows of an observable whose noise resolved numerically keep their own
    # sigma; on the log-normal side of a v2 export it is the ln-scale one.
    one <- unname(has_one[as.character(df$name)])
    noi_par[!is.na(one) & one] <- ""
    own <- unname(numeric_sigma[as.character(df$name)])
    own[is.na(own)] <- FALSE
    if (any(own)) {
      sg <- df$sigma[own]
      sg[unname(.rescale_v2[as.character(df$name)])[own]] <-
        sg[unname(.rescale_v2[as.character(df$name)])[own]] * ln10
      noi_par[own] <- format(sg, digits = 17, scientific = FALSE, trim = TRUE)
    }
    base <- data.frame(
      observableId         = df$name,
      # A row parked at the equilibration horizon is a steady-state
      # measurement and goes back out as `inf`.
      time                 = ifelse(df$time == .petab_ss_time, "inf",
                                    format(df$time, digits = 17,
                                           scientific = FALSE, trim = TRUE)),
      measurement          = unscale_obs(df$value, om$obs_trafo[df$name]),
      observableParameters = obs_par,
      noiseParameters      = noi_par,
      stringsAsFactors = FALSE)
    if (major == 1L) {
      base$simulationConditionId       <- as_id(sim)
      base$preequilibrationConditionId <- if (nzchar(pre)) as_id(pre) else pre
    } else {
      pair_ix <- which(uniq_pairs$simulationConditionId == sim &
                       uniq_pairs$preequilibrationConditionId == pre)
      base$experimentId <- uniq_pairs$experimentId[pair_ix]
    }
    base
  }))
  if (major == 1L) {
    base_cols <- c("observableId", "simulationConditionId", "time",
                   "measurement", "preequilibrationConditionId",
                   "observableParameters", "noiseParameters")
    meas_rows <- meas_rows[, base_cols, drop = FALSE]
  } else {
    base_cols <- c("observableId", "experimentId", "time", "measurement",
                   "observableParameters", "noiseParameters")
    meas_rows <- meas_rows[, base_cols, drop = FALSE]
  }
  drop_optional <- c("preequilibrationConditionId", "observableParameters",
                     "noiseParameters", "experimentId")
  for (col in intersect(drop_optional, colnames(meas_rows))) {
    if (all(meas_rows[[col]] == "")) meas_rows[[col]] <- NULL
  }
  utils::write.table(meas_rows, paths$measurements, sep = "\t", quote = FALSE,
                     row.names = FALSE, na = "")

  ## --- SBML --------------------------------------------------------------
  # SBML stores linear values; a log-scaled nominal 0 would otherwise reach
  # the writer as -Inf, which libsbml rejects.
  all_pars <- c(apply_inv_scale(bestfit, names(bestfit)),
                apply_inv_scale(fixed, names(fixed)), sbml_only)
  inits <- inits_meta %||%
           setNames(rep(0, length(reactions$states)), reactions$states)
  if (length(cond_inits)) {
    inits <- as.list(inits)
    inits[names(cond_inits)] <- cond_inits
  }
  # Event targets promoted to states on import go back to being parameters or
  # compartments; they appear in no reaction, so dropping the column is enough.
  promoted <- meta$promoted %||% petab$promoted %||% character(0)
  promoted <- intersect(promoted, reactions$states)
  if (length(promoted)) {
    keep <- setdiff(reactions$states, promoted)
    smat <- reactions$smatrix
    if (!is.null(smat)) smat <- smat[, match(keep, reactions$states), drop = FALSE]
    reactions <- eqnlist(smatrix = smat, states = keep, rates = reactions$rates,
                         description = reactions$description,
                         compartments = reactions$compartments,
                         compartmentOf = reactions$compartmentOf[keep],
                         reactionCompartment = reactions$reactionCompartment,
                         amountStates = intersect(reactions$amountStates, keep))
    # The promoted state's initial value is the symbol's own value; without it
    # a compartment would export with the default size.
    restore <- setdiff(promoted, names(all_pars))
    vals <- suppressWarnings(as.numeric(inits[restore]))
    restore <- restore[!is.na(vals)]
    if (length(restore))
      all_pars <- c(all_pars, setNames(vals[!is.na(vals)], restore))
    inits <- inits[setdiff(names(inits), promoted)]
  }

  # A condition target nothing else declares (a parameter read only by an
  # initial assignment) still needs its SBML id, and a compartment sized by
  # its own symbol its size; both take the model default.
  sbml_defaults <- meta$sbml_pars %||% numeric(0)
  comps <- names(reactions$compartments)
  own_size <- comps[vapply(comps, function(cid)
    identical(reactions$compartments[[cid]]$volume, cid), logical(1))]
  undeclared <- c(
    setdiff(setdiff(colnames(cond_grid), c("conditionId", "conditionName", "condition")),
            c(names(all_pars), reactions$states, names(reactions$compartments))),
    setdiff(own_size, names(all_pars)))
  if (length(undeclared))
    all_pars <- c(all_pars, vapply(undeclared, function(x)
      if (x %in% names(sbml_defaults)) as.numeric(sbml_defaults[[x]]) else 1,
      numeric(1)))

  .petab_check_self_refs(inits, cond_grid, conditions = major == 1L)
  exportSbml(reactions, parameters = all_pars, inits = inits,
              filepath = paths$sbml, modelID = modelID,
              events = meta$eventsSource %||% petab$eventsSource)

  ## --- YAML manifest -----------------------------------------------------
  if (major == 1L) {
    manifest <- list(
      format_version = 1L,
      parameter_file = basename(paths$parameters),
      problems = list(list(
        sbml_files        = list(basename(paths$sbml)),
        condition_files   = list(basename(paths$conditions)),
        measurement_files = list(basename(paths$measurements)),
        observable_files  = list(basename(paths$observables))
      ))
    )
  } else {
    manifest <- list(
      format_version    = "2.0.0",
      parameter_files   = list(basename(paths$parameters)),
      model_files       = setNames(
        list(list(location = basename(paths$sbml), language = "sbml")),
        modelID),
      observable_files  = list(basename(paths$observables)),
      measurement_files = list(basename(paths$measurements)),
      condition_files   = list(basename(paths$conditions)),
      experiment_files  = list(basename(paths$experiments))
    )
  }
  yaml::write_yaml(manifest, paths$yaml)

  invisible(paths$yaml)
}


# `%||%` is base R since 4.4; relied on per dMod's existing usage.
