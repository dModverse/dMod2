#' Alignment model of replicate measurements
#'
#' @description
#' Builds the model of an alignment: readings of several experiments (gels,
#' plates, runs) that measure the same biological values up to experiment
#' specific factors, written as a free expression per reading as in blotIt.
#' [align()] fits it.
#'
#' Every symbol of the expressions is a linear predictor over the columns of
#' `data`, set by the formulas of the roles in the coding of
#' [stats::model.matrix()] with treatment contrasts, as in [insert()] with
#' effect formulas: `yi ~ target * time` gives one value per combination of
#' target and time, `sj ~ gel + target` additive effects. A column gives one
#' level per value, numeric or not; `I(dose)` enters a numeric column as a
#' covariate with one slope. Random terms, `(1 | gel)`, `(1 + dose | p | gel)`
#' and `gr(gel, dist = spikeSlab())`, add random effects per group as in
#' [randomEffects()]. `trafo` sets the scale of a symbol, `"sj <- 10^sj"`: the
#' predictor is then on the log10 scale.
#'
#' Values and scales share a common factor. As in blotIt, the experiments that
#' share biological values form blocks, and in each block the mean of the
#' predictors of the first biological symbol over the values of that symbol is
#' held at 0, by a
#' residual of weight 1000: the geometric mean of the values is 1 for
#' `"yi <- exp(yi)"` or `"yi <- 2^yi"`, the mean of a symbol without `trafo`,
#' a log scale or an offset, is 0.
#'
#' Absolute readings, on the scale the values are wanted on, such as
#' concentrations in nM by targeted mass spectrometry, are marked in a logical
#' column `absolute`. Their prediction is the expression `absolute` in the
#' biological symbols and data columns, so they depend on no scaling
#' parameter; the scaling formulas are coded over the other readings, and
#' their columns may be `NA` in absolute rows. A block that holds absolute
#' readings is not normalized: the values in it and the experiments connected to
#' them come out in the units of `absolute`.
#'
#' @param data data.frame with a column `value`, the readings, and the columns
#'   named in the formulas and expressions. An optional column `lloq` holds the
#'   lower limit of quantification of each reading (`NA` for none); a reading
#'   with a limit and `value <= lloq` or `value` `NA` is censored (M3). A
#'   reading without a limit whose value is not finite is dropped with a
#'   warning. An optional logical column `absolute` marks the absolute readings.
#' @param model Expression of the prediction of a relative reading in the
#'   symbols of the formulas and data columns (covariates), e.g.
#'   `"yi * sj + bj"`.
#' @param errmodel Expression of the standard deviation of a relative reading;
#'   `prediction` stands for
#'   the prediction of the reading, e.g. `"sqrt(sigmaA^2 + (sigmaR * prediction)^2)"`.
#' @param biological,scaling,error Formulas, or lists of formulas, with one or
#'   more symbols on the left, `sj + bj ~ gel`, and their predictor on the
#'   right; all three are required, `error = NULL` for an error model without
#'   symbols, a known standard deviation. `biological` names the values that
#'   are reported. Random terms group by a column or a combination of
#'   columns, `(1 | gel:lane)`; terms of several groupings, nested or crossed,
#'   are normal and integrated per connected component, see [em()].
#' @param absolute Expression of the prediction of an absolute reading in the
#'   biological symbols, on the scale their `trafo` gives them, and data
#'   columns, e.g. `"yi + log2(1000)"` for readings in pM of values in nM.
#'   `NULL` (default) is the first biological symbol. A scaling symbol in it
#'   is an error.
#' @param trafo Character vector of `"symbol <- expression"`, the scale of a
#'   symbol in the predictor of that symbol, e.g. `c("yi <- 10^yi", "sj <- 10^sj")`. Symbols
#'   not named are their predictor.
#' @param technical Names of the columns of `data` that tell technical
#'   replicates apart, such as the lanes of one sample on a gel. Readings that
#'   agree in all other columns but `value` and `lloq` are one sample and are
#'   replaced by their mean before the alignment. The mean is geometric,
#'   taken on the log scale, for an error model in `prediction` (relative
#'   errors), and arithmetic otherwise, as for readings already on a log
#'   scale. A sample with censored readings
#'   is censored at the mean of the sample readings with the limits in place of the
#'   censored ones. `NULL` (default): every reading is a sample.
#' @param distribution `"normal"` (default) or `"t"` (Student-t, [normT()]).
#' @param df Degrees of freedom for `distribution = "t"`.
#' @param prior Optional prior of parameters on their predictor scale, a
#'   [constraintL2()] that [align()] adds to the objective. The log standard
#'   deviation of the normal random effect of a symbol `s` without correlation is
#'   `omega_s_s`, e.g. `constraintL2(mu = c(omega_bj_bj = 0), sigma = 0.3)`.
#'
#' @return Object of class `alignmentmodel`: the data, one row per sample
#'   with `technical`, the design of every symbol, the compiled expressions,
#'   the prior and the `normalization`: the block of every reading and per
#'   block the coefficients of the mean.
#' @seealso [align()]
#' @references Kemmer, S., Bang, S., Rosenblatt, M., Timmer, J. and Kaschek, D. (2022).
#'   BlotIt: Optimal alignment of Western blot and qPCR experiments.
#'   \emph{PLOS ONE} 17(8), e0264295. \doi{10.1371/journal.pone.0264295}
#' @example inst/examples/alignment.R
#' @export
alignmentModel <- function(data, model, errmodel, biological, scaling, error, absolute = NULL, trafo = NULL,
                           technical = NULL, distribution = c("normal", "t"), df = NULL, prior = NULL) {
  distribution <- match.arg(distribution)
  if (missing(model) || missing(errmodel))
    stop("alignmentModel: `model` and `errmodel` have no default; write both.", call. = FALSE)
  if (missing(biological) || missing(scaling) || missing(error))
    stop("alignmentModel: give the formulas of all three roles, `biological`, `scaling` and `error` ",
         "(`error = NULL` for an error model without symbols).", call. = FALSE)
  if (distribution == "t" && (is.null(df) || !is.finite(df) || df <= 0))
    stop("alignmentModel: distribution = \"t\" needs a positive `df`.", call. = FALSE)
  if (!is.null(prior) && !inherits(prior, "objfn"))
    stop("alignmentModel: `prior` is a constraintL2().", call. = FALSE)
  data <- as.data.frame(data)
  if (!"value" %in% names(data)) stop("alignmentModel: `data` needs a column `value`.", call. = FALSE)
  lloq <- if ("lloq" %in% names(data)) as.numeric(data$lloq) else rep(NA_real_, nrow(data))
  isAbs <- if ("absolute" %in% names(data)) data$absolute else rep(FALSE, nrow(data))
  if (!is.logical(isAbs) || anyNA(isAbs))
    stop("alignmentModel: the column `absolute` is TRUE or FALSE in every row.", call. = FALSE)
  if (!is.null(absolute) && !any(isAbs))
    stop("alignmentModel: `absolute` is the prediction of the rows marked in a column `absolute`, ",
         "which the data do not have.", call. = FALSE)
  censored <- !is.na(lloq) & (is.na(data$value) | data$value <= lloq)
  bad <- !censored & !is.finite(data$value)
  if (any(bad)) {
    warning(sprintf("alignmentModel: %d readings not finite dropped.", sum(bad)), call. = FALSE)
    data <- data[!bad, , drop = FALSE]; lloq <- lloq[!bad]; censored <- censored[!bad]; isAbs <- isAbs[!bad]
  }
  rownames(data) <- NULL

  # the predictor of every symbol over the readings
  # the symbols of one formula share the hyperparameters of the non-normal random terms of the formula
  if (is.null(biological) || is.null(scaling))
    stop("alignmentModel: `biological` and `scaling` need formulas.", call. = FALSE)
  roles <- list(biological = biological, scaling = scaling, error = error)
  roles <- roles[!vapply(roles, is.null, TRUE)]
  forms <- unlist(lapply(names(roles), function(r) {
    f <- roles[[r]]
    lapply(if (inherits(f, "formula")) list(f) else f, function(x) list(formula = x, role = r))
  }), recursive = FALSE)
  parsed <- unlist(lapply(seq_along(forms), function(k)
    lapply(.parseEffectFormula(forms[[k]]$formula), function(pf)
      c(pf[names(pf) != "origin"], origin = k, role = forms[[k]]$role))), recursive = FALSE)
  symbols <- vapply(parsed, `[[`, "", "symbol")
  if (anyDuplicated(symbols))
    stop("alignmentModel: a symbol is assigned twice: ", symbols[duplicated(symbols)][1], call. = FALSE)
  tech <- NULL
  if (length(technical)) {
    if (!is.character(technical) || !all(technical %in% names(data)) ||
        any(technical %in% c("value", "lloq", "absolute")))
      stop("alignmentModel: `technical` names columns of `data` other than value, lloq and absolute.",
           call. = FALSE)
    read <- intersect(technical, c(unlist(lapply(forms, function(f) all.vars(f$formula))),
                                   getSymbols(c(model, errmodel, absolute))))
    if (length(read))
      stop("alignmentModel: the technical column ", read[1], " is in a formula or expression; ",
           "technical replicates are averaged before the alignment.", call. = FALSE)
    geometric <- "prediction" %in% getSymbols(errmodel)
    tech <- list(columns = technical, readings = nrow(data), mean = if (geometric) "geometric" else "arithmetic")
    data <- .technicalMeans(data, technical, lloq, censored, geometric)
    lloq <- if ("lloq" %in% names(data)) data$lloq else rep(NA_real_, nrow(data))
    censored <- !is.na(lloq) & (is.na(data$value) | data$value <= lloq)
    isAbs <- if ("absolute" %in% names(data)) data$absolute else rep(FALSE, nrow(data))
  }
  # a bare column is a factor, one level per value; I(x) keeps it a covariate
  coded <- data
  for (v in intersect(unlist(lapply(parsed, function(pf) .bareVars(pf$fixed))), names(data)))
    if (is.numeric(data[[v]])) coded[[v]] <- factor(data[[v]], levels = sort(unique(data[[v]])))
  design <- setNames(lapply(parsed, function(pf) {
    em <- if (pf$role == "scaling" && any(isAbs)) .relativeDesign(pf, coded, isAbs)
          else .effectMatrices(pf, coded, sparse = TRUE)
    c(em, origin = pf$origin)
  }), symbols)
  role <- setNames(vapply(parsed, `[[`, "", "role"), symbols)

  # the scale of every symbol, substituted into the expressions
  tr <- setNames(symbols, symbols)
  for (t in trafo) {
    s <- .decolonize(.parseSubstitution(t))
    if (!s[1] %in% symbols)
      stop("alignmentModel: `trafo` names ", s[1], ", not a symbol of the formulas.", call. = FALSE)
    tr[s[1]] <- s[2]
  }
  sub <- function(e) replaceSymbols(symbols, paste0("(", tr, ")"), e)
  errmodel2 <- replaceSymbols("prediction", paste0("(", model, ")"), errmodel)
  bio <- symbols[role == "biological"]
  if (any(isAbs)) {
    absolute <- absolute %||% bio[1]
    off <- intersect(getSymbols(absolute), setdiff(symbols, bio))
    if (length(off))
      stop("alignmentModel: `absolute` is the prediction in the biological symbols; ", paste(off, collapse = ", "),
           if (all(role[off] == "scaling")) " scaling." else " not biological.", call. = FALSE)
    errmodelAbs <- replaceSymbols("prediction", paste0("(", absolute, ")"), errmodel)
  }
  used <- unique(c(getSymbols(model), getSymbols(errmodel2), if (any(isAbs)) getSymbols(errmodelAbs)))
  miss <- setdiff(used, c(symbols, names(data)))
  if (length(miss))
    stop("alignmentModel: the expressions use ", paste(miss, collapse = ", "),
         ", neither a symbol of the formulas nor a data column.", call. = FALSE)
  unused <- setdiff(symbols, used)
  if (length(unused))
    stop("alignmentModel: the symbols ", paste(unused, collapse = ", "), " are not in the expressions.",
         call. = FALSE)
  covariates <- setdiff(used, symbols)

  random <- unlist(lapply(design, `[[`, "random"), recursive = FALSE)
  groupings <- unique(vapply(random, `[[`, "", "group"))
  if (length(groupings) > 1L && any(vapply(random, `[[`, "", "dist") != "normal"))
    stop("alignmentModel: random terms over several groupings are normal.", call. = FALSE)

  m <- list(data = data, y = ifelse(censored, lloq, data$value), lloq = lloq, censored = censored, absolute = isAbs,
            symbols = data.frame(symbol = symbols, role = role, stringsAsFactors = FALSE),
            design = design, trafo = tr, covariates = covariates, model = model, errmodel = errmodel,
            absoluteModel = if (any(isAbs)) absolute,
            fn = .alignmentFn(sub(model), sub(errmodel2), c(symbols, covariates)),
            fnAbsolute = if (any(isAbs)) .alignmentFn(sub(absolute), sub(errmodelAbs), c(symbols, covariates)),
            grouping = if (length(groupings)) groupings, distribution = distribution, df = df,
            prior = prior, technical = tech)
  m$coef <- .alignmentCoefficients(m)
  m$normalization <- .alignmentNormalizationDesign(m)
  ap <- .alignmentPars(m)
  unknown <- setdiff(attr(prior, "parameters"), ap$fixed)
  unknown <- unknown[!(length(groupings) > 0 & grepl("^omega_", unknown))]
  if (length(unknown))
    stop("alignmentModel: the prior names ", paste(unknown, collapse = ", "),
         ", no parameter of the model.", call. = FALSE)
  class(m) <- "alignmentmodel"
  m
}

## The design of a scaling formula over the relative readings alone, columns NA or not in
## the absolute ones: those rows have no coefficient and no group of the random terms.
.relativeDesign <- function(pf, coded, isAbs) {
  rel <- which(!isAbs)
  sub <- .effectMatrices(pf, droplevels(coded[rel, , drop = FALSE]), sparse = TRUE)
  T <- methods::as(methods::as(sub$X, "generalMatrix"), "TsparseMatrix")
  X <- Matrix::sparseMatrix(i = rel[T@i + 1L], j = T@j + 1L, x = T@x, dims = c(nrow(coded), ncol(sub$X)),
                            dimnames = list(NULL, colnames(sub$X)))
  random <- lapply(sub$random, function(r) {
    z <- numeric(nrow(coded)); groups <- etas <- rep(NA_character_, nrow(coded))
    z[rel] <- r$z; groups[rel] <- r$groups; etas[rel] <- r$etas
    r$z <- z; r$groups <- groups; r$etas <- etas
    r
  })
  list(X = X, random = random, columns = sub$columns)
}

## One row per sample, the readings that agree in all columns but `value`, `lloq` and `technical`:
## their mean, on the log scale with `geometric`. A sample with censored readings is censored at
## the mean of the sample readings with the limits in place of the censored ones, a bound of the sample mean.
.technicalMeans <- function(data, technical, lloq, censored, geometric) {
  key <- setdiff(names(data), c("value", "lloq", technical))
  id <- if (length(key)) do.call(paste, c(lapply(data[key], function(x) ifelse(is.na(x), "\rNA", as.character(x))),
                                          sep = "\r"))
        else rep("1", nrow(data))
  sample <- factor(id, unique(id))
  bound <- ifelse(censored, lloq, data$value)
  if (geometric && any(bound <= 0))
    stop("alignmentModel: the geometric mean of technical replicates, for an error model in `prediction`, ",
         "needs positive readings and positive limits of the censored ones.", call. = FALSE)
  to <- if (geometric) function(x) log(pmax(x, 0)) else identity
  from <- if (geometric) exp else identity
  avg <- function(x) from(as.vector(tapply(to(x), sample, mean)))
  cens <- as.vector(tapply(censored, sample, any))
  out <- data[!duplicated(id), setdiff(names(data), technical), drop = FALSE]
  out$value <- ifelse(cens, NA_real_, avg(bound))
  if ("lloq" %in% names(data)) out$lloq <- ifelse(cens, avg(bound), avg(lloq))
  rownames(out) <- NULL
  out
}

## The variables of a formula that stand outside a function call such as factor().
.bareVars <- function(e) {
  if (is.name(e)) return(as.character(e))
  if (!is.call(e) || !as.character(e[[1]]) %in% c("~", "+", "-", "*", ":", "/", "("))
    return(character(0))
  unique(unlist(lapply(as.list(e)[-1], .bareVars)))
}

#' @export
#' @rdname alignmentModel
#' @param x An `alignmentmodel`.
#' @param ... Ignored.
print.alignmentmodel <- function(x, ...) {
  cat(sprintf("alignmentModel: %d readings, model %s, error %s\n", nrow(x$data), x$model, x$errmodel))
  for (i in seq_len(nrow(x$symbols))) {
    s <- x$symbols$symbol[i]; d <- x$design[[s]]
    cat(sprintf("  %-10s %-10s %4d %s%s%s\n", s, x$symbols$role[i], ncol(d$X),
                if (ncol(d$X) == 1L) "parameter" else "parameters",
                if (x$trafo[[s]] != s) paste0(", ", s, " = ", x$trafo[[s]]) else "",
                if (length(d$random)) paste0(", random over ", d$random[[1]]$group) else ""))
  }
  if (!is.null(x$technical))
    cat(sprintf("  technical replicates over %s: %d readings, their %s means\n",
                paste(x$technical$columns, collapse = ", "), x$technical$readings, x$technical$mean))
  if (any(x$censored)) cat(sprintf("  %d readings censored at their lloq\n", sum(x$censored)))
  if (length(x$covariates)) cat("  covariates:", paste(x$covariates, collapse = ", "), "\n")
  if (!is.null(x$prior)) cat("  prior on", paste(attr(x$prior, "parameters"), collapse = ", "), "\n")
  nb <- nrow(x$normalization$A)
  if (nb) cat(sprintf("  mean of %s normalized in %d %s of experiments\n", x$normalization$symbol, nb,
                      if (nb == 1L) "block" else "blocks"))
  if (any(x$absolute)) {
    nk <- length(unique(x$normalization$block[x$absolute]))
    nb <- length(unique(x$normalization$block))
    cat(sprintf("  %d absolute readings, model %s: %d of %d %s on their scale\n", sum(x$absolute),
                x$absoluteModel, nk, nb, if (nb == 1L) "block" else "blocks"))
  }
  invisible(x)
}


# Compiled expressions, one per model, error model and variables in the session, built in the
# output directory of the session, in a folder per host and process: forks that build the same
# expressions at once share no file.
.alignmentCache <- new.env(parent = emptyenv())
.alignmentFn <- function(model, errmodel, variables) {
  key <- substr(digest::digest(list(model, errmodel, variables)), 1, 10)
  if (!is.null(.alignmentCache[[key]])) return(.alignmentCache[[key]])
  outdir <- file.path(.dmodOutdir(), "align", paste0(Sys.info()[["nodename"]], "_", Sys.getpid()))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  .alignmentCache[[key]] <- cppDE::cppFUN(c(prediction = model, sigma = errmodel), variables = variables,
                                          modelname = paste0("align_", key), outdir = outdir,
                                          compile = TRUE, deriv2 = TRUE, convenient = FALSE)
}

## The fixed and random parameters of the model.
.alignmentPars <- function(m) {
  fixed <- unique(unlist(lapply(m$design, function(d) colnames(d$X)), use.names = FALSE))
  eta <- unique(unlist(lapply(m$design, function(d) lapply(d$random, `[[`, "etas")), use.names = FALSE))
  list(fixed = fixed, eta = eta)
}

## Per symbol the coefficients of every reading in every parameter, a sparse Matrix: the
## predictor is coef %*% pars. Scaling symbols have none on absolute readings.
.alignmentCoefficients <- function(m) {
  ap <- .alignmentPars(m)
  all <- c(ap$fixed, ap$eta)
  n <- nrow(m$data)
  scaling <- m$symbols$symbol[m$symbols$role == "scaling"]
  Map(function(d, s) {
    X <- methods::as(methods::as(Matrix::Matrix(d$X, sparse = TRUE), "generalMatrix"), "TsparseMatrix")
    i <- X@i + 1L; j <- match(colnames(d$X), all)[X@j + 1L]; x <- X@x
    for (r in d$random) {
      k <- which(r$z != 0)
      i <- c(i, k); j <- c(j, match(r$etas[k], all)); x <- c(x, r$z[k])
    }
    keep <- !(s %in% scaling & m$absolute[i])
    Matrix::sparseMatrix(i = i[keep], j = j[keep], x = x[keep], dims = c(n, length(all)),
                         dimnames = list(NULL, all))
  }, m$design, names(m$design))
}

## The predictor of every symbol at the readings `rows`.
.alignmentPredictors <- function(m, rows, pars, coef = .rowCoefficients(m, rows)) {
  matrix(vapply(coef, function(C) {
    p <- pars[intersect(colnames(C), names(pars))]
    as.numeric(C[, names(p), drop = FALSE] %*% p)
  }, numeric(length(rows))), length(rows), dimnames = list(NULL, m$symbols$symbol))
}

## The coefficients of every symbol at the readings `rows`.
.rowCoefficients <- function(m, rows) lapply(m$coef, function(C) C[rows, , drop = FALSE])

## The readings `rows` for their evaluation: the parameters they depend on, per symbol the
## columns of the parameters of the symbol among them with the coefficients there as a dense matrix, and the
## covariates. `coef` holds the coefficients of the rows.
.alignmentLocal <- function(m, rows, coef = .rowCoefficients(m, rows)) {
  pn <- colnames(coef[[1]])[Matrix::colSums(Reduce(`+`, lapply(coef, abs))) > 0]
  coef <- lapply(coef, function(C) C[, pn, drop = FALSE])
  cols <- lapply(coef, function(C) which(Matrix::colSums(C != 0) > 0))
  list(pars = pn, cols = cols, dense = Map(function(C, k) as.matrix(C[, k, drop = FALSE]), coef, cols),
       covariates = as.matrix(m$data[rows, m$covariates, drop = FALSE]), absolute = m$absolute[rows])
}

## The local descriptions of the readings of every condition, as .alignmentLocal(), from the
## nonzero coefficients split once by condition.
.alignmentLocals <- function(m, cond, rowsOf) {
  pos <- integer(nrow(m$data))
  pos[unlist(rowsOf, use.names = FALSE)] <- unlist(lapply(rowsOf, seq_along), use.names = FALSE)
  trip <- lapply(m$coef, function(C) {
    T <- methods::as(methods::as(C, "generalMatrix"), "TsparseMatrix")
    list(i = T@i + 1L, j = T@j + 1L, x = T@x)
  })
  byCondition <- lapply(trip, function(t) split(seq_along(t$i), factor(cond[t$i], names(rowsOf))))
  pnames <- colnames(m$coef[[1]])
  covariates <- as.matrix(m$data[, m$covariates, drop = FALSE])
  lapply(setNames(nm = names(rowsOf)), function(cn) {
    r <- rowsOf[[cn]]
    js <- sort(unique(unlist(lapply(seq_along(trip), function(a) trip[[a]]$j[byCondition[[a]][[cn]]]))))
    cols <- dense <- vector("list", length(trip))
    for (a in seq_along(trip)) {
      k <- byCondition[[a]][[cn]]
      ja <- sort(unique(trip[[a]]$j[k]))
      D <- matrix(0, length(r), length(ja))
      D[cbind(pos[trip[[a]]$i[k]], match(trip[[a]]$j[k], ja))] <- trip[[a]]$x[k]
      cols[[a]] <- match(ja, js)
      dense[[a]] <- D
    }
    list(pars = pnames[js], cols = setNames(cols, names(m$coef)), dense = setNames(dense, names(m$coef)),
         covariates = covariates[r, , drop = FALSE], absolute = m$absolute[r])
  })
}

## The readings `k` of a local description.
.localRows <- function(local, k) {
  local$dense <- lapply(local$dense, function(D) D[k, , drop = FALSE])
  local$covariates <- local$covariates[k, , drop = FALSE]
  local$absolute <- local$absolute[k]
  local
}

## Prediction and sigma of the readings `rows` with first and second derivatives in the parameters they
## depend on, from the compiled derivatives in the predictors times their coefficients; `local` describes the
## rows (.alignmentLocal()), without it the values alone come from the sparse coefficients.
.alignmentEval <- function(m, rows, pars, deriv = TRUE, deriv2 = FALSE, local = NULL) {
  sym <- m$symbols$symbol
  n <- length(rows)
  if (!deriv && !deriv2 && is.null(local)) {
    vars <- cbind(.alignmentPredictors(m, rows, pars), as.matrix(m$data[rows, m$covariates, drop = FALSE]))
    colnames(vars) <- c(sym, m$covariates)
    return(list(y = .alignmentRun(m, m$absolute[rows], function(fn, k)
      fn$func(vars = .rowsOf(vars, k), params = numeric(0)))))
  }
  local <- local %||% .alignmentLocal(m, rows)
  ev <- .alignmentSymbolEval(m, local, pars, n, deriv || deriv2, deriv2)
  if (!deriv && !deriv2) return(list(y = ev))
  used <- local$pars
  P <- length(used)
  J <- array(0, c(n, 2L, P), dimnames = list(NULL, c("prediction", "sigma"), used))
  H <- if (deriv2) array(0, c(n, 2L, P, P), dimnames = list(NULL, c("prediction", "sigma"), used, used))
  for (a in seq_along(sym)) {
    ia <- local$cols[[a]]
    if (!length(ia)) next
    Da <- local$dense[[a]]
    for (o in 1:2) J[, o, ia] <- J[, o, ia] + ev$tangent[, o, a] * Da
    if (!deriv2) next
    for (b in seq_along(sym)) {
      ib <- local$cols[[b]]
      if (!length(ib)) next
      DD <- Da[, rep(seq_along(ia), length(ib)), drop = FALSE] * local$dense[[b]][, rep(seq_along(ib), each = length(ia)), drop = FALSE]
      for (o in 1:2) {
        cur <- H[, o, ia, ib, drop = FALSE]
        H[, o, ia, ib] <- cur + array(ev$hessian[, o, a, b] * DD, dim(cur))
      }
    }
  }
  list(y = ev$y, J = J, H = H)
}

## The compiled model at the `n` readings of `local`: the values, with `deriv` also the
## derivatives in the symbols, `tangent` [n, 2, symbol] and with `deriv2` `hessian`
## [n, 2, symbol, symbol].
.alignmentSymbolEval <- function(m, local, pars, n, deriv, deriv2 = FALSE) {
  sym <- m$symbols$symbol
  i <- match(local$pars, names(pars))
  p <- pars[i]
  p[is.na(i)] <- 0
  lin <- vapply(seq_along(sym), function(a) as.numeric(local$dense[[a]] %*% p[local$cols[[a]]]), numeric(n))
  vars <- cbind(matrix(lin, n), local$covariates)
  colnames(vars) <- c(sym, m$covariates)
  if (!deriv) return(.alignmentRun(m, local$absolute, function(fn, k) fn$func(vars = .rowsOf(vars, k), params = numeric(0))))
  tX <- array(0, c(n, ncol(vars), length(sym)), dimnames = list(NULL, colnames(vars), sym))
  for (a in seq_along(sym)) tX[, a, a] <- 1
  .alignmentRun(m, local$absolute, function(fn, k)
    fn$evaluate(vars = .rowsOf(vars, k), params = numeric(0), tangentX = .rowsOf(tX, k), deriv2 = deriv2))
}

## The compiled expressions at readings of which `isAbs` marks the absolute ones: `f(fn, rows)`
## with the model on the relative rows and the expression of `absolute` on the others, the
## results stacked in the order of the readings.
.alignmentRun <- function(m, isAbs, f) {
  if (!any(isAbs)) return(f(m$fn, NULL))
  if (all(isAbs)) return(f(m$fnAbsolute, NULL))
  rows <- list(which(!isAbs), which(isAbs))
  .joinRows(list(f(m$fn, rows[[1]]), f(m$fnAbsolute, rows[[2]])), rows, length(isAbs))
}

## The rows `k` of a matrix or array, all for NULL.
.rowsOf <- function(x, k) {
  if (is.null(k)) return(x)
  if (length(dim(x)) == 3L) x[k, , , drop = FALSE] else x[k, , drop = FALSE]
}

## Arrays, or lists of them, with rows `rows[[i]]` of `n` in their first dimension, joined.
.joinRows <- function(parts, rows, n) {
  first <- parts[[1]]
  if (is.list(first) && is.null(dim(first)))
    return(setNames(lapply(names(first), function(nm) .joinRows(lapply(parts, `[[`, nm), rows, n)), names(first)))
  if (is.null(first)) return(NULL)
  d <- dim(first) %||% length(first)
  out <- matrix(0, n, prod(d[-1]))
  for (i in seq_along(parts)) out[rows[[i]], ] <- matrix(parts[[i]], length(rows[[i]]))
  if (length(d) == 1L) return(setNames(as.vector(out), NULL))
  array(out, c(n, d[-1]), dimnames = c(list(NULL), dimnames(first)[-1]))
}


## ---- scales and normalization ------------------------------------------------------------

## A symbol's value from the predictor and the derivative in it, by the trafo expression of the symbol.
.trafoValue <- function(m, s, lin) {
  e <- str2lang(m$trafo[[s]])
  env <- setNames(list(lin), s)
  slope <- eval(stats::D(e, s), env)
  list(value = eval(e, env), slope = rep_len(slope, length(lin)))
}

## The predictor at which a symbol has the value `target`, by root finding on the trafo; NA
## where none is found.
.trafoInverse <- function(m, s, target) {
  f <- function(v) .trafoValue(m, s, v)$value - target
  for (b in c(1, 10, 100)) {
    lo <- f(-b); hi <- f(b)
    if (is.finite(lo) && is.finite(hi) && sign(lo) != sign(hi))
      return(stats::uniroot(f, c(-b, b), tol = 1e-12)$root)
  }
  tryCatch(stats::uniroot(f, c(-1, 1), extendInt = "yes", tol = 1e-12, maxiter = 5000L)$root,
           error = function(e) NA_real_, warning = function(w) NA_real_)
}

## The blocks of experiments connected by shared values of the first biological symbol, as in
## blotIt, and per block without absolute readings the coefficients of the mean of the block
## predictors over those values.
.alignmentNormalizationDesign <- function(m) {
  s <- m$symbols$symbol[m$symbols$role == "biological"][1]
  key <- function(cols) if (length(cols)) do.call(paste, c(lapply(m$data[cols], as.character), sep = "\r"))
                        else rep("1", nrow(m$data))
  bio <- key(.alignmentColumns(m, s))
  # the fixed parts of the scaling formulas: random effects have mean 0 and shift no block
  scaling <- m$symbols$symbol[m$symbols$role == "scaling"]
  cols <- unique(unlist(lapply(scaling, function(x) m$design[[x]]$columns)))
  block <- .connectedBlocks(bio, key(cols))
  X <- m$design[[s]]$X
  relative <- setdiff(sort(unique(block)), unique(block[m$absolute]))
  A <- C <- matrix(0, 0, ncol(X), dimnames = list(NULL, colnames(X)))
  if (length(relative)) {
    first <- lapply(relative, function(b) {
      rows <- which(block == b)
      rows[!duplicated(bio[rows])]
    })
    A <- do.call(rbind, lapply(first, function(r) Matrix::colMeans(X[r, , drop = FALSE])))
    C <- do.call(rbind, lapply(first, function(r) as.numeric(X[r[1], ])))
    colnames(C) <- colnames(X)
  }
  list(symbol = s, block = block, A = A, C = C)
}

## Connected components of the bipartite graph of values and experiments, the component of every
## reading.
.connectedBlocks <- function(values, experiments) {
  v <- match(values, unique(values)); e <- match(experiments, unique(experiments))
  label <- seq_len(max(v))
  repeat {
    byExperiment <- tapply(label[v], e, min)
    new <- pmin(label, as.vector(tapply(byExperiment[e], v, min)))
    if (identical(new, label)) break
    label <- new
  }
  match(label[v], unique(label[v]))
}

## The normalization as a term of the objective: per block the residual 1000 times the mean of the
## predictors; asked for hessian = "sparse", the Hessian w a a' is w c c' over the predictor c of
## one value of the block plus the low-rank part w (a a' - c c').
.alignmentNormalization <- function(m) {
  A <- Matrix::Matrix(m$normalization$A, sparse = TRUE)
  C <- Matrix::Matrix(m$normalization$C %||% .alignmentNormalizationDesign(m)$C, sparse = TRUE)
  pn <- colnames(A)
  fn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, hessian = NULL, conditions = NULL,
                 env = NULL, cores = 1L) {
    p <- list(...)[[match.fnargs(list(...), "pars")]]
    r <- 1e3 * as.numeric(A %*% .localPars(pn, p, fixed)[pn])
    if (!deriv) return(objlist(value = sum(r^2), gradient = NULL, hessian = NULL))
    free <- intersect(pn, names(p))
    J <- 1e3 * A[, free, drop = FALSE]
    gradient <- setNames(numeric(length(p)), names(p))
    gradient[free] <- 2 * as.numeric(Matrix::crossprod(J, r))
    cv <- .resolveCurvature(deriv, deriv2, hessian, "forward")
    H <- NULL
    if (cv$sparse) {
      Jc <- 1e3 * C[, free, drop = FALSE]
      split <- which(Matrix::rowSums(abs(J - Jc)) > 0)
      out <- objlist(value = sum(r^2), gradient = gradient, hessian = .sparseOver(2 * Matrix::crossprod(Jc), names(p)))
      if (length(split))
        out$lowrank <- list(U = t(as.matrix(rbind(J[split, , drop = FALSE], Jc[split, , drop = FALSE]))),
                            d = rep(c(2, -2), each = length(split)))
      return(out)
    }
    if (cv$hessian) {
      Hf <- 2 * Matrix::crossprod(J)
      if (cv$sparse) H <- .sparseOver(Hf, names(p))
      else {
        H <- matrix(0, length(p), length(p), dimnames = list(names(p), names(p)))
        H[free, free] <- as.matrix(Hf)
      }
    }
    objlist(value = sum(r^2), gradient = gradient, hessian = H)
  }
  class(fn) <- c("objfn", "fn")
  attr(fn, "parameters") <- pn
  fn
}


## ---- objective --------------------------------------------------------------------

## Prediction and error model of the readings of one condition, the model in `shared`, an environment of all
## conditions that serialize() writes once. The prediction frame holds the error-model frame, which the error
## model returns when asked at the same parameters, split the same way, for no higher derivatives.
.alignmentOutputs <- function(shared, rows, condition, times, pars, fixed, deriv, deriv2) {
  m <- shared$m
  local <- shared$local[[condition]]
  at <- match(local$pars, names(pars))
  allPars <- .localPars(local$pars, pars, fixed, at)
  k <- match(times, rows); ok <- which(!is.na(k))
  if (length(ok) < length(rows) || any(k[ok] != seq_along(ok))) local <- .localRows(local, k[ok])
  ev <- .alignmentEval(m, rows[k[ok]], allPars, deriv || deriv2, deriv2, local)
  free <- local$pars[!is.na(at)]
  # tangents of the parameters: the derivatives along them, and the parameters keep them
  V <- attr(pars, "deriv", exact = TRUE)
  if (!is.null(V)) {
    V <- .localTangents(V, free)
    allPars <- as.parvec(allPars, deriv = rbind(V, matrix(0, length(allPars) - nrow(V), ncol(V),
                                                          dimnames = list(setdiff(names(allPars), free), NULL)))[names(allPars), , drop = FALSE])
  }
  frame <- function(o, empty) {
    y <- rep(empty, length(times)); y[ok] <- ev$y[, o]
    jac <- if (deriv || deriv2) .alignmentJacobianRows(ev$J[, o, , drop = FALSE], ok, length(times), free)
    jac2 <- if (deriv2) .alignmentHessianRows(ev$H[, o, , , drop = FALSE], ok, length(times), free)
    if (!is.null(V) && !is.null(jac)) {
      jac2 <- if (!is.null(jac2)) .chainSecond(jac2, V)
      jac <- array(matrix(jac, length(times)) %*% V, c(length(times), 1L, ncol(V)),
                   dimnames = list(NULL, "y", colnames(V)))
    }
    prdframe(cbind(time = times, y = y), deriv = jac, deriv2 = jac2, parameters = allPars)
  }
  out <- frame("prediction", 0)
  attr(out, "sigma") <- list(pars = unclass(allPars)[names(allPars)], cols = free, deriv = deriv || deriv2,
                             deriv2 = deriv2, tangents = V, frame = frame("sigma", 1))
  out
}

## The tangents `V` (parameters by directions) on the parameters `free`, zero where `V` has no row.
.localTangents <- function(V, free) {
  out <- matrix(0, length(free), ncol(V), dimnames = list(free, colnames(V)))
  hit <- intersect(free, rownames(V))
  out[hit, ] <- V[hit, , drop = FALSE]
  out
}

## Second derivatives [time, 1, p, p] chained with the tangents V: [time, 1, D, D] = V' H V.
.chainSecond <- function(H, V) {
  n <- dim(H)[1]; D <- ncol(V)
  out <- array(0, c(n, 1L, D, D), dimnames = list(NULL, "y", colnames(V), colnames(V)))
  for (r in seq_len(n)) out[r, 1, , ] <- crossprod(V, H[r, 1, , ] %*% V)
  out
}

## The vjp of an output (1 prediction, 2 sigma) of the readings `rows` at `times`: per seed J' w
## and sum_r w_r Hess(f_r) v_k + J' wdot_k with the tangents v_k of `pars`, contracted per pair
## of symbols with the coefficients of their predictors.
.alignmentVjp <- function(shared, rows, condition, times, pars, fixed, cotangent, o) {
  m <- shared$m
  local <- shared$local[[condition]]
  at <- match(local$pars, names(pars))
  allPars <- .localPars(local$pars, pars, fixed, at)
  k <- match(times, rows); ok <- which(!is.na(k))
  if (length(ok) < length(rows) || any(k[ok] != seq_along(ok))) local <- .localRows(local, k[ok])
  W <- .dropTime(cotangent)
  K <- dim(W)[3]
  w <- matrix(W[ok, 1L, , drop = FALSE], length(ok), K)
  bl <- .ctBlocks(K)
  out <- matrix(0, length(pars), K, dimnames = list(names(pars), NULL))
  if (o == 2L && !shared$sigmaVjp) return(out)
  # the vjps of prediction and sigma of a sweep meet the same point: the first evaluates and
  # holds the evaluation, the second takes it out
  key <- list(pars = allPars, rows = k[ok], deriv2 = bl$D > 0L)
  hold <- shared$vjpEval[[condition]]
  if (identical(hold$key, key)) {
    ev <- hold$ev
    shared$vjpEval[[condition]] <- NULL
  } else {
    ev <- .alignmentSymbolEval(m, local, allPars, length(ok), TRUE, bl$D > 0L)
    if (shared$sigmaVjp) shared$vjpEval[[condition]] <- list(key = key, ev = ev)
  }
  free <- local$pars[!is.na(at)]
  V <- if (bl$D > 0L) .localTangents(attr(pars, "deriv", exact = TRUE) %||% matrix(0, 0, bl$D),
                                     local$pars)[, seq_len(bl$D), drop = FALSE]
  ct <- matrix(0, length(local$pars), K)
  for (a in seq_along(m$symbols$symbol)) {
    ia <- local$cols[[a]]
    if (!length(ia)) next
    Da <- local$dense[[a]]
    ct[ia, ] <- ct[ia, ] + crossprod(Da, ev$tangent[, o, a] * w)
    if (bl$D == 0L) next
    for (b in seq_along(m$symbols$symbol)) {
      ib <- local$cols[[b]]
      if (!length(ib)) next
      Zb <- local$dense[[b]] %*% V[ib, , drop = FALSE]
      for (s in seq_len(bl$S))
        ct[ia, bl$rest[, s]] <- ct[ia, bl$rest[, s]] +
          crossprod(Da, (w[, bl$first[s]] * ev$hessian[, o, a, b]) * Zb)
    }
  }
  out[free, ] <- ct[match(free, local$pars), , drop = FALSE]
  out
}

## The values of the parameters `names` among `pars` and `fixed`, without copying either; `at`
## is the position of `names` in `pars`.
.localPars <- function(names, pars, fixed, at = match(names, names(pars))) {
  v <- .subset(pars, at)
  if (length(fixed)) {
    j <- match(names, names(fixed))
    take <- is.na(v) & !is.na(j)
    v[take] <- .subset(fixed, j[take])
  }
  keep <- !is.na(v)
  setNames(as.numeric(v[keep]), names[keep])
}

## The prediction of the readings of all conditions, which every call names; the times are rows.
## The batch entries evaluate the compiled model once at the rows of all conditions asked for.
.alignmentP2X <- function(shared, parNames) {
  P2X <- function(times, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE, condition = NULL) {
    cn <- .alignmentCondition(shared, condition, times)
    .alignmentOutputs(shared, shared$rowsOf[[cn]], cn, times, pars, fixed, deriv, deriv2)
  }
  attr(P2X, "batchfn") <- function(times, parsList, fixedList, deriv = TRUE, deriv2 = FALSE, cores = 1L,
                                   conditions = NULL)
    .alignmentOutputsMany(shared, times, parsList, fixedList, deriv, deriv2, cores, conditions)
  attr(P2X, "vjpfn") <- function(times, pars, fixed = NULL, cotangent, condition = NULL) {
    cn <- .alignmentCondition(shared, condition, times)
    .alignmentVjp(shared, shared$rowsOf[[cn]], cn, times, pars, fixed, cotangent, 1L)
  }
  attr(P2X, "vjpbatchfn") <- function(times, parsList, fixedList, cotangentList, conditions = NULL, cores = 1L,
                                      seeds = FALSE)
    .alignmentVjpMany(shared, times, parsList, fixedList, cotangentList, 1L, cores, seeds, conditions)
  # the compiled expressions, for loadDLL() in a new process
  fns <- Filter(Negate(is.null), list(shared$m$fn, shared$m$fnAbsolute))
  attr(P2X, "modelname") <- vapply(fns, attr, "", "modelname")
  attr(P2X, "compileInfo") <- do.call(.collectCompileInfo, fns)
  prdfn(P2X, parNames, NULL)
}

## The error model of the readings of all conditions: sigma depends on the parameters alone, the error
## model has the prediction substituted, so the vjp has no half on the prediction frame.
.alignmentX2E <- function(shared, parNames) {
  X2Y <- function(out, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE, condition = NULL) {
    times <- out[, "time"]
    cn <- .alignmentCondition(shared, condition, times)
    s <- attr(out, "sigma", exact = TRUE)
    if (!is.null(s) && (s$deriv || !(deriv || deriv2)) && (s$deriv2 || !deriv2)) {
      local <- shared$local[[cn]]
      allPars <- .localPars(local$pars, pars, fixed)
      V <- attr(pars, "deriv", exact = TRUE)
      if (identical(allPars, s$pars) && identical(intersect(local$pars, names(pars)), s$cols) &&
          identical(if (!is.null(V)) .localTangents(V, s$cols), s$tangents)) return(s$frame)
    }
    attr(.alignmentOutputs(shared, shared$rowsOf[[cn]], cn, times, pars, fixed, deriv, deriv2), "sigma")$frame
  }
  attr(X2Y, "batchfn") <- function(outList, parsList, fixedList, deriv = TRUE, deriv2 = FALSE, cores = 1L,
                                   conditions = NULL)
    Map(function(o, p, f, cn) X2Y(o, p, f, deriv, deriv2, cn), outList, parsList, fixedList,
        conditions %||% rep(list(NULL), length(outList)))
  attr(X2Y, "vjpfn") <- function(out, pars, fixed = NULL, cotangent, condition = NULL) {
    times <- out[, "time"]
    cn <- .alignmentCondition(shared, condition, times)
    .ct(pars = .alignmentVjp(shared, shared$rowsOf[[cn]], cn, times, pars, fixed, cotangent, 2L))
  }
  attr(X2Y, "vjpbatchfn") <- function(outList, parsList, fixedList, cotangentList, conditions = NULL, cores = 1L)
    lapply(.alignmentVjpMany(shared, lapply(outList, function(o) o[, "time"]), parsList, fixedList,
                             cotangentList, 2L, cores, FALSE, conditions), function(ct) .ct(pars = ct))
  obsfn(X2Y, parNames, NULL)
}

## The condition of a call: the one it names, else the one of the first row.
.alignmentCondition <- function(shared, condition, times)
  if (length(condition)) condition[[1L]] else shared$cond[[times[1L]]]

## The compiled model at the rows `rows` with the predictors `lin` [rows x symbol], with the
## derivatives in the symbols: `y` [n, 2], `tangent` [n, 2, symbol] and with `deriv2` `hessian`
## [n, 2, symbol, symbol].
.alignmentEvalRows <- function(m, rows, lin, deriv2) {
  sym <- m$symbols$symbol
  vars <- cbind(lin, as.matrix(m$data[rows, m$covariates, drop = FALSE]))
  colnames(vars) <- c(sym, m$covariates)
  tX <- array(0, c(length(rows), ncol(vars), length(sym)), dimnames = list(NULL, colnames(vars), sym))
  for (a in seq_along(sym)) tX[, a, a] <- 1
  .alignmentRun(m, m$absolute[rows], function(fn, k)
    fn$evaluate(vars = .rowsOf(vars, k), params = numeric(0), tangentX = .rowsOf(tX, k), deriv2 = deriv2))
}

## The conditions of the calls, NULL unless the times of every condition are the rows of the condition.
.alignmentRowConditions <- function(shared, times, conditions = NULL) {
  if (!length(times) || !is.list(times)) return(NULL)
  conds <- if (length(conditions) == length(times)) vapply(conditions, function(cn) cn[[1L]], "")
           else vapply(times, function(t) shared$cond[[t[1L]]], "")
  same <- vapply(seq_along(times), function(j) identical(as.integer(times[[j]]), shared$rowsOf[[conds[j]]]), TRUE)
  if (all(same)) conds
}

## The blocks of the conditions `conds` for the kernels, their rows laid out one after the other.
## Kept for the last set of conditions.
.alignmentBlocks <- function(shared, conds) {
  if (identical(shared$blocks$conds, conds)) return(shared$blocks$blocks)
  start <- c(0L, cumsum(vapply(conds, function(cn) nrow(shared$local[[cn]]$covariates), 0L)))
  blocks <- lapply(seq_along(conds), function(j) {
    local <- shared$local[[conds[j]]]
    list(start = start[j], n = start[j + 1L] - start[j], P = length(local$pars), cols = unname(local$cols),
         dense = unname(local$dense))
  })
  shared$blocks <- list(conds = conds, blocks = blocks)
  blocks
}

## Per condition the positions of the condition parameters among the names `pn` of the parameters and `fn` of
## the fixed ones of a request. Kept for the last two sets of names, as the prediction and the
## error model of a sweep ask with their names each.
.alignmentPositions <- function(shared, conds, pn, fn) {
  key <- list(conds, pn, fn)
  for (hold in shared$positions) if (identical(hold$key, key)) return(hold$pos)
  lp <- lapply(conds, function(cn) shared$local[[cn]]$pars)
  flat <- unlist(lp, use.names = FALSE)
  grp <- factor(rep(seq_along(conds), lengths(lp)), seq_along(conds))
  pos <- Map(function(ip, ifx) list(ip = ip, ifx = ifx), split(match(flat, pn), grp), split(match(flat, fn), grp))
  shared$positions <- c(list(list(key = key, pos = pos)), shared$positions)[seq_len(min(2L, length(shared$positions) + 1L))]
  pos
}

## Per condition of a batch the condition parameters in the request: `theta`, the values of all
## parameters of the condition with 0 where the request has none, `keep` and `free`, the free ones, `pars`, the
## values it has, from the parameters before the fixed ones, and `V`, the tangents of the free ones.
.alignmentRequests <- function(shared, conds, parsList, fixedList) {
  p0 <- parsList[[1L]]; f0 <- fixedList[[1L]]
  same <- all(vapply(parsList, function(p) identical(names(p), names(p0)), TRUE)) &&
          all(vapply(fixedList, function(f) identical(names(f), names(f0)), TRUE))
  pos <- if (same) .alignmentPositions(shared, conds, names(p0), names(f0))
  # tangents shared by all requests: their rows found once for all conditions
  V0 <- attr(p0, "deriv", exact = TRUE)
  oneV <- !is.null(V0) && all(vapply(parsList, function(p) identical(attr(p, "deriv", exact = TRUE), V0), TRUE))
  rowsV <- if (oneV && same) lapply(pos, function(ps) match(names(p0)[ps$ip[!is.na(ps$ip)]], rownames(V0)))
  lapply(seq_along(conds), function(j) {
    lp <- shared$local[[conds[j]]]$pars
    p <- parsList[[j]]; f <- fixedList[[j]]
    ps <- if (same) pos[[j]] else list(ip = match(lp, names(p)), ifx = match(lp, names(f)))
    v <- .subset(unclass(p), ps$ip)
    fromFixed <- is.na(ps$ip) & !is.na(ps$ifx)
    v[fromFixed] <- .subset(unclass(f), ps$ifx[fromFixed])
    has <- !is.na(v)
    keep <- which(!is.na(ps$ip))
    V <- attr(p, "deriv", exact = TRUE)
    Vl <- if (!is.null(rowsV)) {
      out <- matrix(0, length(keep), ncol(V0), dimnames = list(lp[keep], colnames(V0)))
      hit <- !is.na(rowsV[[j]])
      out[hit, ] <- V0[rowsV[[j]][hit], , drop = FALSE]
      out
    } else if (!is.null(V)) .localTangents(V, lp[keep])
    list(theta = ifelse(has, v, 0), keep = keep, free = lp[keep], pars = setNames(v[has], lp[has]), V = Vl)
  })
}

## The prediction frames of many conditions, by one evaluation at all their rows where the times of
## every condition are the rows of the condition, else one by one.
.alignmentOutputsMany <- function(shared, times, parsList, fixedList, deriv, deriv2, cores, conditions = NULL) {
  conds <- .alignmentRowConditions(shared, times, conditions)
  if (is.null(conds))
    return(lapply(seq_along(times), function(j) {
      cn <- .alignmentCondition(shared, conditions[[j]], times[[j]])
      .alignmentOutputs(shared, shared$rowsOf[[cn]], cn, times[[j]], parsList[[j]], fixedList[[j]], deriv, deriv2)
    }))
  m <- shared$m
  d <- deriv || deriv2
  blocks <- .alignmentBlocks(shared, conds)
  req <- .alignmentRequests(shared, conds, parsList, fixedList)
  rows <- as.integer(unlist(times, use.names = FALSE))
  lin <- alignmentPredictors(blocks, lapply(req, `[[`, "theta"), nrow(m$symbols), as.integer(cores))
  ev <- .alignmentEvalRows(m, rows, lin, deriv2)
  if (!d) ev <- list(y = ev$y)
  for (j in seq_along(blocks)) blocks[[j]]$keep <- req[[j]]$keep
  jac <- if (d) alignmentForwardBlocks(ev$tangent, if (deriv2) ev$hessian, blocks,
                                       lapply(req, `[[`, "free"), as.integer(cores))
  lapply(seq_along(conds), function(j) {
    b <- blocks[[j]]; r <- req[[j]]
    allPars <- r$pars
    at <- b$start + seq_len(b$n)
    tm <- times[[j]]
    # tangents of the parameters: the derivatives along them, as in .alignmentOutputs()
    Vl <- r$V
    if (!is.null(Vl)) {
      others <- setdiff(names(allPars), r$free)
      allPars <- as.parvec(allPars, deriv = rbind(Vl, matrix(0, length(others), ncol(Vl), dimnames = list(others, NULL)))[
        names(allPars), , drop = FALSE])
    }
    frame <- function(o, nm) {
      jac1 <- if (d) jac[[j]][[nm]]
      jac2 <- if (deriv2) jac[[j]][[paste0(nm, "2")]]
      if (!is.null(Vl) && !is.null(jac1)) {
        jac2 <- if (!is.null(jac2)) .chainSecond(jac2, Vl)
        jac1 <- array(matrix(jac1, b$n) %*% Vl, c(b$n, 1L, ncol(Vl)), dimnames = list(NULL, "y", colnames(Vl)))
      }
      prdframe(cbind(time = tm, y = ev$y[at, o]), deriv = jac1, deriv2 = jac2, parameters = allPars)
    }
    out <- frame(1L, "prediction")
    attr(out, "sigma") <- list(pars = unclass(allPars)[names(allPars)], cols = r$free, deriv = d, deriv2 = deriv2,
                               tangents = Vl, frame = frame(2L, "sigma"))
    out
  })
}

## The vjps of output `o` (1 prediction, 2 sigma) of many conditions, [free x K] each: one evaluation at all
## rows where the times of a condition are the condition rows, shared by both outputs as in .alignmentVjp(), else one by
## one. With `seeds` every column is a first-order seed.
.alignmentVjpMany <- function(shared, times, parsList, fixedList, cotangentList, o, cores, seeds = FALSE,
                              conditions = NULL) {
  conds <- .alignmentRowConditions(shared, times, conditions)
  one <- function(j, w) {
    cn <- .alignmentCondition(shared, conditions[[j]], times[[j]])
    .alignmentVjp(shared, shared$rowsOf[[cn]], cn, times[[j]], parsList[[j]], fixedList[[j]], w, o)
  }
  if (is.null(conds))
    return(lapply(seq_along(times), function(j) {
      if (!seeds) return(one(j, cotangentList[[j]]))
      w <- cotangentList[[j]]
      do.call(cbind, lapply(seq_len(dim(w)[3]), function(k) one(j, w[, , k, drop = FALSE])))
    }))
  m <- shared$m
  K <- dim(cotangentList[[1L]])[3]
  bl <- if (seeds) list(S = K, D = 0L, first = seq_len(K), rest = matrix(0L, 0L, K)) else .ctBlocks(K)
  blocks <- .alignmentBlocks(shared, conds)
  req <- .alignmentRequests(shared, conds, parsList, fixedList)
  if (o == 2L && !shared$sigmaVjp)
    return(lapply(req, function(r) matrix(0, length(r$free), K, dimnames = list(r$free, NULL))))
  rows <- as.integer(unlist(times, use.names = FALSE))
  theta <- lapply(req, `[[`, "theta")
  key <- list(theta = theta, rows = rows, deriv2 = bl$D > 0L)
  hold <- shared$vjpEvalMany
  if (identical(hold$key, key)) {
    ev <- hold$ev
    shared$vjpEvalMany <- NULL
  } else {
    lin <- alignmentPredictors(blocks, theta, nrow(m$symbols), as.integer(cores))
    ev <- .alignmentEvalRows(m, rows, lin, bl$D > 0L)
    if (shared$sigmaVjp) shared$vjpEvalMany <- list(key = key, ev = ev)
  }
  for (j in seq_along(blocks)) blocks[[j]]$keep <- seq_len(blocks[[j]]$P)
  W <- lapply(seq_along(conds), function(j) matrix(.dropTime(cotangentList[[j]]), blocks[[j]]$n, K))
  # the tangents of all parameters of a condition, 0 for those it does not have free
  Vl <- if (bl$D > 0L) lapply(seq_along(conds), function(j) {
    lp <- shared$local[[conds[j]]]$pars
    out <- matrix(0, length(lp), bl$D)
    if (!is.null(req[[j]]$V)) out[req[[j]]$keep, ] <- req[[j]]$V[, seq_len(bl$D)]
    out
  })
  ct <- alignmentVjpBlocks(ev$tangent, if (bl$D > 0L) ev$hessian, as.integer(o), blocks, W, Vl %||% list(),
                           as.integer(bl$first), matrix(as.integer(bl$rest), nrow(bl$rest)), as.integer(cores))
  lapply(seq_along(req), function(j) {
    r <- req[[j]]
    matrix(ct[[j]][r$keep, , drop = FALSE], length(r$keep), K, dimnames = list(r$free, NULL))
  })
}

## [time, "y", parameter] over the parameters `free`, rows `ok` of `n` filled.
.alignmentJacobianRows <- function(J, ok, n, free) {
  out <- array(0, c(n, 1L, length(free)), dimnames = list(NULL, "y", free))
  out[ok, 1, ] <- J[, 1, free, drop = FALSE]
  out
}

.alignmentHessianRows <- function(H, ok, n, free) {
  out <- array(0, c(n, 1L, length(free), length(free)), dimnames = list(NULL, "y", free, free))
  out[ok, 1, , ] <- H[, 1, free, free, drop = FALSE]
  out
}

## The condition of every reading: the experiment, the columns of the fixed parts of the scaling formulas and
## the component of the random-term groups, cut so that each condition depends on few parameters. The name
## joins the parts set and not repeated, "absolute" for a reading without any.
.alignmentConditions <- function(m) {
  scaling <- m$symbols$symbol[m$symbols$role == "scaling"]
  cols <- unique(unlist(lapply(scaling, function(s) m$design[[s]]$columns)))
  keys <- c(list(.alignmentComponents(m)), lapply(m$data[cols], as.character))
  keys <- Filter(Negate(is.null), keys)
  if (!length(keys)) return(rep("all", nrow(m$data)))
  full <- do.call(paste, c(keys, sep = "_"))
  short <- character(nrow(m$data))
  for (j in seq_along(keys)) {
    k <- keys[[j]]
    add <- !is.na(k)
    for (p in keys[seq_len(j - 1L)]) add <- add & (is.na(p) | p != k)
    short[add] <- ifelse(nzchar(short[add]), paste(short[add], k[add], sep = "_"), k[add])
  }
  short[!nzchar(short)] <- "absolute"
  distinct <- all(lengths(lapply(split(full, short), unique)) == 1L)
  .alignmentChunks(m, .effectName(if (distinct) short else full))
}

## The conditions `cond` split into parts of at most `size` values of the first biological
## symbol each: the derivatives of a condition are dense over the parameters of the condition, so the memory of an
## evaluation grows with the readings times the square of the parameters per condition.
.alignmentChunks <- function(m, cond, size = 16L) {
  s <- m$symbols$symbol[m$symbols$role == "biological"][1]
  cols <- .alignmentColumns(m, s)
  if (!length(cols)) return(cond)
  bio <- do.call(paste, c(lapply(m$data[cols], as.character), sep = "\r"))
  part <- integer(length(cond))
  for (rows in split(seq_along(cond), factor(cond, unique(cond)))) {
    v <- match(bio[rows], unique(bio[rows]))
    part[rows] <- (v - 1L) %/% size + 1L
  }
  many <- cond %in% unique(cond[part > 1L])
  if (!any(many)) return(cond)
  out <- ifelse(many, paste0(cond, "_part", part), cond)
  # distinct parts keep distinct names
  key <- paste(cond, part, sep = "\r")
  first <- !duplicated(key)
  named <- setNames(make.unique(out[first], sep = "_"), key[first])
  unname(named[key])
}

## The connected component of every reading in the groups of the random terms, named by the
## group of the first grouping; NA for a reading in no group, NULL without random terms.
.alignmentComponents <- function(m) {
  r <- unlist(lapply(m$design, `[[`, "random"), recursive = FALSE)
  if (!length(r)) return(NULL)
  byGrouping <- lapply(split(r, vapply(r, `[[`, "", "group")), function(rs) rs[[1]]$groups)
  first <- r[[1]]$groups
  if (length(byGrouping) == 1L) return(first)
  label <- seq_along(first)
  repeat {
    old <- label
    for (g in byGrouping) {
      ok <- !is.na(g)
      low <- tapply(label[ok], g[ok], min)
      label[ok] <- pmin(label[ok], low[g[ok]])
    }
    if (identical(label, old)) break
  }
  out <- ifelse(is.na(first[label]), paste0("C", label), first[label])
  out[Reduce(`&`, lapply(byGrouping, is.na))] <- NA_character_
  out
}

## The objective of the readings, one condition per experiment and group, answered by one
## prediction and one error-model leaf for all of them.
.alignmentObjective <- function(m) {
  n <- nrow(m$data)
  cond <- .alignmentConditions(m)
  conds <- unique(cond)
  rowsOf <- split(seq_len(n), factor(cond, conds))
  shared <- new.env(parent = emptyenv())
  shared$m <- m
  shared$local <- .alignmentLocals(m, cond, rowsOf)
  # sigma has a vjp to answer where it depends on a symbol
  sigma <- vapply(c(m$model, m$absoluteModel), function(e)
    replaceSymbols("prediction", paste0("(", e, ")"), m$errmodel), "")
  shared$sigmaVjp <- any(getSymbols(sigma) %in% m$symbols$symbol)
  shared$cond <- cond
  shared$rowsOf <- rowsOf
  parNames <- unique(unlist(lapply(shared$local, `[[`, "pars"), use.names = FALSE))
  prd <- .fnWithConditions(.alignmentP2X(shared, parNames), conds)
  ef <- .fnWithConditions(.alignmentX2E(shared, parNames), conds)
  # sigma is evaluated with the prediction and travels in the prediction frame
  attr(ef, "sigmaInPrediction") <- TRUE
  dl <- as.datalist(data.frame(name = "y", time = seq_len(n), value = m$y, sigma = NA_real_,
                               lloq = ifelse(is.na(m$lloq), -Inf, m$lloq), condition = cond,
                               stringsAsFactors = FALSE))
  times <- lapply(rowsOf, as.numeric)
  if (m$distribution == "t") normT(dl, prd, errmodel = ef, df = m$df, times = times)
  else normL2(dl, prd, errmodel = ef, times = times)
}


## The random effects of the model as the random effects of em(): per grouping the groups, the
## conditions that hold readings of a group and the etas of a group; a condition holds the
## readings of several groups of a grouping crossed with the experiment of the condition.
.alignmentRandomEffects <- function(m) {
  cond <- .alignmentConditions(m)
  conds <- unique(cond)
  rows <- unlist(lapply(m$design, `[[`, "random"), recursive = FALSE)
  if (!length(rows)) return(NULL)
  origin <- unlist(lapply(m$design, function(d) rep(d$origin, length(d$random))))
  terms <- do.call(rbind, lapply(seq_along(rows), function(i) data.frame(
    symbol = rows[[i]]$symbol, origin = origin[i], term = rows[[i]]$term, column = rows[[i]]$column,
    group = rows[[i]]$group, id = rows[[i]]$id, dist = rows[[i]]$dist,
    prob = rows[[i]]$prob, concentration = rows[[i]]$concentration,
    tau = rows[[i]]$tau,
    stem = rows[[i]]$stem, stringsAsFactors = FALSE)))
  # per condition the value of a reading column where all readings of the condition share it, else NA
  single <- function(x) vapply(split(x, factor(cond, conds)), function(v) {
    u <- unique(v)
    if (length(u) == 1L) u else NA_character_
  }, "")
  etas <- as.data.frame(lapply(setNames(rows, terms$stem), function(r) single(r$etas)),
                        row.names = conds, stringsAsFactors = FALSE, optional = TRUE)
  groups <- lapply(setNames(rows, terms$stem), function(r) setNames(.effectName(single(r$groups)), conds))
  subjects <- etaNames <- list()
  for (gr in unique(terms$group)) {
    k <- which(terms$group == gr)
    g <- .effectName(rows[[k[1]]]$groups)
    use <- which(!is.na(g))
    at <- use[!duplicated(g[use])]
    subjects[[gr]] <- lapply(split(cond[use], factor(g[use], g[at])), unique)
    etaNames[[gr]] <- matrix(vapply(rows[k], function(r) r$etas[at], character(length(at))), length(at),
                             dimnames = list(g[at], terms$stem[k]))
  }
  .newRandomEffects(terms, etas, groups, character(0), subjects = subjects, etaNames = etaNames)
}


## ---- fit ----------------------------------------------------------------------------

#' Fit an alignment model
#'
#' @description Estimates the parameters of an [alignmentModel()] by maximum
#'   likelihood, the marginal one with random terms ([em()]), with the
#'   normalization of the values per block of experiments (see
#'   [alignmentModel()]) in the objective; the error parameters and the
#'   hyperparameters of normal random terms by default by restricted
#'   likelihood. The Hessian is sparse over the experiments, with the
#'   normalization as a low-rank part (see [objlist()]); [trust()] and [em()]
#'   take the sparse path as `sparse = "auto"` decides.
#'
#' @param model An `alignmentmodel`.
#' @param fits Number of fits from perturbed starts. Default 10.
#' @param cores Number of fits run in parallel; every fit runs on one core.
#'   Default 1.
#' @param method With random terms: `"focei"` (default) or `"quadrature"`, see
#'   [em()].
#' @param sd Standard deviation of the start perturbation on the parameter
#'   scale. Default 0.5.
#' @param iterlim Iteration limit of the optimizer; with spike-and-slab or
#'   Dirichlet-process terms the outer iterations of [em()]. Default 1000.
#' @param init Optional named start of some parameters; the intercepts of the
#'   biological symbols start where the symbol is the median reading, the other
#'   intercepts where their symbol is 1, the other parameters at 0. The
#'   coefficients of the first scaling symbol then start where the predictor
#'   of that symbol
#'   is, per experiment, the median of the predictors at which the prediction
#'   meets the relative readings.
#' @param reml `TRUE` (default): the error parameters, and the hyperparameters
#'   of normal random terms, by restricted likelihood, the parameters of the
#'   symbols of the prediction (`model` and `absolute` of the
#'   [alignmentModel()]) integrated out under a flat prior (Lindstrom and
#'   Bates 1990), in at most 20 rounds from the best maximum likelihood fit
#'   until the error parameters and hyperparameters move by less than 1e-4.
#'   Without random terms or with normal ones and `method = "focei"`; else
#'   maximum likelihood. Where the information of the values and scales is
#'   singular, at the maximum likelihood fit or at a round, the restricted
#'   likelihood is not defined, and the fit is by maximum likelihood with a
#'   warning. `FALSE`: maximum likelihood.
#'
#' @details The objective of an alignment is exact to rounding, so the fits stop
#'   at changes of the objective, predicted decreases and gradients below
#'   `sqrt(.Machine$double.eps)`, the inner fits of the random effects too.
#'
#' @return Object of class `alignment`: `aligned` (one row per combination of
#'   the columns of the biological formulas: the columns, `symbol`, `estimate`
#'   and `se` of the predictor, `value` and `sigma` on the scale of the
#'   symbol; with REML `se` from the information of the mean parameters alone,
#'   the error parameters held at their estimates), `scaled` (the readings with `prediction`, `sigma`, `residual`
#'   and `bloq`), `parameters` (with attribute `estimation`, `"REML"` or
#'   `"ML"`), `value` (-2 log L, of the restricted likelihood with REML, with
#'   the normalization and the prior), `npar` (one less per block without
#'   absolute readings), `AIC`
#'   (with REML comparable only between fits with the same parameters of the
#'   prediction), `values` (-2 log L of all fits, sorted), `multistart`
#'   (`value`, `converged` and `iterations` per fit, see [plotValues()];
#'   with REML the maximum likelihood fits the rounds start from),
#'   `remlRounds` (the rounds of the restricted likelihood, 0 without),
#'   `seSource` (`"Hessian"`, or `"Gauss-Newton"` where the reported Hessian
#'   of a fit with random terms is not positive definite: then the
#'   Gauss-Newton Schur complement of the expected information over the random effects),
#'   `fit`, `objfn`, `omega`, `random` and `randomEffects` with random terms,
#'   and the `model`. `scaled` also has `scaled`, the value of the biological symbol at which the prediction meets
#'   the reading with the other symbols at the fit (`NaN` where none does), and
#'   `scaledSigma`, the error model over the slope of the prediction there;
#'   `aligned` has `lloq`, the smallest such value of the censored readings of
#'   a combination.
#' @seealso [alignmentModel()], [plotAligned()]; the guide chapter
#'   [Alignment of relative data measurements](https://dmodverse.github.io/dMod2/articles/guide-alignment.html).
#' @references Kemmer, S., Bang, S., Rosenblatt, M., Timmer, J. and Kaschek, D. (2022).
#'   BlotIt: Optimal alignment of Western blot and qPCR experiments.
#'   \emph{PLOS ONE} 17(8), e0264295. \doi{10.1371/journal.pone.0264295}
#'
#'   Lindstrom, M. J. and Bates, D. M. (1990). Nonlinear mixed effects models
#'   for repeated measures data. \emph{Biometrics} 46(3), 673-687.
#'   \doi{10.2307/2532087}
#' @example inst/examples/alignment.R
#' @export
align <- function(model, fits = 10L, cores = 1L, method = c("focei", "quadrature"),
                  sd = 0.5, iterlim = 1000L, init = NULL, reml = TRUE) {
  method <- match.arg(method)
  m <- model
  if (!inherits(m, "alignmentmodel")) stop("align: `model` is an alignmentModel().", call. = FALSE)
  if (length(cores) != 1L || is.na(cores) || cores < 1)
    stop("align: `cores` is the number of fits run in parallel.", call. = FALSE)
  cores <- as.integer(cores)
  # more cores inside one fit gain little, so every fit runs on one core
  old <- options(dMod2.cores = 1L)
  on.exit(options(old), add = TRUE)
  tol <- sqrt(.Machine$double.eps)
  tc <- list(ftol = tol, mtol = tol, gtol = tol)
  obj <- .alignmentObjective(m) + .alignmentNormalization(m) + m$prior
  re <- .alignmentRandomEffects(m)
  free <- .alignmentPars(m)$fixed
  # the restricted likelihood without random terms and with normal ones by FOCEI
  reml <- isTRUE(reml) && (is.null(re) || (method == "focei" && all(re$terms$dist == "normal")))
  mean <- if (reml) intersect(.alignmentMean(m), free)
  # without error parameters and random terms the restricted likelihood is the likelihood
  if (reml && is.null(re) && !length(setdiff(free, mean))) reml <- FALSE
  st <- .alignmentStart(m, free)
  st[intersect(names(init), names(st))] <- init[intersect(names(init), names(st))]
  if (anyNA(st))
    stop("align: no start value for ", paste(names(st)[is.na(st)], collapse = ", "),
         ", whose symbols never reach the readings on their scale; give `init`.", call. = FALSE)
  if (!is.null(re)) {
    ctrl <- list(sd = 1, trustControl = c(list(iterlim = iterlim), tc), innerControl = tc["mtol"],
                 quadrature = list(cm1Control = list(tolControl = tc)),
                 sparseFit = list(maxOuter = iterlim, trustControl = list(tolControl = tc)))
    best <- if (fits > 1L) msem(obj, re, st, method = method, control = ctrl, fits = fits, cores = cores, sd = sd)
            else em(obj, re, st, method = method, control = ctrl, cores = 1L, verbose = FALSE)
    multistart <- if (!is.null(best$multistart)) best$multistart[c("value", "converged", "iterations")]
                  else data.frame(value = best$value, converged = isTRUE(best$converged),
                                  iterations = as.numeric(best$iterations %||% NA)[1])
    # the random effects the fit kept, a variance it set to 0 with the random effects of that variance fixed at 0
    re <- best$design$randomEffects
    pars <- c(best$argument, .etaVector(best), best$fixed)
    ctrl$boundary <- FALSE
  } else {
    fit <- if (fits > 1L)
      mstrust(obj, st, rinit = 1, rmax = 10, fits = fits, cores = cores, sd = sd, iterlim = iterlim,
              tolControl = tc, startFromCenter = TRUE, resultPath = tempdir(), name = "align")
    else list(trust(obj, st, rinit = 1, rmax = 10, iterlim = iterlim, tolControl = tc))
    multistart <- data.frame(value = vapply(fit, function(f) if (is.null(f$value)) NA_real_ else f$value, 0),
                             converged = vapply(fit, function(f) isTRUE(f$converged), TRUE),
                             iterations = vapply(fit, function(f) as.numeric(f$iterations %||% NA)[1], 0))
    if (all(is.na(multistart$value))) stop("align: every fit failed.", call. = FALSE)
    best <- fit[[which.min(multistart$value)]]
    pars <- best$argument
  }
  rounds <- 0L
  if (reml) {
    # rounds of the restricted likelihood from the best fit, the REML term linearized at the last round
    design <- if (!is.null(re)) .emDesign(obj, re) else .remlDesign(obj)
    theta <- setdiff(names(best$argument), mean)
    ml <- list(best = best, pars = pars)
    repeat {
      term <- .remlAnchor(design, pars, mean, theta)
      # where the information of the mean parameters is singular, the restricted likelihood is not
      # defined, as where the rounds run off to a spread of random effects without bound
      if (!is.finite(term(best$argument)$value)) {
        warning("align: the values and scales are not identifiable at ",
                if (rounds) sprintf("round %d of the restricted likelihood", rounds) else "the maximum likelihood fit",
                ", where the restricted likelihood is not defined; the fit is by maximum likelihood.", call. = FALSE)
        best <- ml$best; pars <- ml$pars
        reml <- FALSE; rounds <- 0L
        break
      }
      rounds <- rounds + 1L
      nxt <- if (!is.null(re)) em(obj + term, re, best$argument, fixed = best$fixed, method = method,
                                  control = ctrl, cores = 1L, verbose = FALSE)
             else trust(obj + term, best$argument, rinit = 1, rmax = 10, iterlim = iterlim, tolControl = tc)
      # the mean parameters are at their optimum for the error parameters of every round
      moved <- max(abs(nxt$argument[theta] - best$argument[theta]))
      nxt$boundary <- best$boundary
      best <- nxt
      pars <- if (!is.null(re)) c(best$argument, .etaVector(best), best$fixed) else best$argument
      if (moved < 1e-4) break
      if (rounds == 20L) {
        warning("align: the restricted likelihood moved the error parameters by ", signif(moved, 2),
                " in the last of 20 rounds.", call. = FALSE)
        break
      }
    }
  }
  H <- best$hessian
  if (is.null(H) && !is.null(best$information)) H <- 2 * best$information
  # without random terms the exact Hessian: the one of the fit leaves out the second derivatives
  # of prediction and sigma, and can be indefinite at a minimum with large residuals
  lowrank <- NULL
  if (is.null(re)) {
    exact <- obj(pars, deriv2 = TRUE, hessian = if (is.null(H) || .isSparse(H)) "sparse")
    H <- exact$hessian
    lowrank <- exact$lowrank
  } else lowrank <- best$lowrank
  allPars <- pars
  on <- intersect(names(best$argument), rownames(H))
  # where the reported Hessian of a fit with random terms is not positive definite, the
  # Gauss-Newton Schur complement over the random effects
  seSource <- "Hessian"
  if (!is.null(re) && !is.null(H) && !.positiveDefinite(H[on, on, drop = FALSE], .lowrankOver(lowrank, on))) {
    # a singular eta block, as at a fit that does not identify the values, leaves no standard errors
    gn <- tryCatch(.alignmentGaussNewton(m, obj, re, best), error = function(e) NULL)
    H <- gn$hessian
    lowrank <- gn$lowrank
    seSource <- "Gauss-Newton"
  }
  # with REML the standard errors of the mean parameters at the estimated error parameters
  # (Lindstrom and Bates 1990): the inverse of their block of the information
  variance <- if (!is.null(H)) .alignmentVariance(H, if (reml) intersect(on, mean) else on, lowrank)
  aligned <- .alignmentValues(m, allPars, variance)
  scaled <- .alignmentScaled(m, allPars)
  if (any(m$censored)) aligned$lloq <- .alignmentLimits(m, aligned, scaled)
  npar <- length(best$argument) - nrow(m$normalization$A)
  attr(allPars, "estimation") <- if (reml) "REML" else "ML"
  structure(list(aligned = aligned, scaled = scaled, parameters = allPars,
                 value = best$value, npar = npar, AIC = best$value + 2 * npar,
                 values = sort(multistart$value, na.last = TRUE), multistart = multistart, remlRounds = rounds,
                 seSource = seSource,
                 fit = best, objfn = obj,
                 omega = if (!is.null(re)) best$omega, random = if (!is.null(re)) best$random,
                 randomEffects = re, model = m), class = "alignment")
}

## The fixed parameters of the symbols of the prediction, of any role, integrated out by the
## restricted likelihood.
.alignmentMean <- function(m) {
  s <- intersect(m$symbols$symbol, getSymbols(c(m$model, m$absoluteModel)))
  unique(unlist(lapply(m$design[s], function(d) colnames(d$X)), use.names = FALSE))
}

## Whether a symmetric Hessian, dense or sparse plus the low-rank part `lowrank`, is positive definite.
.positiveDefinite <- function(H, lowrank = NULL) {
  if (.isSparse(H)) {
    if (!is.null(.lowrankCholesky(H, lowrank))) return(TRUE)
    # with H positive definite the factor decides by the inertia; else the sum is factorized as one
    if (is.null(lowrank) || !is.null(.lowrankCholesky(H))) return(FALSE)
    return(!is.null(.lowrankCholesky(.lowrankDense(H, lowrank))))
  }
  H <- as.matrix(H)
  !is.null(tryCatch(chol((H + t(H)) / 2), error = function(e) NULL))
}

## The Gauss-Newton Schur complement over the random effects at the modes of a fit: of the
## expected information of the readings and the normal random effects plus the Hessian of the
## normalization and the prior, positive semidefinite by construction; `hessian` and `lowrank`.
.alignmentGaussNewton <- function(m, obj, re, fit) {
  design <- .emDesign(obj, re)
  eta <- .etaVector(fit)
  on <- names(fit$argument)
  e <- names(eta)
  free <- c(on, e)
  pars <- c(setNames(as.numeric(fit$argument), on), eta, fit$fixed)
  extra <- .alignmentNormalization(m)
  if (!is.null(m$prior)) extra <- extra + m$prior
  ex <- .hessianParts(extra(pars = pars[free], fixed = fit$fixed, deriv = TRUE, hessian = "sparse"), free, on)
  H <- .alignmentFisher(m, pars, free) + .normalFisher(design, pars, free) + ex$hessian
  list(hessian = Matrix::forceSymmetric(H[on, on] - Matrix::crossprod(H[e, on], Matrix::solve(H[e, e], H[e, on]))),
       lowrank = ex$lowrank)
}

## The expected information of the readings in -2 log L over the parameters `free` at `pars`, from
## the first derivatives of prediction and sigma: of the value for a reading, of lying below the
## limit for a censored one.
.alignmentFisher <- function(m, pars, free) {
  cond <- .alignmentConditions(m)
  rowsOf <- split(seq_len(nrow(m$data)), factor(cond, unique(cond)))
  locals <- .alignmentLocals(m, cond, rowsOf)
  nu <- if (identical(m$distribution, "t")) m$df else Inf
  # location and scale information of a normal or Student-t reading
  iLoc <- if (is.finite(nu)) 2 * (nu + 1) / (nu + 3) else 2
  iScale <- if (is.finite(nu)) 4 * nu / (nu + 3) else 4
  logF <- function(x, lower) if (is.finite(nu)) stats::pt(x, nu, lower.tail = lower, log.p = TRUE)
                             else stats::pnorm(x, lower.tail = lower, log.p = TRUE)
  logf <- function(x) if (is.finite(nu)) stats::dt(x, nu, log = TRUE) else stats::dnorm(x, log = TRUE)
  parts <- lapply(names(rowsOf), function(cn) {
    r <- rowsOf[[cn]]
    k <- match(locals[[cn]]$pars, free)
    use <- which(!is.na(k))
    if (!length(use)) return(NULL)
    ev <- .alignmentEval(m, r, pars, local = locals[[cn]])
    s <- ev$y[, "sigma"]
    Jp <- matrix(ev$J[, "prediction", use], length(r)) / s
    Js <- matrix(ev$J[, "sigma", use], length(r)) / s
    cens <- m$censored[r]
    # a censored reading is a Bernoulli one with p = F(z), z = (lloq - prediction) / sigma
    z <- ifelse(cens, (m$lloq[r] - ev$y[, "prediction"]) / s, 0)
    wz <- ifelse(cens, 2 * exp(2 * logf(z) - logF(z, TRUE) - logF(z, FALSE)), 0)
    B <- crossprod(Jp * sqrt(iLoc * !cens)) + crossprod(Js * sqrt(iScale * !cens)) +
      crossprod((Jp + z * Js) * sqrt(wz))
    list(i = rep(k[use], length(use)), j = rep(k[use], each = length(use)), x = as.vector(B))
  })
  Matrix::sparseMatrix(i = unlist(lapply(parts, `[[`, "i")), j = unlist(lapply(parts, `[[`, "j")),
                       x = unlist(lapply(parts, `[[`, "x")), dims = rep(length(free), 2L),
                       dimnames = list(free, free))
}

## The expected information of the normal random effects in -2 log L over `free`: 2 Omega^-1 per
## group for the etas, G tr(Omega^-1 Omega_a Omega^-1 Omega_b) for the Cholesky parameters.
.normalFisher <- function(design, pars, free) {
  H <- .sparseOver(NULL, free)
  for (sd in design$groupings %||% list(design)) {
    lay <- sd$normal
    if (is.null(lay) || !lay$K) next
    E <- sd$etaNames[, lay$stems, drop = FALSE]
    L <- .choleskyFactor(lay, pars)
    cd <- .choleskyDerivatives(lay, L)
    M <- lapply(cd$dOinv, function(d) -d %*% tcrossprod(L))
    # tr(M_a M_b) of M_a = Omega^-1 Omega_a
    Fc <- crossprod(vapply(M, as.vector, numeric(lay$K^2)), vapply(M, function(x) as.vector(t(x)), numeric(lay$K^2)))
    nm <- c(as.vector(E), lay$cholPars)
    blk <- Matrix::bdiag(2 * Matrix::kronecker(cd$Oinv, Matrix::Diagonal(nrow(E))), nrow(E) * Fc)
    dimnames(blk) <- list(nm, nm)
    H <- H + .sparseOver(blk, free)
  }
  H
}

## The columns of the data the fixed effects of a symbol read.
.alignmentColumns <- function(m, s) m$design[[s]]$columns

## The values of the biological symbols per combination of the columns of their formulas, on the
## scale of the symbol, with standard errors by the delta method on `variance`, from
## .alignmentVariance() (NA without).
.alignmentValues <- function(m, allPars, variance) {
  bio <- m$symbols$symbol[m$symbols$role == "biological"]
  .bindFill(lapply(bio, function(s) {
    X <- m$design[[s]]$X
    cols <- .alignmentColumns(m, s)
    key <- if (length(cols)) do.call(paste, c(lapply(m$data[cols], as.character), sep = "\r")) else rep("1", nrow(X))
    first <- !duplicated(key)
    Xc <- X[first, , drop = FALSE]
    lin <- as.numeric(Xc %*% allPars[colnames(X)])
    se <- if (is.null(variance)) rep(NA_real_, nrow(Xc)) else sqrt(variance(Xc))
    tv <- .trafoValue(m, s, lin)
    out <- m$data[first, cols, drop = FALSE]
    rownames(out) <- NULL
    cbind(out, symbol = s, estimate = lin, se = se, value = tv$value, sigma = abs(tv$slope) * se,
          stringsAsFactors = FALSE)
  }))
}

## Data frames bound by rows over the union of their columns, NA where one has none.
.bindFill <- function(frames) {
  cols <- unique(unlist(lapply(frames, names)))
  do.call(rbind, lapply(frames, function(f) {
    f[setdiff(cols, names(f))] <- NA
    f[cols]
  }))
}

## The readings with prediction, sigma, residual and their values on the scale of the first
## biological symbol.
.alignmentScaled <- function(m, allPars) {
  ev <- .alignmentEval(m, seq_len(nrow(m$data)), allPars, deriv = FALSE)
  out <- m$data
  out$bloq <- m$censored
  out$prediction <- ev$y[, "prediction"]
  out$sigma <- ev$y[, "sigma"]
  out$residual <- ifelse(m$censored, NA_real_, (m$y - out$prediction) / out$sigma)
  inv <- .alignmentInvert(m, allPars, m$y)
  out$scaled <- inv$value
  out$scaledSigma <- out$sigma / abs(inv$slope)
  out
}

## Per combination of the first biological symbol the smallest scaled value of the censored
## readings.
.alignmentLimits <- function(m, aligned, scaled) {
  s <- m$symbols$symbol[m$symbols$role == "biological"][1]
  cols <- .alignmentColumns(m, s)
  key <- function(d) if (length(cols)) do.call(paste, c(lapply(d[cols], as.character), sep = "\r")) else rep("1", nrow(d))
  lim <- tapply(ifelse(scaled$bloq, scaled$scaled, NA), key(scaled),
                function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE))
  ifelse(aligned$symbol == s, as.numeric(lim[key(aligned)]), NA_real_)
}

## Readings mapped onto the first biological symbol: per reading the value of the symbol at which
## the prediction meets the target with the other symbols at the fit; NaN where none does. Also
## d prediction / d value.
.alignmentInvert <- function(m, allPars, target) {
  bio <- m$symbols$symbol[m$symbols$role == "biological"][1]
  sol <- .alignmentSolve(m, allPars, target, bio)
  tv <- .trafoValue(m, bio, sol$predictor)
  v <- tv$value; v[is.nan(sol$predictor)] <- NaN
  list(value = v, slope = sol$slope / tv$slope)
}

## Per reading the predictor of the symbol `s` at which the prediction meets the target, the other
## symbols at `allPars`, by Newton steps with the exact derivative; NaN where none does. Also
## d prediction / d predictor there.
.alignmentSolve <- function(m, allPars, target, s) {
  n <- nrow(m$data)
  vars <- cbind(.alignmentPredictors(m, seq_len(n), allPars), as.matrix(m$data[, m$covariates, drop = FALSE]))
  colnames(vars) <- c(m$symbols$symbol, m$covariates)
  tX <- array(0, c(n, ncol(vars), 1L), dimnames = list(NULL, colnames(vars), s))
  tX[, s, 1] <- 1
  ok <- is.finite(target)
  evaluate <- function() .alignmentRun(m, m$absolute, function(fn, k)
    fn$evaluate(vars = .rowsOf(vars, k), params = numeric(0), tangentX = .rowsOf(tX, k)))
  for (it in 1:50) {
    ev <- evaluate()
    dl <- ev$tangent[, "prediction", 1]
    step <- ifelse(ok & dl != 0, (ev$y[, "prediction"] - target) / dl, 0)
    ok <- ok & is.finite(step)
    step <- pmax(pmin(replace(step, !ok, 0), 2), -2)
    vars[, s] <- vars[, s] - step
    if (all(abs(step) < 1e-12)) break
  }
  ev <- evaluate()
  miss <- !ok | !(abs(ev$y[, "prediction"] - target) <= 1e-8 * pmax(1, abs(target)))
  v <- vars[, s]; v[miss] <- NaN
  list(predictor = v, slope = ev$tangent[, "prediction", 1])
}


## Start values: biological intercepts at the median reading (the median positive one where the scale cannot
## reach it), other intercepts where their symbol is 1, the rest 0; the coefficients of the first scaling
## symbol where the predictor per experiment is the median of those meeting the relative readings.
.alignmentStart <- function(m, free) {
  y <- m$y[is.finite(m$y)]
  st <- setNames(numeric(length(free)), free)
  for (s in intersect(m$symbols$symbol, free)) {
    if (m$symbols$role[m$symbols$symbol == s] != "biological") {
      st[s] <- .trafoInverse(m, s, 1)
      next
    }
    st[s] <- .trafoInverse(m, s, stats::median(y))
    if (is.na(st[s])) st[s] <- .trafoInverse(m, s, stats::median(y[y > 0]))
  }
  if (anyNA(st)) return(st)
  s <- m$symbols$symbol[m$symbols$role == "scaling"][1]
  X <- m$design[[s]]$X
  cols <- intersect(colnames(X), free)
  if (!length(cols)) return(st)
  v <- .alignmentSolve(m, st, ifelse(m$censored | m$absolute, NA_real_, m$y), s)$predictor
  cs <- .alignmentColumns(m, s)
  key <- if (length(cs)) do.call(paste, c(lapply(m$data[cs], as.character), sep = "\r")) else rep("1", nrow(m$data))
  target <- tapply(v, key, function(x) if (any(is.finite(x))) stats::median(x[is.finite(x)]) else NA_real_)
  rows <- match(names(target), key)
  ok <- is.finite(target)
  if (!any(ok)) return(st)
  Xe <- X[rows[ok], cols, drop = FALSE]
  use <- cols[Matrix::colSums(Xe != 0) > 0]
  coef <- tryCatch(as.numeric(Matrix::qr.coef(Matrix::qr(Xe[, use, drop = FALSE]), as.numeric(target[ok]))),
                   error = function(e) NULL)
  if (length(coef) == length(use) && all(is.finite(coef))) st[use] <- coef
  st
}


#' @export
print.alignment <- function(x, ...) {
  reml <- identical(attr(x$parameters, "estimation"), "REML")
  cat(sprintf("alignment: %d readings, %d biological values, %d parameters\n", nrow(x$scaled), nrow(x$aligned),
              x$npar))
  cat(sprintf("  %s: -2 log L %.2f, AIC %.2f\n",
              if (!reml) "maximum likelihood"
              else sprintf("restricted likelihood (REML, %d %s)", x$remlRounds, if (x$remlRounds == 1L) "round" else "rounds"),
              x$value, x$AIC))
  # with REML the multistart is by maximum likelihood, the best multistart fit the start of the rounds
  nf <- length(x$values)
  if (nf > 1L || reml)
    cat(sprintf("  %s%s: -2 log L %s\n",
                if (!reml) "fits"
                else if (nf == 1L) "maximum likelihood fit, the start of REML"
                else "maximum likelihood fits, REML from the best",
                if (nf > 10L) sprintf(" (best 10 of %d)", nf) else "",
                paste(sprintf("%.2f", utils::head(x$values, 10L)), collapse = " ")))
  nid <- sum(is.infinite(x$aligned$se))
  if (nid) cat(sprintf("  %d biological values not identifiable (se Inf)\n", nid))
  if (identical(x$seSource, "Gauss-Newton"))
    cat("  standard errors from the Gauss-Newton information, the Hessian of the marginal is",
        "not positive definite\n")
  if (all(is.na(x$aligned$se))) cat("  standard errors NA: the random effects give no Hessian of the marginal\n")
  invisible(x)
}

#' @export
as.data.frame.alignment <- function(x, ...) x$aligned


#' Aligned values as dMod2 data
#'
#' @param x An `alignment` from [align()].
#' @param time,name,condition Columns of the aligned values that become
#'   `time`, `name` and `condition`; by default `"time"`, `"symbol"` (the
#'   biological symbol) and the remaining columns pasted.
#' @param ... Ignored.
#' @return A [datalist] with `value` and `sigma` on the scale of the
#'   biological symbol.
#' @export
as.datalist.alignment <- function(x, time = "time", name = "symbol", condition = NULL, ...) {
  a <- x$aligned
  if (!time %in% names(a)) stop("as.datalist: no column ", time, " in the aligned values.", call. = FALSE)
  vars <- setdiff(.alignedColumns(a), c(time, name))
  if (is.null(condition)) condition <- vars
  cn <- if (length(condition)) do.call(paste, c(lapply(a[condition], as.character), sep = "_")) else "C1"
  as.datalist(data.frame(name = as.character(a[[name]]),
                         time = as.numeric(as.character(a[[time]])),
                         value = a$value, sigma = a$sigma, condition = cn, stringsAsFactors = FALSE))
}

.alignedColumns <- function(a) setdiff(names(a), c("symbol", "estimate", "se", "value", "sigma", "lloq"))

## The levels of character columns on one axis: in the order of their values where all are numbers,
## else sorted.
.axisLevels <- function(...) {
  u <- unique(unlist(lapply(list(...), as.character)))
  num <- suppressWarnings(as.numeric(u))
  if (anyNA(num)) .naturalLevels(u) else u[order(num)]
}


#' Diagnostic plots of an alignment
#'
#' @description `plotAligned()` shows the aligned values with their standard
#'   errors, `plotScaled()` the measured readings mapped onto the values next to
#'   them (a censored reading has no value and is left out), `plotResiduals()` the
#'   standardized residuals against the prediction.
#'
#' @param x An `alignment` from [align()].
#' @param ... Ignored.
#' @return A ggplot.
#' @importFrom rlang .data
#' @export
plotAligned <- function(x, ...) {
  a <- x$aligned
  cols <- .alignedColumns(a)
  xv <- if ("time" %in% cols) "time" else cols[1]
  grp <- setdiff(cols, xv)
  a$group <- if (length(grp)) do.call(paste, c(lapply(a[grp], as.character), sep = " ")) else "all"
  if (is.character(a[[xv]])) a[[xv]] <- factor(a[[xv]], levels = .axisLevels(a[[xv]]))
  ggplot2::ggplot(a, ggplot2::aes(x = .data[[xv]], y = .data$value, color = .data$group)) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$value - .data$sigma, ymax = .data$value + .data$sigma),
                           width = 0) +
    ggplot2::labs(x = xv, y = "aligned value", color = NULL) +
    scale_color_dMod2(guide = if (length(grp)) "legend" else "none") + theme_dMod2()
}

#' @rdname plotValues
#' @export
plotValues.alignment <- function(x, ..., tol = 1, showSteps = FALSE, pseudoLog = FALSE, relative = FALSE)
  .plotWaterfall(list(.valueFrame(x, ...)), tol, showSteps, pseudoLog, relative)

#' @rdname plotAligned
#' @export
plotScaled <- function(x, ...) {
  s <- x$scaled
  bio <- x$aligned
  cols <- .alignedColumns(bio)
  xv <- if ("time" %in% cols) "time" else cols[1]
  s <- s[is.finite(s$scaled) & !s$bloq, ]
  if (is.character(bio[[xv]])) {
    lev <- .axisLevels(bio[[xv]], s[[xv]])
    bio[[xv]] <- factor(bio[[xv]], levels = lev)
    s[[xv]] <- factor(s[[xv]], levels = lev)
  }
  ggplot2::ggplot(s, ggplot2::aes(x = .data[[xv]], y = .data$scaled)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_point(data = bio, ggplot2::aes(x = .data[[xv]], y = .data$value), color = dMod2_colors[2]) +
    ggplot2::labs(x = xv, y = "reading on the scale of the values") +
    theme_dMod2()
}

#' @rdname plotAligned
#' @param parframe An `alignment`, for `plotResiduals()`.
#' @export
plotResiduals.alignment <- function(parframe, ...) {
  s <- parframe$scaled[!parframe$scaled$bloq, ]
  ggplot2::ggplot(s, ggplot2::aes(x = .data$prediction, y = .data$residual)) +
    ggplot2::geom_point() + ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::labs(x = "prediction", y = "standardized residual") +
    theme_dMod2()
}


## The variances x' V x of the rows x of a matrix over the parameters, V the inverse of half the Hessian of
## -2 log L plus `lowrank`: sparse Cholesky and Woodbury for a sparse Hessian, else the eigendecomposition;
## Inf along the null space of a singular one.
.alignmentVariance <- function(H, free, lowrank = NULL) {
  if (is.null(H) || !length(free)) return(NULL)
  H <- H[free, free, drop = FALSE] / 2
  lr <- .lowrankOver(lowrank, free)
  if (!is.null(lr)) lr$d <- lr$d / 2
  if (.isSparse(H)) {
    f <- .lowrankCholesky(H, lr)
    if (!is.null(f))
      return(function(X) {
        ok <- intersect(colnames(X), free)
        T <- methods::as(methods::as(Matrix::Matrix(X[, ok, drop = FALSE], sparse = TRUE), "generalMatrix"),
                         "TsparseMatrix")
        E <- Matrix::sparseMatrix(i = match(ok, free)[T@j + 1L], j = T@i + 1L, x = T@x,
                                  dims = c(length(free), nrow(X)))
        v <- numeric(nrow(X))
        for (start in seq(1L, nrow(X), by = 512L)) {
          j <- start:min(nrow(X), start + 511L)
          v[j] <- .lowrankQuadratic(f, .lowrankZ(f, E[, j, drop = FALSE]))
        }
        pmax(v, 0)
      })
    if (length(free) > 5000L) return(NULL)
    H <- as.matrix(.lowrankDense(H, lr))
  }
  e <- eigen((H + t(H)) / 2, symmetric = TRUE)
  null <- e$values < 1e-10 * max(e$values)
  U <- e$vectors[, !null, drop = FALSE]
  bad <- if (any(null)) free[rowSums(e$vectors[, null, drop = FALSE]^2) > 1e-3] else character(0)
  function(X) {
    ok <- intersect(colnames(X), free)
    Y <- matrix(0, nrow(X), length(free))
    Y[, match(ok, free)] <- as.matrix(X[, ok, drop = FALSE])
    v <- rowSums(sweep(Y %*% U, 2, sqrt(e$values[!null]), "/")^2)
    v[Matrix::rowSums(X[, intersect(bad, colnames(X)), drop = FALSE] != 0) > 0] <- Inf
    v
  }
}


#' @rdname anova.em
#' @order 2
#' @export
anova.alignment <- function(object, ...) {
  fits <- c(list(object), list(...))
  # restricted likelihoods differ by the log determinant of their mean parameters
  reml <- vapply(fits, function(f) identical(attr(f$parameters, "estimation"), "REML"), TRUE)
  mean <- lapply(fits, function(f) sort(.alignmentMean(f$model)))
  if (any(reml) && (!all(reml) || length(unique(mean)) > 1L))
    stop("anova: restricted likelihoods compare only between fits with the same parameters of the prediction; ",
         "fit with reml = FALSE.", call. = FALSE)
  .anovaTable(fits)
}


#' Simulate readings from an alignment
#'
#' @description Draws the random effects from their fitted distributions,
#' predicts every reading, adds the noise of the error model (normal or
#' Student-t) and censors at the lower limits.
#'
#' @param object An `alignment` from [align()].
#' @param nsim Number of data sets. Default 1.
#' @param seed Seed of the random number generator. Default `NULL`, the
#'   current state.
#' @param ... Ignored.
#' @return A list of `nsim` data.frames, the data of the model with simulated
#'   `value`.
#' @export
simulate.alignment <- function(object, nsim = 1, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  m <- object$model
  pars <- object$parameters
  lapply(seq_len(nsim), function(i) {
    if (!is.null(object$randomEffects))
      pars[names(.etaVector(object$fit))] <- .drawEtas(object$fit$design, c(object$fit$argument, object$fit$fixed))[names(.etaVector(object$fit))]
    y <- .alignmentEval(m, seq_len(nrow(m$data)), pars, deriv = FALSE)$y
    noise <- if (m$distribution == "t") stats::rt(nrow(y), m$df) else stats::rnorm(nrow(y))
    d <- m$data
    d$value <- y[, "prediction"] + y[, "sigma"] * noise
    if (!is.null(d$lloq)) d$value <- ifelse(!is.na(d$lloq) & d$value <= d$lloq, NA_real_, d$value)
    d
  })
}
