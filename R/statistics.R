# Intermediate profile per parameter and direction, written under cautiousMode.
.profileDump <- function(out, name, pars, obj.attributes, folder, side) {
  frame <- as.data.frame(out)
  frame$whichPar <- name
  dput(parframe(frame,
                parameters     = names(pars),
                metanames      = c("value", "constraint", "stepsize", "gamma", "whichPar"),
                objAttributes  = names(obj.attributes)),
       file = file.path(folder, paste0(name, "-", side, ".R")))
}


#' Profile Likelihood
#'
#' Computes the profile likelihood of one or more parameters around an
#' optimum. `method = "integrate"` follows the profile path by Euler steps of
#' the path equation of the Lagrangian; `method = "optimize"` re-optimises the
#' other parameters at every step. A method of [stats::profile()] for dMod
#' objective functions and plain R functions.
#'
#' @param fitted Objective function \code{fitted(pars, fixed, ...)} returning a
#'   list with \code{value}, \code{gradient} and \code{hessian}. Numeric
#'   attributes of the result, such as \code{"data"} and \code{"prior"}, are
#'   kept as columns of the profile.
#' @param pars Named numeric parameter vector at the optimum.
#' @param whichPar Numeric or character vector, the parameters to profile.
#' @param alpha Significance level; the profile stops where the objective has
#'   risen by \code{qchisq(1 - alpha, 1)}. Default \code{0.05}.
#' @param delta Rise in the objective at which a profile stops. Default
#'   \code{NULL}, which uses \code{qchisq(1 - alpha, 1)}. Pass
#'   [profileThreshold()] with the arguments later given to
#'   [confint.parframe()], so the profile reaches that threshold.
#' @param limits Numeric of length 2, the lower and upper deviation from
#'   \code{pars[whichPar]} at which a profile stops. Default
#'   \code{c(lower = -Inf, upper = Inf)}.
#' @param method \code{"integrate"} (default) or \code{"optimize"}. Selects the
#'   defaults of \code{stepControl}, \code{algoControl} and \code{optControl}.
#' @param stepControl List of step-size settings, merged into the defaults.
#'   Defaults are given as \code{"integrate"} / \code{"optimize"}.
#'   \describe{
#'     \item{\code{stepsize}}{Initial step in the profiled parameter.
#'       \code{1e-4} / \code{1e-2}.}
#'     \item{\code{min}, \code{max}}{Bounds on the step size. \code{1e-4} and
#'       \code{Inf} for both methods.}
#'     \item{\code{atol}}{The step shrinks by a factor 1.5 when actual and
#'       predicted change of the objective differ by more than \code{atol}.
#'       \code{1e-2} / \code{1e-1}.}
#'     \item{\code{rtol}}{The step doubles when the difference is below
#'       \code{0.3 * atol}, or below \code{0.3 * rtol} relative to the actual
#'       change. \code{1e-2} / \code{1e-1}.}
#'     \item{\code{limit}}{Maximum number of steps per side. \code{500} /
#'       \code{100}.}
#'     \item{\code{stop}}{Component of the objective whose rise is compared
#'       with \code{delta}: \code{"value"} or the name of a numeric attribute
#'       such as \code{"data"}. \code{"value"} for both methods.}
#'   }
#' @param algoControl List of path settings, merged into the defaults.
#'   Defaults are given as \code{"integrate"} / \code{"optimize"}.
#'   \describe{
#'     \item{\code{gamma}}{Initial factor of the correction towards the
#'       optimum of the other parameters. \code{1} / \code{0}.}
#'     \item{\code{W}}{Matrix in the path equation, \code{"hessian"} or
#'       \code{"identity"}. \code{"hessian"} / \code{"identity"}.}
#'     \item{\code{reoptimize}}{Whether every step is followed by a fit of the
#'       other parameters with [trust()]. \code{FALSE} / \code{TRUE}.}
#'     \item{\code{correction}}{Relative correction above which \code{gamma}
#'       is halved; below half of it \code{gamma} doubles up to its initial
#'       value. \code{1} for both methods.}
#'     \item{\code{reg}}{Singular values below \code{reg} are dropped when the
#'       Hessian is inverted. \code{.Machine$double.eps} / \code{0}.}
#'   }
#' @param optControl List of arguments for [trust()] in the re-optimisation,
#'   merged into \code{list(rinit = 0.1, rmax = 10, iterlim = 10)} for
#'   \code{"integrate"} and \code{list(rinit = 0.1, rmax = 10, iterlim = 100)}
#'   for \code{"optimize"}. Any other argument of [trust()] can be added.
#' @param verbose Logical, print progress. Default \code{FALSE}.
#' @param cores Number of parallel workers over parameters, or
#'   `c(pars = , conditions = )` to also give each worker a condition axis,
#'   which selects a PSOCK backend. Multiplies with the threads each objective
#'   evaluation uses; keep the product below the number of cores. Default
#'   \code{1}.
#' @param cautiousMode Logical. If \code{TRUE}, every step is written to
#'   \code{profiles-interRes/<parameter>-<side>.R} in the working directory,
#'   and the files are kept. Default \code{FALSE}.
#' @param side \code{"both"} (default), \code{"left"} or \code{"right"}, the
#'   side of the optimum to profile.
#' @param ... Further arguments passed to \code{fitted}, for example
#'   \code{fixed}.
#'
#' @return A [parframe()] with one row per profile point of all profiled
#'   parameters, with metanames \code{value} (the objective), \code{constraint}
#'   (deviation of the profiled parameter from its value in \code{pars}),
#'   \code{stepsize}, \code{gamma} and \code{whichPar} (the profiled
#'   parameter), columns for the numeric attributes of the objective, and one
#'   column per parameter. Profiles that fail are dropped with a message.
#'
#' @seealso [confint.parframe()] for confidence intervals, [plotProfile()] and
#'   [plotPaths()] for plots, [profileThreshold()].
#' @example inst/examples/profiles.R
#' @importFrom stats profile
#' @rdname profile
#' @export
profile.objfn <- function(fitted, pars, whichPar, alpha = 0.05,
                          limits = c(lower = -Inf, upper = Inf),
                          method = c("integrate", "optimize"),
                          stepControl = NULL,
                          algoControl = NULL,
                          optControl  = NULL,
                          verbose = FALSE,
                          cores = 1,
                          cautiousMode = FALSE,
                          side = c("both", "left", "right"),
                          delta = NULL,
                          ...) {

  dotArgs <- list(...)
  objfun <- fitted

  sanePars <- sanitizePars(pars, dotArgs$fixed)
  pars <- sanePars$pars
  fixed <- sanePars$fixed

  method  <- match.arg(method)
  side <- match.arg(side)
  if (is.null(delta)) delta <- qchisq(1 - alpha, 1)

  sControl <- switch(method,
    integrate = list(stepsize = 1e-4, min = 1e-4, max = Inf, atol = 1e-2, rtol = 1e-2,
                     limit = 500, stop = "value"),
    optimize  = list(stepsize = 1e-2, min = 1e-4, max = Inf, atol = 1e-1, rtol = 1e-1,
                     limit = 100, stop = "value"))
  aControl <- switch(method,
    integrate = list(gamma = 1, W = "hessian", reoptimize = FALSE, correction = 1,
                     reg = .Machine$double.eps),
    optimize  = list(gamma = 0, W = "identity", reoptimize = TRUE, correction = 1, reg = 0))
  oControl <- switch(method,
    integrate = list(rinit = .1, rmax = 10, iterlim = 10),
    optimize  = list(rinit = .1, rmax = 10, iterlim = 100))
  sControl <- .mergeControl(sControl, stepControl, "stepControl")
  aControl <- .mergeControl(aControl, algoControl, "algoControl")
  aControl$W <- match.arg(aControl$W[1], c("hessian", "identity"))
  oControl <- .trustControl(oControl, optControl, label = "optControl")

  if (is.character(whichPar)) {
    lost <- setdiff(whichPar, names(pars))
    if (length(lost))
      stop("profile: ", paste(lost, collapse = ", "),
           if (all(lost %in% names(fixed))) " fixed" else " not among the parameters",
           call. = FALSE)
    whichPar <- match(whichPar, names(pars))
  }
  if (!length(whichPar) || anyNA(whichPar) || any(whichPar < 1 | whichPar > length(pars)))
    stop("profile: whichPar must name or index parameters of pars.", call. = FALSE)

  # cores = c(pars = , conditions = ) splits the two axes; see mstrust().
  .cc <- .splitCores(cores, "pars")
  coresConditions <- if (!is.null(.cc$conditions)) .sanitizeCores(.cc$conditions)
  cores <- .sanitizeCores(min(length(whichPar), .cc$outer))

  folder <- "profiles-interRes"
  if (cautiousMode) dir.create(folder, showWarnings = FALSE)

  out <- .parallelLapply(whichPar, function(whichIndex) tryCatch(
    .profileOne(whichIndex, objfun, pars, fixed, dotArgs, sControl, aControl, oControl,
                limits, delta, side, verbose, if (cautiousMode) folder),
    error = function(e) e), cores = cores, coresConditions = coresConditions)

  names(out) <- names(pars)[whichPar]
  failed <- vapply(out, inherits, TRUE, "error")
  reasons <- vapply(out[failed], conditionMessage, "")
  if (all(failed))
    stop("profile: ", paste(unique(reasons), collapse = "; "), call. = FALSE)
  if (any(failed))
    message("Profiles failed: ", paste0(names(reasons), " (", reasons, ")", collapse = ", "))
  do.call(rbind, out[!failed])
}


# The profile of parameter whichIndex. Each side follows the path equation of the
# Lagrangian by steps from the optimum, until the objective has risen by delta, the
# limit or the step limit is reached, or the objective fails.
.profileOne <- function(whichIndex, objfun, pars, fixed, dotArgs, sControl, aControl,
                        oControl, limits, delta, side, verbose, folder) {

  if (inherits(objfun, "fn")) loadDLL(objfun)
  name <- names(pars)[whichIndex]
  n <- length(pars)
  # A profile step starts next to an optimum, where single shooting is enough.
  obj.opt <- .singleShooting(objfun)
  nSingular <- 0L
  nFailedOpt <- 0L

  # With W = "identity" the steps go along whichPar only.
  Id <- NULL
  evalObj <- function(p) {
    out <- do.call(objfun, c(list(p), dotArgs))
    if (aControl$W == "identity") {
      if (is.null(Id)) {
        Id <<- diag(1/.Machine$double.eps, length(out$gradient))
        Id[whichIndex, whichIndex] <<- 1
        dimnames(Id) <<- list(names(out$gradient), names(out$gradient))
      }
      out$hessian <- Id
    }
    out
  }

  # Solutions of M x = rhs; a matrix too ill-conditioned for solve() takes the
  # pseudoinverse, whose singular values below reg are flagged invalid.
  solveM <- function(M, rhs) {
    if (rcond(M) > sqrt(.Machine$double.eps)) {
      sol <- tryCatch(solve(M, rhs), error = function(e) NULL)
      if (!is.null(sol)) return(list(sol = sol, valid = rep(TRUE, nrow(M))))
    }
    s <- svd(M)
    idx <- which(abs(s$d) > max(dim(M)) * max(s$d) * aControl$reg)
    sol <- if (length(idx))
      s$v[, idx, drop = FALSE] %*% (crossprod(s$u[, idx, drop = FALSE], rhs) / s$d[idx])
    else matrix(0, ncol(M), ncol(rhs))
    list(sol = sol, valid = seq_along(s$d) %in% idx)
  }

  lagrange <- function(y, direction, gamma, out = evalObj(y)) {
    gdot <- direction * (seq_len(n) == whichIndex)
    M <- rbind(cbind(out$hessian, gdot), c(gdot, 0))
    st <- solveM(M, cbind(c(-gamma * out$gradient, 1), c(rep(0, n), 1)))
    dy <- st$sol[seq_len(n), 1]
    dy0 <- st$sol[seq_len(n), 2]
    valid <- st$valid
    if (!all(valid)) {
      dy[!valid[seq_len(n)]] <- dy0[!valid[seq_len(n)]] <- 0
      dy[whichIndex] <- dy0[whichIndex] <- direction
      nSingular <<- nSingular + 1L
    }
    names(dy) <- names(y)
    # Numeric attributes become profile modes that add up to the total; chi2
    # decomposes one of them and is not a mode of its own.
    attrs <- Filter(is.numeric, attributes(out))
    attrs <- attrs[!grepl("^chi2($|_)", names(attrs))]
    c(list(dy = dy, value = out$value, gradient = out$gradient,
           correction = sqrt(sum((dy - dy0)^2)) / sqrt(sum(dy^2)),
           valid = valid, attributes = names(attrs), obj = out),
      attrs)
  }

  iterate <- function(y, dy, lag, direction, stepsize) {
    optimize <- aControl$reoptimize
    if (anyNA(dy)) {
      optimize <- TRUE
      y.try <- y
      y.try[whichIndex] <- y[whichIndex] + direction * stepsize
      rinit <- oControl$rinit
    } else {
      rinit <- min(oControl$rinit, 3 * sqrt(sum(dy^2)))
      y.try <- y + dy
      if (!all(lag$valid)) optimize <- TRUE
    }
    if (!optimize) return(y.try)
    arglist <- c(list(objfun = obj.opt, parinit = y.try[-whichIndex],
                      fixed = c(fixed, y.try[whichIndex]), rinit = rinit),
                 oControl[names(oControl) != "rinit"],
                 dotArgs[names(dotArgs) != "fixed"])
    fit <- try(do.call(trust, arglist), silent = TRUE)
    if (inherits(fit, "try-error")) nFailedOpt <<- nFailedOpt + 1L
    else y.try[names(fit$argument)] <- as.vector(fit$argument)
    y.try
  }

  # Accepts y.try if the objective changed as predicted, and adapts step and gamma.
  adapt <- function(y, y.try, lag, gamma, stepsize, direction) {
    lag.try <- lagrange(y.try, direction, gamma)
    valid <- TRUE
    dobj.pred <- sum(lag$gradient * (y.try - y))
    dobj.fact <- lag.try$value - lag$value
    if (lag.try$correction > aControl$correction) gamma <- gamma / 2
    if (lag.try$correction < 0.5 * aControl$correction) gamma <- min(aControl$gamma, gamma * 2)
    mismatch <- abs(dobj.fact - dobj.pred)
    if (mismatch > sControl$atol && stepsize > sControl$min) {
      stepsize <- max(stepsize / 1.5, sControl$min)
      valid <- FALSE
    }
    if (mismatch < .3 * sControl$atol || abs(mismatch / dobj.fact) < .3 * sControl$rtol)
      stepsize <- min(stepsize * 2, sControl$max)
    list(y = y.try, lag = lag.try, stepsize = stepsize, gamma = gamma, valid = valid)
  }

  profileRow <- function(lag, con, stepsize, gamma, y)
    c(value = lag$value, constraint = as.vector(con), stepsize = stepsize, gamma = gamma,
      whichPar = whichIndex, unlist(lag[lag$attributes]), y)

  origin <- lagrange(pars, 1, aControl$gamma)
  if (!is.numeric(origin[[sControl$stop]]))
    stop("stepControl$stop = '", sControl$stop, "' is no element of the objective's value.",
         call. = FALSE)
  threshold <- origin[[sControl$stop]] + delta
  originRow <- profileRow(origin, 0, sControl$stepsize, aControl$gamma, pars)

  dump <- function(rows, label) {
    if (is.null(folder)) return(invisible(NULL))
    frame <- do.call(rbind, if (label == "right") c(list(originRow), rows)
                            else c(rev(rows), list(originRow)))
    .profileDump(frame, name, pars, origin[origin$attributes], folder, label)
  }

  # Rows away from the optimum and the reason the side ended
  walk <- function(direction, lag) {
    label <- if (direction > 0) "right" else "left"
    if (verbose) cat("\nCompute", label, "profile\n")
    limit <- if (direction > 0) limits[2] else limits[1]
    gamma <- aControl$gamma
    stepsize <- sControl$stepsize
    y <- pars
    con <- 0
    rows <- list()
    reason <- "step limit"
    for (i in seq_len(sControl$limit)) {
      # The last step lands on the limit instead of overshooting.
      stepsize <- min(stepsize, direction * (limit - con))
      sufficient <- FALSE
      retry <- 0
      while (!sufficient && retry < 5) {
        step <- tryCatch(adapt(y, iterate(y, stepsize * lag$dy, lag, direction, stepsize), lag,
                               gamma, stepsize, direction),
                         error = function(e) e)
        if (inherits(step, "error")) {
          stepsize <- stepsize / 1.5
          retry <- retry + 1
        } else {
          sufficient <- step$valid
          stepsize <- step$stepsize
        }
      }
      if (inherits(step, "error")) {
        reason <- conditionMessage(step)
        break
      }
      y <- step$y
      lag <- step$lag
      con <- y[[whichIndex]] - pars[[whichIndex]]
      gamma <- step$gamma
      rows[[i]] <- profileRow(lag, con, stepsize, gamma, y)
      dump(rows, label)
      if (verbose) {
        progress <- max(c(if (is.finite(threshold)) 1 - max(0, min(1, (threshold - lag$value) / delta)),
                          i / sControl$limit,
                          if (is.finite(limit) && limit != 0) con / limit), na.rm = TRUE)
        .progressBar(100 * progress)
        cat("\tvalue:", formatC(lag$value, digits = 4, format = "g", width = 10),
            "constraint:", formatC(con, digits = 4, format = "g", width = 10),
            "gamma:", formatC(gamma, digits = 4, format = "g", width = 10),
            "valid:", all(lag$valid))
      }
      if (lag[[sControl$stop]] > threshold) { reason <- "threshold"; break }
      if (direction * (con - limit) >= -sqrt(.Machine$double.eps)) { reason <- "limit"; break }
    }
    list(rows = rows, reason = reason, label = label)
  }

  right <- if (side %in% c("right", "both")) walk(1, origin)
  left <- if (side %in% c("left", "both"))
    walk(-1, lagrange(pars, -1, aControl$gamma, out = origin$obj))
  sides <- list(left, right)
  if (verbose) cat("\n")

  for (sd in Filter(Negate(is.null), sides))
    if (!sd$reason %in% c("threshold", "limit", "step limit"))
      message("profile ", name, ", ", sd$label, " side: stopped by an error of the objective: ",
              sd$reason)
  if (nSingular)
    warning("profile ", name, ": the Hessian was singular at ", nSingular,
            " steps, which were reoptimized.", call. = FALSE)
  if (nFailedOpt)
    warning("profile ", name, ": ", nFailedOpt,
            " reoptimizations failed; the profile may be erroneous.", call. = FALSE)

  out <- as.data.frame(do.call(rbind, c(rev(sides[[1]]$rows), list(originRow), sides[[2]]$rows)))
  out$whichPar <- name
  parframe(out,
           parameters = names(pars),
           metanames = c("value", "constraint", "stepsize", "gamma", "whichPar"),
           objAttributes = names(unlist(origin[origin$attributes])))
}

#' @rdname profile
#' @export
profile.function <- function(fitted, ...) profile.objfn(fitted, ...)



# A progress bar on the current console line, percentage in 0 to 100
.progressBar <- function(percentage, size = 50) {
  filled <- round(size * min(max(percentage, 0), 100) / 100)
  cat("\r|", strrep("=", filled), strrep(" ", size - filled), "| ",
      format(paste0(round(percentage), "%"), width = 4), sep = "")
  invisible(NULL)
}

#' Threshold for a Profile Likelihood
#'
#' @description
#' How far a profile has to rise for a parameter value to leave the confidence
#' region. `"chisq"` is the asymptotic \eqn{\chi^2_1} quantile. `"F"` is the
#' finite-sample threshold \eqn{n\log(1 + F_{1,\nu,1-\alpha}/\nu)},
#' \eqn{\nu = n - p}, exact for Gaussian errors with a single sigma profiled
#' out and the generalisation of the classical t interval.
#'
#' @param level the confidence level required.
#' @param method `"chisq"` (default) or `"F"`.
#' @param n number of data points. Only for `method = "F"`.
#' @param p number of estimated mean parameters, that is without the
#' error-model parameters. Only for `method = "F"`. [remlLeverage] reports the
#' effective count, which is what to pass when directions are non-identifiable.
#'
#' @return A single number.
#' @seealso [profile()][profile.objfn], [confint.parframe()]
#' @importFrom stats qf
#' @export
profileThreshold <- function(level = 0.95, method = c("chisq", "F"),
                             n = NULL, p = NULL) {

  method <- match.arg(method)
  if (method == "chisq") return(qchisq(level, df = 1))

  if (is.null(n) || is.null(p))
    stop("profileThreshold: method = \"F\" needs n and p.", call. = FALSE)
  nu <- n - p
  if (nu <= 0)
    stop("profileThreshold: n - p must be positive, got ", nu, ".", call. = FALSE)
  n * log(1 + qf(level, 1, nu) / nu)
}


#' Confidence Intervals from Profiles
#'
#' @description Reads confidence intervals off profile likelihoods. Where a
#'   profile does not reach the threshold, the last part of it is extrapolated
#'   linearly.
#' @param object A `parframe` returned by [profile()][profile.objfn].
#' @param parm Names of the parameters to give intervals for. Default `NULL`,
#'   all profiled parameters.
#' @param level The confidence level. Default `0.95`.
#' @param ... Not used.
#' @param val.column Column of the profile the threshold applies to, for
#'   example `"value"` or `"data"`. Default `"data"`, or `"value"` for a profile
#'   without a `data` column.
#' @param method Calibration of the threshold, `"chisq"` (default) or `"F"`,
#'   see [profileThreshold()].
#' @inheritParams profileThreshold
#' @return A data frame with one row per parameter and columns `name`,
#'   `value` (the value at the optimum), `lower` and `upper`. A bound is `NA`
#'   or infinite when the profile on that side is too short or does not rise
#'   towards the threshold.
#' @seealso [profile()][profile.objfn], [profileThreshold()]
#' @importFrom stats qf
#' @export
confint.parframe <- function(object, parm = NULL, level = 0.95, ...,
                             val.column = "data",
                             method = c("chisq", "F"), n = NULL, p = NULL) {

  method <- match.arg(method)
  profile <- object
  if (missing(val.column) && !val.column %in% names(profile)) val.column <- "value"
  if (!val.column %in% names(profile))
    stop("confint: the profile has no column '", val.column, "'.", call. = FALSE)

  if (is.null(parm))
    parm <- unique(profile[["whichPar"]])

  threshold <- profileThreshold(level, method, n, p)
  
  # Reduce to profiles for parm
  profile <- profile[profile[["whichPar"]] %in% parm,]
  
  # Evaluate confidence intervals per parameter
  CIs <- lapply(split(profile, profile[["whichPar"]]), function(d) {
    
    # Get origin of profile
    origin <- which.min(abs(d[["constraint"]]))
    whichPar <- d[["whichPar"]][origin]
    
    # Define function to return constraint value where threshold is passed
    get_xThreshold <- function(branch) {
      
      y <- branch[[val.column]] - d[[val.column]][origin]
      x <- branch[["constraint"]]
      
      # If less than 3 points, return NA
      if (length(x) < 3)
        return(NA)
      
      # If threshold exceeded, take closest points below and above threshold
      # and interpolate
      if (any(y > threshold)) {
        i.above <- utils::head(which(y > threshold), 1)
        i.below <- utils::tail(which(y < threshold), 1)
        if (i.below > i.above) {
          return(NA)
        } else {
          slope <- (y[i.above] - y[i.below])/(x[i.above] - x[i.below])
          dy <- threshold - y[i.below]
          dx <- dy/slope
          x_threshold <- x[i.below] + dx
          return(x_threshold)
        }
      }
      
      # If threshold not exceeded,
      # take the last 20% of points (at least 3) an perform linear fit
      n_last20 <- max(3, length(which(x - x[1] > 0.8 * (max(x) - x[1]))))
      x <- tail(x, n_last20)
      y <- tail(y, n_last20)
      slope <- sum((x - mean(x))*(y - mean(y)))/sum((x - mean(x))^2)
      
      # If slope < 0, return Inf
      if (slope < 0)
        return(Inf)
      
      # Extrapolate until threshold is passed
      dy <- threshold - tail(y, 1)
      dx <- dy/slope
      x_threshold <- tail(x, 1) + dx
      
      # Test if extrapolation takes the point of passage very far
      # Set to Inf in that case
      if (x_threshold > 10*(max(x) - min(x)))
        x_threshold <- Inf
      
      return(x_threshold)
      
      
    } 
    
    # Right profiles
    right <- d[d[["constraint"]] >= 0,]
    upper <- d[[whichPar]][origin] + get_xThreshold(right)
    
    # Left profile
    left <- d[d[["constraint"]] <= 0,]
    left[["constraint"]] <- - left[["constraint"]]
    left <- left[order(left[["constraint"]]), ]
    lower <- d[[whichPar]][origin] - get_xThreshold(left)
    
    
    data.frame(name = whichPar, 
               value = d[[whichPar]][origin], 
               lower = lower, 
               upper = upper)
    
    
  })
  
  do.call(rbind, CIs)
  
}


#' Variance-Covariance Matrix of a Fit
#'
#' Methods of [stats::vcov()] for an evaluation of an objective function and
#' for the result of [trust()].
#'
#' @param object An evaluation of an objective function, class `objlist`, such
#'   as `obj(pars)`, or a result of [trust()], class `trustfit`.
#' @param fixed Names of parameters held fixed, for example at a bound.
#'   Default `NULL`.
#' @param parupper,parlower Named upper and lower bounds. Parameters at or
#'   beyond a bound count as fixed, as do those that `object$atBound` marks.
#'   Default `NULL`.
#' @param ... Not used.
#'
#' @return Square matrix with the parameter names as dimnames: the inverse of
#'   half the Hessian of `object` (pseudo-inverse when singular), with rows
#'   and columns of fixed parameters set to zero. For a `trustfit`, `NULL`
#'   when it has no `argument`, and an error when it has no Hessian; the fits
#'   of [mstrust()] have none, use `vcov(obj(fit$argument))` for them.
#' @seealso [trust()], [mstrust()]
#' @examples
#' obj <- function(x, ...) {
#'   list(value = sum(x^2 / c(1, 4)),
#'        gradient = 2 * x / c(1, 4),
#'        hessian = diag(2 / c(1, 4)))
#' }
#' fit <- trust(obj, c(a = 1, b = 1))
#' vcov(fit)
#'
#' # The same from an evaluation of the objective
#' ev <- obj(fit$argument)
#' vcov(objlist(ev$value, ev$gradient, ev$hessian))
#' @importFrom stats vcov
#' @export
vcov.objlist <- function(object, fixed = NULL, ...) {

  H <- object[["hessian"]]
  if (is.null(H)) stop("The objective was evaluated without a Hessian.")
  if (is.null(colnames(H))) {
    nm <- names(object[["gradient"]])
    dimnames(H) <- list(nm, nm)
  }

  V <- 0*H
  isFixed <- colnames(H) %in% fixed
  Hfree <- H[!isFixed, !isFixed, drop = FALSE]
  Vfree <- try(solve(0.5*Hfree), silent = TRUE)
  if (inherits(Vfree, "try-error")) Vfree <- MASS::ginv(0.5*Hfree)
  V[!isFixed, !isFixed] <- Vfree

  V

}

#' @rdname vcov.objlist
#' @export
vcov.trustfit <- function(object, parupper = NULL, parlower = NULL, ...) {

  argument <- object[["argument"]]
  if (is.null(argument)) return()
  if (is.null(object[["hessian"]]))
    stop("The fit has no Hessian. Evaluate the objective at the fit: ",
         "vcov(obj(fit$argument)).")

  # Parameters held by a bound rather than by the data. A reflective boundary
  # keeps iterates strictly inside and reports `atBound` itself; otherwise the
  # bounds are compared with a relative tolerance.
  fixed <- NULL
  atBound <- object[["atBound"]]
  if (!is.null(atBound) && !is.null(names(atBound)))
    fixed <- union(fixed, names(atBound)[which(as.logical(atBound))])
  reltol <- 1e-6
  if (!is.null(parupper)) {
    value <- argument[names(parupper)]
    fixed <- union(fixed, names(value)[
      which(value >= parupper - reltol * pmax(1, abs(parupper)))])
  }
  if (!is.null(parlower)) {
    value <- argument[names(parlower)]
    fixed <- union(fixed, names(value)[
      which(value <= parlower + reltol * pmax(1, abs(parlower)))])
  }

  vcov.objlist(objlist(object[["value"]], object[["gradient"]], object[["hessian"]]),
               fixed = fixed)

}




#' Non-Linear Optimisation, Multi-Start
#'
#' @description Runs [trust()] from several starting points: `center` plus
#'   random offsets drawn by `samplefun`, or the rows of a parframe.
#'
#' @param objfun Objective function, see [trust()].
#' @param center Named numeric parameter centre, or a parframe whose rows are
#'   used as starting points; `fits` is then the number of rows. The random
#'   starts, those of retries included, are drawn before the fits are spread
#'   over `cores`, so [set.seed()] reproduces a run.
#' @param rinit,rmax Initial and maximum trust-region radius, see [trust()].
#'   Default `0.1` and `10`.
#' @param fits Number of fits. Default `20`.
#' @param cores Number of parallel fits, or `c(fits = , conditions = )` to
#'   also give each fit a condition axis, which selects a PSOCK backend. Keep
#'   the product below the number of cores. Default `1`.
#' @param optmethod Optimiser, a function or its name, called via
#'   `do.call(optmethod, ...)`. Default `"trust"`. Arguments in `...` are
#'   routed by the formals of this optimiser.
#' @param samplefun Sampler for the random offsets, a function or its name,
#'   with an argument `n`, the number of parameters. Default `"rnorm"`.
#' @param resultPath Folder under which files are written. Default `"."`, the
#'   working directory.
#' @param name Folder name under `resultPath`. Default `"mstrust"`.
#' @param stats Logical, print the summary to the console. Default `FALSE`.
#' @param output Logical, write every fit and the final parlist to disk.
#'   Default `FALSE`.
#' @param cautiousMode Logical, also write every fit deparsed to an `.R` file
#'   in `interRes/`, which is kept after the run while its `.Rda` files are
#'   removed. Default `FALSE`.
#' @param startFromCenter Logical. If `TRUE`, the first fit starts at
#'   `center` without offset. Default `FALSE`.
#' @param retry Logical. If `TRUE` (default), a fit that fails with an error
#'   is repeated from a new random start, up to `nTries` attempts. Not used
#'   when `center` is a parframe.
#' @param nTries Maximum number of attempts per fit, including the first.
#'   Default `10`.
#' @param ... Arguments passed by name to `optmethod`, `samplefun` or
#'   `objfun`. Names matching neither the optimiser nor the sampler go to
#'   `objfun`.
#'
#' @details Nothing is written to disk unless `output`, `cautiousMode` or a
#'   `traceFile` argument for [trust()] asks for it. Files then go to
#'   `<resultPath>/<name>/trial-<k>-<date>/`: the per-fit files in its
#'   subfolder `interRes/`, which is removed on completion unless
#'   `cautiousMode = TRUE`, the final `parameterList.Rda` (with `output`) and
#'   the log `mstrust.log`. [load.parlist()] reads the fits from `interRes/`
#'   of a run that is still going or was aborted.
#'
#' @return A parlist with one entry per fit, failed fits included. Each entry
#'   is the result of the optimiser without `hessian`, plus `parinit`, its
#'   starting point, and `elapsed`, its wall time in seconds including
#'   retries; a fit that failed is a list with element `error`. For the
#'   covariance at a fit, use `vcov(objfun(fit$argument))`. The attribute
#'   `"run"` holds the start time, total wall time, number of fits, cores,
#'   optimiser and sampler with their arguments, and the versions of R, dMod2
#'   and cppDE. [as.parframe()] turns the parlist into a table.
#'
#' @seealso [trust()], [msParframe()], [as.parframe()], [load.parlist()]
#'
#' @examples
#' rosenbrock <- function(x, ...) {
#'   a <- x[["a"]]
#'   b <- x[["b"]]
#'   list(value = 100 * (b - a^2)^2 + (1 - a)^2,
#'        gradient = c(a = -400 * a * (b - a^2) - 2 * (1 - a),
#'                     b = 200 * (b - a^2)),
#'        hessian = matrix(c(1200 * a^2 - 400 * b + 2, -400 * a, -400 * a, 200),
#'                         2, 2, dimnames = list(c("a", "b"), c("a", "b"))))
#' }
#'
#' # Random starts around center, sd goes to rnorm()
#' fits <- mstrust(rosenbrock, center = c(a = 0, b = 0), fits = 5, sd = 2,
#'                 iterlim = 200)
#' as.parframe(fits)
#'
#' # Reproducible starts from a parframe
#' starts <- msParframe(c(a = 0, b = 0), n = 5, sd = 2)
#' fits <- mstrust(rosenbrock, center = starts, iterlim = 200)
#' as.parframe(fits)
#'
#' @export
mstrust <- function(objfun, center, rinit = .1, rmax = 10, fits = 20, cores = 1, optmethod = "trust",
                    samplefun = "rnorm", resultPath = ".", name = "mstrust",
                    stats = FALSE, output = FALSE, cautiousMode = FALSE,
                    startFromCenter = FALSE, retry = TRUE, nTries = 10L,
                    ...) {

  samplefunName <- if (is.character(samplefun)) samplefun else
    paste(deparse(substitute(samplefun)), collapse = "")
  samplefun <- match.fun(samplefun)
  optmethodName <- if (is.character(optmethod)) optmethod else
    paste(deparse(substitute(optmethod)), collapse = "")
  optfun <- if (is.function(optmethod)) optmethod else match.fun(optmethod)

  given <- c(mget(c("objfun", "rinit", "rmax", .mstrustOwn), envir = environment()), list(...))
  given$n <- length(center)
  route <- .mstrustRoute(given, optfun, samplefun, samplefunName)
  traceFile <- given[["traceFile"]]

  if (is.parframe(center)) fits <- nrow(center)
  if (length(fits) != 1L || is.na(fits) || fits < 1)
    stop("mstrust: fits must be a positive number.", call. = FALSE)

  # Nothing is written unless asked for: `output` and `cautiousMode` keep fits,
  # `traceFile` puts the per-fit solver traces into the result folder.
  writeFiles <- isTRUE(output) || isTRUE(cautiousMode) || !is.null(traceFile)

  started <- Sys.time()
  paths <- .mstrustPaths(resultPath, name, started)

  # The log is written by this process only, workers cannot share its connection.
  writeLog <- function(msg) if (writeFiles) cat(msg, file = paths$log, append = TRUE, sep = "\n")
  if (writeFiles) dir.create(path = paths$fits, showWarnings = FALSE, recursive = TRUE)
  assigned <- .mstrustColon(c(mstrust = toString(.mstrustOwn),
                              setNames(toString(names(route$optimiser)), optmethodName),
                              setNames(toString(names(route$sampler)), samplefunName)))
  writeLog(paste0("Parameter assignment information\n", paste(assigned, collapse = "\n"), "\n"))

  # cores = c(fits = , conditions = ) splits the two axes. cppDE's batch runs
  # serially inside a fork, so an inner axis > 1 selects PSOCK.
  .cc <- .splitCores(cores, "fits")
  coresConditions <- if (!is.null(.cc$conditions)) .sanitizeCores(.cc$conditions)
  cores <- .sanitizeCores(min(fits, .cc$outer))
  if (!is.null(coresConditions)) {
    oldCores <- options(dMod2.cores = coresConditions)
    on.exit(options(oldCores), add = TRUE)
  }

  # Every start, those of retries included, is drawn here before the fits are
  # spread over workers, so set.seed() reproduces a run whatever `cores`.
  maxTries <- if (isTRUE(retry) && !is.parframe(center)) as.integer(nTries) else 1L
  draw <- function() center + do.call(samplefun, route$sampler)
  starts <- lapply(seq_len(fits), function(i) {
    if (is.parframe(center)) as.parvec(center, i)
    else if (i == 1 && startFromCenter) center
    else draw()
  })
  restarts <- lapply(seq_len(fits), function(i) lapply(seq_len(maxTries - 1L), function(k) draw()))

  fitOne <- function(i) {
    # A plain R function has no shared objects to load.
    if (inherits(objfun, "fn")) suppressMessages(loadDLL(objfun))
    optArgs <- route$optimiser
    if (!is.null(traceFile)) {
      optArgs$traceFile <- file.path(paths$run, paste0(
        formatC(i, digits = floor(log10(fits)), flag = "0"), "_", traceFile))
    }

    # Each fit starts without the warm-start roots of the previous fit
    # (cores = 1) or of the parent process (cores > 1, forked).
    try(resetWarmStarts(objfun, verbose = FALSE), silent = TRUE)

    # Up to maxTries attempts, each from the next drawn start after a try-error
    # or fit$error, with warm starts reset. Rows of a parframe are never retried.
    t0 <- proc.time()[["elapsed"]]
    start <- starts[[i]]
    for (attempt in seq_len(maxTries)) {
      optArgs$parinit <- start
      fit <- try(do.call(optmethod, c(optArgs, route$objective)), silent = TRUE)
      failed <- inherits(fit, "try-error") || !is.list(fit) || "error" %in% names(fit)
      if (!failed || attempt == maxTries) break
      start <- restarts[[i]][[attempt]]
      try(resetWarmStarts(objfun, verbose = FALSE), silent = TRUE)
    }

    # Only the numeric attributes of the optimiser's result are kept.
    attr.fit <- attributes(fit)
    keep.attr <- vapply(attr.fit, is.numeric, TRUE)
    fit <- fit[seq_along(fit)]
    if (any(keep.attr)) attributes(fit) <- c(attributes(fit), attr.fit[keep.attr])
    if ("trustfit" %in% attr.fit$class) fit <- .trustFit(fit)
    # A crash can return a try-error, which is no list.
    if (!is.list(fit)) fit <- list(error = fit)

    fit[["parinit"]] <- start
    fit$hessian <- NULL
    fit$elapsed <- proc.time()[["elapsed"]] - t0

    if (output) {
      saveRDS(fit, file = file.path(paths$fits, paste0("fit-", i, ".Rda")))
      if (cautiousMode) dput(fit[c("value", "argument", "iterations", "converged", "elapsed")],
                           file = file.path(paths$fits, paste0("fit-", i, ".R")),
                           control = c("keepNA", "keepInteger", "niceNames", "showAttributes", "digits17"))
    }
    fit
  }

  parlist <- as.parlist(.parallelLapply(seq_len(fits), fitOne, cores = cores,
                                        coresConditions = coresConditions))

  if (output) for (i in seq_along(parlist)) writeLog(.mstrustFitEntry(i, parlist[[i]]))
  # objfun is a closure and parinit varies per fit, so neither is kept.
  attr(parlist, "run") <- list(
    started = started,
    elapsed = as.numeric(difftime(Sys.time(), started, units = "secs")),
    fits = fits,
    cores = c(fits = cores, conditions = coresConditions),
    optimiser = optmethodName,
    optimiserArgs = route$optimiser[setdiff(names(route$optimiser), c("objfun", "parinit"))],
    sampler = samplefunName,
    samplerArgs = route$sampler,
    versions = c(R = as.character(getRversion()),
                 dMod2 = as.character(utils::packageVersion("dMod2")),
                 cppDE = as.character(utils::packageVersion("cppDE"))))

  if (output) saveRDS(parlist, file = paths$parlist)
  if (writeFiles) {
    if (!cautiousMode) unlink(paths$fits, recursive = TRUE)
    else unlink(list.files(paths$fits, "\\.Rda$", full.names = TRUE))
  }

  summaryText <- .mstrustSummary(parlist, fits)
  writeLog(summaryText)
  if (stats) cat(summaryText)

  parlist
}

# Arguments mstrust() keeps for itself, never passed on.
.mstrustOwn <- c("name", "center", "fits", "cores", "optmethod", "samplefun",
                 "resultPath", "stats", "output", "cautiousMode", "retry", "nTries")

# Splits the named arguments of mstrust() by destination: the formals of the
# optimiser, then those of the sampler, which must be disjoint, and every
# remaining name that is not mstrust()'s own goes to the objective.
.mstrustRoute <- function(given, optfun, samplefun, samplerName) {
  accepts <- function(f) intersect(setdiff(names(formals(f)), "..."), names(given))
  toOptimiser <- accepts(optfun)
  .trustRejectMoved(setdiff(names(given), toOptimiser), "mstrust")
  toSampler <- accepts(samplefun)
  if (length(intersect(toOptimiser, toSampler)))
    stop("Argument names of the optimiser and ", samplerName, "() clash.")
  toObjective <- setdiff(names(given), c(toOptimiser, toSampler, .mstrustOwn))
  list(optimiser = given[toOptimiser], sampler = given[toSampler],
       objective = given[toObjective])
}

# Folder of one run, `<resultPath>/<name>/trial-<k>-<date>`, numbered after the
# runs already there, and the files inside it.
.mstrustPaths <- function(resultPath, name, started) {
  base <- file.path(resultPath, name)
  k <- length(list.files(base, pattern = "trial")) + 1L
  run <- file.path(base, sprintf("trial-%d-%s", k, format(started, "%d-%m-%Y-%H%M%S")))
  list(run = run, fits = file.path(run, "interRes"), log = file.path(run, "mstrust.log"),
       parlist = file.path(run, "parameterList.Rda"))
}

# "label       : value" lines, labels padded to 12 characters.
.mstrustColon <- function(x) paste0(formatC(names(x), width = -12), ": ", x)

# Log entry of fit i: its outcome between two rules.
.mstrustFitEntry <- function(i, fit) {
  body <- if (!is.null(fit$error)) {
    c(sprintf("Fit %d failed with error", i),
      paste("-->", trimws(paste(as.character(fit$error), collapse = " "))))
  } else {
    shown <- list(iterations = fit$iterations, converged = fit$converged,
                  `obj. value` = round(fit$value, digits = 2))
    c(sprintf("Fit %d completed", i),
      sprintf("--> %10s : %s", names(shown), vapply(shown, paste, "", collapse = " ")))
  }
  paste(c("-------", body, "-------"), collapse = "\n")
}

# Count of fits per outcome and, where the optimiser reports one, per
# termination reason.
.mstrustSummary <- function(parlist, fits) {
  status <- .statParlist(parlist)
  counts <- c(Outcome = "Occurrence",
              Error = sum(status == "error"),
              Unconverged = sum(status == "notconverged"),
              Converged = sum(status == "converged"))
  out <- c("Multi start summary", .mstrustColon(counts),
           paste0(strrep(" ", 11), strrep("-", 11)), .mstrustColon(c(Total = fits)))
  reasons <- .stopReasonTable(parlist)
  if (!is.null(reasons))
    out <- c(out, "", "Termination reason",
             .mstrustColon(setNames(as.integer(reasons), names(reasons))))
  paste0(paste(out, collapse = "\n"), "\n")
}


#' Reproducible Random Starting Points
#'
#' Draws a parframe of starting points for the `center` argument of
#' [mstrust()] from the seed `seed`. The state of the global random number
#' generator is left unchanged.
#'
#' @param pars Named numeric vector. If `samplefun` has an argument `mean`,
#'   the draws are centred on `pars`.
#' @param n Integer, number of rows. Default `20`.
#' @param seed Seed for the random number generator. Default `12345`.
#' @param samplefun Random number generator such as [rnorm()] or [runif()].
#'   Default [stats::rnorm()].
#' @param keepfirst Logical. If `TRUE` (default), the first row is `pars`.
#' @param ... Arguments passed to `samplefun`.
#'
#' @return A parframe without metanames, one column per parameter.
#' @export
#'
#' @seealso [mstrust()], [parframe()]
#'
#' @examples
#' # Five starts around pars, the first one at pars
#' msParframe(c(k1 = 1, k2 = 10), n = 5)
#'
#' # A standard deviation per parameter
#' msParframe(c(k1 = 1, k2 = 10), n = 5, sd = c(0.1, 2))
msParframe <- function(pars, n = 20, seed = 12345, samplefun = stats::rnorm,
                       keepfirst = TRUE, ...) {

  oldSeed <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (is.null(oldSeed)) rm(".Random.seed", envir = globalenv())
          else assign(".Random.seed", oldSeed, envir = globalenv()))
  set.seed(seed)

  if (keepfirst && n == 1) return(parframe(as.data.frame(t(pars))))

  # `keepfirst` spends one row on `pars` itself.
  ndraw <- if (keepfirst) n - 1L else n
  draws <- matrix(samplefun(ndraw * length(pars), ...), nrow = ndraw, byrow = TRUE)

  # A sampler with a `mean` argument is centred on `pars`.
  if ("mean" %in% names(formals(samplefun)))
    draws <- draws + t(matrix(pars, nrow = length(pars), ncol = ndraw))

  if (keepfirst) draws <- rbind(t(pars), draws)

  parframe(`names<-`(as.data.frame(draws), names(pars)))
}


#' Load Fits Written by mstrust
#'
#' @description Reads the per-fit files that [mstrust()] writes with
#'   `output = TRUE`, for example after an aborted run.
#'
#' @param folder The folder `<resultPath>/<name>/trial-<k>-<date>/interRes`
#'   of the run, which holds one `fit-<i>.Rda` per completed fit.
#'
#' @return A parlist.
#'
#' @seealso [mstrust()]
#'
#' @export
load.parlist <- function(folder) {
  files <- dir(folder, pattern = "\\.Rda$")
  number <- suppressWarnings(as.integer(gsub("\\D", "", files)))
  files <- files[order(number, files)]
  as.parlist(lapply(file.path(folder, files), readRDS))
}


#' Reduce Replicates to Mean and Standard Error
#'
#' @description
#' Replaces the replicates of each condition by their mean and the standard
#' error of the mean.
#'
#' @param data A data frame with the columns `name` (observable), `time`,
#'   `value` and `condition`, plus any further columns that define conditions.
#' @param select Names of the columns whose values, together with `name`,
#'   `time` and `condition`, define a condition. Default `"condition"`.
#' @param datatrans Character, an expression in `x` applied to `value` before
#'   the reduction, for example `"log(x)"`. Default `NULL`.
#' @param keep Names of columns that are kept although their values differ
#'   within a condition; the first value is used. Default `NULL`.
#' @param weighted Logical. If `TRUE`, mean and standard error are weighted by
#'   `1 / sigma^2` from the column `sigma` of `data`. Default `FALSE`.
#'
#' @details
#' A column not in `select` is kept when its value is the same for all
#' replicates of every condition, otherwise it is dropped with a message.
#'
#' @return
#' A data frame with the columns
#' \describe{
#'  \item{time}{Time of the measurement.}
#'  \item{value}{Mean of the replicates.}
#'  \item{sigma}{Standard error of the mean, `NA` for a single measurement.}
#'  \item{n}{Number of replicates.}
#'  \item{name}{Observable.}
#'  \item{condition}{The values of the `select` columns other than `name` and
#'    `time`, joined by `"_"`.}
#' }
#' followed by the kept columns.
#'
#' @seealso [fitErrorModel()]
#'
#' @examples
#' data <- data.frame(name = "y", time = rep(c(0, 1, 2), each = 3),
#'                    value = c(1.0, 1.2, 0.9, 0.6, 0.5, 0.7, 0.3, 0.2, 0.25),
#'                    condition = "ctrl", replicate = rep(1:3, 3))
#' reduceReplicates(data)
#'
#' @export
reduceReplicates <- function(data, select = "condition", datatrans = NULL, keep = NULL, weighted = FALSE) {
  UseMethod("reduceReplicates")
}

#' @rdname reduceReplicates
#' @export
reduceReplicates.data.frame <- function(data, select = "condition", datatrans = NULL, keep = NULL, weighted = FALSE) {

  required <- c("name", "time", "value", "condition")
  if (!all(required %in% names(data)))
    stop("Mandatory column names are: ", paste(required, collapse = ", "))
  if (weighted && !"sigma" %in% names(data))
    stop("Column 'sigma' is required for weighted = TRUE but was not found in the data.")

  if (is.character(datatrans))
    data$value <- eval(parse(text = datatrans), list(x = data$value), parent.frame())

  # A group is one combination of observable, time and the `select` columns,
  # in the order of its first row.
  select <- unique(c("name", "time", "condition", select))
  key <- do.call(paste, c(unname(as.list(data[select])), sep = "\r"))
  group <- factor(key, levels = unique(key))
  first <- !duplicated(key)

  isConstant <- function(column)
    all(vapply(split(column, group), function(v) length(unique(v)) == 1L, TRUE))
  candidates <- setdiff(names(data), c("value", "sigma", "n", select))
  extra <- union(Filter(function(col) isConstant(data[[col]]), candidates),
                 intersect(keep, names(data)))

  moments <- vapply(seq_len(nlevels(group)), function(g) {
    rows <- which(as.integer(group) == g)
    v <- data$value[rows]
    k <- length(v)
    if (k == 1L) return(c(v, NA_real_, 1))
    if (!weighted) return(c(mean(v), stats::sd(v) / sqrt(k), k))
    w <- 1 / data$sigma[rows]^2
    m <- sum(w * v) / sum(w)
    spread <- sum(w * (v - m)^2) / (sum(w) - sum(w^2) / sum(w))
    c(m, sqrt(spread) / sqrt(k), k)
  }, numeric(3))

  label <- setdiff(select, c("name", "time"))
  out <- data.frame(
    time = data$time[first],
    value = moments[1, ],
    sigma = moments[2, ],
    n = as.integer(moments[3, ]),
    name = data$name[first],
    condition = do.call(paste, c(lapply(data[first, label, drop = FALSE], as.character), sep = "_")),
    data[first, extra, drop = FALSE])

  message("Dropped columns: ", paste(setdiff(names(data), names(out)), collapse = ", "))
  out
}


# A bound vector in the order of `par`: unnamed bounds are positional, named
# ones are matched by name and leave the other parameters unbounded.
.parBound <- function(bound, par, default) {
  if (is.null(bound)) return(rep(default, length(par)))
  if (is.null(names(bound))) return(rep_len(bound, length(par)))
  out <- setNames(rep(default, length(par)), names(par))
  hit <- intersect(names(bound), names(par))
  out[hit] <- bound[hit]
  unname(out)
}


#' Fit an Error Model to Replicate Data
#'
#' @description Fits the variance of replicate measurements as a function of
#' their mean by maximum likelihood, assuming that the sample variance follows
#' a scaled \eqn{\chi^2} distribution with \eqn{n - 1} degrees of freedom.
#' Requires the package \pkg{optimx}.
#'
#' @param data A data frame as returned by [reduceReplicates()], with columns
#'   `value` (mean of the replicates), `sigma` (standard error of the mean) and
#'   `n` (number of replicates).
#' @param factors Character vector, the columns of `data` whose value
#'   combinations define groups. The model is fitted separately per group.
#' @param errorModel Character, the variance of a single measurement as an
#'   expression in the mean `x` and the parameters. Default
#'   `"exp(s0)+exp(srel)*x^2"`.
#' @param par Named numeric vector of starting values for the parameters of
#'   `errorModel`. Default `c(s0 = 1, srel = 0.1)`.
#' @param lower,upper Bounds on the parameters, named or in the order of
#'   `par`. Default `NULL`, unbounded.
#' @param plotting Logical. If `TRUE`, plot the variances, the fitted model
#'   and its 68 and 95 percent bands per group. Default `FALSE`.
#' @param blather Logical. If `TRUE`, return the extended data frame described
#'   under Value. Default `FALSE`.
#' @param ... Further arguments passed to `optimx::optimr()`, which runs
#'   `"L-BFGS-B"`.
#'
#' @return `data` with `sigma` replaced by the standard error of the mean that
#'   the fitted model gives.
#'
#'   With `blather = TRUE`, a data frame with the columns of `data`, `sigma`
#'   replaced as above, one column per parameter of `errorModel` with its
#'   fitted value, the variance bands `cbLower68`, `cbUpper68`, `cbLower95` and
#'   `cbUpper95`, the group label `condidnt` and the input `sigma` as
#'   `sigmaLS`. Its attribute `"errorModel"` holds `errorModel`.
#'
#' @seealso [reduceReplicates()]
#'
#' @examplesIf requireNamespace("optimx", quietly = TRUE)
#' set.seed(1)
#' mu <- rep(c(0.1, 1, 10, 100), each = 4)
#' data <- data.frame(name = "y", time = rep(1:4, each = 4), condition = "ctrl",
#'                    value = mu + rnorm(16, sd = sqrt(0.01 + 0.04 * mu^2)))
#' reduced <- reduceReplicates(data)
#' fitErrorModel(reduced, factors = "name", par = c(s0 = -4, srel = -3))
#'
#' @export
#' @importFrom stats qchisq
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_ribbon ylab facet_wrap scale_y_log10 theme
fitErrorModel <- function(data, factors, errorModel = "exp(s0)+exp(srel)*x^2",
                          par = c(s0 = 1, srel = .1), lower = NULL, upper = NULL,
                          plotting = FALSE, blather = FALSE, ...) {

  .require_ns("optimx", "fitErrorModel()")
  model <- str2expression(errorModel)
  callEnv <- parent.frame()
  variance <- function(p, x) eval(model[[length(model)]], c(as.list(p), list(x = x)), callEnv)

  # A single factor column labels the groups as it is, several are pasted.
  label <- if (length(factors) == 1L) data[[factors]] else
    do.call(paste, unname(as.list(data[factors])))
  sampleSd <- data$sigma * sqrt(data$n)
  dof <- data$n - 1
  bounds <- list(lower = .parBound(lower, par, -Inf), upper = .parBound(upper, par, Inf))

  # The sample variance times dof / variance is chi^2(dof), whose likelihood
  # depends on the model variance only through dof * (log v + s^2 / v).
  estimates <- matrix(NA_real_, nrow(data), length(par), dimnames = list(NULL, names(par)))
  single <- rep(NA_real_, nrow(data))
  for (g in unique(label)) {
    rows <- which(label == g)
    nll <- function(p) {
      v <- variance(p, data$value[rows])
      sum(dof[rows] * (log(v) + sampleSd[rows]^2 / v), na.rm = TRUE)
    }
    best <- optimx::optimr(par, nll, method = "L-BFGS-B",
                           lower = bounds$lower, upper = bounds$upper, ...)$par
    estimates[rows, ] <- rep(best, each = length(rows))
    single[rows] <- variance(best, data$value[rows])
  }

  result <- data
  result$sigma <- sqrt(single)
  band <- function(level, upper)
    result$sigma^2 * stats::qchisq((1 - level) / 2, dof, lower.tail = !upper) / dof
  result <- cbind(result, as.data.frame(estimates),
                  cbLower68 = band(.683, FALSE), cbUpper68 = band(.683, TRUE),
                  cbLower95 = band(.955, FALSE), cbUpper95 = band(.955, TRUE))
  result$condidnt <- label
  result$sigmaLS <- data$sigma
  attr(result, "errorModel") <- errorModel

  if (plotting) {
    shown <- transform(result, observed = sampleSd^2, predicted = single)
    print(ggplot(shown, aes(x = value)) +
            geom_point(aes(y = observed)) +
            geom_line(aes(y = predicted)) +
            geom_ribbon(aes(ymin = cbLower95, ymax = cbUpper95), alpha = 0.3) +
            geom_ribbon(aes(ymin = cbLower68, ymax = cbUpper68), alpha = 0.3) +
            facet_wrap("condidnt", scales = "free") +
            scale_y_log10() + ylab("variance") + theme_dMod())
  }

  result$sigma <- result$sigma / sqrt(data$n)
  if (blather) return(result)
  data$sigma <- result$sigma
  data
}

