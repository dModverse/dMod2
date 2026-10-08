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
#'   \code{fixed}. `objfun` is deprecated, use `fitted`.
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

  dotArgs <- .renameArgs(list(...), c(objfun = "fitted"), "profile")
  # Bound in this frame, so the parallel workers receive it.
  objfun <- fitted

  # Guarantee that pars is named numeric without deriv attribute
  sanePars <- sanitizePars(pars, dotArgs$fixed)
  pars <- sanePars$pars
  fixed <- sanePars$fixed
  
  
  # Initialize control parameters depending on method
  method  <- match.arg(method)
  side <- match.arg(side)
  
  
  if (method == "integrate") {
    sControl <- list(stepsize = 1e-4, min = 1e-4, max = Inf, atol = 1e-2, rtol = 1e-2, limit = 500, stop = "value")
    aControl <- list(gamma = 1, W = "hessian", reoptimize = FALSE, correction = 1, reg = .Machine$double.eps)
    oControl <- list(rinit = .1, rmax = 10, iterlim = 10)
  }
  if (method == "optimize") {
    sControl <- list(stepsize = 1e-2, min = 1e-4, max = Inf, atol = 1e-1, rtol = 1e-1, limit = 100, stop = "value")
    aControl <- list(gamma = 0, W = "identity", reoptimize = TRUE, correction = 1, reg = 0)
    oControl <- list(rinit = .1, rmax = 10, iterlim = 100)
  }
  
  # cores = c(pars = , conditions = ) splits the two axes; see mstrust().
  .cc <- .splitCores(cores, "pars")
  coresConditions <- .cc$conditions
  cores <- .sanitizeCores(min(length(whichPar), .cc$outer))
  if (!is.null(coresConditions) && cores == 1L)
    options(dMod.cores = coresConditions)
  
  # Substitute user-set control parameters
  if (!is.null(stepControl)) sControl[match(names(stepControl), names(sControl))] <- stepControl
  if (!is.null(algoControl)) aControl[match(names(algoControl), names(aControl))] <- algoControl
  oControl <- .trustControl(oControl, optControl, label = "optControl")
  
  
  # Create interRes folder for cautiousMode
  if (cautiousMode){
    interResFolder <- "profiles-interRes"
    dir.create(interResFolder,showWarnings = FALSE)
  }
  
  
  # Start cluster if on windows
  if (cores > 1) {
    
    if (Sys.info()[['sysname']] == "Windows" ||
        (!is.null(coresConditions) && coresConditions > 1L)) {
      
      cluster <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl = cluster)
      
      parallel::clusterCall(cl = cluster, function(x) .libPaths(x), .libPaths())
      
      varlist <- ls()
      # Exclude things like "missing argument"
      varlist <- c("objfun", "whichPar", "alpha", "delta", "limits", "method", "verbose", "cores",
                   "pars", "fixed", "dotArgs", "coresConditions",
                   "sControl", "aControl", "oControl")
      parallel::clusterExport(cluster, envir = environment(), varlist = varlist)
      
    } else {
      
      doParallel::registerDoParallel(cores = cores)
      
    }
    
    "%mydo%" <- foreach::"%dopar%"
    
    
  } else {
    
    "%mydo%" <- foreach::"%do%"
    
  }
  
  
  # Convert whichPar to index vector
  if (is.character(whichPar)) whichPar <- which(names(pars) %in% whichPar)
  
  loaded_packages <- .packages()  
  out <- foreach::foreach(whichIndex = whichPar, 
                          .packages = loaded_packages, 
                          .inorder = TRUE,
                          .options.multicore = list(preschedule = FALSE)) %mydo% {
                            
                            if (inherits(objfun, "fn")) loadDLL(objfun)
                            if (!is.null(coresConditions))
                              options(dMod.cores = coresConditions)
                            
                            
                            whichPar.name <- names(pars)[whichIndex]
                            
                            
                            
                            ## Functions needed during profile computation -----------------------
                            # A profile step starts next to an optimum, where
                            # single shooting is enough and cheaper.
                            obj.opt <- .singleShooting(objfun)
                            obj.prof <- function(p, ...) {
                              out <- objfun(p, ...)
                              # If "identity", substitute hessian such that steps are in whichPar-direction.
                              Id <- diag(1/.Machine$double.eps, length(out$gradient))
                              Id[whichIndex, whichIndex] <- 1
                              colnames(Id) <- rownames(Id) <- names(out$gradient)
                              
                              W <- match.arg(aControl$W[1], c("hessian", "identity"))
                              out$hessian <- switch(W,
                                                    "hessian" = out$hessian,
                                                    "identity" = Id)
                              return(out)    
                            }
                            
                            pseudoinverse <- function(m, tol) {
                              msvd <- svd(m)
                              index <- which(abs(msvd$d) > max(dim(m))*max(msvd$d)*tol) 
                              if (length(index) == 0) {
                                out <- array(0, dim(m)[2:1])
                              }
                              else {
                                out <- msvd$u[,index] %*% (1/msvd$d[index] * t(msvd$v)[index,])
                              }
                              attr(out, "valid") <- 1:length(msvd$d) %in% index
                              return(out)
                            }
                            
                            constraint <- function(p) {
                              value <- p[whichIndex] - pars[whichIndex]
                              gradient <- rep(0, length(p))
                              gradient[whichIndex] <- 1
                              return(list(value = value, gradient = gradient))
                            }
                            lagrange <- function(y) {
                              
                              # initialize values
                              p <- y
                              lambda <- 0
                              out <- do.call(obj.prof, c(list(p = p), dotArgs))
                              g.original <- constraint(p)
                              
                              # evaluate derivatives and constraints
                              g     <- direction * g.original$value
                              gdot  <- direction * g.original$gradient
                              ldot  <- out$gradient
                              lddot <- out$hessian 
                              
                              # compute rhs of profile ODE
                              M <- rbind(cbind(lddot, gdot), 
                                         matrix(c(gdot, 0), nrow=1))
                              
                              v <- c(-rep(gamma, length(p))*(ldot + lambda*gdot), 1)
                              v0 <- c(-rep(0, length(p))*(ldot + lambda*gdot), 1)
                              
                              W <- pseudoinverse(M, tol = aControl$reg)
                              valid <- attr(W, "valid")
                              if(any(!valid)) {
                                dy <- try(as.vector(W%*%v)[1:length(p)], silent=FALSE)
                                dy0 <- try(as.vector(W%*%v0)[1:length(p)], silent=FALSE)
                                dy[!valid[1:length(p)]] <- dy0[!valid[1:length(p)]] <- 0
                                dy[whichIndex] <- dy0[whichIndex] <- direction
                                warning(paste0("Iteration ", i, ": Some singular values of the Hessian are below the threshold. Optimization will be performed."))
                              } else {
                                dy <- try(as.vector(W%*%v)[1:length(p)], silent=FALSE)
                                dy0 <- try(as.vector(W%*%v0)[1:length(p)], silent=FALSE)
                              }
                              
                              
                              if(!inherits(dy, "try-error")) {
                                names(dy) <- names(y) 
                                correction <- sqrt(sum((dy-dy0)^2))/sqrt(sum(dy^2))
                              } else {
                                dy <- NA
                                correction <- 0
                                warning(paste0("Iteration ", i, ": Impossible to invert Hessian. Trying to optimize instead."))
                              }
                              
                              # Numeric attributes of the objective value become
                              # profile modes, and those have to add up to the
                              # total. chi2 is a decomposition of one of them,
                              # not a term next to it.
                              out.attributes <- attributes(out)[sapply(attributes(out), is.numeric)]
                              out.attributes <- out.attributes[
                                !grepl("^chi2($|_)", names(out.attributes))]
                              out.attributes.names <- names(out.attributes)
                              
                              
                              return(c(list(dy = dy, 
                                            value = out$value, 
                                            gradient = out$gradient, 
                                            correction = correction, 
                                            valid = valid, 
                                            attributes = out.attributes.names),
                                       out.attributes))
                              
                              
                              
                              
                            }
                            doIteration <- function() {
                              
                              optimize <- aControl$reoptimize
                              # Check for error in evaluation of lagrange()
                              if(is.na(dy[1])) {
                                #cat("Evaluation of lagrange() not successful. Will optimize instead.\n")
                                optimize <- TRUE
                                y.try <- y
                                y.try[whichIndex] <- y[whichIndex] + direction*stepsize
                                rinit <- oControl$rinit
                              } else {
                                dy.norm <- sqrt(sum(dy^2))
                                rinit <- min(c(oControl$rinit, 3*dy.norm))
                                y.try <- y + dy
                                if(any(!lagrange.out$valid)) optimize <- TRUE
                              }
                              
                              # Do reoptimization if requested or necessary
                              if(optimize) {      
                                parinit.opt <- y.try[-whichIndex]
                                fixed.opt <- c(fixed, y.try[whichIndex])
                                
                                arglist <- c(list(objfun = obj.opt, parinit = parinit.opt, fixed = fixed.opt, rinit = rinit), 
                                             oControl[names(oControl)!="rinit"],
                                             dotArgs[names(dotArgs) != "fixed"])
                                
                                
                                myfit <- try(do.call(trust, arglist), silent=FALSE)
                                if(!inherits(myfit, "try-error")) {
                                  y.try[names(myfit$argument)] <- as.vector(myfit$argument)  
                                } else {
                                  warning("Optimization not successful. Profile may be erroneous.")
                                }
                                
                              }
                              
                              return(y.try)
                              
                            }
                            doAdaption <- function() {
                              
                              lagrange.out.try <- lagrange(y.try)
                              valid <- TRUE
                              
                              # Predicted change of the objective value
                              dobj.pred <- sum(lagrange.out$gradient*(y.try - y))
                              dobj.fact <- lagrange.out.try$value - lagrange.out$value
                              correction <- lagrange.out.try$correction
                              
                              # Gamma adaption based on amount of actual correction
                              if (correction > aControl$correction) gamma <- gamma/2
                              if (correction < 0.5*aControl$correction) gamma <- min(c(aControl$gamma, gamma*2))
                              
                              # Stepsize adaption based on difference in predicted change of objective value
                              if (abs(dobj.fact - dobj.pred) > sControl$atol & stepsize > sControl$min) {
                                stepsize <- max(c(stepsize/1.5, sControl$min))
                                valid <- FALSE
                              }
                              if (abs(dobj.fact - dobj.pred) < .3*sControl$atol | abs((dobj.fact - dobj.pred)/dobj.fact) < .3*sControl$rtol) {
                                stepsize <- min(c(stepsize*2, sControl$max))
                              }
                              
                              ## Verbose
                              if (verbose) {
                                # Compute progres
                                diff.thres <- diff.steps <- diff.limit <- 0
                                if (threshold < Inf)
                                  diff.thres <- 1 - max(c(0, min(c(1, (threshold - lagrange.out.try$value)/delta))))
                                if (sControl$limit < Inf)
                                  diff.steps <- i/sControl$limit
                                diff.limit <- switch(as.character(sign(constraint.out$value)),
                                                     "1"  = 1 - (limits[2] - constraint.out$value)/limits[2],
                                                     "-1" = diff.limit <- 1 - (limits[1] - constraint.out$value)/limits[1],
                                                     "0"  = 0)
                                
                                percentage <- max(c(diff.thres, diff.steps, diff.limit), na.rm = TRUE)*100
                                progressBar(percentage)
                                
                                
                                #cat("diff.thres:", diff.thres, "diff.steps:", diff.steps, "diff.limit:", diff.limit)
                                myvalue <- format(substr(lagrange.out$value  , 0, 8), width = 8)
                                myconst <- format(substr(constraint.out$value, 0, 8), width = 8)
                                mygamma <- format(substr(gamma               , 0, 8), width = 8)
                                myvalid <- all(lagrange.out$valid)
                                cat("\tvalue:", myvalue, "constraint:", myconst, "gamma:", mygamma, "valid:", myvalid) 
                              }
                              
                              
                              
                              return(list(lagrange = lagrange.out.try, stepsize = stepsize, gamma = gamma, valid = valid))
                              
                              
                            }
                            
                            ## Compute profile -------------------------------------------------
                            
                            # Initialize profile
                            i <- 0 
                            direction <- 1
                            gamma <- aControl$gamma
                            stepsize <- sControl$stepsize
                            ini <- pars
                            
                            lagrange.out <- lagrange(ini)
                            constraint.out <- constraint(pars)
                            
                            if (is.null(delta)) delta <- qchisq(1 - alpha, 1)
                            
                            
                            
                            threshold <- lagrange.out[[sControl$stop]] + delta
                            out.attributes <- unlist(lagrange.out[lagrange.out$attributes])
                            
                            out <- c(value = lagrange.out$value, 
                                     constraint = as.vector(constraint.out$value), 
                                     stepsize = stepsize, 
                                     gamma = gamma, 
                                     whichPar = whichIndex,
                                     out.attributes, ini)
                            if (side %in% c("right", "both")) {
                              # Compute right profile
                              if (verbose) {
                                cat("Compute right profile\n")
                              }
                              direction <- 1
                              gamma <- aControl$gamma
                              stepsize <- sControl$stepsize
                              y <- ini
                              
                              lagrange.out <- lagrange.out
                              constraint.out <- constraint.out
                              
                              while (i < sControl$limit) {

                                ## Iteration step
                                # Clamp so the last step lands on the limit instead of overshooting
                                stepsize <- min(stepsize, limits[2] - constraint.out$value)
                                sufficient <- FALSE
                                retry <- 0
                                while (!sufficient & retry < 5) {
                                  dy <- stepsize*lagrange.out$dy
                                  y.try <- try(doIteration(), silent = TRUE)
                                  out.try <- try(doAdaption(), silent = TRUE)
                                  if (inherits(y.try, "try-error") | inherits(out.try, "try-error")) {
                                    sufficient <- FALSE
                                    stepsize <- stepsize/1.5
                                    retry <- retry + 1
                                  } else {
                                    sufficient <- out.try$valid
                                    stepsize <- out.try$stepsize
                                  }

                                }
                                if (inherits(y.try, "try-error") | inherits(out.try, "try-error")) break
                                
                                
                                ## Set values
                                y <- y.try
                                lagrange.out <- out.try$lagrange
                                constraint.out <- constraint(y.try)
                                stepsize <- out.try$stepsize
                                gamma <- out.try$gamma
                                out.attributes <- unlist(lagrange.out[lagrange.out$attributes])
                                
                                ## Return values 
                                out <- rbind(out, 
                                             c(value = lagrange.out$value, 
                                               constraint = as.vector(constraint.out$value), 
                                               stepsize = stepsize, 
                                               gamma = gamma, 
                                               whichPar = whichIndex,
                                               out.attributes, 
                                               y))
                                
                                if (cautiousMode)
                                  .profileDump(out, whichPar.name, pars, out.attributes,
                                               interResFolder, "right")
                                
                                value <- lagrange.out[[sControl$stop]]
                                if (value > threshold | constraint.out$value >= limits[2] - sqrt(.Machine$double.eps)) break
                                
                                i <- i + 1
                                
                              }
                            }
                            
                            if (side %in% c("left", "both")) {
                              # Compute left profile
                              if (verbose) {
                                cat("\nCompute left profile\n")
                              }
                              i <- 0
                              direction <- -1
                              gamma <- aControl$gamma
                              stepsize <- sControl$stepsize
                              y <- ini
                              
                              lagrange.out <- lagrange(ini)
                              constraint.out <- constraint(pars)
                              
                              while (i < sControl$limit) {

                                ## Iteration step
                                # Clamp so the last step lands on the limit instead of overshooting
                                stepsize <- min(stepsize, constraint.out$value - limits[1])
                                sufficient <- FALSE
                                retry <- 0
                                while (!sufficient & retry < 5) {
                                  dy <- stepsize*lagrange.out$dy
                                  y.try <- try(doIteration(), silent = TRUE)
                                  out.try <- try(doAdaption(), silent = TRUE)
                                  if (inherits(y.try, "try-error") | inherits(out.try, "try-error")) {
                                    sufficient <- FALSE
                                    stepsize <- stepsize/1.5
                                    retry <- retry + 1
                                  } else {
                                    sufficient <- out.try$valid
                                    stepsize <- out.try$stepsize
                                  }

                                }
                                if (inherits(y.try, "try-error") | inherits(out.try, "try-error")) break
                                
                                ## Set values
                                y <- y.try
                                lagrange.out <- out.try$lagrange
                                constraint.out <- constraint(y.try)
                                stepsize <- out.try$stepsize
                                gamma <- out.try$gamma
                                out.attributes <- unlist(lagrange.out[lagrange.out$attributes])
                                
                                
                                ## Return values
                                out <- rbind(c(value = lagrange.out$value, 
                                               constraint = as.vector(constraint.out$value), 
                                               stepsize = stepsize, 
                                               gamma = gamma,
                                               whichPar = whichIndex,
                                               out.attributes,
                                               y), 
                                             out)
                                
                                if (cautiousMode)
                                  .profileDump(out, whichPar.name, pars, out.attributes,
                                               interResFolder, "left")
                                
                                
                                value <- lagrange.out[[sControl$stop]]
                                if (value > threshold | constraint.out$value <= limits[1] + sqrt(.Machine$double.eps)) break
                                
                                i <- i + 1
                                
                              }
                            }
                            # Output
                            out <- as.data.frame(out)
                            out$whichPar <- whichPar.name
                            parframe(
                              out,
                              parameters = names(pars),
                              metanames = c("value", "constraint", "stepsize", "gamma", "whichPar"),
                              objAttributes = names(out.attributes)
                            )
                            
                          }
  
  
  
  if (Sys.info()[['sysname']] == "Windows" & cores > 1) {
    
    parallel::stopCluster(cluster)
    doParallel::stopImplicitCluster()
    
  }
  
  # .. Prepare output -----
  outncol <- vapply(out, ncol, 1)
  failed <- outncol != max(outncol)
  if (any(failed)) message("The following profiles failed: ", paste0(whichPar[failed], collapse = ", "))
  out <- out[!failed]  
  
  do.call(rbind, out)
  
  
}

#' @rdname profile
#' @export
profile.function <- function(fitted, ...) profile.objfn(fitted, ...)



#' Progress Bar
#'
#' Prints a progress bar on the current console line.
#'
#' @param percentage Numeric between 0 and 100.
#' @param size Integer, width of the bar in characters. Default \code{50}.
#' @param number Logical, whether the percentage is printed after the bar.
#'   Default \code{TRUE}.
#' @return \code{NULL}, invisibly; called for its output.
#' @keywords internal
#' @export
progressBar <- function(percentage, size = 50, number = TRUE) {
  
  if(percentage < 0) percentage <- 0
  if(percentage > 100) percentage <- 100
  
  out <- paste("\r|", paste(rep("=", round(size*percentage/100)), collapse=""), paste(rep(" ", size-round(size*percentage/100)), collapse=""), "|", sep="")
  cat(out)
  if(number) cat(format(paste(" ", round(percentage), "%", sep=""), width=5))
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
#'   example `"value"` or `"data"`. Default `"data"`.
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
  obj.attributes <- attr(profile, "obj.attributes")

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
#' A method of [stats::vcov()] for the result of [trust()].
#'
#' @param object A result of [trust()], class `trustfit`, or one fit of a
#'   parlist.
#' @param parupper,parlower Named upper and lower bounds. Parameters at or
#'   beyond a bound count as fixed, as do those that `object$atBound` marks.
#'   Default `NULL`.
#' @param ... `fit` is deprecated, use `object`.
#'
#' @return Square matrix with the parameter names as dimnames: the inverse of
#'   half the Hessian of `object` (pseudo-inverse when singular), with rows
#'   and columns of fixed parameters set to zero. All `NA` when `object` has
#'   no Hessian, `NULL` when it has no `argument`.
#' @seealso [trust()]
#' @examples
#' obj <- function(x, ...) {
#'   list(value = sum(x^2 / c(1, 4)),
#'        gradient = 2 * x / c(1, 4),
#'        hessian = diag(2 / c(1, 4)))
#' }
#' fit <- trust(obj, c(a = 1, b = 1))
#' vcov(fit)
#' @importFrom stats vcov
#' @export
vcov.trustfit <- function(object, parupper = NULL, parlower = NULL, ...) {

  .renameArgs(list(...), c(fit = "object"), "vcov", strict = TRUE)
  fit <- object
  hessian__ <- fit[["hessian"]]
  arg__ <- fit[["argument"]]
  
  if (is.null(arg__)) return()

  if (is.null(hessian__)) {
    vcov__ <- matrix(NA, nrow = length(arg__), ncol = length(arg__), 
                     dimnames = list(names(arg__), names(arg__)))
    return(vcov__)
  }
  
  # Which parameters are held by a bound rather than determined by the data.
  # `stepControl$boundary = "reflective"` keeps iterates strictly inside the
  # box, so an exact comparison against the bound never fires; it reports the
  # activity itself. The comparison below is the fallback for fits without that
  # field, and is relaxed to a relative tolerance for the same reason.
  fixed <- NULL
  atBound__ <- fit[["atBound"]]
  if (!is.null(atBound__) && !is.null(names(atBound__)))
    fixed <- union(fixed, names(atBound__)[which(as.logical(atBound__))])
  reltol__ <- 1e-6
  if (!is.null(parupper)) {
    myarg__ <- arg__[names(parupper)]
    fixed <- union(fixed, names(myarg__)[
      which(myarg__ >= parupper - reltol__ * pmax(1, abs(parupper)))])
  }
  if (!is.null(parlower)) {
    myarg__ <- arg__[names(parlower)]
    fixed <- union(fixed, names(myarg__)[
      which(myarg__ <= parlower + reltol__ * pmax(1, abs(parlower)))])
  }

  vcov__ <- 0*hessian__
  is_fixed__ <- colnames(hessian__) %in% fixed
  subhessian__ <- hessian__[!is_fixed__, !is_fixed__, drop = FALSE]
  subvcov__ <- try(solve(0.5*subhessian__), silent = TRUE)
  if (inherits(subvcov__, "try-error")) subvcov__ <- MASS::ginv(0.5*subhessian__)
  vcov__[!is_fixed__, !is_fixed__] <- subvcov__
  
  # This part should not be necessary due to regularization usually done
  # Perform identifiability check based on
  
  return(vcov__) 
  
}




#' Non-Linear Optimisation, Multi-Start
#'
#' @description Runs [trust()] from several starting points: `center` plus
#'   random offsets drawn by `samplefun`, or the rows of a parframe.
#'
#' @param objfun Objective function, see [trust()].
#' @param center Named numeric parameter centre, or a parframe whose rows are
#'   used as starting points; `fits` is then the number of rows. [msParframe()]
#'   gives reproducible starts.
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
#'   and keep the intermediate files. Default `FALSE`.
#' @param startFromCenter Logical. If `TRUE`, the first fit starts at
#'   `center` without offset. Default `FALSE`.
#' @param retry Logical. If `TRUE` (default), a fit that fails with an error
#'   is repeated from a new random start, up to `nTries` attempts. Not used
#'   when `center` is a parframe.
#' @param nTries Maximum number of attempts per fit, including the first.
#'   Default `10`.
#' @param ... Arguments passed by name to `optmethod`, `samplefun` or
#'   `objfun`. Names matching neither the optimiser nor the sampler go to
#'   `objfun`. `start1stfromCenter` is deprecated, use `startFromCenter`.
#'
#' @details Nothing is written to disk unless `output`, `cautiousMode` or a
#'   `traceFile` argument for [trust()] asks for it. Files then go to
#'   `<resultPath>/<name>/trial-<k>-<date>/`: the per-fit files in its
#'   subfolder `interRes/`, which is removed on completion unless
#'   `cautiousMode = TRUE`, the final `parameterList.Rda` (with `output`) and
#'   the log `mstrust.log`. [load.parlist()] reads the fits from `interRes/`.
#'
#' @return A parlist with one entry per fit, failed fits included. Each entry
#'   is the result of the optimiser plus `parinit`, its starting point, or a
#'   list with element `error` for a fit that failed. [as.parframe()] turns
#'   it into a table.
#'
#' @seealso [trust()], [msParframe()], [as.parframe()], [load.parlist()]
#'
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
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
#' @import parallel
mstrust <- function(objfun, center, rinit = .1, rmax = 10, fits = 20, cores = 1, optmethod = "trust",
                    samplefun = "rnorm", resultPath = ".", name = "mstrust",
                    stats = FALSE, output = FALSE, cautiousMode = FALSE,
                    startFromCenter = FALSE, retry = TRUE, nTries = 10L,
                    ...) {

  narrowing <- NULL
  varargslist <- .renameArgs(list(...),
                             c(start1stfromCenter = "startFromCenter"),
                             "mstrust")
  samplefunName <- if (is.character(samplefun)) samplefun else
    paste(deparse(substitute(samplefun)), collapse = "")
  samplefun <- match.fun(samplefun)

  # Check if on Windows
  cores <- .sanitizeCores(cores)
  
  argslist <- list(
    objfun = objfun, center = center, name = name,
    rinit = rinit, rmax = rmax, fits = fits,
    cores = cores, optmethod = optmethod, samplefun = samplefun,
    resultPath = resultPath, stats = stats, output = output,
    cautiousMode = cautiousMode, retry = retry, nTries = nTries)
  argslist <- c(argslist, varargslist)
  
  # Add extra arguments
  argslist$n <- length(center) # How many inital values do we need?
  
  # Determine target function for each function argument.
  # First, define argument names used locally in mstrust().
  # Second, check what trust() and samplefun() accept and check for name clashes.
  # Third, whatever is unused is passed to the objective function objfun().
  nameslocal <- c("name", "center", "fits", "cores", "optmethod", "samplefun",
                  "resultPath", "stats", "narrowing", "output", "cautiousMode",
                  "retry", "nTries")
  # A name that moved into one of trust()'s control lists is no longer a
  # formal, so without this it would silently be routed to objfun instead.
  .trustRejectMoved(names(argslist), "mstrust")
  # Routed by the optimiser actually called: reading trust()'s formals sent
  # every argument another optimiser has and trust() lacks to the objective.
  optfun <- if (is.function(optmethod)) optmethod else match.fun(optmethod)
  namestrust <- intersect(setdiff(names(formals(optfun)), "..."), names(argslist))
  namessample <- intersect(names(formals(samplefun)), names(argslist))
  if (length(intersect(namestrust, namessample)) != 0) {
    stop("Argument names of the optimiser and ", samplefunName, "() clash.")
  }
  
  # Default optimizer
  namesobj <- setdiff(names(argslist), c(namestrust, namessample, nameslocal))
  
  # Assemble argument lists common to all calls in mclapply
  # Sample function
  argssample <- structure(vector("list", length = length(namessample)), names = namessample)
  for (name in namessample) {
    argssample[[name]] <- argslist[[name]]
  }
  
  # Objective function
  argsobj <- structure(vector("list", length = length(namesobj)), names = namesobj)
  for (name in namesobj) {
    argsobj[[name]] <- argslist[[name]]
  }
  
  # Trust optimizer, except for initial values
  argstrust <- structure(vector("list", length = length(namestrust)), names = namestrust)
  for (name in namestrust) {
    argstrust[[name]] <- argslist[[name]]
  }
  
  # Assemble output filenames and folders. Nothing is created unless something
  # was asked for: `output` and `cautiousMode` persist fits, `traceFile` puts
  # the per-fit solver traces here.
  writeFiles <- isTRUE(output) || isTRUE(cautiousMode) ||
    !is.null(argslist[["traceFile"]])

  m_timeStamp <- paste0(format(Sys.time(), "%d-%m-%Y-%H%M%S"))

  resultFolderBase <- file.path(argslist$resultPath, argslist$name)
  m_trial <- paste0("trial-", length(dir(resultFolderBase, pattern = "trial*")) + 1)
  resultFolder <- file.path(resultFolderBase, paste0(m_trial, "-", m_timeStamp))
  interResultFolder <- file.path(resultFolder, "interRes")
  fileLog <- file.path(resultFolder, "mstrust.log")
  fileParList <- file.path(resultFolder, "parameterList.Rda")

  logfile <- NULL
  if (writeFiles) {
    dir.create(path = interResultFolder, showWarnings = FALSE, recursive = TRUE)
    # The error checking leverages that mclapply runs each job in a try().
    logfile <- file(fileLog, open = "w")
  }
  
  # Parameter assignment information
  if (is.null(narrowing) || narrowing[1] == 1) {
    msg <- paste0("Parameter assignment information\n",
                  strpad("mstrust", 12),                        ": ", paste0(nameslocal, collapse = ", "), "\n",
                  strpad("trust", 12),                          ": ", paste0(namestrust, collapse = ", "), "\n",
                  strpad(samplefunName, 12), ": ", paste0(namessample, collapse = ", "), "\n\n")
    #strpad(as.character(argslist$objfun), 12),    ": ", paste0(namesobj, collapse = ", "), "\n\n")
    if (!is.null(logfile)) { writeLines(msg, logfile); flush(logfile) }
  }
  
  # Write narrowing status information to file
  if (!is.null(narrowing)) {
    msg <- paste0("--> Narrowing, run ", narrowing[1], " of ", narrowing[2], "\n",
                  "--> " , fits, " fits to run\n")
    if (!is.null(logfile)) writeLines(msg, logfile)
    flush(logfile)
  }
  
  if(is.parframe(center)) {
    fits <- nrow(center)
  }
  
  
  # cores = c(fits = , conditions = ) splits the two axes. A forked outer axis
  # cannot nest an inner one, cppDE's batch runs serially inside a fork --
  # so an inner axis > 1 selects PSOCK.
  .cc <- .splitCores(cores, "fits")
  coresConditions <- .cc$conditions
  cores <- min(fits, .cc$outer)
  if (!is.null(coresConditions) && cores == 1L)
    options(dMod.cores = coresConditions)
  if (cores > 1) {
    
    # Start cluster if on windows
    if (Sys.info()[['sysname']] == "Windows" ||
        (!is.null(coresConditions) && coresConditions > 1L)) {
      
      cluster <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cluster)
      parallel::clusterCall(cl = cluster, function(x) .libPaths(x), .libPaths())
      
      varlist <- ls()
      # Exclude things like "missing argument"
      varlist <- c("objfun", "center", "argstrust", 
                   "samplefun", "argssample", "argsobj", 
                   "output", "interResultFolder", "logfile", "coresConditions")
      parallel::clusterExport(cluster, envir = environment(), varlist = varlist)
      
    } else {
      
      doParallel::registerDoParallel(cores = cores)
      
    }
    
    "%mydo%" <- foreach::"%dopar%"
    
  } else {
    
    "%mydo%" <- foreach::"%do%"
    
  }
  
  
  loaded_packages <- .packages()  
  m_parlist <- as.parlist(foreach::foreach(i = 1:fits, 
                                           .packages = loaded_packages, 
                                           .inorder = TRUE,
                                           .options.multicore = list(preschedule = FALSE)
  ) %mydo% {
    
    # A plain R function has no shared objects to load.
    if (inherits(objfun, "fn")) suppressMessages(loadDLL(objfun))
    # PSOCK workers do not inherit options; forks do, harmlessly.
    if (!is.null(coresConditions)) options(dMod.cores = coresConditions)
    
    if(is.parframe(center)) {
      argstrust$parinit <- as.parvec(center, i)
    } else {
      if (i == 1 && startFromCenter) {
        # First fit always starts from center
        argstrust$parinit <- center
      } else {
        # All other fits start from random positions
        argstrust$parinit <- center + do.call(samplefun, argssample)
      }
    }

    # Check if traceFile is requested. In that case combine tracefile, with logfolder and fit number
    if (!is.null(argstrust[["traceFile"]])) {
      digits <- floor(log10(fits))
      argstrust[["traceFile"]] <- file.path(resultFolder, paste0(formatC(i, digits = digits, flag = "0"), "_", argslist[["traceFile"]]))
    }

    # Each fit starts without the warm-start roots of the previous fit
    # (cores = 1) or of the parent process (cores > 1, forked).
    try(resetWarmStarts(objfun, verbose = FALSE), silent = TRUE)

    # Retry loop: a try-error or fit$error triggers re-sampling parinit and
    # re-running optmethod, up to nTries times. parframe-supplied centers
    # are skipped (rows are taken as given). Warm-start caches are reset
    # between attempts.
    max_tries <- if (isTRUE(retry) && !is.parframe(center)) as.integer(nTries) else 1L
    fit <- NULL
    for (try_i in seq_len(max_tries)) {
      fit <- try(do.call(optmethod, c(argstrust, argsobj)), silent = !output)
      ok  <- !inherits(fit, "try-error") && is.list(fit) && !any(names(fit) == "error")
      if (ok || try_i == max_tries) break
      argstrust$parinit <- center + do.call(samplefun, argssample)
      try(resetWarmStarts(objfun, verbose = FALSE), silent = TRUE)
    }
    
    # Keep only numeric attributes of object returned by trust()
    attr.fit <- attributes(fit)
    keep.attr <- sapply(attr.fit, is.numeric)
    fit <- fit[1:length(fit)] # deletes attributes
    if (any(keep.attr)) attributes(fit) <- c(attributes(fit), attr.fit[keep.attr]) # attach numeric attributes
    if ("trustfit" %in% attr.fit$class) fit <- .trustFit(fit)
    
    
    # In some crashes a try-error object is returned which is not a list. Since
    # each element in the parlist is assumed to be a list, these cases are wrapped.
    if (!is.list(fit)) {
      f <- list()
      f$error <- fit
      fit <- f
    }
    
    fit$parinit <- argstrust$parinit
    
    # Write current fit to disk
    if (output) {
      saveRDS(fit, file = file.path(interResultFolder, paste0("fit-", i, ".Rda")))
      if (cautiousMode) dput(fit[c("value", "argument", "iterations", "converged")], file = file.path(interResultFolder, paste0("fit-", i, ".R")))
      # Reporting
      # With concurent jobs and everyone reporting, this is a classic race
      # condition. Assembling the message beforhand lowers the risk of interleaved
      # output to the log.
      msgSep <- "-------"
      if (any(names(fit) == "error")) {
        msg <- paste0(msgSep, "\n",
                      "Fit ", i, " failed after ", fit$iterations, " iterations with error\n",
                      "--> ", fit$error,
                      msgSep, "\n")
        
        writeLines(msg, logfile)
        flush(logfile)
      } else {
        msg <- paste0(msgSep, "\n",
                      "Fit ", i, " completed\n",
                      "--> iterations : ", fit$iterations, "\n",
                      "-->  converged : ", fit$converged, "\n",
                      "--> obj. value : ", round(fit$value, digits = 2), "\n",
                      msgSep)
        
        writeLines(msg, logfile)
        flush(logfile)
      }
    }
    return(fit)
  })
  
  if (!is.null(logfile)) close(logfile)

  if (Sys.info()[['sysname']] == "Windows" & cores > 1) {

    parallel::stopCluster(cluster)
    doParallel::stopImplicitCluster()

  }
  
  
  
  # Cull failed and completed fits Two kinds of errors occure. The first returns
  # an object of class "try-error". The reason for these failures are unknown to
  # me. The second returns a list of results from trust(), where one name of the
  # list is error holding an object of class "try-error". These abortions are 
  # due to errors which are captured within trust(). Completed fits return with 
  # a valid result list from trust(), with "error" not part of its names. These
  # fits, can still be unconverged, if the maximim number of iterations was the
  # reason for the return of trust(). Be also aware of fits which converge due
  # to the trust radius hitting rmin. Such fits are reported as converged but
  # are not in truth.
  m_trustFlags.converged = 0
  m_trustFlags.unconverged = 1
  m_trustFlags.error = 2
  m_trustFlags.fatal = 3
  idxStatus <- sapply(m_parlist, function(fit) {
    if (inherits(fit, "try-error") || any(names(fit) == "error")) {
      return(m_trustFlags.error)
    } else if (!any(names(fit) == "converged")) {
      return(m_trustFlags.fatal)
    } else if (fit$converged) {
      return(m_trustFlags.converged)
    } else {
      return(m_trustFlags.unconverged)
    }
  })
  
  
  # Wrap up
  # Write out results
  if (output) saveRDS(m_parlist, file = fileParList)
  
  # Remove temporary files
  if (writeFiles) {
    if (!cautiousMode) {
      unlink(interResultFolder, recursive = TRUE)
    } else {
      unlink(list.files(interResultFolder, "\\.Rda$", full.names = TRUE))
    }
  }
  
  # Show summary
  sum.error <- sum(idxStatus == m_trustFlags.error)
  sum.fatal <- sum(idxStatus == m_trustFlags.fatal)
  sum.unconverged <- sum(idxStatus == m_trustFlags.unconverged)
  sum.converged <- sum(idxStatus == m_trustFlags.converged)
  reasons <- .stopReasonTable(m_parlist)
  msg <- paste0("Multi start trust summary\n",
                "Outcome     : Occurrence\n",
                "Error       : ", sum.error, "\n",
                "Fatal       : ", sum.fatal, " must be 0\n",
                "Unconverged : ", sum.unconverged, "\n",
                "Converged   : ", sum.converged, "\n",
                "           -----------\n",
                "Total       : ", sum.error + sum.fatal + sum.unconverged + sum.converged, paste0("[", fits, "]"), "\n",
                if (!is.null(reasons))
                  paste0("\nTermination reason\n",
                         paste0(formatC(names(reasons), width = -12), ": ",
                                as.integer(reasons), collapse = "\n"), "\n"))
  if (writeFiles) {
    logfile <- file(fileLog, open = "a")
    writeLines(msg, logfile)
    flush(logfile)
    close(logfile)
  }
  
  if (stats) {
    cat(msg)
  }
  
  return(m_parlist)
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
#' msParframe(c(a = 0, b = 100000), 5)
#' 
#' # Parameter specific sigma
#' msParframe(c(a = 0, b = 100000), 5, samplefun = rnorm, sd = c(100, 0.5))
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
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#'
#' @export
load.parlist <- function(folder) {
  # Read in all fits
  m_fileList <- dir(folder, pattern = "\\.Rda$")
  m_parVec <- lapply(m_fileList, function(file) {
    return(readRDS(file.path(folder, file)))
  })
  
  return(as.parlist(m_parVec))
}


#' Reduce Replicates to Mean and Standard Error
#'
#' @description
#' Replaces the replicates of each condition by their mean and the standard
#' error of the mean.
#'
#' @param data A data frame with the columns `name` (observable), `time`,
#'   `value` and `condition`, plus any further columns that define conditions.
#'   Alternatively, the path of a `.csv`, `.xls` or `.xlsx` file with these
#'   columns; Excel files need the package \pkg{openxlsx}.
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
#'  \item{time}{Measurement time point.}
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
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#' @author Simon Beyer, \email{simon.beyer@@fdm.uni-freiburg.de}
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
  # File format definition
  fmtnames <- c("name", "time", "value", "condition")
  if (length(intersect(names(data), fmtnames)) != length(fmtnames)) {
    stop(paste("Mandatory column names are:", paste(fmtnames, collapse = ", ")))
  }
  
  # Check if sigma column is present if weighted = TRUE
  if (weighted && !"sigma" %in% names(data)) {
    stop("Column 'sigma' is required for weighted = TRUE but was not found in the data.")
  }
  
  # Transform data if requested
  if (!is.null(datatrans) && is.character(datatrans)) {
    x <- data$value
    data$value <- eval(parse(text = datatrans))
  }
  
  # Define grouping conditions
  select <- unique(c("name", "time", "condition", select))
  condidnt <- apply(data[select], 1, paste, collapse = "_")
  conditions <- unique(condidnt)
  
  # Identify columns that are consistent across replicates
  potential_cols <- setdiff(names(data), c("value", "sigma", "n", select))
  if (length(potential_cols) == 0) {
    stable_cols <- character(0)
  } else {
    stable_cols <- potential_cols[which(sapply(potential_cols, function(col) {
      res <- tapply(data[[col]], condidnt, function(x) length(unique(x)) == 1)
      all(as.logical(res))
    }))]
  }
  
  # Add columns from 'keep' (if any), even if unstable
  if (!is.null(keep)) {
    keep <- intersect(keep, names(data))  # Make sure the columns exist
    stable_cols <- union(stable_cols, keep)
  }
  
  dropped_cols <- setdiff(names(data), c("time", "value", "sigma", "n", stable_cols))
  
  # Reduce data
  reduct <- do.call(rbind, lapply(conditions, function(cond) {
    conddata <- data[condidnt == cond, ]
    mergecond <- paste(unique(conddata[setdiff(select, c("name", "time"))]), collapse = "_")
    
    # Determine value and sigma
    if (weighted && nrow(conddata) > 1) {
      weights <- 1 / conddata$sigma^2
      mean_val <- sum(weights * conddata$value) / sum(weights)
      sigma_val <- sqrt(sum(weights * (conddata$value - mean_val)^2) / 
                          (sum(weights) - sum(weights^2) / sum(weights))) / sqrt(nrow(conddata))
    } else if (nrow(conddata) > 1) {
      mean_val <- mean(conddata$value)
      sigma_val <- sd(conddata$value) / sqrt(nrow(conddata))
    } else {
      mean_val <- conddata$value
      sigma_val <- NA
    }
    
    data.frame(
      time = conddata[1, "time"],
      value = mean_val,
      sigma = sigma_val,
      n = nrow(conddata),
      name = conddata[1, "name"],
      condition = mergecond,
      conddata[1, stable_cols, drop = FALSE]
    )
  }))
  
  message("Dropped columns: ", paste(setdiff(dropped_cols, names(reduct)), collapse = ", "))
  return(reduct)
}


#' @rdname reduceReplicates
#' @export
reduceReplicates.character <- function(data, select = "condition", datatrans = NULL, keep = NULL, weighted = FALSE) {
  # Ensure the file exists
  if (!file.exists(data)) {
    stop("The specified file does not exist.")
  }
  
  # Determine the file type based on extension
  ext <- tools::file_ext(data)
  if (tolower(ext) == "csv") {
    data <- read.csv(data, stringsAsFactors = FALSE)
  } else if (tolower(ext) %in% c("xls", "xlsx")) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("The 'openxlsx' package is required to read Excel files. Please install it.")
    }
    data <- openxlsx::read.xlsx(data)
  } else {
    stop("Unsupported file format. Only .csv, .xls, and .xlsx are supported.")
  }
  
  # Call the data.frame method
  reduceReplicates(as.data.frame(data), select = select, datatrans = datatrans,
                   keep = keep, weighted = weighted)
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
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#' @author Simon Beyer, \email{simon.beyer@@fdm.uni-freiburg.de}
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
                          par = c(s0 = 1, srel = .1),
                          lower = NULL, upper = NULL,  # Optional: Parametergrenzen
                          plotting = FALSE, blather = FALSE, ...) {

  .require_ns("optimx", "fitErrorModel()")
  # Assemble conditions
  condidnt <- Reduce(paste, subset(data, select = factors))
  conditions <- unique(condidnt)
  
  # Fit error model
  nColData <- ncol(data)
  dataErrorModel <- cbind(data, as.list(par))
  
  for (cond in conditions) {
    subdata <- dataErrorModel[condidnt == cond,]
    x <- subdata$value
    n <- subdata$n
    y <- subdata$sigma * sqrt(n)
    
    # Zielfunktion mit der analytischen Maximum-Likelihood
    obj <- function(par) {
      with(as.list(par), {
        sigma2 <- eval(parse(text = errorModel)) 
        negLogLik <- sum((n - 1) * (log(sigma2) + (y^2 / sigma2)), na.rm = TRUE)
        return(negLogLik)
      })
    }
    
    fit <- optimx::optimr(par, obj, method = "L-BFGS-B",
                          lower = .parBound(lower, par, -Inf),
                          upper = .parBound(upper, par, Inf), ...)
    
    sigma <- sqrt(with(as.list(fit$par), eval(parse(text = errorModel))))
    dataErrorModel[condidnt == cond, ]$sigma <- sigma 
    dataErrorModel[condidnt == cond, -(nColData:1)] <- data.frame(as.list(fit$par))
  }
  
  # Calculate confidence bounds for sigma
  p68 <- (1 - .683) / 2
  p95 <- (1 - .955) / 2
  dataErrorModel$cbLower68 <- dataErrorModel$sigma^2 * qchisq(p = p68, df = dataErrorModel$n - 1) / (dataErrorModel$n - 1)
  dataErrorModel$cbUpper68 <- dataErrorModel$sigma^2 * qchisq(p = p68, df = dataErrorModel$n - 1, lower.tail = FALSE) / (dataErrorModel$n - 1)
  dataErrorModel$cbLower95 <- dataErrorModel$sigma^2 * qchisq(p = p95, df = dataErrorModel$n - 1) / (dataErrorModel$n - 1)
  dataErrorModel$cbUpper95 <- dataErrorModel$sigma^2 * qchisq(p = p95, df = dataErrorModel$n - 1, lower.tail = FALSE) / (dataErrorModel$n - 1)
  
  # Assemble result
  dataErrorModel <- cbind(dataErrorModel, condidnt, sigmaLS = data$sigma)
  attr(dataErrorModel, "errorModel") <- errorModel
  
  # Plot if requested
  if (plotting) {
    print(ggplot(dataErrorModel, aes(x = value)) +
            geom_point(aes(y = sigmaLS^2 * (n))) +
            geom_line(aes(y = sigma^2)) +
            geom_ribbon(aes(ymin = cbLower95, ymax = cbUpper95), alpha = .3) +
            geom_ribbon(aes(ymin = cbLower68, ymax = cbUpper68), alpha = .3) +
            ylab("variance") +
            facet_wrap(~condidnt, scales = "free") +
            scale_y_log10() +
            theme_dMod()
    )}
  
  # Return standard error of the mean
  dataErrorModel$sigma <- dataErrorModel$sigma / sqrt(dataErrorModel$n)
  data$sigma <- dataErrorModel$sigma
  if (blather)
    return(dataErrorModel)
  else 
    return(data)
}
