#' dMod2: Dynamic Modeling and Parameter Estimation in ODE Models
#'
#' Models are built from reactions or equations, compiled to C++ with
#' sensitivities, and composed from observation, prediction and parameter
#' transformation functions as `g * x * p`. Objective functions over several
#' experimental conditions are minimized by trust-region optimization, and
#' parameter uncertainty is assessed by profile likelihoods.
#'
#' @seealso
#' Models: [eqnlist()], [addReaction()], [odemodel()], [Xs()], [Y()], [P()].
#' Data and objectives: [as.datalist()], [normL2()], [constraintL2()].
#' Fitting and uncertainty: [trust()], [mstrust()], [profile()][profile.objfn],
#' [plotProfile()]. Exchange formats: [importPEtab()], [exportPEtab()].
#' Identifiability: [symmetryDetection()].
#'
"_PACKAGE"

#' @useDynLib dMod2, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom cOde getSymbols replaceSymbols prodSymb sensitivitiesSymb
#' @importFrom cOde setForcings odeC
#' @importFrom tidyselect all_of
#' @importFrom stats D approx rnorm runif sd
#' @importFrom utils capture.output glob2rx head modifyList
NULL

utils::globalVariables(c(
  # plotting / dplyr / data.table NSE column names
  "value", "sigma", "condition", "x", "y", "name", "proflist", "delta", "loq", "bloq",
  "sigmaLS", "cbLower95", "cbUpper95", "cbLower68", "cbUpper68",
  "iteration", "idx", "is.zero", "index", "converged", "iterations",
  ".runbgOutput", "terminal", "parvalue", "dataErrorModel",
  "Rate", "weighted.residual", "i", "whichIndex",
  # additional NSE column names flagged by R CMD check
  ".", "..required", ".calls", "ID", "ParValue",
  "combination", "constraint", "data", "e", "fits", "g",
  "hessian", "hyp_step", "keys", "label", "max.dev", "maxres",
  "myeqnlist", "myprofiles", "p", "partable", "partner",
  "termcd_label", "time", "whichPar"
))
