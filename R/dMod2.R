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
#' @section Options:
#' \describe{
#'   \item{`dMod2.cores`}{Cores over the experimental conditions inside one
#'     evaluation of a prediction, observation or objective, the default of
#'     their `cores` argument. Default 1.}
#'   \item{`dMod2.outdir`}{Directory the model constructors such as
#'     [odemodel()], [Y()] and [P()] write their sources and shared objects
#'     to, and [distributedComputing()] takes the model files from. Default
#'     the working directory.}
#'   \item{`dMod2.batch.check`}{`TRUE` evaluates every batch of conditions a
#'     second time condition by condition and stops where the two differ.
#'     Default `FALSE`.}
#'   \item{`dMod2.compile.cmdlimit`}{Longest command line in characters that
#'     [compile()] passes to the compiler and archiver. Default 96000, 8000 on
#'     Windows.}
#' }
#'
#' @section Environment variables:
#' \describe{
#'   \item{`DMOD_LIBSBML_PYTHON`}{Python interpreter with libsbml for the SBML
#'     interface.}
#'   \item{`SYMIDENT_<NAME>`}{Switch `<name>` of symident for
#'     [symmetryDetection()] and [symmetryReduction()], see the Details there.}
#' }
#'
"_PACKAGE"

#' @useDynLib dMod2, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom cOde sensitivitiesSymb
#' @importFrom cOde setForcings odeC
#' @importFrom tidyselect all_of
#' @importFrom stats D approx rnorm runif sd
#' @importFrom utils capture.output glob2rx head modifyList
NULL

utils::globalVariables(c(
  # plotting / dplyr / data.table NSE column names
  "value", "sigma", "condition", "x", "y", "name", "proflist", "delta", "loq", "bloq",
  "sigmaLS", "cbLower95", "cbUpper95", "cbLower68", "cbUpper68",
  "iteration", "path", "idx", "is.zero", "index", "step", "converged", "iterations",
  "terminal", "parvalue", "dataErrorModel",
  "Rate", "weighted.residual", "i", "whichIndex",
  # additional NSE column names flagged by R CMD check
  "..required", ".calls", "ID", "ParValue", "par", "set",
  "combination", "constraint", "data", "e", "fits", "g",
  "hessian", "hyp_step", "keys", "label", "max.dev", "maxres",
  "myeqnlist", "p", "partner",
  "termcd_label", "time", "whichPar"
))
