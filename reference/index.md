# Package index

## Package

- [`dMod2`](https://dmodverse.github.io/dMod2/reference/dMod2-package.md)
  [`dMod2-package`](https://dmodverse.github.io/dMod2/reference/dMod2-package.md)
  : dMod2: Parameter Estimation in Dynamical Models

## Reaction networks and equations

Reactions, compartments, equation vectors, conserved quantities and
steady states.

- [`as.eqnlist()`](https://dmodverse.github.io/dMod2/reference/eqnlist.md)
  [`is.eqnlist()`](https://dmodverse.github.io/dMod2/reference/eqnlist.md)
  [`eqnlist()`](https://dmodverse.github.io/dMod2/reference/eqnlist.md)
  : Generate or Coerce to an Equation List
- [`addReaction()`](https://dmodverse.github.io/dMod2/reference/addReaction.md)
  : Add Reaction to an Equation List
- [`assignCompartment()`](https://dmodverse.github.io/dMod2/reference/assignCompartment.md)
  : Assign States to Compartments
- [`setCompartmentVolume()`](https://dmodverse.github.io/dMod2/reference/setCompartmentVolume.md)
  : Set Compartment Volumes
- [`getReactions()`](https://dmodverse.github.io/dMod2/reference/getReactions.md)
  : Generate a Table of Reactions from an Equation List
- [`getFluxes()`](https://dmodverse.github.io/dMod2/reference/getFluxes.md)
  : Generate List of Fluxes from Equation List
- [`plotFluxes()`](https://dmodverse.github.io/dMod2/reference/plotFluxes.md)
  : Plot Fluxes Given a List of Flux Equations
- [`print(`*`<eqnlist>`*`)`](https://dmodverse.github.io/dMod2/reference/print.eqnlist.md)
  : Print or Pander Equation List
- [`as.data.frame(`*`<eqnlist>`*`)`](https://dmodverse.github.io/dMod2/reference/as.data.frame.eqnlist.md)
  : Coerce Equation List into a Data Frame
- [`subset(`*`<eqnlist>`*`)`](https://dmodverse.github.io/dMod2/reference/subset.eqnlist.md)
  : Subset of an Equation List
- [`write.eqnlist()`](https://dmodverse.github.io/dMod2/reference/write.eqnlist.md)
  : Write Equation List into a CSV File
- [`is.eqnvec()`](https://dmodverse.github.io/dMod2/reference/eqnvec.md)
  [`eqnvec()`](https://dmodverse.github.io/dMod2/reference/eqnvec.md) :
  Generate Equation Vector Object
- [`as.eqnvec()`](https://dmodverse.github.io/dMod2/reference/as.eqnvec.md)
  : Coerce to an Equation Vector
- [`format(`*`<eqnvec>`*`)`](https://dmodverse.github.io/dMod2/reference/format.eqnvec.md)
  : Format Equation Vector with Spaces
- [`print(`*`<eqnvec>`*`)`](https://dmodverse.github.io/dMod2/reference/print.eqnvec.md)
  : Print Equation Vector
- [`summary(`*`<eqnvec>`*`)`](https://dmodverse.github.io/dMod2/reference/summary.eqnvec.md)
  : Summary of an Equation Vector
- [`conservedQuantities()`](https://dmodverse.github.io/dMod2/reference/conservedQuantities.md)
  : Conserved Quantities of an Equation List
- [`customTotals()`](https://dmodverse.github.io/dMod2/reference/customTotals.md)
  : Set or Reset User-Defined Conserved Quantities
- [`steadyStates()`](https://dmodverse.github.io/dMod2/reference/steadyStates.md)
  : Analytical Steady States
- [`log10Transform()`](https://dmodverse.github.io/dMod2/reference/log10Transform.md)
  : Transform an ODE to log10 Coordinates
- [`eventlist()`](https://dmodverse.github.io/dMod2/reference/eventlist.md)
  [`addEvent()`](https://dmodverse.github.io/dMod2/reference/eventlist.md)
  : Eventlist
- [`as.eventlist()`](https://dmodverse.github.io/dMod2/reference/as.eventlist.md)
  : Coerce to Eventlist
- [`getSymbols()`](https://dmodverse.github.io/dMod2/reference/getSymbols.md)
  : Symbols in Expressions
- [`replaceSymbols()`](https://dmodverse.github.io/dMod2/reference/replaceSymbols.md)
  : Replace Symbols in Expressions
- [`resolveRecurrence()`](https://dmodverse.github.io/dMod2/reference/resolveRecurrence.md)
  : Substitute Earlier Entries into Later Ones

## Prediction

Compiled ODE models and the prediction functions built on them.

- [`odemodel()`](https://dmodverse.github.io/dMod2/reference/odemodel.md)
  : Compiled ODE Model
- [`Xs()`](https://dmodverse.github.io/dMod2/reference/Xs.md) :
  Prediction Function of an ODE Model with Sensitivities
- [`Xf()`](https://dmodverse.github.io/dMod2/reference/Xf.md) :
  Prediction Function of an ODE Model without Sensitivities
- [`Xd()`](https://dmodverse.github.io/dMod2/reference/Xd.md) :
  Prediction Function Interpolating Parameters
- [`Xt()`](https://dmodverse.github.io/dMod2/reference/Xt.md) :
  Prediction Function Returning Time
- [`prdfn()`](https://dmodverse.github.io/dMod2/reference/prdfn.md) :
  Prediction Function
- [`predict(`*`<prdfn>`*`)`](https://dmodverse.github.io/dMod2/reference/predict.prdfn.md)
  : Model Predictions
- [`prdframe()`](https://dmodverse.github.io/dMod2/reference/prdframe.md)
  : Prediction Frame
- [`as.prdlist()`](https://dmodverse.github.io/dMod2/reference/prdlist.md)
  [`prdlist()`](https://dmodverse.github.io/dMod2/reference/prdlist.md)
  : Prediction List
- [`compile()`](https://dmodverse.github.io/dMod2/reference/compile.md)
  : Compile Model-Related C/C++ Code
- [`loadDLL()`](https://dmodverse.github.io/dMod2/reference/loadDLL.md)
  : Load Shared Objects for dMod Objects
- [`getLocalDLLs()`](https://dmodverse.github.io/dMod2/reference/getLocalDLLs.md)
  : Loaded Shared Objects in the Working Directory
- [`modelname()`](https://dmodverse.github.io/dMod2/reference/modelname.md)
  [`` `modelname<-`() ``](https://dmodverse.github.io/dMod2/reference/modelname.md)
  : Get and Set Modelname
- [`controls()`](https://dmodverse.github.io/dMod2/reference/controls.md)
  [`` `controls<-`() ``](https://dmodverse.github.io/dMod2/reference/controls.md)
  : List, Get and Set Controls

## Observation and parameter functions

The function objects composed as `g * x * p`: observation functions,
parameter transformations, concatenation and direct sum.

- [`Y()`](https://dmodverse.github.io/dMod2/reference/Y.md) :
  Observation Function
- [`obsfn()`](https://dmodverse.github.io/dMod2/reference/obsfn.md) :
  Observation Function
- [`P()`](https://dmodverse.github.io/dMod2/reference/P.md) : Parameter
  Transformation
- [`Pexpl()`](https://dmodverse.github.io/dMod2/reference/Pexpl.md) :
  Explicit Parameter Transformation
- [`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md) :
  Implicit Parameter Transformation
- [`parfn()`](https://dmodverse.github.io/dMod2/reference/parfn.md) :
  Parameter Transformation Function
- [`define()`](https://dmodverse.github.io/dMod2/reference/define.md)
  [`insert()`](https://dmodverse.github.io/dMod2/reference/define.md)
  [`branch()`](https://dmodverse.github.io/dMod2/reference/define.md) :
  Construct and Modify Parameter Transformations
- [`repar()`](https://dmodverse.github.io/dMod2/reference/repar.md) :
  Reparameterization
- [`print(`*`<parfn>`*`)`](https://dmodverse.github.io/dMod2/reference/print.parfn.md)
  : Print Parameter Transformations
- [`Id()`](https://dmodverse.github.io/dMod2/reference/Id.md) : Identity
  Function
- [`` `*`( ``*`<fn>`*`)`](https://dmodverse.github.io/dMod2/reference/times-.fn.md)
  : Concatenation of Functions
- [`` `+`( ``*`<fn>`*`)`](https://dmodverse.github.io/dMod2/reference/plus-.fn.md)
  : Direct Sum of Functions
- [`dModfn`](https://dmodverse.github.io/dMod2/reference/dModfn.md) :
  Calling dMod Functions

## Data

- [`as.datalist()`](https://dmodverse.github.io/dMod2/reference/datalist.md)
  [`` `names<-`( ``*`<datalist>`*`)`](https://dmodverse.github.io/dMod2/reference/datalist.md)
  [`is.datalist()`](https://dmodverse.github.io/dMod2/reference/datalist.md)
  [`c(`*`<datalist>`*`)`](https://dmodverse.github.io/dMod2/reference/datalist.md)
  [`datalist()`](https://dmodverse.github.io/dMod2/reference/datalist.md)
  : Generate a Datalist Object
- [`` `+`( ``*`<datalist>`*`)`](https://dmodverse.github.io/dMod2/reference/plus-.datalist.md)
  : Direct Sum of Datasets
- [`as.data.frame(`*`<datalist>`*`)`](https://dmodverse.github.io/dMod2/reference/as.data.frame.dMod.md)
  [`as.data.frame(`*`<prdlist>`*`)`](https://dmodverse.github.io/dMod2/reference/as.data.frame.dMod.md)
  : Coerce to a Data Frame
- [`covariates()`](https://dmodverse.github.io/dMod2/reference/covariates.md)
  : Access the Covariates in the Data
- [`lbind()`](https://dmodverse.github.io/dMod2/reference/lbind.md) :
  Bind a Named List of Data Frames into One Data Frame
- [`combine()`](https://dmodverse.github.io/dMod2/reference/combine.md)
  : Combine Data Frames or Matrices by Rows
- [`reduceReplicates()`](https://dmodverse.github.io/dMod2/reference/reduceReplicates.md)
  : Reduce Replicates to Mean and Standard Error
- [`coordTransform()`](https://dmodverse.github.io/dMod2/reference/coordTransform.md)
  : Coordinate Transformation for Data Frames

## Structural identifiability

Symmetries of ODE models by the Python package symident.

- [`symmetryDetection()`](https://dmodverse.github.io/dMod2/reference/symmetryDetection.md)
  : Structural Non-Identifiabilities of an ODE Model
- [`symmetryReduction()`](https://dmodverse.github.io/dMod2/reference/symmetryReduction.md)
  : Remove Detected Symmetries by Reparametrization
- [`discreteSymmetries()`](https://dmodverse.github.io/dMod2/reference/discreteSymmetries.md)
  : Discrete Symmetries
- [`reconstControl()`](https://dmodverse.github.io/dMod2/reference/reconstControl.md)
  : Settings of the Observability Engine

## Likelihood and prior

Likelihoods of the data, priors and constraints, sums of objectives.

- [`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md) :
  L2 Norm Between Data and Model Prediction
- [`normT()`](https://dmodverse.github.io/dMod2/reference/normT.md) :
  Student-t likelihood between data and model prediction
- [`res()`](https://dmodverse.github.io/dMod2/reference/res.md) :
  Residuals Between Data and Model Prediction
- [`datapointL2()`](https://dmodverse.github.io/dMod2/reference/datapointL2.md)
  : L2 Objective of a Validation Data Point
- [`constraintL2()`](https://dmodverse.github.io/dMod2/reference/constraintL2.md)
  : Soft L2 Constraint on Parameters
- [`constraintStability()`](https://dmodverse.github.io/dMod2/reference/constraintStability.md)
  : Stability of the steady state of Pimpl() as a constraint
- [`constraintL1()`](https://dmodverse.github.io/dMod2/reference/constraintL1.md)
  : Soft L1 Constraint on Parameters
- [`constraintCauchy()`](https://dmodverse.github.io/dMod2/reference/constraintCauchy.md)
  : Soft Cauchy Constraint on Parameters
- [`constraintChisq()`](https://dmodverse.github.io/dMod2/reference/constraintChisq.md)
  : Soft Chi-Squared Constraint on Positive Parameters
- [`constraintExponential()`](https://dmodverse.github.io/dMod2/reference/constraintExponential.md)
  : Soft Exponential Constraint on Positive Parameters
- [`constraintGamma()`](https://dmodverse.github.io/dMod2/reference/constraintGamma.md)
  : Soft Gamma Constraint on Positive Parameters
- [`constraintRayleigh()`](https://dmodverse.github.io/dMod2/reference/constraintRayleigh.md)
  : Soft Rayleigh Constraint on Positive Parameters
- [`objlist()`](https://dmodverse.github.io/dMod2/reference/objlist.md)
  : Objective List
- [`as.objlist()`](https://dmodverse.github.io/dMod2/reference/as.objlist.md)
  : Zero Objective List
- [`` `+`( ``*`<objfn>`*`)`](https://dmodverse.github.io/dMod2/reference/plus-.objfn.md)
  : Direct Sum of Objective Functions
- [`` `+`( ``*`<objlist>`*`)`](https://dmodverse.github.io/dMod2/reference/plus-.objlist.md)
  : Add Two Objective Lists
- [`` `%.*%` ``](https://dmodverse.github.io/dMod2/reference/grapes-.-times-grapes.md)
  : Multiplication of Objective Functions with Scalars
- [`objframe()`](https://dmodverse.github.io/dMod2/reference/objframe.md)
  : Objective Frame
- [`vcov(`*`<objlist>`*`)`](https://dmodverse.github.io/dMod2/reference/vcov.objlist.md)
  [`vcov(`*`<trustfit>`*`)`](https://dmodverse.github.io/dMod2/reference/vcov.objlist.md)
  : Variance-Covariance Matrix of a Fit
- [`reml()`](https://dmodverse.github.io/dMod2/reference/reml.md) :
  Estimate the Error Model by Restricted Maximum Likelihood
- [`remlLeverage()`](https://dmodverse.github.io/dMod2/reference/remlLeverage.md)
  : Leverage of the Data Points in the Mean Model
- [`fitErrorModel()`](https://dmodverse.github.io/dMod2/reference/fitErrorModel.md)
  : Fit an Error Model to Replicate Data

## Optimization and profiles

- [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md)
  [`trustL1()`](https://dmodverse.github.io/dMod2/reference/trust.md) :
  Non-Linear Optimization by a Trust-Region Method
- [`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md)
  : Non-Linear Optimization, Multi-Start
- [`msParframe()`](https://dmodverse.github.io/dMod2/reference/msParframe.md)
  : Reproducible Random Starting Points
- [`as.parlist()`](https://dmodverse.github.io/dMod2/reference/parlist.md)
  [`print(`*`<parlist>`*`)`](https://dmodverse.github.io/dMod2/reference/parlist.md)
  [`summary(`*`<parlist>`*`)`](https://dmodverse.github.io/dMod2/reference/parlist.md)
  [`c(`*`<parlist>`*`)`](https://dmodverse.github.io/dMod2/reference/parlist.md)
  [`parlist()`](https://dmodverse.github.io/dMod2/reference/parlist.md)
  : Parameter List
- [`load.parlist()`](https://dmodverse.github.io/dMod2/reference/load.parlist.md)
  : Load Fits Written by mstrust
- [`plot(`*`<parlist>`*`)`](https://dmodverse.github.io/dMod2/reference/plot.parlist.md)
  : Plot a Parameter List
- [`is.parframe()`](https://dmodverse.github.io/dMod2/reference/parframe.md)
  [`` `[`( ``*`<parframe>`*`)`](https://dmodverse.github.io/dMod2/reference/parframe.md)
  [`subset(`*`<parframe>`*`)`](https://dmodverse.github.io/dMod2/reference/parframe.md)
  [`rbind(`*`<parframe>`*`)`](https://dmodverse.github.io/dMod2/reference/parframe.md)
  [`parframe()`](https://dmodverse.github.io/dMod2/reference/parframe.md)
  : Generate a Parameter Frame
- [`as.parframe()`](https://dmodverse.github.io/dMod2/reference/as.parframe.md)
  : Coerce to a Parameter Frame
- [`as.parvec(`*`<parframe>`*`)`](https://dmodverse.github.io/dMod2/reference/as.parvec.parframe.md)
  [`as.parvec(`*`<parlist>`*`)`](https://dmodverse.github.io/dMod2/reference/as.parvec.parframe.md)
  : Select a Parameter Vector from a Parameter Frame
- [`unique(`*`<parframe>`*`)`](https://dmodverse.github.io/dMod2/reference/unique.parframe.md)
  : Rows of a Parameter Frame with Unique Values
- [`as.parvec()`](https://dmodverse.github.io/dMod2/reference/parvec.md)
  [`parvec()`](https://dmodverse.github.io/dMod2/reference/parvec.md) :
  Parameter Vector
- [`print(`*`<parvec>`*`)`](https://dmodverse.github.io/dMod2/reference/print.parvec.md)
  : Print Parameter Vectors
- [`c(`*`<parvec>`*`)`](https://dmodverse.github.io/dMod2/reference/c.parvec.md)
  : Concatenate Parameter Vectors
- [`` `[`( ``*`<parvec>`*`)`](https://dmodverse.github.io/dMod2/reference/sub-.parvec.md)
  : Subset a Parameter Vector
- [`profile(`*`<objfn>`*`)`](https://dmodverse.github.io/dMod2/reference/profile.md)
  [`profile(`*`<function>`*`)`](https://dmodverse.github.io/dMod2/reference/profile.md)
  : Profile Likelihood
- [`profileThreshold()`](https://dmodverse.github.io/dMod2/reference/profileThreshold.md)
  : Threshold for a Profile Likelihood
- [`confint(`*`<parframe>`*`)`](https://dmodverse.github.io/dMod2/reference/confint.parframe.md)
  : Confidence Intervals from Profiles
- [`plotProfile()`](https://dmodverse.github.io/dMod2/reference/plotProfile.md)
  : Profile Likelihood Plot
- [`plotPaths()`](https://dmodverse.github.io/dMod2/reference/plotPaths.md)
  : Profile Likelihood: Plot of the Parameter Paths
- [`plotPathsMulti()`](https://dmodverse.github.io/dMod2/reference/plotPathsMulti.md)
  : Profile Likelihood: Plot All Parameter Paths of One Profile Together
- [`plotProfilesAndPaths()`](https://dmodverse.github.io/dMod2/reference/plotProfilesAndPaths.md)
  : Profile Likelihood: Plot Profiles Along with Their Parameter Paths
- [`plotArray()`](https://dmodverse.github.io/dMod2/reference/plotArray.md)
  : Plot an Array of Trajectories Along the Profile of a Parameter
- [`plotValues()`](https://dmodverse.github.io/dMod2/reference/plotValues.md)
  : Plot Objective Values of a Collection of Fits
- [`plotPars()`](https://dmodverse.github.io/dMod2/reference/plotPars.md)
  : Plot Parameter Values of a Collection of Fits
- [`plotPairs()`](https://dmodverse.github.io/dMod2/reference/plotPairs.md)
  : Pair Plot of Parameter Samples
- [`plotResiduals()`](https://dmodverse.github.io/dMod2/reference/plotResiduals.md)
  : Plot Residuals of a Collection of Fits
- [`plotTrace()`](https://dmodverse.github.io/dMod2/reference/plotTrace.md)
  : Convergence Trace of an Iterative Fit

## Mixed effects

Random effects on parameter transformations, fits by marginal likelihood
and their diagnostics.

- [`randomEffects()`](https://dmodverse.github.io/dMod2/reference/randomEffects.md)
  [`print(`*`<randomeffects>`*`)`](https://dmodverse.github.io/dMod2/reference/randomEffects.md)
  : Random effects of a parameter transformation
- [`normal()`](https://dmodverse.github.io/dMod2/reference/normal.md)
  [`spikeSlab()`](https://dmodverse.github.io/dMod2/reference/normal.md)
  [`dirichletProcess()`](https://dmodverse.github.io/dMod2/reference/normal.md)
  : Distributions of random effects
- [`em()`](https://dmodverse.github.io/dMod2/reference/em.md) :
  Mixed-effects fit by marginal maximum likelihood
- [`msem()`](https://dmodverse.github.io/dMod2/reference/msem.md) :
  Multistart mixed-effects fit
- [`print(`*`<em>`*`)`](https://dmodverse.github.io/dMod2/reference/print.em.md)
  : Print a mixed-effects fit
- [`summary(`*`<em>`*`)`](https://dmodverse.github.io/dMod2/reference/summary.em.md)
  [`print(`*`<summary.em>`*`)`](https://dmodverse.github.io/dMod2/reference/summary.em.md)
  : Summarize a mixed-effects fit
- [`coef(`*`<em>`*`)`](https://dmodverse.github.io/dMod2/reference/coef.em.md)
  : Estimates of a mixed-effects fit
- [`vcov(`*`<em>`*`)`](https://dmodverse.github.io/dMod2/reference/vcov.em.md)
  : Covariance of the estimates of a mixed-effects fit
- [`confint(`*`<em>`*`)`](https://dmodverse.github.io/dMod2/reference/confint.em.md)
  : Confidence intervals of a mixed-effects fit
- [`logLik(`*`<em>`*`)`](https://dmodverse.github.io/dMod2/reference/logLik.em.md)
  : Log-likelihood of a mixed-effects fit
- [`anova(`*`<em>`*`)`](https://dmodverse.github.io/dMod2/reference/anova.em.md)
  [`anova(`*`<alignment>`*`)`](https://dmodverse.github.io/dMod2/reference/anova.em.md)
  : Compare nested mixed-effects fits
- [`predict(`*`<em>`*`)`](https://dmodverse.github.io/dMod2/reference/predict.em.md)
  : Predictions of a mixed-effects fit
- [`residuals(`*`<em>`*`)`](https://dmodverse.github.io/dMod2/reference/residuals.em.md)
  : Residuals of a mixed-effects fit
- [`simulate(`*`<em>`*`)`](https://dmodverse.github.io/dMod2/reference/simulate.em.md)
  : Simulate data from a mixed-effects fit
- [`plot(`*`<em>`*`)`](https://dmodverse.github.io/dMod2/reference/plot.em.md)
  : Observed vs predicted scatter (DV vs IPRED and DV vs PRED)
- [`plotCombined(`*`<em>`*`)`](https://dmodverse.github.io/dMod2/reference/plotCombined.em.md)
  : Data and predictions of a mixed-effects fit
- [`plotResiduals(`*`<em>`*`)`](https://dmodverse.github.io/dMod2/reference/plotResiduals.em.md)
  : Residual diagnostics of a mixed-effects fit
- [`plotRandomEffects()`](https://dmodverse.github.io/dMod2/reference/plotRandomEffects.md)
  : Random-effect distribution diagnostics
- [`plotVPC()`](https://dmodverse.github.io/dMod2/reference/plotVPC.md)
  : Visual predictive check

## Alignment

Measurements on unknown scales aligned to common values.

- [`align()`](https://dmodverse.github.io/dMod2/reference/align.md) :
  Fit an alignment model
- [`alignmentModel()`](https://dmodverse.github.io/dMod2/reference/alignmentModel.md)
  [`print(`*`<alignmentmodel>`*`)`](https://dmodverse.github.io/dMod2/reference/alignmentModel.md)
  : Alignment model of replicate measurements
- [`as.datalist(`*`<alignment>`*`)`](https://dmodverse.github.io/dMod2/reference/as.datalist.alignment.md)
  : Aligned values as dMod2 data
- [`simulate(`*`<alignment>`*`)`](https://dmodverse.github.io/dMod2/reference/simulate.alignment.md)
  : Simulate readings from an alignment
- [`plotAligned()`](https://dmodverse.github.io/dMod2/reference/plotAligned.md)
  [`plotScaled()`](https://dmodverse.github.io/dMod2/reference/plotAligned.md)
  [`plotResiduals(`*`<alignment>`*`)`](https://dmodverse.github.io/dMod2/reference/plotAligned.md)
  : Diagnostic plots of an alignment

## PEtab and SBML

- [`importPEtab()`](https://dmodverse.github.io/dMod2/reference/importPEtab.md)
  : Import a PEtab Problem
- [`exportPEtab()`](https://dmodverse.github.io/dMod2/reference/exportPEtab.md)
  : Export a dMod Problem to PEtab
- [`exportPEtabObject()`](https://dmodverse.github.io/dMod2/reference/exportPEtabObject.md)
  : Export an Imported PEtab Problem
- [`readPEtabTables()`](https://dmodverse.github.io/dMod2/reference/readPEtabTables.md)
  : Read the Tables of a PEtab Problem
- [`readPEtabYaml()`](https://dmodverse.github.io/dMod2/reference/readPEtabYaml.md)
  : Read a PEtab YAML File
- [`print(`*`<petabproblem>`*`)`](https://dmodverse.github.io/dMod2/reference/print.petabproblem.md)
  : Print a PEtab Problem
- [`importSbml()`](https://dmodverse.github.io/dMod2/reference/importSbml.md)
  : Import an SBML Model
- [`exportSbml()`](https://dmodverse.github.io/dMod2/reference/exportSbml.md)
  : Export an Equation List to SBML

## Distributed computing

- [`distributedComputing()`](https://dmodverse.github.io/dMod2/reference/distributedComputing.md)
  : Run R Code as Tasks on Other Processes, Machines or a SLURM Cluster
- [`localNode()`](https://dmodverse.github.io/dMod2/reference/localNode.md)
  [`sshNodes()`](https://dmodverse.github.io/dMod2/reference/localNode.md)
  [`slurmCluster()`](https://dmodverse.github.io/dMod2/reference/localNode.md)
  : Where distributedComputing() Runs Tasks
- [`slurmTemplate()`](https://dmodverse.github.io/dMod2/reference/slurmTemplate.md)
  : Job script templates of SLURM clusters
- [`dMod2Jobs()`](https://dmodverse.github.io/dMod2/reference/dMod2Jobs.md)
  : Jobs of distributedComputing()
- [`detectFreeCores()`](https://dmodverse.github.io/dMod2/reference/detectFreeCores.md)
  : Detect Number of Free Cores

## Plots and colors

- [`plot(`*`<datalist>`*`)`](https://dmodverse.github.io/dMod2/reference/plotData.md)
  [`plotData()`](https://dmodverse.github.io/dMod2/reference/plotData.md)
  : Plot Observed Data
- [`plotPrediction()`](https://dmodverse.github.io/dMod2/reference/plotPrediction.md)
  : Plot a List of Model Predictions
- [`plotCombined()`](https://dmodverse.github.io/dMod2/reference/plotCombined.md)
  [`plot(`*`<prdlist>`*`)`](https://dmodverse.github.io/dMod2/reference/plotCombined.md)
  [`plot(`*`<prdframe>`*`)`](https://dmodverse.github.io/dMod2/reference/plotCombined.md)
  : Plot Model Predictions and Data in One Plot
- [`theme_dMod()`](https://dmodverse.github.io/dMod2/reference/theme_dMod.md)
  : Standard Plotting Theme of dMod
- [`scale_color_dMod()`](https://dmodverse.github.io/dMod2/reference/scale_color_dMod.md)
  [`scale_fill_dMod()`](https://dmodverse.github.io/dMod2/reference/scale_color_dMod.md)
  : Discrete dMod Color Scales
- [`scale_color_dMod_c()`](https://dmodverse.github.io/dMod2/reference/scale_color_dMod_c.md)
  [`scale_fill_dMod_c()`](https://dmodverse.github.io/dMod2/reference/scale_color_dMod_c.md)
  [`scale_color_dMod_div()`](https://dmodverse.github.io/dMod2/reference/scale_color_dMod_c.md)
  [`scale_fill_dMod_div()`](https://dmodverse.github.io/dMod2/reference/scale_color_dMod_c.md)
  : Continuous dMod Color Scales
- [`dMod_colors`](https://dmodverse.github.io/dMod2/reference/dMod_colors.md)
  : Seed Colors of the dMod Palette
- [`dMod_colors_cb`](https://dmodverse.github.io/dMod2/reference/dMod_colors_cb.md)
  : Colorblind-Safe Seed Colors
- [`dMod_palette()`](https://dmodverse.github.io/dMod2/reference/dMod_palette.md)
  : Generate Distinct Colors from a dMod Palette
- [`dMod_palettes()`](https://dmodverse.github.io/dMod2/reference/dMod_palettes.md)
  : Overview of the dMod Palettes
- [`dMod_gradient`](https://dmodverse.github.io/dMod2/reference/dMod_gradient.md)
  [`dMod_gradient_cb`](https://dmodverse.github.io/dMod2/reference/dMod_gradient.md)
  [`dMod_divergent`](https://dmodverse.github.io/dMod2/reference/dMod_gradient.md)
  [`dMod_divergent_cb`](https://dmodverse.github.io/dMod2/reference/dMod_gradient.md)
  : Continuous dMod Color Ramps

## Data sets

- [`bachmann`](https://dmodverse.github.io/dMod2/reference/bachmann.md)
  : Time-Course Data for JAK2-STAT5 Signaling in CFU-E Cells
- [`badata`](https://dmodverse.github.io/dMod2/reference/badata.md) :
  Time-Course Data for the Bile Acid Demonstration Model
- [`boehm`](https://dmodverse.github.io/dMod2/reference/boehm.md) :
  Time-Course Data for STAT5 Dimerization in BaF3-EpoR Cells
- [`jakstat`](https://dmodverse.github.io/dMod2/reference/jakstat.md) :
  Time-Course Data for the JAK-STAT Signaling Pathway
- [`kok`](https://dmodverse.github.io/dMod2/reference/kok.md) :
  Time-Course Data for IFN-alpha Signaling in Huh7.5 Cells
- [`warfarin`](https://dmodverse.github.io/dMod2/reference/warfarin.md)
  : Warfarin Pharmacokinetics and Pharmacodynamics
- [`mavoglurant`](https://dmodverse.github.io/dMod2/reference/mavoglurant.md)
  : Mavoglurant Pharmacokinetics

## Accessors and linear algebra

- [`getConditions()`](https://dmodverse.github.io/dMod2/reference/getConditions.md)
  : Extract the Conditions of an Object
- [`getDerivs()`](https://dmodverse.github.io/dMod2/reference/getDerivs.md)
  : Extract the First Derivatives of an Object
- [`getDerivs2()`](https://dmodverse.github.io/dMod2/reference/getDerivs2.md)
  : Extract the Second Derivatives of an Object
- [`getEquations()`](https://dmodverse.github.io/dMod2/reference/getEquations.md)
  : Extract the Equations of an Object
- [`getObservables()`](https://dmodverse.github.io/dMod2/reference/getObservables.md)
  : Extract the Observables of an Object
- [`getParameters()`](https://dmodverse.github.io/dMod2/reference/getParameters.md)
  : Extract the Parameters of an Object
- [`attrs()`](https://dmodverse.github.io/dMod2/reference/attrs.md) :
  Select Attributes
