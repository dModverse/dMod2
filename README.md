dMod2 - Dynamic Modeling and Parameter Estimation in R
================

<!-- badges: start -->

[![R-CMD-check](https://github.com/dModverse/dMod2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dModverse/dMod2/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

dMod2 builds ODE models of reaction networks, generates and compiles C++
code for them, and estimates their parameters from data. The framework
follows the paradigm that derivative information should be used whenever
possible: every object it produces carries its symbolic derivatives, so
prediction, observation, parameter transformation and objective
functions all deliver gradient and Hessian along with their value.

Models are assembled from four kinds of function objects that compose
with `*` (concatenation, chain rule) and `+` (direct sum over
experimental conditions):

| Object | Constructor | Role |
|----|----|----|
| `odemodel`, `prdfn` | `odemodel()`, `Xs()` | the compiled ODE and its sensitivities |
| `obsfn` | `Y()` | observation and error functions |
| `parfn` | `P()` | parameter transformations |
| `objfn` | `normL2()`, `constraintL2()` | objective functions |

## System requirements

Model code is generated and compiled at runtime, so **C and C++
compilers** are required. On Linux they are usually present, Windows
users need [Rtools](https://cran.r-project.org/bin/windows/Rtools/). The
default backend is [cppDE](https://github.com/dModverse/cppDE), which
writes C++ with first- and second-order sensitivities;
[cOde](https://github.com/dkaschek/cOde) with deSolve is available as an
alternative through `odemodel(backend = "deSolve")`.

## Installation

``` r
remotes::install_github("dModverse/dMod2")
```

The backend [cppDE](https://github.com/dModverse/cppDE) comes along as a
`Remotes:` dependency and does not have to be installed separately.

### As an RStudio project

Working on the package itself is easier from a checkout.

1.  **File → New Project → Version Control → Git**, repository URL
    `https://github.com/dModverse/dMod2`.

2.  Install the dependencies, including the ones declared under
    `Remotes:`:

    ``` r
    remotes::install_deps(dependencies = TRUE)
    ```

3.  **Build → Install and Restart**.

`cppDE` installs the same way and has to be built first: dMod2 compiles
the generated model sources against its headers.

## Example: STAT5 dimerisation after Epo stimulation

The model of [Boehm et al. (2014)](https://doi.org/10.1021/pr5006923)
describes what happens in BaF3-EpoR cells after stimulation with
erythropoietin: STAT5A and STAT5B are phosphorylated, form the three
dimers ApA, ApB and BpB, are imported into the nucleus and return to the
cytoplasm. The data are relative quantities from mass spectrometry, and
they ship with the package.

``` r
library(dMod2)
library(ggplot2)
set.seed(5555)

# generated sources and the shared object go here, not into the working directory
outdir <- tempfile("boehm")
dir.create(outdir)
```

### Reactions

Reactions are added one at a time and carry their compartment. The Epo
stimulus decays exponentially and enters the phosphorylation rates in
closed form, so `time` appears directly in a rate expression. Nuclear
species are assigned before they are produced, so that the reactions
importing them do not decide which compartment they belong to.

``` r
epo <- "1.25e-7*exp(-Epo_degradation_BaF3*time)"

reactions <- eqnlist() |>
  assignCompartment(nucpApA = "nuc", nucpApB = "nuc", nucpBpB = "nuc", volume = "0.45") |>
  addReaction("2*STAT5A", "pApA", paste0("k_phos*", epo, "*STAT5A^2"),
              "STAT5A-STAT5A phosphorylation", compartment = "cyt") |>
  addReaction("STAT5A + STAT5B", "pApB", paste0("k_phos*", epo, "*STAT5A*STAT5B"),
              "STAT5A-STAT5B phosphorylation", compartment = "cyt") |>
  addReaction("2*STAT5B", "pBpB", paste0("k_phos*", epo, "*STAT5B^2"),
              "STAT5B-STAT5B phosphorylation", compartment = "cyt") |>
  addReaction("pApA", "nucpApA", "k_imp_homo*pApA", "pApA import", compartment = "cyt") |>
  addReaction("pApB", "nucpApB", "k_imp_hetero*pApB", "pApB import", compartment = "cyt") |>
  addReaction("pBpB", "nucpBpB", "k_imp_homo*pBpB", "pBpB import", compartment = "cyt") |>
  addReaction("nucpApA", "2*STAT5A", "k_exp_homo*nucpApA", "pApA export", compartment = "nuc") |>
  addReaction("nucpApB", "STAT5A + STAT5B", "k_exp_hetero*nucpApB", "pApB export", compartment = "nuc") |>
  addReaction("nucpBpB", "2*STAT5B", "k_exp_homo*nucpBpB", "pBpB export", compartment = "nuc") |>
  setCompartmentVolume(cyt = "1.4")
```

`as.eqnvec()` turns the table into the ODE right-hand side. The
compartment volumes appear as the ratios that convert a flux from the
frame it was written in into the frame of the species it acts on.

``` r
as.eqnvec(reactions)
```

    ## Idx   Inner <- Outer
    ##   1 nucpApA <- 1*(k_imp_homo*pApA)*(1.4/0.45)-1*(k_exp_homo*nucpApA)
    ##   2 nucpApB <- 1*(k_imp_hetero*pApB)*(1.4/0.45)-1*(k_exp_hetero*nucpApB)
    ##   3 nucpBpB <- 1*(k_imp_homo*pBpB)*(1.4/0.45)-1*(k_exp_homo*nucpBpB)
    ##   4    pApA <- 1*(k_phos*1.25e-7*exp(-Epo_degradation_BaF3*time)*STAT5A^2)-1*(k_imp_homo*pApA)
    ##   5    pApB <- 1*(k_phos*1.25e-7*exp(-Epo_degradation_BaF3*time)*STAT5A*STAT5B)-1*(k_imp_hetero*pApB)
    ##   6    pBpB <- 1*(k_phos*1.25e-7*exp(-Epo_degradation_BaF3*time)*STAT5B^2)-1*(k_imp_homo*pBpB)
    ##   7  STAT5A <- -2*(k_phos*1.25e-7*exp(-Epo_degradation_BaF3*time)*STAT5A^2)-1*(k_phos*1.25e-7*exp(-Epo_degradation_BaF3*time)*STAT5A*STAT5B)+2*(k_exp_homo*nucpApA)*(0.45/1.4)+1*(k_exp_hetero*nucpApB)*(0.45/1.4)
    ##   8  STAT5B <- -1*(k_phos*1.25e-7*exp(-Epo_degradation_BaF3*time)*STAT5A*STAT5B)-2*(k_phos*1.25e-7*exp(-Epo_degradation_BaF3*time)*STAT5B^2)+1*(k_exp_hetero*nucpApB)*(0.45/1.4)+2*(k_exp_homo*nucpBpB)*(0.45/1.4)

### Prediction, observation and error model

Only relative quantities were measured, mixed by the isotope ratio
`specC17`. The error model contributes one standard deviation per
observable, estimated alongside the dynamic parameters.

``` r
model <- odemodel(reactions, modelname = "boehm_ode", compile = FALSE, outdir = outdir)
x <- Xs(model, 
        optionsOde = list(atol = 1e-8, rtol = 1e-6, maxattemps = 100L, maxsteps = 1e6), 
        optionsSens = list(atol = 1e-8, rtol = 1e-6, maxattemps = 100L, maxsteps = 1e6))

observables <- eqnvec(
  pSTAT5A_rel = "(100*pApB + 200*pApA*specC17)/(pApB + STAT5A*specC17 + 2*pApA*specC17)",
  pSTAT5B_rel = "-(100*pApB - 200*pBpB*(specC17 - 1))/((STAT5B*(specC17 - 1) - pApB) + 2*pBpB*(specC17 - 1))",
  rSTAT5A_rel = "(100*pApB + 100*STAT5A*specC17 + 200*pApA*specC17)/(2*pApB + STAT5A*specC17 + 2*pApA*specC17 - STAT5B*(specC17 - 1) - 2*pBpB*(specC17 - 1))"
)

errors <- eqnvec(
  pSTAT5A_rel = "sd_pSTAT5A_rel",
  pSTAT5B_rel = "sd_pSTAT5B_rel",
  rSTAT5A_rel = "sd_rSTAT5A_rel"
)

g <- Y(observables, x, modelname = "boehm_obs", attach.input = FALSE,
       compile = FALSE, outdir = outdir)
e <- Y(errors, g, modelname = "boehm_err", attach.input = FALSE,
       compile = FALSE, outdir = outdir)
```

### Parameter transformation

The transformation fixes the initial values, splits the total STAT5 pool
by the measured ratio, and puts every estimated parameter on a log10
scale. What comes out is the set of nine parameters the fit sees.

``` r
innerpars <- getParameters(model, g, e)
estimated <- c("Epo_degradation_BaF3", "k_exp_hetero", "k_exp_homo", "k_imp_hetero",
               "k_imp_homo", "k_phos", "sd_pSTAT5A_rel", "sd_pSTAT5B_rel", "sd_rSTAT5A_rel")

p <- eqnvec() |>
  define("x~x", x = innerpars) |>
  define("x~0", x = c("pApA", "pApB", "pBpB", "nucpApA", "nucpApB", "nucpBpB")) |>
  insert("STAT5A ~ 207.6*ratio") |>
  insert("STAT5B ~ 207.6 - 207.6*ratio") |>
  insert("ratio ~ 0.693") |>
  insert("specC17 ~ 0.107") |>
  insert("x~exp10(x)", x = estimated) |>
  P(modelname = "boehm_trafo", condition = "Boehm2014",
    compile = FALSE, outdir = outdir)

getParameters(p)
```

    ## [1] "k_phos"               "Epo_degradation_BaF3" "k_exp_homo"          
    ## [4] "k_exp_hetero"         "k_imp_homo"           "k_imp_hetero"        
    ## [7] "sd_pSTAT5A_rel"       "sd_pSTAT5B_rel"       "sd_rSTAT5A_rel"

### Compilation

Every object was built with `compile = FALSE`, so one call compiles them
together into a single shared object. Compile flags are per file: the
ODE and its sensitivities are built with OpenMP, the observation and
transformation sources without, and the report says so.

``` r
compile(g, x, p, e, output = "boehm", cores = 4)
```

    ## using C++ compiler: g++
    ##   every source : -O2 -DNDEBUG -w -fPIC -DKLU -I/home/simon/.cache/R/cppDE/deps-sundials-7.4.0-suitesparse-7.10.0/include/suitesparse
    ##   2 of 5 also  : -DKLUBTF=0 -DKLUAMD=0 -fopenmp

### Data

``` r
data(boehm)
head(boehm)
```

    ##   time        name     value sigma condition
    ## 1  0.0 pSTAT5A_rel  7.901073    NA Boehm2014
    ## 2  2.5 pSTAT5A_rel 66.363494    NA Boehm2014
    ## 3  5.0 pSTAT5A_rel 81.171324    NA Boehm2014
    ## 4 10.0 pSTAT5A_rel 94.730308    NA Boehm2014
    ## 5 15.0 pSTAT5A_rel 95.116483    NA Boehm2014
    ## 6 20.0 pSTAT5A_rel 91.441717    NA Boehm2014

### Objective and multi-start fit

`normL2()` combines data, the composed prediction `g*x*p` and the error
model into an objective returning value, gradient and Hessian. Two
directions of this model are practically non-identifiable, so a weak
prior on the dynamic parameters keeps them at finite values. It
deliberately stays off the `sd_*`: `reml()` below solves the
stationarity condition of the data term for those and would not see a
prior on them.

``` r
pouter <- structure(rep(-1, length(getParameters(p))), names = getParameters(p))
dyn <- setdiff(names(pouter), grep("^sd_", names(pouter), value = TRUE))
data <- as.datalist(boehm)
obj <- normL2(data, g*x*p, e) + constraintL2(pouter[dyn], sigma = 4)
```

A single trust-region run finds a local optimum. Multi-start search
shows how the objective is shaped: fifty fits from randomly drawn
starting points, sorted by their value, give the usual staircase with a
plateau of runs that reached the best optimum.

``` r
fits <- mstrust(obj, center = pouter, sd = 3, fits = 100, cores = 8,
                rinit = 0.1, rmax = 10, iterlim = 5000)
outframe <- as.parframe(fits)

plotValues(outframe, tol = 0.1, value < 1e4)
```

<img src="man/figures/README-mstrust-1.svg" alt="" width="100%" />

``` r
best <- as.parvec(outframe)
```

### Error model by restricted maximum likelihood

Estimating the standard deviations by maximum likelihood leaves them too
small: the residuals have already spent parameters on the mean. `reml()`
corrects for that by charging every data point its own leverage, which
is what the restricted likelihood asks for, rather than giving each the
same share of the parameter budget.

The restricted likelihood is the likelihood of error contrasts,
combinations of the data that carry no information about the mean
parameters and therefore leave sigma unbiased ([Patterson and Thompson
(1971)](https://doi.org/10.1093/biomet/58.3.545)).

``` r
remlfit <- reml(obj, best)

rbind(ML = best[remlfit$errpars], REML = remlfit$argument[remlfit$errpars])
```

    ##      sd_pSTAT5A_rel sd_pSTAT5B_rel sd_rSTAT5A_rel
    ## ML        0.5880246      0.8200799      0.4987346
    ## REML      0.6254604      0.8442266      0.5248558

``` r
remlfit$dof     # n_g minus the leverage each observable spends
```

    ## pSTAT5A_rel pSTAT5B_rel rSTAT5A_rel 
    ##    13.67822    14.22151    14.10027

``` r
remlfit$rank    # effective number of mean parameters
```

    ## [1] 6

`remlLeverage()` shows where the budget goes. The rank is what the model
actually spends on the mean, not the nominal parameter count.

``` r
lev <- remlLeverage(obj, best)
tapply(lev$leverage, lev$name, sum)
```

    ## pSTAT5A_rel pSTAT5B_rel rSTAT5A_rel 
    ##    2.368112    1.743765    1.888123

### Fit and uncertainty band

The band is one standard deviation and therefore belongs to the REML
estimate. The dynamic parameters move a little as well: the three
correction factors differ, which reweights the observables against each
other.

``` r
pars <- remlfit$argument
times <- seq(0, 240, length.out = 200)

prdout <- (g*x*p)(times, pars, deriv = FALSE) |> as.data.frame(errfn = e)

plot((g*x*p)(times, pars, deriv = FALSE), data) +
  geom_ribbon(data = prdout, aes(x = time, ymin = value - sigma, ymax = value + sigma),
              linetype = "dashed", alpha = 0.3)
```

<img src="man/figures/README-bestfit-1.svg" alt="" width="100%" />

### Profile likelihood

Profiles re-optimise every other parameter while one is held away from
its optimum. They are the identifiability statement the Hessian at the
optimum cannot give ([Raue et
al. (2009)](https://doi.org/10.1093/bioinformatics/btp358)).

Profiling a nuisance parameter out leaves the profile biased as a
likelihood. The general correction is Barndorff-Nielsen’s adjusted
profile likelihood ([1983](https://doi.org/10.1093/biomet/70.2.343),
[1994](https://doi.org/10.1111/j.2517-6161.1994.tb01965.x)), with [Cox
and Reid (1987)](https://doi.org/10.1111/j.2517-6161.1987.tb01422.x) as
the tractable approximation. For Gaussian errors with the error model
profiled out it reduces to the restricted likelihood above, which is
what `reml()` computes. dMod2 implements that case only, and the
effective rank it reports is the `p` the threshold below is calibrated
on.

Three details decide whether the result is read off or guessed.

The profile must reach the threshold the interval is taken at, so
`profileThreshold()` computes it once and `profile()`, `plotProfile()`
and `confint()` all use the same number. `method = "F"` is the
finite-sample calibration, exact for Gaussian errors with a single sigma
profiled out and the generalisation of the classical t interval.

Everything refers to the objective that was actually optimised, prior
included. `stop = "value"` and `val.column = "value"` select the total
rather than the data term, and only for the total is the fit the
minimum, so only there does the line marking the optimum sit where it
belongs.

And the profiles are anchored at that fit, because a profile is measured
against the value at its origin.

``` r
nd     <- sum(vapply(data, nrow, 0L))
thr_95 <- profileThreshold(0.95, method = "F", n = nd, p = remlfit$rank)
thr_90 <- profileThreshold(0.90, method = "F", n = nd, p = remlfit$rank)

profiles <- profile(
  obj, best, names(best), limits = c(-5, 5), cores = length(best), delta = thr_95,
  stepControl = list(stepsize = 1e-4, min = 1e-4, max = Inf,
                     atol = 1e-3, rtol = 1e-3, limit = 200, stop = "data"),
  algoControl = list(reoptimize = TRUE),
  optControl  = list(rinit = 0.1, rmax = 10, iterlim = 20))

plotProfile(profiles, mode %in% c("data", "prior"), threshold = c("90%" = thr_90, "95%" = thr_95))
```

<img src="man/figures/README-profiles-1.svg" alt="" width="100%" />

``` r
confint(profiles, level = 0.95, val.column = "data",
        method = "F", n = nd, p = remlfit$rank)
```

    ##                                      name      value       lower      upper
    ## Epo_degradation_BaF3 Epo_degradation_BaF3 -1.5676092 -1.74038103 -1.3898086
    ## k_exp_hetero                 k_exp_hetero -4.3927603        -Inf -2.9516401
    ## k_exp_homo                     k_exp_homo -2.2018715 -2.55803260 -1.8883996
    ## k_imp_hetero                 k_imp_hetero -1.7876328 -1.91114147 -1.6579472
    ## k_imp_homo                     k_imp_homo  1.4806915  0.09815481        Inf
    ## k_phos                             k_phos  4.1978110  4.10619897  4.2967967
    ## sd_pSTAT5A_rel             sd_pSTAT5A_rel  0.5880246  0.40940895  0.7947085
    ## sd_pSTAT5B_rel             sd_pSTAT5B_rel  0.8200799  0.66445676  1.0133337
    ## sd_rSTAT5A_rel             sd_rSTAT5A_rel  0.4987346  0.34842048  0.6917000

The plot separates the two terms, and that is where the two
non-identifiable directions become visible: for `k_imp_homo` and
`k_exp_hetero` the `data` curve stays level and only `prior` rises,
slowly. Nuclear import of the homodimers is already fast enough to be
rate-limited by phosphorylation, and the heterodimer barely returns to
the cytoplasm within the measured 240 minutes. Neither reaches the
threshold inside the searched range, so those intervals stay open on one
side. The prior keeps the fit at finite values without pretending to
close them.

## A larger model

[`inst/examples/example_BachmannMSB2011.R`](inst/examples/example_BachmannMSB2011.R)
builds the JAK2-STAT5 model of [Bachmann et
al. (2011)](https://doi.org/10.1038/msb.2011.50): 25 states in two
compartments, the three negative feedbacks CIS, SOCS3 and SHP1, and 113
parameters estimated from 541 measurements across thirteen experiments.

## Citation

There is no publication for dMod2 yet. The paper below describes
[dMod](https://github.com/JetiLab/dMod), its predecessor, and is the
reference for the modelling framework:

Kaschek D, Mader W, Fehling-Kaschek M, Rosenblatt M, Timmer J (2019).
Dynamic Modeling, Parameter Estimation, and Uncertainty Analysis in R.
*Journal of Statistical Software*, 88(10), 1-32.
[doi:10.18637/jss.v088.i10](https://doi.org/10.18637/jss.v088.i10)
