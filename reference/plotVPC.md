# Visual predictive check

Per observable and bin of time the 5th, 50th and 95th percentiles of the
readings against the 95 percent prediction intervals of the same
percentiles of `nsim` data sets simulated from the fit
([simulate()](https://dmodverse.github.io/dMod2/reference/simulate.em.md)).
Censored readings and simulations enter at their limit, and a percentile
at or below the limit of the bin is not drawn. Where the data have a
limit of quantification a second panel shows the fraction of censored
readings per bin against the prediction interval of that fraction.

## Usage

``` r
plotVPC(x, ...)

# S3 method for class 'em'
plotVPC(
  x,
  nsim = 200L,
  bins = 8L,
  pc = FALSE,
  logScale = FALSE,
  strata = NULL,
  seed = NULL,
  ...
)
```

## Arguments

- x:

  An `em` fit, see
  [`em()`](https://dmodverse.github.io/dMod2/reference/em.md).

- ...:

  Ignored.

- nsim:

  Number of simulated data sets.

- bins:

  Number of bins, at quantiles of the times of the readings, or the
  breaks of the bins.

- pc:

  `TRUE`: prediction corrected (Bergstrand et al. 2011), every reading
  and simulation times the median population prediction of the bin over
  the population prediction of that reading. Without `logScale`,
  readings whose population prediction is not positive are left out,
  with their simulations.

- logScale:

  With `pc`, `TRUE` for observables on a log scale: the correction adds
  the median population prediction of the bin minus the own one instead.

- strata:

  Optional column of the condition grid of the data of the fit
  ([datalist](https://dmodverse.github.io/dMod2/reference/datalist.md)),
  one panel per value.

- seed:

  Seed of the random number generator.

## Value

A ggplot, with the BLQ panel a list of two ggplots where cowplot is not
available; attribute `vpc` holds the table of the percentiles and
intervals per bin.

## References

Bergstrand, M., Hooker, A. C., Wallin, J. E. and Karlsson, M. O. (2011).
Prediction-corrected visual predictive checks for diagnosing nonlinear
mixed-effects models. *The AAPS Journal* 13(2), 143-151.

## Examples

``` r
# \donttest{
## Theophylline, 12 subjects: one-compartment model with first-order absorption in closed
## form, log parameters with a normal random intercept per subject, additive error
th <- datasets::Theoph
th$Subject <- as.character(th$Subject)
subjects <- unique(th$Subject)
data <- as.datalist(data.frame(name = "y", time = th$Time, value = th$conc, sigma = NA_real_,
                               condition = th$Subject))
## second derivatives forward over reverse, the default of em()
dm <- c("forward", "reverse", "forward-reverse")
x <- Xt()
g <- Y(c(y = "amt * Ka / (V * Ka - Cl) * (exp(-Cl / V * time) - exp(-Ka * time))"),
       states = "time", modelname = "theoph_obs", compile = FALSE, derivMode = dm,
       outdir = tempdir())
#> Overwriting: /tmp/RtmptyY1gm/theoph_obs.cpp
err <- Y(c(y = "sigma_add"), g, attachInput = FALSE, modelname = "theoph_err", compile = FALSE,
         derivMode = dm, outdir = tempdir())
#> Overwriting: /tmp/RtmptyY1gm/theoph_err.cpp
## one condition per subject with the subject dose; the effect formula replaces each log parameter
## by an intercept plus a random effect per subject
doses <- tapply(th$Dose * th$Wt, th$Subject, `[`, 1)[subjects]
trafo <- branch(define(NULL, "x <- x", x = c("Ka", "V", "Cl", "amt", "sigma_add")),
                table = data.frame(ID = subjects, dose = doses, row.names = subjects))
trafo <- insert(trafo, "amt <- dose", dose = dose)
trafo <- insert(trafo, "x <- exp(x)", x = c("Ka", "V", "Cl", "sigma_add"))
p <- P(insert(trafo, Ka + V + Cl ~ 1 + (1 | ID)), method = "explicit", compile = FALSE,
       modelname = "theoph_p", derivMode = dm, outdir = tempdir())
compile(g, err, p, output = tempfile("theoph"))
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC]
#> reusing 14 unchanged object(s)
obj <- normL2(data, g * x * p, errmodel = err)
fit <- em(obj, randomEffects(p), init = c(Ka = 0, V = 3, Cl = 0.7, sigma_add = log(0.5)))

## percentiles of the data in six bins of time against 200 simulated data sets
vpc <- plotVPC(fit, bins = 6, seed = 1)
vpc

head(attr(vpc, "vpc"))
#>   name stratum bin time percentile observed         lo        mid         hi
#> 1    y     all   1 0.00         5%    0.000 -1.4162763 -0.7843804 -0.3069035
#> 2    y     all   1 0.00        50%    0.490  0.2199461  0.7968084  1.3836547
#> 3    y     all   1 0.00        95%    4.837  3.6921855  5.5539821  7.8604911
#> 4    y     all   2 0.59         5%    2.377  1.5481615  2.8228618  4.4320750
#> 5    y     all   2 0.59        50%    6.490  4.7121152  6.3270944  7.6635283
#> 6    y     all   2 0.59        95%    8.994  7.8375444  9.5905519 11.8168105
#>   limit limited blq blqLo blqMid blqHi
#> 1  -Inf   FALSE   0     0      0     0
#> 2  -Inf   FALSE   0     0      0     0
#> 3  -Inf   FALSE   0     0      0     0
#> 4  -Inf   FALSE   0     0      0     0
#> 5  -Inf   FALSE   0     0      0     0
#> 6  -Inf   FALSE   0     0      0     0

## prediction corrected
plotVPC(fit, bins = 6, pc = TRUE, seed = 1)

# }
```
