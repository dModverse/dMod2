# Residuals of a mixed-effects fit

Residuals of a mixed-effects fit

## Usage

``` r
# S3 method for class 'em'
residuals(
  object,
  type = c("IWRES", "IRES", "CWRES", "PRES", "PWRES", "NPDE"),
  nsim = 1000L,
  seed = NULL,
  ...
)
```

## Arguments

- object:

  An `em` object, see
  [`em()`](https://dmodverse.github.io/dMod2/reference/em.md).

- type:

  `"IWRES"` (default; data minus individual prediction over the error
  model), `"IRES"` (data minus individual prediction), `"CWRES"`
  (conditional weighted residuals of the FOCE linearization at the
  modes: per group the residual of the group data against \\f(\eta^\*) -
  G \eta^\*\\, whitened by \\G \Omega G' + \Sigma\\, \\G\\ the
  sensitivities in the normal random effects; without those equal to
  IWRES), `"PRES"` or `"PWRES"` (against the population prediction),
  `"NPDE"` (normalized prediction distribution errors from `nsim`
  simulated data sets, decorrelated per group, readings at or below
  their limit imputed as in Nguyen et al. 2012, drawn given the other
  readings of their group).

- nsim, seed:

  With `"NPDE"`: number of simulated data sets and the seed of the
  random number generator.

- ...:

  Ignored.

## Value

Numeric vector, one entry per data row, in the order of the rows of
`predict(object)` with `source == "obs"`. NPDE have the attribute
`tests`, the tests of mean 0 (t-test), variance 1 (chi-square) and
normality (Shapiro-Wilk) with their p-values and the global p-value,
three times the smallest.

## References

Comets, E., Brendel, K. and Mentré, F. (2008). Computing normalized
prediction distribution errors to evaluate nonlinear mixed-effect
models. *Computer Methods and Programs in Biomedicine* 90(2), 154-166.

Nguyen, T. H. T., Comets, E. and Mentré, F. (2012). Extension of NPDE
for evaluation of nonlinear mixed effect models in presence of data
below the quantification limit with applications to HIV dynamic model.
*Journal of Pharmacokinetics and Pharmacodynamics* 39(5), 499-518.

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
#> Overwriting: /tmp/Rtmp2N4Zxw/theoph_obs.cpp
err <- Y(c(y = "sigma_add"), g, attachInput = FALSE, modelname = "theoph_err", compile = FALSE,
         derivMode = dm, outdir = tempdir())
#> Overwriting: /tmp/Rtmp2N4Zxw/theoph_err.cpp
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

## residuals in the order of the data rows of predict()
pr <- predict(fit)
pr <- pr[pr$source == "obs", ]
head(data.frame(pr[c("condition", "time", "observed", "IPRED", "PRED")],
                IWRES = residuals(fit), CWRES = residuals(fit, type = "CWRES"),
                PWRES = residuals(fit, type = "PWRES")))
#>   condition time observed    IPRED     PRED      IWRES      CWRES      PWRES
#> 1         1 0.00     0.74 0.000000 0.000000  1.0656621  1.0656621  1.0656621
#> 2         1 0.25     2.84 3.845625 3.262278 -1.4481847 -0.1771436 -0.6081158
#> 3         1 0.57     6.57 6.784953 5.829128 -0.3095509  1.0950932  1.0669177
#> 4         1 1.12    10.50 9.043486 7.864023  2.0975030  2.3313980  3.7960283
#> 5         1 2.02     9.66 9.784612 8.504583 -0.1794513 -0.6321836  1.6638977
#> 6         1 3.82     8.58 9.093280 7.620456 -0.7391657 -0.3220713  1.3818233

## NPDE from 200 simulated data sets, with the tests of mean, variance and normality
npde <- residuals(fit, type = "NPDE", nsim = 200, seed = 1)
attr(npde, "tests")
#>        test statistic           p
#> 1      mean 0.8845372 0.378026872
#> 2  variance 0.8560526 0.236255368
#> 3 normality 0.9643513 0.001549586
#> 4    global        NA 0.004648758
# }
```
