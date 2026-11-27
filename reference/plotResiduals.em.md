# Residual diagnostics of a mixed-effects fit

One row of panels per type of residual: against the population
prediction (IWRES against the individual one) and against time, both
with a loess smoother, and a QQ-plot against N(0, 1). See
[residuals()](https://dmodverse.github.io/dMod2/reference/residuals.em.md).

## Usage

``` r
# S3 method for class 'em'
plotResiduals(
  parframe,
  type = c("IWRES", "CWRES", "NPDE"),
  nsim = 1000L,
  seed = NULL,
  ...
)
```

## Arguments

- parframe:

  An `em` fit, see
  [`em()`](https://dmodverse.github.io/dMod2/reference/em.md).

- type:

  Any of `"IWRES"`, `"CWRES"` and `"NPDE"`, by default all three.

- nsim, seed:

  With `"NPDE"`: the number of simulated data sets, default 1000, and
  the seed of the random number generator, default `NULL`.

- ...:

  Ignored.

## Value

A ggplot.

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
#> Overwriting: /tmp/RtmpX1Wc0y/theoph_obs.cpp
err <- Y(c(y = "sigma_add"), g, attachInput = FALSE, modelname = "theoph_err", compile = FALSE,
         derivMode = dm, outdir = tempdir())
#> Overwriting: /tmp/RtmpX1Wc0y/theoph_err.cpp
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

## IWRES, CWRES and NPDE from 200 simulated data sets
plotResiduals(fit, nsim = 200, seed = 1)

plotResiduals(fit, type = "CWRES")

# }
```
