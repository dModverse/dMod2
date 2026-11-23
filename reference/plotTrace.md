# Convergence Trace of an Iterative Fit

Generic for plotting the quantities an iterative fitting method records,
one panel each, with a method for
[`em()`](https://dmodverse.github.io/dMod2/reference/em.md) fits.

## Usage

``` r
plotTrace(x, ...)

# S3 method for class 'em'
plotTrace(x, ...)
```

## Arguments

- x:

  An `em` fit, see
  [`em()`](https://dmodverse.github.io/dMod2/reference/em.md), with a
  stage trace: every method except `"focei"` with normal random effects
  only. The hyperparameters of spike-and-slab and Dirichlet-process
  random effects that the trace records (`prob_<grouping>`,
  `tau_<grouping>`, `concentration_<grouping>`) get one panel each.

- ...:

  Ignored.

## Value

A `ggplot` object, as returned by the method.

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
#> Overwriting: /tmp/RtmpnLJd0N/theoph_obs.cpp
err <- Y(c(y = "sigma_add"), g, attachInput = FALSE, modelname = "theoph_err", compile = FALSE,
         derivMode = dm, outdir = tempdir())
#> Overwriting: /tmp/RtmpnLJd0N/theoph_err.cpp
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
## SAEM: the changes of the parameters per iteration, burn-in and averaging
fit <- em(obj, randomEffects(p), init = c(Ka = 0, V = 3, Cl = 0.7, sigma_add = log(0.5)),
          method = "saem", verbose = FALSE)
plotTrace(fit)

# }
```
