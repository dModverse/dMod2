# Simulate data from a mixed-effects fit

Draws the random effects of every group from their fitted distributions,
predicts the data of every condition, adds the noise of the error model
(normal, or Student-t for
[`normT()`](https://dmodverse.github.io/dMod2/reference/normT.md)) and
censors readings at their lower limit.

## Usage

``` r
# S3 method for class 'em'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- object:

  An `em` object, see
  [`em()`](https://dmodverse.github.io/dMod2/reference/em.md).

- nsim:

  Number of data sets.

- seed:

  Seed of the random number generator, as in
  [`stats::simulate()`](https://rdrr.io/r/stats/simulate.html).

- ...:

  Ignored.

## Value

A list of `nsim`
[datalist](https://dmodverse.github.io/dMod2/reference/datalist.md)s in
the layout of the data of the fit.

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
#> Overwriting: /tmp/Rtmp8qFUqD/theoph_obs.cpp
err <- Y(c(y = "sigma_add"), g, attachInput = FALSE, modelname = "theoph_err", compile = FALSE,
         derivMode = dm, outdir = tempdir())
#> Overwriting: /tmp/Rtmp8qFUqD/theoph_err.cpp
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

## two data sets from the fit: random effects drawn per subject, noise from the error model
sims <- simulate(fit, nsim = 2, seed = 1)
head(as.data.frame(sims[[1]]))
#>   name time      value sigma lloq condition
#> 1    y 0.00 -0.2737965    NA -Inf         1
#> 2    y 0.25  2.4925352    NA -Inf         1
#> 3    y 0.57  5.5951301    NA -Inf         1
#> 4    y 1.12  7.6380937    NA -Inf         1
#> 5    y 2.02  8.2504030    NA -Inf         1
#> 6    y 3.82  7.6247819    NA -Inf         1

## a refit to the first simulated data set, started at the estimates
refit <- em(normL2(sims[[1]], g * x * p, errmodel = err), randomEffects(p), init = fit)
rbind(fit = coef(fit), refit = coef(refit))
#>              Ka        V        Cl  sigma_add omega_Ka_Ka omega_V_V omega_Cl_Cl
#> fit   0.4638825 3.459748 1.0121568 -0.3647014   -0.455224 -1.977427   -1.336115
#> refit 0.7119933 3.484717 0.9898651 -0.4233285   -1.035738 -1.767662   -1.835440
# }
```
