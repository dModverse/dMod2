# Compare nested mixed-effects fits

Likelihood ratio tests of fits ordered by their number of parameters,
each against the one before, with AIC and BIC.

## Usage

``` r
# S3 method for class 'em'
anova(object, ...)

# S3 method for class 'alignment'
anova(object, ...)
```

## Arguments

- object, ...:

  `em` objects of nested models on the same data, or `alignment`
  objects.

## Value

A data.frame with `model` (the position of the fit among the arguments),
`npar`, `value` (-2 log L), `AIC`, `BIC`, `chisq`, `df` and `p` per fit.

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
err <- Y(c(y = "sigma_add"), g, attachInput = FALSE, modelname = "theoph_err", compile = FALSE,
         derivMode = dm, outdir = tempdir())
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
obj <- normL2(data, g * x * p, errmodel = err)
fit <- em(obj, randomEffects(p), init = c(Ka = 0, V = 3, Cl = 0.7, sigma_add = log(0.5)))

## the nested model without a random effect on V
p0 <- P(insert(trafo, Ka + Cl ~ 1 + (1 | ID)), method = "explicit", compile = FALSE,
        modelname = "theoph_p0", derivMode = dm, outdir = tempdir())
compile(p0, output = tempfile("theoph"))
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC]
fit0 <- em(normL2(data, g * x * p0, errmodel = err), randomEffects(p0),
           init = c(Ka = 0, V = 3, Cl = 0.7, sigma_add = log(0.5)))
anova(fit0, fit)
#>   model npar    value      AIC      BIC   chisq df            p
#> 1     1    6 372.8893 384.8893 402.1861      NA NA           NA
#> 2     2    7 359.4032 373.4032 393.5828 13.4861  1 0.0002403376
# }
```
