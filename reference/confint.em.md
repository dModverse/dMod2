# Confidence intervals of a mixed-effects fit

`method = "wald"` from the standard errors of
[vcov()](https://dmodverse.github.io/dMod2/reference/vcov.em.md);
`method = "profile"` refits with each parameter held at trial values and
finds where -2 log L has risen by the \\\chi^2_1\\ quantile of `level`
on either side, by bisection from the Wald limits. A side on which the
rise is not reached within 128 times the distance of the Wald limit from
the estimate gives `-Inf` or `Inf`.

## Usage

``` r
# S3 method for class 'em'
confint(
  object,
  parm = NULL,
  level = 0.95,
  method = c("wald", "profile"),
  control = NULL,
  ...
)
```

## Arguments

- object:

  An `em` object, see
  [`em()`](https://dmodverse.github.io/dMod2/reference/em.md).

- parm:

  Parameter names; default all estimated parameters.

- level:

  Confidence level. Default 0.95.

- method:

  `"wald"` (default) or `"profile"`.

- control:

  With `method = "profile"`: `control` of
  [`em()`](https://dmodverse.github.io/dMod2/reference/em.md) for the
  refits, by default that of the fit.

- ...:

  Ignored.

## Value

A data.frame with columns `name`, `value`, `lower`, `upper`.

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

## estimates and their standard errors
coef(fit)
#>          Ka           V          Cl   sigma_add omega_Ka_Ka   omega_V_V 
#>   0.4638825   3.4597479   1.0121568  -0.3647014  -0.4552240  -1.9774269 
#> omega_Cl_Cl 
#>  -1.3361153 
sqrt(diag(vcov(fit)))
#>          Ka           V          Cl   sigma_add omega_Ka_Ka   omega_V_V 
#>  0.19181010  0.04717589  0.08333122  0.07267625  0.23604232  0.29443384 
#> omega_Cl_Cl 
#>  0.25056516 

## marginal log-likelihood, AIC and BIC
logLik(fit)
#> 'log Lik.' -179.7016 (df=7)
AIC(fit)
#> [1] 373.4032
BIC(fit)
#> [1] 393.5828

## Wald intervals of all parameters, the profile interval of the clearance
confint(fit)
#>          name      value       lower        upper
#> 1          Ka  0.4638825  0.08794164  0.839823427
#> 2           V  3.4597479  3.36728482  3.552210899
#> 3          Cl  1.0121568  0.84883060  1.175482991
#> 4   sigma_add -0.3647014 -0.50714422 -0.222258542
#> 5 omega_Ka_Ka -0.4552240 -0.91785850  0.007410402
#> 6   omega_V_V -1.9774269 -2.55450662 -1.400347157
#> 7 omega_Cl_Cl -1.3361153 -1.82721395 -0.845016587
confint(fit, parm = "Cl", method = "profile")
#>   name    value     lower    upper
#> 1   Cl 1.012157 0.8360508 1.191971
# }
if (FALSE) { # \dontrun{
## profile intervals of all parameters
confint(fit, method = "profile")
} # }
```
