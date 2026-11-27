# Print a mixed-effects fit

Print a mixed-effects fit

## Usage

``` r
# S3 method for class 'em'
print(x, ...)
```

## Arguments

- x:

  An `em` object, see
  [`em()`](https://dmodverse.github.io/dMod2/reference/em.md).

- ...:

  Ignored.

## Value

`x` invisibly.

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
init <- c(Ka = 0, V = 3, Cl = 0.7, sigma_add = log(0.5))

## FOCEI fit, the estimates and the random effects per subject
fit <- em(obj, randomEffects(p), init = init)
fit
#> EM (focei): -2 log L 359.403186, converged TRUE, 10 iterations, 12 groups
#>          Ka           V          Cl   sigma_add omega_Ka_Ka   omega_V_V 
#>   0.4638825   3.4597479   1.0121568  -0.3647014  -0.4552240  -1.9774269 
#> omega_Cl_Cl 
#>  -1.3361153 
#> Random effects:
#>  type eta1 eta2  estimate         se
#>    sd   Ka   Ka 0.6343058 0.14972302
#>    sd    V    V 0.1384250 0.04075699
#>    sd   Cl   Cl 0.2628648 0.06586477
summary(fit)
#> EM summary (focei): -2 log L 359.403186, converged TRUE, 10 iterations
#>   without the 2 pi constant of the data (NONMEM OBJ): 116.803413
#> 
#> Parameters:
#>             estimate      se rsePct
#> Ka            0.4639 0.19180 41.350
#> V             3.4600 0.04718  1.364
#> Cl            1.0120 0.08333  8.233
#> sigma_add    -0.3647 0.07268 19.930
#> omega_Ka_Ka  -0.4552 0.23600 51.850
#> omega_V_V    -1.9770 0.29440 14.890
#> omega_Cl_Cl  -1.3360 0.25060 18.750
#> Standard errors without the curvature of the prediction in the random effects within the volume term (forward over reverse; derivMode = "forward-forward" gives the full one)
#> 
#> Random effects:
#>  type eta1 eta2 estimate      se
#>    sd   Ka   Ka   0.6343 0.14970
#>    sd    V    V   0.1384 0.04076
#>    sd   Cl   Cl   0.2629 0.06586
#> 
#> Shrinkage: eta_Ka = 0.018, eta_V = 0.103, eta_Cl = 0.039 
fit$omega
#>   type eta1 eta2  estimate         se
#> 1   sd   Ka   Ka 0.6343058 0.14972302
#> 2   sd    V    V 0.1384250 0.04075699
#> 3   sd   Cl   Cl 0.2628648 0.06586477
head(fit$random)
#>   group eta parameter        mode        se
#> 1     1  Ka  eta_Ka_1  0.08621976 0.1515854
#> 2     2  Ka  eta_Ka_2  0.19818408 0.1725709
#> 3     3  Ka  eta_Ka_3  0.37186008 0.1810720
#> 4     4  Ka  eta_Ka_4 -0.28173371 0.1640986
#> 5     5  Ka  eta_Ka_5 -0.04271272 0.1461175
#> 6     6  Ka  eta_Ka_6 -0.38568935 0.2043190

## quadrature started from the FOCEI estimates, and SAEM
fitQ <- em(obj, randomEffects(p), init = fit, method = "quadrature", verbose = FALSE)
fitS <- em(obj, randomEffects(p), init = init, method = "saem", verbose = FALSE)
c(focei = fit$value, quadrature = fitQ$value, saem = fitS$value)
#>      focei quadrature       saem 
#>   359.4032   360.7089   360.7384 
# }
if (FALSE) { # \dontrun{
## Ka and Cl with correlated random effects, V with a spike-and-slab one: a second transformation
p2 <- P(insert(trafo, Ka + Cl ~ 1 + (1 | p | ID), V ~ 1 + (1 | gr(ID, dist = spikeSlab()))),
        method = "explicit", compile = FALSE, modelname = "theoph_p2", derivMode = dm,
        outdir = tempdir())
compile(p2, output = tempfile("theoph"))
fit2 <- em(normL2(data, g * x * p2, errmodel = err), randomEffects(p2), init = init)
fit2$omega
fit2$slab
AIC(fit, fit2)
} # }
```
