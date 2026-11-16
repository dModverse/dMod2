# Random-effect distribution diagnostics

Per random effect a histogram against the estimated `N(0, Omega_kk)`
density and a QQ-plot against the estimated normal, of draws from the
conditional distributions of the random effects or of their modes; with
several groupings one panel per grouping and random effect.

## Usage

``` r
plotRandomEffects(x, ...)

# S3 method for class 'em'
plotRandomEffects(
  x,
  type = c("conditional", "modes"),
  draws = 20L,
  bins = 30,
  covariates = NULL,
  ...
)
```

## Arguments

- x:

  An `em` fit, see
  [`em()`](https://dmodverse.github.io/dMod2/reference/em.md).

- ...:

  Ignored.

- type:

  `"conditional"` (default): per group `draws` values from the
  conditional distribution of the random effects of the group, normal
  around the mode with the conditional standard error; under the model
  they follow `N(0, Omega_kk)` whatever the shrinkage. `"modes"`: the
  modes, narrower than `N(0, Omega_kk)` by the shrinkage.

- draws:

  Number of draws per group with `type = "conditional"`. Default 20.

- bins:

  Number of bars of every histogram. Default 30.

- covariates:

  Columns of the condition grid of the data of the fit
  ([datalist](https://dmodverse.github.io/dMod2/reference/datalist.md)),
  or a data.frame with one row per condition (row names or a column
  `condition`): instead the modes of the random effects of every group
  against the covariates of the first condition of the group, a scatter
  with a loess trend for a numeric covariate of more than five values,
  boxes for the others.

## Value

A ggplot (or a list of ggplots if `cowplot` is unavailable).

A ggplot, or without cowplot a list with `hist` and `qq`. With
`covariates`, a ggplot for one covariate; for several a cowplot grid, or
without cowplot a list of ggplots named by covariate.

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
#> Overwriting: /tmp/RtmpriQ23f/theoph_obs.cpp
err <- Y(c(y = "sigma_add"), g, attachInput = FALSE, modelname = "theoph_err", compile = FALSE,
         derivMode = dm, outdir = tempdir())
#> Overwriting: /tmp/RtmpriQ23f/theoph_err.cpp
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

## draws from the conditional distributions, then the modes, against N(0, Omega_kk)
plotRandomEffects(fit)

plotRandomEffects(fit, type = "modes")


## the modes against the body weight of the subjects
weight <- tapply(th$Wt, th$Subject, `[`, 1)[subjects]
plotRandomEffects(fit, covariates = data.frame(weight = weight, row.names = subjects))

# }
```
