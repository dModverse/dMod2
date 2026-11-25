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
  type = c("conditional", "modes", "inclusion"),
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
  conditional distribution of the random effects of the group; under the
  model they follow the random-effect distribution whatever the
  shrinkage. A normal random effect is drawn normal around the mode with
  the conditional standard error, against `N(0, Omega_kk)`. A
  spike-and-slab random effect is 0 with the posterior probability `p0`
  and otherwise normal, with the reported value, or where the reported
  structure has the spike the mode of the slab, as mean and the variance
  that gives the conditional standard error; drawn against the prior,
  the slab density `prob N(0, tau^2)` and a cap at the height of a bar
  of area `1 - prob` at 0. `"modes"`: the modes, or the reported values
  of spike-and-slab random effects, narrower than the distribution by
  the shrinkage. `"inclusion"`: a heatmap of the posterior probability
  `1 - p0` of a deviation per group and sparse random effect, labeled
  with the reported value where the reported structure has a deviation.

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

A ggplot, or without cowplot a list with `hist` and `qq`; the QQ plot
has the normal random effects only, a fit with sparse random effects
alone gives the histograms as one ggplot. With `type = "inclusion"` a
ggplot. With `covariates`, a ggplot for one covariate; for several a
cowplot grid, or without cowplot a list of ggplots named by covariate.
Spike-and-slab random effects that reach no data of the group are left
out; Dirichlet-process random effects, which couple the groups, appear
in the heatmap only.

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
#> Overwriting: /tmp/RtmpJDzFHt/theoph_obs.cpp
err <- Y(c(y = "sigma_add"), g, attachInput = FALSE, modelname = "theoph_err", compile = FALSE,
         derivMode = dm, outdir = tempdir())
#> Overwriting: /tmp/RtmpJDzFHt/theoph_err.cpp
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
