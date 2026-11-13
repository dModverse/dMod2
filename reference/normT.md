# Student-t likelihood between data and model prediction

Builds the objective of a fit whose residuals follow a Student-t
distribution with `df` degrees of freedom and scale given by the error
model or the data's `sigma`: the `-2 log` likelihood of the data, summed
over the conditions. The tails are heavier than those of
[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md), so
single outliers move the fit less; `df = Inf` would be
[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md).

## Usage

``` r
normT(
  data,
  x,
  errmodel = NULL,
  df,
  times = NULL,
  attrName = "data",
  optBLOQ = c("M3", "M1")
)
```

## Arguments

- data:

  Object of class
  [datalist](https://dmodverse.github.io/dMod2/reference/datalist.md).
  Each condition of the data has to be a condition of `x`.

- x:

  Object of class
  [prdfn](https://dmodverse.github.io/dMod2/reference/prdfn.md).

- errmodel:

  Optional object of class
  [obsfn](https://dmodverse.github.io/dMod2/reference/obsfn.md). The
  error model may be defined only for a subset of conditions.

- df:

  Degrees of freedom, a positive number, fixed during the fit.

- times:

  Optional time points at which the prediction is evaluated in addition
  to the data times: a numeric vector for all conditions, which then
  share one grid, or a list named by condition, which gives each
  condition a separate grid. The prediction of a condition starts at the
  first time of the grid, where the initial values take effect, so a
  single time per condition sets the start of that condition. A model
  built with `odemodel(..., includeTimeZero = TRUE)`, the default, has 0
  in every grid and so starts at 0 at the latest. Fixed-time events
  before the start do not fire. Event times should be included if the
  model uses events.

- attrName:

  Character string. The objective value is additionally returned as an
  attribute of this name, and the sum of squares behind it under `chi2`.
  Adding objectives pools the terms sharing an `attrName`, so `chi2`
  stays one number; where two `attrName`s meet it splits into
  `chi2_<attrName>` per contribution.

- optBLOQ:

  `"M3"` (default, the probability of lying below the limit under the t
  distribution) or `"M1"` (censored rows dropped).

## Value

An objective function of class `objfn`, as
[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md)
returns it, usable in
[`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md),
[`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md),
[`profile()`](https://rdrr.io/r/stats/profile.html) and
[`em()`](https://dmodverse.github.io/dMod2/reference/em.md).

## See also

[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md)

## Examples

``` r
## A straight line with two gross outliers: the Student-t fit stays on the line
x <- Xd(data.frame(name = "y", time = 0:9, row.names = paste0("a", 0:9)), condition = "C1")
set.seed(1)
y <- 1 + 0.5 * (0:9) + rnorm(10, 0, 0.1)
y[c(3, 8)] <- y[c(3, 8)] + c(3, -4)
data <- as.datalist(data.frame(name = "y", time = 0:9, value = y, sigma = 0.1, condition = "C1"))
pars <- setNames(1 + 0.5 * (0:9), paste0("a", 0:9))
obj <- normT(data, x, df = 3)
obj(pars)$value
#> [1] 26.33673
```
