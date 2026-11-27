# Plot Observed Data

Plots data with error bars and marks values below the limit of
quantification (BLoQ).

## Usage

``` r
# S3 method for class 'datalist'
plotData(
  data,
  ...,
  scales = "free",
  facet = c("wrap", "grid", "wrap_plain"),
  transform = NULL
)

# S3 method for class 'datalist'
plot(
  x,
  ...,
  scales = "free",
  facet = c("wrap", "grid", "wrap_plain"),
  transform = NULL
)

plotData(data, ...)

# S3 method for class 'data.frame'
plotData(data, ...)
```

## Arguments

- data:

  A [datalist](https://dmodverse.github.io/dMod2/reference/datalist.md),
  i.e. a named list of data frames with columns `name`, `time`, `value`
  and `sigma`, or such a data frame, split into conditions by
  [`as.datalist()`](https://dmodverse.github.io/dMod2/reference/datalist.md).

- ...:

  Logical expressions to subset the plotted data. They can refer to the
  columns of the data and of the condition grid.

- scales:

  The `scales` argument of
  [`ggplot2::facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html)
  or
  [`ggplot2::facet_grid()`](https://ggplot2.tidyverse.org/reference/facet_grid.html):
  `"free"` (default), `"fixed"`, `"free_x"` or `"free_y"`.

- facet:

  `"wrap"` (default): one panel per name, color by condition. `"grid"`:
  names as rows, conditions as columns. `"wrap_plain"`: one panel per
  combination of name and condition.

- transform:

  Named list of transformations of the states, see
  [`coordTransform()`](https://dmodverse.github.io/dMod2/reference/coordTransform.md).

- x:

  object of class `datalist`.

## Value

A `ggplot` object. The attribute `"data"` holds the plotted data frame.

## Details

The plotted data frame has columns `time`, `value`, `sigma`, `name` and
`condition`.

## Functions

- `plot(datalist)`: S3 plot method for `datalist` objects.

## See also

[`plotCombined()`](https://dmodverse.github.io/dMod2/reference/plotCombined.md),
[`plotPrediction()`](https://dmodverse.github.io/dMod2/reference/plotPrediction.md)

[`plotPrediction()`](https://dmodverse.github.io/dMod2/reference/plotPrediction.md),
[`plotCombined()`](https://dmodverse.github.io/dMod2/reference/plotCombined.md)

## Examples

``` r
data <- datalist(
  C1 = data.frame(name = "A", time = 0:5, value = 0:5, sigma = 0.1),
  C2 = data.frame(name = "A", time = 0:5, value = sin(0:5), sigma = 0.1)
)
plotData(data, time < 4)

plotData(data, facet = "grid")


# \donttest{
## Observation function and one parameter transformation per condition
fn <- eqnvec(
  sine = "1 + sin(6.28*omega*time)",
  cosine = "cos(6.28*omega*time)"
)
g <- Y(fn, parameters = "omega", outdir = tempdir())
#> Overwriting: /tmp/RtmpX1Wc0y/obsfn.cpp
x <- Xt()
trafo <- lapply(1:3, function(i) eqnvec(omega = paste0("omega_", i)))
names(trafo) <- paste0("frequency_", 1:3)
p <- P(trafo, outdir = tempdir())
compile(g, p, output = "plotting_example")
#> Warning: A shared library named 'plotting_example' is already loaded; overwriting it is not portable. Using 'plotting_example_2' instead.
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC]
#> reusing 4 unchanged object(s)

## Evaluate prediction
times <- seq(0, 1, .01)
pars <- structure(seq(1, 2, length.out = 3), names = getParameters(p))
prediction <- (g*x*p)(times, pars)

## Plot prediction
plotPrediction(prediction)

plotPrediction(prediction, scales = "fixed")

plotPrediction(prediction, facet = "grid")

plotPrediction(prediction, scales = "fixed",
               transform = list(sine = "x^2", cosine = "x - 1"))


## Simulate data
set.seed(1)
dataset <- as.data.frame(prediction)
dataset <- dataset[seq(1, nrow(dataset), 5), ]
dataset$value <- dataset$value + rnorm(nrow(dataset), 0, .1)
dataset$sigma <- 0.1
data <- as.datalist(dataset, splitBy = "condition")

## Plot data, and data with prediction
plotData(data)

plotCombined(prediction, data)

plotCombined(prediction, data, time <= 0.5 & condition == "frequency_1")

plotCombined(prediction, data, time <= 0.5 & condition != "frequency_1",
             facet = "grid")

plotCombined(prediction, data, aesthetics = list(linetype = "condition"))

# }
```
