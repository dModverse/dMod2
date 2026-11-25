# Direct Sum of Functions

Adds prediction functions, parameter transformations or observation
functions of the same class.

## Usage

``` r
# S3 method for class 'fn'
e1 + e2
```

## Arguments

- e1:

  function of class `obsfn`, `prdfn` or `parfn`.

- e2:

  function of the same class as `e1`.

## Value

Object of the same class as `e1` and `e2`, defined for the union of
conditions.

## Details

Each function is defined for a set of conditions. The sum is defined for
their union. A condition present in both is taken from `e2`, with a
warning. A function defined for no specific condition
(`conditions = NULL`) cannot be added and gives an error.

## See also

[`P()`](https://dmodverse.github.io/dMod2/reference/P.md),
[`Y()`](https://dmodverse.github.io/dMod2/reference/Y.md),
[`Xs()`](https://dmodverse.github.io/dMod2/reference/Xs.md),
[\*.fn](https://dmodverse.github.io/dMod2/reference/times-.fn.md)

## Examples

``` r
## A prediction function interpolating parameter values piecewise linearly
times <- 0:5
grid <- data.frame(name = "A", time = times, row.names = paste0("p", times))
x <- Xd(grid)

## An observation function
g <- Y(c(Aobs = "s*A"), states = "A", parameters = "s",
       modelname = "prediction_obs", outdir = tempdir())
#> Overwriting: /tmp/RtmpJDzFHt/prediction_obs.cpp

## One parameter transformation per condition
innerpars <- c(getParameters(x), getParameters(g))
trafo <- setNames(innerpars, innerpars)
p <- P(list(C1 = replaceSymbols(innerpars, paste0(innerpars, "_C1"), trafo),
            C2 = replaceSymbols(innerpars, paste0(innerpars, "_C2"), trafo)),
       modelname = "prediction_p", outdir = tempdir())

# \donttest{
## Compile the observation function and the transformations together
compile(g, p, output = "prediction_ex")
#> Warning: A shared library named 'prediction_ex' is already loaded; overwriting it is not portable. Using 'prediction_ex_4' instead.
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC]
#> reusing 3 unchanged object(s)

outerpars <- getParameters(p)
pars <- setNames(seq(0.1, 1, length.out = length(outerpars)), outerpars)

## Unobserved states, then observables
plot((x * p)(times, pars))

plot((g * x * p)(times, pars))

# }
```
