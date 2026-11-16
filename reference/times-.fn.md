# Concatenation of Functions

Concatenates observation functions, prediction functions and parameter
transformations: `(e1 * e2)(times, pars)` evaluates `e2` first and
passes the output of `e2` to `e1`.

## Usage

``` r
# S3 method for class 'fn'
e1 * e2
```

## Arguments

- e1:

  function of class `objfn`, `obsfn`, `prdfn`, `parfn` or `idfn`.

- e2:

  function of class `obsfn`, `prdfn`, `parfn` or `idfn`.

## Value

`obsfn * obsfn` and `obsfn * parfn` return an `obsfn`, `obsfn * prdfn`
and `prdfn * parfn` a `prdfn`, `parfn * parfn` a `parfn` and
`objfn * parfn` an `objfn`. `Id() * f` and `f * Id()` return `f`.

## Details

A function defined for one condition combines only with one for the same
condition or for no specific condition (`conditions = NULL`); with one
for another condition or for several conditions it is an error.

## See also

[+.fn](https://dmodverse.github.io/dMod2/reference/plus-.fn.md),
[`Id()`](https://dmodverse.github.io/dMod2/reference/Id.md)

## Examples

``` r
## A prediction function interpolating parameter values piecewise linearly
times <- 0:5
grid <- data.frame(name = "A", time = times, row.names = paste0("p", times))
x <- Xd(grid)

## An observation function
g <- Y(c(Aobs = "s*A"), states = "A", parameters = "s",
       modelname = "prediction_obs", outdir = tempdir())
#> Overwriting: /tmp/Rtmpr3VxZ8/prediction_obs.cpp

## One parameter transformation per condition
innerpars <- c(getParameters(x), getParameters(g))
trafo <- setNames(innerpars, innerpars)
p <- P(list(C1 = replaceSymbols(innerpars, paste0(innerpars, "_C1"), trafo),
            C2 = replaceSymbols(innerpars, paste0(innerpars, "_C2"), trafo)),
       modelname = "prediction_p", outdir = tempdir())

# \donttest{
## Compile the observation function and the transformations together
compile(g, p, output = "prediction_ex")
#> Warning: A shared library named 'prediction_ex' is already loaded; overwriting it is not portable. Using 'prediction_ex_6' instead.
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC]
#> reusing 3 unchanged object(s)

outerpars <- getParameters(p)
pars <- setNames(seq(0.1, 1, length.out = length(outerpars)), outerpars)

## Unobserved states, then observables
plot((x * p)(times, pars))

plot((g * x * p)(times, pars))

# }
```
