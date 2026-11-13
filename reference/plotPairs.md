# Pair Plot of Parameter Samples

Generic for a corner plot of sampled parameters: scatter plots of pairs
in the lower triangle, marginal densities on the diagonal. dMod2 defines
no method for it; packages and extensions that return samples do.

## Usage

``` r
plotPairs(x, ...)
```

## Arguments

- x:

  Object holding parameter samples.

- ...:

  Method-specific arguments.

## Value

A `ggplot` object, as returned by the method.
