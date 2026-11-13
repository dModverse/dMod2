# Simulate readings from an alignment

Draws the random effects from their fitted distributions, predicts every
reading, adds the noise of the error model (normal or Student-t) and
censors at the lower limits.

## Usage

``` r
# S3 method for class 'alignment'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- object:

  An `alignment` from
  [`align()`](https://dmodverse.github.io/dMod2/reference/align.md).

- nsim:

  Number of data sets.

- seed:

  Seed of the random number generator.

- ...:

  Ignored.

## Value

A list of `nsim` data.frames, the data of the model with simulated
`value`.
