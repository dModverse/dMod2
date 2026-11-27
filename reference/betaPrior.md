# Hyperpriors of spike-and-slab random effects

Priors on the hyperparameters of
[`spikeSlab()`](https://dmodverse.github.io/dMod2/reference/normal.md),
in place of a given value or an estimate, for blocks with few groups,
where an estimate of `prob` tends to 0 or 1 and an estimate of `tau`
follows the modes.

- `betaPrior(a, b)`: `prob` ~ Beta(a, b), integrated out. The prior
  weight of a subset of `k` slab deviations among `K` is \\B(a + k, b +
  K - k) / B(a, b)\\; `betaPrior(1, 1)` makes every number of slab
  deviations equally likely.

- `halfCauchyPrior(scale)`, `halfNormalPrior(scale)`: on `tau`.

- `invGammaPrior(shape, rate)`: on `tau^2`.

[`em()`](https://dmodverse.github.io/dMod2/reference/em.md) estimates
`tau` under a scale prior as the posterior mode, or with
`integrate = TRUE` integrates the phase `"exact"` over `log tau`; see
section Sparse random effects of
[`em()`](https://dmodverse.github.io/dMod2/reference/em.md).

## Usage

``` r
betaPrior(a = 1, b = 1)

halfCauchyPrior(scale = 1, integrate = FALSE)

halfNormalPrior(scale = 1, integrate = FALSE)

invGammaPrior(shape = 1, rate = 1, integrate = FALSE)
```

## Arguments

- a, b:

  Shape parameters of the Beta prior, positive.

- scale:

  Scale of the half-Cauchy or half-normal prior, positive.

- integrate:

  Logical. `FALSE` (default) takes the posterior mode of `tau`, `TRUE`
  integrates over `log tau` in the phase `"exact"`.

- shape, rate:

  Shape and rate of the inverse-gamma prior, positive.

## Value

An object of class `hyperprior`.

## See also

[`spikeSlab()`](https://dmodverse.github.io/dMod2/reference/normal.md),
[`em()`](https://dmodverse.github.io/dMod2/reference/em.md)

## Examples

``` r
spikeSlab(prob = betaPrior(1, 1), tau = halfCauchyPrior(1), mean = -1, spike = -Inf)
#> spikeSlab(prob = betaPrior(1, 1), tau = halfCauchyPrior(1), mean = -1, spike = -Inf) 
```
