# Distributions of random effects

The distribution of the random effects of a term, set in an effect
formula by `gr(group, dist = )`, e.g.
`cl ~ 1 + (1 | gr(cellline, dist = spikeSlab()))`. Their parameters are
estimated by [`em()`](https://dmodverse.github.io/dMod2/reference/em.md)
unless given.

- `normal()`: multivariate normal with mean 0 and a covariance estimated
  per correlation block. The default.

- `spikeSlab()`: the deviation of a group is exactly 0 (spike) or drawn
  from \\N(0, \tau^2)\\ (slab); every group deviates with probability
  `prob`, independently. The spike-and-slab random effects of a group
  are integrated jointly, see section Sparse random effects of
  [`em()`](https://dmodverse.github.io/dMod2/reference/em.md).

- `dirichletProcess()`: per symbol, the deviations follow a Dirichlet
  process with base measure \\N(0, \tau^2)\\ and `concentration`, the
  fixed effect seated first at deviation 0: the groups and the fixed
  effect form clusters under the Ewens prior of the Chinese restaurant
  process, the cluster of the fixed effect at exactly 0 and every other
  cluster one common deviation from the base measure. The prior depends
  on the sizes of the clusters only, not on the order of the groups.
  [`em()`](https://dmodverse.github.io/dMod2/reference/em.md) sums
  exactly over all partitions of at most 14 groups and reports the most
  probable clusters; `method = "saem"` draws the clusters of any number
  of groups.

[`em()`](https://dmodverse.github.io/dMod2/reference/em.md) estimates
`prob`, `concentration` and `tau` unless given.

## Usage

``` r
normal()

spikeSlab(prob = NULL, tau = NULL)

dirichletProcess(concentration = NULL, tau = NULL)
```

## Arguments

- prob:

  Probability that a group deviates, in (0, 1); `NULL` (default)
  estimates it.

- tau:

  Standard deviation of the slab; `NULL` (default) estimates it.

- concentration:

  Concentration of the Ewens prior, positive; `NULL` (default) estimates
  it.

## Value

An object of class `redist`.

## See also

[`randomEffects()`](https://dmodverse.github.io/dMod2/reference/randomEffects.md),
[`insert()`](https://dmodverse.github.io/dMod2/reference/define.md)

## Examples

``` r
## The distributions of a random term, written into gr() of an effect formula
normal()
#> normal() 
spikeSlab()
#> spikeSlab(prob = NULL, tau = NULL) 
spikeSlab(prob = 0.2, tau = 0.5)
#> spikeSlab(prob = 0.2, tau = 0.5) 
dirichletProcess()
#> dirichletProcess(concentration = NULL, tau = NULL) 
dirichletProcess(concentration = 1)
#> dirichletProcess(concentration = 1, tau = NULL) 

## cell lines whose rates share the population value or a common deviation in clusters
cov <- data.frame(line = paste0("L", 1:4), row.names = paste0("C", 1:4))
trafo <- branch(define(NULL, "x <- x", x = c("kp", "kd")), table = cov)
trafo <- insert(trafo, "x <- exp(x)", x = c("kp", "kd"))
trafo <- insert(trafo, x ~ 1 + (1 | gr(line, dist = dirichletProcess())), x = c("kp", "kd"))
randomEffects(trafo)
#> Random effects: 2 terms over line
#>   block 1 (dirichletProcess, concentration estimated, tau estimated, line): eta_kp, eta_kd
#>   fixed effects: kp, kd 
```
