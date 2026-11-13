# Stability of the steady state of Pimpl() as a constraint

A term of the objective that keeps the steady state of a
[`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md)
stable. For every eigenvalue \\\lambda\\ of the Jacobian of the
right-hand side at the steady state, on the manifold of the conserved
quantities, it adds \\w\\ r^2\\ with \\r = \max(0,
\mathrm{Re}\\\lambda - \lambda_0) / (\lambda\_{tol} - \lambda_0)\\:
nothing while every real part lies below the onset \\\lambda_0\\, and
the weight \\w\\ for each eigenvalue, both of a complex pair, whose real
part reaches \\\lambda\_{tol}\\, the `stabilityTol` of `pSS`. Steady
states that agree in all inputs count once. The gradient is exact, the
Hessian Gauss-Newton.

## Usage

``` r
constraintStability(
  pSS,
  p = NULL,
  onset = NULL,
  duration = NULL,
  weight = 100,
  attrName = "stability"
)
```

## Arguments

- pSS:

  A parfn of
  [`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md),
  built with `deriv2 = TRUE`.

- p:

  The parfn whose output `pSS` reads, or `NULL` where `pSS` reads the
  parameters of the objective.

- onset:

  \\\lambda_0\\, below `stabilityTol`. `NULL` (default) for
  `stabilityTol - 1 / duration`: the term begins where a perturbation of
  the steady state decays more slowly than the data last.

- duration:

  The time the data span, for the default `onset`.

- weight:

  \\w\\, the value of the term where the real part of an eigenvalue
  reaches `stabilityTol`, in units of \\-2 \log L\\. Default 100.

- attrName:

  Name of the attribute that holds the value of the term, default
  `"stability"`.

## Value

An objfn of the parameters of `p`, or of `pSS` without `p`, to add to
the objective.

## Details

The term is an assumption of the model, that the system rests where the
steady state is taken, not data:
[`profile()`](https://rdrr.io/r/stats/profile.html) with `stop = "data"`
and likelihood ratios read the data term without it. Parameters of `p`
that the evaluation does not give, such as random effects, enter at 0.

Where `pSS` requires a stable root, an evaluation beyond `stabilityTol`
fails. With `Pimpl(controlsPTC = list(stability = FALSE))` it returns
the unstable root there, and the term goes on rising.

## See also

[`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md),
[`constraintL2()`](https://dmodverse.github.io/dMod2/reference/constraintL2.md)

## Examples

``` r
# \donttest{
# The Brusselator rests at x = a, y = b / a, stable while b < 1 + a^2
f <- c(x = "a - (b + 1) * x + x^2 * y", y = "b * x - x^2 * y")
old <- setwd(tempdir())
pSS <- Pimpl(f, deriv2 = TRUE, modelname = "brusselator_ss")
p <- P(c(a = "exp(la)", b = "exp(lb)"), modelname = "brusselator_p")
compile(pSS, p, output = "brusselator")
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC]
stab <- constraintStability(pSS, p, onset = -0.5, weight = 10)
# far from the Hopf point nothing, near it the term with gradient, beyond it the term on the
# unstable steady state
stab(c(la = 0, lb = log(0.5)))$value
#> [1] 0
stab(c(la = 0, lb = log(1.8)))
#> value
#> ==================
#> 12.8
#> gradient[1:2] (full length = 2)
#> ==================
#>    la    lb 
#> -64.0  57.6 
#> 
#> hessian[1:2,1:2]
#> ==================
#>      la     lb
#> la  160 -144.0
#> lb -144  129.6
#> 
#> 
#> attributes
#> ==================
#>  stability  12.8
stab(c(la = 0, lb = log(2.5)))$value
#> [1] 45
setwd(old)
# }
```
