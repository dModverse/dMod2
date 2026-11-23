# Stability of the steady state of Pimpl() as a constraint

A term of the objective that keeps the steady state of a
[`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md)
stable. For every eigenvalue \\\lambda\\ of the Jacobian of the
right-hand side at the steady state, on the manifold of the conserved
quantities, it adds the square of the residual \\\sqrt{w}\\ r^{3/2}\\
for \\r \le 1\\ and \\\sqrt{w}\\ (1 + 3 (r - 1) / 2)\\ beyond, with \\r
= \max(0, \mathrm{Re}\\\lambda - \lambda_0) / (\lambda\_{tol} -
\lambda_0)\\. With \\\lambda\_{tol}\\ the `stabilityTol` of the
`controlsPTC` of `pSS`, the ramp begins at \\\lambda_0 =
\lambda\_{tol} - \max(d, m s)\\, with \\d\\ the `width`, \\m\\ the
`margin` and \\s\\ the largest eigenvalue modulus on the manifold. The
term is zero on a steady state whose leading real part lies below
\\\lambda_0\\, reaches the weight \\w\\ for each eigenvalue, both of a
complex pair, whose real part reaches \\\lambda\_{tol}\\, and grows
quadratically beyond. Value, gradient and Gauss-Newton Hessian are
continuous. Steady states that agree in the parameters, resting states
and start values they read count once. The gradient is exact;
eigenvalues that coincide within the precision of the eigenvalues take
the mean of their derivatives, exact where they stay equal.

## Usage

``` r
constraintStability(
  pSS,
  p = NULL,
  weight = 100,
  width = 0,
  margin = 1e-06,
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

- weight:

  \\w\\, the value of the term where the real part of an eigenvalue
  reaches the stability tolerance, in units of \\-2 \log L\\. Default
  100.

- width:

  \\d\\, the least width of the ramp in the units of the rates, default
  `0`. A ramp wide against the change of the real parts over a step
  shows the boundary to the model of a trust-region step before the step
  crosses it; a resting state whose slowest mode decays faster than
  `width` is not affected.

- margin:

  \\m\\, the least width of the ramp relative to the largest eigenvalue
  modulus, the precision of the eigenvalues; default `1e-6`.

- attrName:

  Name of the attribute that holds the value of the term, default
  `"stability"`.

## Value

An objfn of the parameters of `p`, or of `pSS` without `p`, to add to
the objective.

## Details

The term is an assumption of the model, that the system rests where the
steady state is taken, not data:
[`profile()`](https://dmodverse.github.io/dMod2/reference/profile.md)
with `stop = "data"` and likelihood ratios read the data term without
it. Parameters of `p` that the evaluation does not give, such as random
effects, enter at 0.

Where `pSS` requires a stable root, an evaluation beyond the stability
tolerance fails. With `Pimpl(controlsPTC = list(stability = FALSE))` it
returns the unstable root there, and the term goes on rising.

The attribute `"eigenvalues"` of the result lists the eigenvalues above
the onset, one row each: the conditions of the steady state, real and
imaginary part, the half-life or doubling time \\\log 2 /
\|\mathrm{Re}\\\lambda\|\\, the onset, the contribution to the term, and
the states with the largest entries of the right eigenvector, scaled to
a largest entry of 1. The sum with other terms keeps the value only, so
the attribute is read from a call of the term itself.

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
stab <- constraintStability(pSS, p, weight = 10)
# nothing on a stable steady state away from the Hopf point, the ramp just below it, the weight
# per eigenvalue at it, beyond it the term on the unstable steady state
stab(c(la = 0, lb = log(1.99)))$value
#> [1] 0
stab(c(la = 0, lb = log(2 - 1e-6)))$value
#> [1] 2.499998
out <- stab(c(la = 0, lb = log(2)))
out$value
#> [1] 19.99999
attr(out, "eigenvalues")
#>   condition            re im     halfLife  onset  penalty             states
#> 1           -1.096938e-13  1 6.318928e+12 -1e-06 9.999997 y (1.00), x (0.71)
#> 2           -1.096938e-13 -1 6.318928e+12 -1e-06 9.999997 y (1.00), x (0.71)
stab(c(la = 0, lb = log(2.5)))$value
#> [1] 2.812515e+12
setwd(old)
# }
```
