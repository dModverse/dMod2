# Transform an ODE to log10 Coordinates

Takes `f` as the right-hand side of \\\dot x = f(x)\\ and returns the
same system with the states in `states` in \\x\_{l10} = \log\_{10} x\\,
\$\$\dot x\_{l10} = \frac{f(x)}{\log(10)\\ x}, \qquad x =
10^{x\_{l10}}.\$\$ A state in log10 coordinates stays positive whatever
the error of the solver. Read results back with `10^x_l10`.

## Usage

``` r
log10Transform(f, suffix = "_l10", simplify = TRUE, states = names(f))
```

## Arguments

- f:

  `eqnvec` or named character vector, the right-hand sides.

- suffix:

  character appended to the name of each transformed state. Default
  `"_l10"`.

- simplify:

  logical, cancel the division symbolically with the Python module
  `sympy`. Default `TRUE`. With `FALSE`, or without `sympy`, the plain
  quotient is returned.

- states:

  character, the states to transform, all of `f` by default; the others
  stay linear and read the transformed ones as \\10^{x\_{l10}}\\.

## Value

An `eqnvec` in the order of `f`, the transformed states renamed.

## See also

[`odemodel()`](https://dmodverse.github.io/dMod2/reference/odemodel.md)

## Examples

``` r
log10Transform(eqnvec(x = "-k*x"))
#> Idx Inner <- Outer
#>   1 x_l10 <- -k/log(10)
log10Transform(eqnvec(A = "-k1*A + k2*B", B = "k1*A - k2*B"))
#> Idx Inner <- Outer
#>   1 A_l10 <- (-k1+(10^(B_l10))*k2/(10^(A_l10)))/log(10)
#>   2 B_l10 <- ((10^(A_l10))*k1/(10^(B_l10))-k2)/log(10)
log10Transform(eqnvec(A = "-k1*A + k2*B", B = "k1*A - k2*B"), states = "A")
#> Idx Inner <- Outer
#>   1 A_l10 <- (-k1+B*k2/(10^(A_l10)))/log(10)
#>   2     B <- k1*(10^(A_l10))-k2*B
```
