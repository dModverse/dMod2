# Coerce to an Equation Vector

Coerce to an Equation Vector

## Usage

``` r
as.eqnvec(x, ...)

# S3 method for class 'character'
as.eqnvec(x = NULL, names = NULL, ...)

# S3 method for class 'eqnlist'
as.eqnvec(x, ...)
```

## Arguments

- x:

  object of class `character` or `eqnlist`.

- ...:

  for the `eqnlist` method, arguments going to
  [`getFluxes()`](https://dmodverse.github.io/dMod2/reference/getFluxes.md),
  e.g. `type`; otherwise not used.

- names:

  character, the left-hand sides of the equations. Default `NULL` uses
  `names(x)`.

## Value

Object of class
[eqnvec](https://dmodverse.github.io/dMod2/reference/eqnvec.md).

## Details

For an `eqnlist`, the right-hand sides of the ODE are assembled from
[`getFluxes()`](https://dmodverse.github.io/dMod2/reference/getFluxes.md).

## Examples

``` r
# Generate eqnlist from the constructor
S <- matrix(c(-1, 1, 1, -1),
            nrow = 2, ncol = 2,
            dimnames = list(NULL, c("A", "B")))

rates <- c("k1*A", "k2*B")
description <- c("forward", "backward")

f <- eqnlist(smatrix = S, rates = rates, description = description)
print(f)
#> Conserved quantities:
#>   total_1 = 1*A+1*B
#> 
#>   Check Educt -> Product Rate Description
#> 1           A ->       B k1*A     forward
#> 2           B ->       A k2*B    backward

# Convert to data.frame
fdata <- as.data.frame(f)
print(fdata)
#>   Description Rate  A  B
#> 1     forward k1*A -1  1
#> 2    backward k2*B  1 -1

# Flat `volumes`: each distinct expression becomes a compartment c1, c2, ...
f <- as.eqnlist(fdata, volumes = c(A = "Vcyt", B = "Vnuc"))
print(f)
#> Conserved quantities:
#>   total_1 = 1*Vcyt*A+1*Vnuc*B
#> 
#> Compartments:
#>   c1 (V=Vcyt)
#>   c2 (V=Vnuc)
#> States by compartment:
#>   c1: A
#>   c2: B
#> 
#>   Check Educt -> Product Rate Description
#> 1           A ->       B k1*A     forward
#> 2           B ->       A k2*B    backward
print(as.eqnvec(f))
#> Idx Inner <- Outer
#>   1     A <- -1*(k1*A)+1*(k2*B)*(Vnuc/Vcyt)
#>   2     B <- 1*(k1*A)*(Vcyt/Vnuc)-1*(k2*B)
print(as.eqnvec(f, type = "amount"))
#> Idx Inner <- Outer
#>   1     A <- -1*(k1*A)*(Vcyt)+1*(k2*B)*(Vnuc)
#>   2     B <- 1*(k1*A)*(Vcyt)-1*(k2*B)*(Vnuc)

# Named compartments with explicit state assignment
f2 <- eqnlist(
  smatrix = S, rates = rates, description = description,
  compartments = list(cyt = "Vcyt", nuc = "Vnuc"),
  compartmentOf = c(A = "cyt", B = "nuc")
)
print(f2)
#> Conserved quantities:
#>   total_1 = 1*Vcyt*A+1*Vnuc*B
#> 
#> Compartments:
#>   cyt (V=Vcyt)
#>   nuc (V=Vnuc)
#> States by compartment:
#>   cyt: A
#>   nuc: B
#> 
#>   Check Educt -> Product Rate Description
#> 1           A ->       B k1*A     forward
#> 2           B ->       A k2*B    backward
print(as.eqnvec(f2))
#> Idx Inner <- Outer
#>   1     A <- -1*(k1*A)+1*(k2*B)*(Vnuc/Vcyt)
#>   2     B <- 1*(k1*A)*(Vcyt/Vnuc)-1*(k2*B)
```
