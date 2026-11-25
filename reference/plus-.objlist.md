# Add Two Objective Lists

Add Two Objective Lists

## Usage

``` r
# S3 method for class 'objlist'
e1 + e2
```

## Arguments

- e1, e2:

  [objlist](https://dmodverse.github.io/dMod2/reference/objlist.md)s, or
  `NULL`, which returns the other.

## Value

An [objlist](https://dmodverse.github.io/dMod2/reference/objlist.md).

## Details

Gradients and Hessians are matched by parameter name; a parameter
missing in one operand counts as zero there, and a `NULL` Hessian is
left out. Numeric attributes are added, an absent one counting as zero.
The `chi2` contributions are summed per `attrName` of their objective:
one `attrName` gives one attribute `chi2`, several give
`chi2_<attrName>` each.

## Examples

``` r
## Two priors on overlapping parameters
obj1 <- constraintL2(mu = c(a = 0, b = 0), sigma = 1)
obj2 <- constraintL2(mu = c(b = 1, c = 1), sigma = 1)

## Gradients and Hessians are matched by parameter name
total <- obj1(c(a = 1, b = 2)) + obj2(c(b = 2, c = 3))
total$value
#> [1] 10
total$gradient
#> a b c 
#> 2 6 4 
```
