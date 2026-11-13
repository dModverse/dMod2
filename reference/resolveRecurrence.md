# Substitute Earlier Entries into Later Ones

Replaces the name of each entry in all later entries by the entry in
parentheses, so that no entry depends on an earlier name, e.g. for
steady state expressions.

## Usage

``` r
resolveRecurrence(variables)
```

## Arguments

- variables:

  Named character vector.

## Value

Named character vector of the same length as `variables`.

## Examples

``` r
resolveRecurrence(c(A = "k1*B/k2", C = "A*k3+k4", D="A*C*k5"))
#>                                A                                C 
#>                        "k1*B/k2"                "(k1*B/k2)*k3+k4" 
#>                                D 
#> "(k1*B/k2)*((k1*B/k2)*k3+k4)*k5" 
```
