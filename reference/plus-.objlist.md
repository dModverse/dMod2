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
