# Generate a Parameter Frame

A parameter frame is a data frame with one parameter set per row. The
columns fall into three groups: meta-information (e.g. index, value,
converged), attributes of an objective function (e.g. data and prior
contribution) and the parameters.

## Usage

``` r
is.parframe(x)

# S3 method for class 'parframe'
x[i, j, drop = FALSE]

# S3 method for class 'parframe'
subset(x, ...)

# S3 method for class 'parframe'
rbind(..., deparse.level = 1)

parframe(
  x = NULL,
  parameters = setdiff(colnames(x), c(metanames, objAttributes)),
  metanames = NULL,
  objAttributes = NULL
)
```

## Arguments

- x:

  data frame. Default `NULL` gives an empty parameter frame.

- i:

  row index in any format; missing keeps all rows. Used alone, as in
  `x[j]`, it selects columns, as for a data frame.

- j:

  column index in any format; missing keeps all columns.

- drop:

  logical. Default `FALSE`. `TRUE` returns a plain data frame or vector
  instead of a `parframe`.

- ...:

  for [`subset()`](https://rdrr.io/r/base/subset.html), logical
  expressions in the columns of `x`; a row is kept where all of them are
  `TRUE`. For [`rbind()`](https://rdrr.io/r/base/cbind.html), parameter
  frames with the same columns.

- deparse.level:

  not used.

- parameters:

  character, the names of the parameter columns. Default all columns
  that are neither in `metanames` nor in `objAttributes`.

- metanames:

  character, the names of the meta-information columns. Default `NULL`.

- objAttributes:

  character, the names of the objective function attributes. Default
  `NULL`.

## Value

Object of class `parframe`, a data frame with the attributes
`parameters`, `metanames` and `obj.attributes`. `is.parframe()` returns
a logical.

## Details

[`rbind()`](https://rdrr.io/r/base/cbind.html) keeps the attributes of
the first frame. Profiles of one parameter computed in parts, such as
the left and right side, are joined into one: the repeated optimum is
dropped and the points are sorted by `constraint`.

Parameter frames can be subsetted either by `[ , ]` or by `subset`. If
`[ , index]` is used, the names of the removed columns will also be
removed from the corresponding attributes, i.e. metanames,
obj.attributes and parameters.

## See also

[`as.parframe()`](https://dmodverse.github.io/dMod2/reference/as.parframe.md),
[profile()](https://dmodverse.github.io/dMod2/reference/profile.md),
[`mstrust()`](https://dmodverse.github.io/dMod2/reference/mstrust.md)

## Examples

``` r
pf <- parframe(data.frame(value = c(3, 1, 2), a = c(1, 2, 3), b = c(0, 1, 0)),
               parameters = c("a", "b"), metanames = "value")
pf[pf$value < 3, ]
#>   value a b
#> 1     1 2 1
#> 2     2 3 0
as.parvec(pf)
#> Parameter vector:
#>   a :       2
#>   b :       1
#> 
#> Attributes:
#>   deriv  : <none>
#>   deriv2 : <none>
#>   fixed  : <none>
```
