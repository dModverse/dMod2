# Prediction List

A prediction list holds model predictions, one
[prdframe](https://dmodverse.github.io/dMod2/reference/prdframe.md) per
condition or parameter set. Prediction functions return prediction
lists.

## Usage

``` r
as.prdlist(x, ...)

# S3 method for class 'list'
as.prdlist(x = NULL, names = NULL, ...)

prdlist(...)
```

## Arguments

- x:

  list of prediction frames. Default `NULL` gives an empty list.

- ...:

  for `prdlist()`, objects of class
  [prdframe](https://dmodverse.github.io/dMod2/reference/prdframe.md);
  unnamed entries are named `"1"`, `"2"`, ... For `as.prdlist()`, not
  used.

- names:

  character, the list names, e.g. the condition names. Default `NULL`
  keeps `names(x)`.

## Value

Object of class `prdlist`.

## See also

[`prdfn()`](https://dmodverse.github.io/dMod2/reference/prdfn.md),
[`plotCombined()`](https://dmodverse.github.io/dMod2/reference/plotCombined.md)
