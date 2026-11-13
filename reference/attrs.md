# Select Attributes

Keeps or removes attributes of an object.

## Usage

``` r
attrs(x, which = NULL, keep = TRUE)
```

## Arguments

- x:

  An object.

- which:

  Character vector of attribute names. `NULL` (default) stands for
  `"class"`, `"dim"`, `"dimnames"`, `"names"`, `"col.names"` and
  `"row.names"`.

- keep:

  Logical. `TRUE` (default) keeps the attributes in `which` and removes
  all others; `FALSE` removes those in `which`.

## Value

`x` with the selected attributes.

## Examples

``` r
x <- structure(1:3, names = c("a", "b", "c"), note = "temporary")
attrs(x)
#> a b c 
#> 1 2 3 
attrs(x, "note", keep = FALSE)
#> a b c 
#> 1 2 3 
  
  
```
