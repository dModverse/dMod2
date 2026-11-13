# Extract the Conditions of an Object

Extract the Conditions of an Object

## Usage

``` r
getConditions(x, ...)

# S3 method for class 'list'
getConditions(x, ...)

# S3 method for class 'fn'
getConditions(x, ...)
```

## Arguments

- x:

  object of class `fn` or a named list such as a `datalist`.

- ...:

  not used.

## Value

Character vector of condition names, `NULL` for a function defined for
every condition.
