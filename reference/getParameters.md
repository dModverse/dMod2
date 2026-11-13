# Extract the Parameters of an Object

Extract the Parameters of an Object

## Usage

``` r
getParameters(x, ..., conditions = NULL)

# S3 method for class 'odemodel'
getParameters(x, ..., conditions = NULL)

# S3 method for class 'fn'
getParameters(x, ..., conditions = NULL)

# S3 method for class 'parvec'
getParameters(x, ..., conditions = NULL)

# S3 method for class 'prdframe'
getParameters(x, ..., conditions = NULL)

# S3 method for class 'prdlist'
getParameters(x, ..., conditions = NULL)

# S3 method for class 'eqnlist'
getParameters(x, ..., conditions = NULL)

# S3 method for class 'eventlist'
getParameters(x, ..., conditions = NULL)

# S3 method for class 'eqnvec'
getParameters(x, ..., conditions = NULL)
```

## Arguments

- x:

  object from which the parameters are extracted.

- ...:

  further objects; their parameters are added to the result.

- conditions:

  character, restrict the result to these conditions. Default `NULL`.
  Used by the methods for `fn` and `prdlist`.

## Value

Character vector of parameter names. For an `odemodel`, the states and
parameters; for a `prdlist`, a list with one entry per condition.
