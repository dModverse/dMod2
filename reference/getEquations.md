# Extract the Equations of an Object

Extract the Equations of an Object

## Usage

``` r
getEquations(x, conditions = NULL)

# S3 method for class 'odemodel'
getEquations(x, conditions = NULL)

# S3 method for class 'fn'
getEquations(x, conditions = NULL)
```

## Arguments

- x:

  object of class `odemodel` or `fn`.

- conditions:

  character or numeric, restrict the result to these conditions. Default
  `NULL`.

## Value

A list of `eqnvec` objects, one per condition. If `conditions` has
length one, the `eqnvec` itself. For an `odemodel`, the `eqnvec` of the
model.
