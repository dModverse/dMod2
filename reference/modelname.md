# Get and Set Modelname

The model name is the base name of the generated sources and shared
objects of a function such as a prediction function, parameter
transformation or objective function.

## Usage

``` r
modelname(x = NULL, ..., conditions = NULL)

# S3 method for class '`NULL`'
modelname(x = NULL, ..., conditions = NULL)

# S3 method for class 'character'
modelname(x = NULL, ..., conditions = NULL)

# S3 method for class 'objfn'
modelname(x = NULL, ..., conditions = NULL)

# S3 method for class 'fn'
modelname(x = NULL, ..., conditions = NULL)

modelname(x, ...) <- value

# S3 method for class 'fn'
modelname(x, conditions = NULL, ...) <- value

# S3 method for class 'objfn'
modelname(x, conditions = NULL, ...) <- value
```

## Arguments

- x:

  object of class `prdfn`, `parfn`, `obsfn` or `objfn`, or a character
  naming such an object in the global environment.

- ...:

  further objects; their model names are added to the result.

- conditions:

  character, restrict the result to these conditions. Default `NULL`.

- value:

  character, the new model name. Files are not renamed.

## Value

Character vector of model names. The setter returns `x`.

## See also

[`compile()`](https://dmodverse.github.io/dMod2/reference/compile.md),
[`loadDLL()`](https://dmodverse.github.io/dMod2/reference/loadDLL.md)
