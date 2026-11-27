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

## Examples

``` r
## Accessors read the generated objects, no compiled code is needed
g <- Y(c(Aobs = "s*A"), states = "A", parameters = "s", compile = FALSE,
       modelname = "accessors_obs", outdir = tempdir())
#> Overwriting: /tmp/RtmpX1Wc0y/accessors_obs.cpp
p <- P(list(C1 = eqnvec(A = "A0", s = "s"), C2 = eqnvec(A = "2*A0", s = "s")),
       compile = FALSE, modelname = "accessors_p", outdir = tempdir())

getParameters(g)
#> [1] "s"
getParameters(p)
#> [1] "A0" "s" 
getConditions(p)
#> [1] "C1" "C2"
getEquations(p, conditions = "C2")
#> Idx Inner <- Outer
#>   1     A <- 2*A0
#>   2     s <- s

modelname(g)
#> [1] "accessors_obs"
modelname(g, p)
#> [1] "accessors_obs"  "accessors_p_C1" "accessors_p_C2"
modelname(p) <- "accessors_trafo"
modelname(p)
#> [1] "accessors_trafo"
```
