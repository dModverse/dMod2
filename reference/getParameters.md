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

## Examples

``` r
## Accessors read the generated objects, no compiled code is needed
g <- Y(c(Aobs = "s*A"), states = "A", parameters = "s", compile = FALSE,
       modelname = "accessors_obs", outdir = tempdir())
#> Overwriting: /tmp/RtmpEh2W6Q/accessors_obs.cpp
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
