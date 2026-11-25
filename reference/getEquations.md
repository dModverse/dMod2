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
