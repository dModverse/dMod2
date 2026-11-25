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

## Examples

``` r
## Accessors read the generated objects, no compiled code is needed
g <- Y(c(Aobs = "s*A"), states = "A", parameters = "s", compile = FALSE,
       modelname = "accessors_obs", outdir = tempdir())
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
