# Print Parameter Transformations

Print Parameter Transformations

## Usage

``` r
# S3 method for class 'parfn'
print(x, ...)
```

## Arguments

- x:

  object of class `parfn`.

- ...:

  not used.

## Value

`x`, invisibly.

## Examples

``` r
# \donttest{
## One transformation per condition, all compiled into one shared object
trafo <- list(
  ctrl  = eqnvec(k1 = "exp(log_k1)", k2 = "exp(log_k2)", A = "A0"),
  treat = eqnvec(k1 = "exp(log_k1 + log_fold)", k2 = "exp(log_k2)", A = "A0"))
p <- P(trafo, compile = TRUE, modelname = "P_example", outdir = tempdir())
#> Warning: A shared library named 'P_example' is already loaded; overwriting it is not portable. Using 'P_example_2' instead.
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC]
#> reusing 2 unchanged object(s)
print(p)
#> Parameter transformation:
#> function (..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, hessian = TRUE, 
#>     conditions = st$default_conditions, env = NULL, cores = getOption("dMod2.cores", 
#>         1L), sweep = "forward")  
#> 
#> ... conditions: ctrl, treat 
#> ... parameters: log_k1, log_k2, A0, log_fold 
getParameters(p)
#> [1] "log_k1"   "log_k2"   "A0"       "log_fold"
p(c(log_k1 = 0, log_k2 = -1, A0 = 2, log_fold = log(3)))
#> $ctrl
#> Parameter vector:
#>   k1 :       1
#>   k2 :  0.367879
#>    A :       2
#> 
#> Attributes:
#>   deriv  : 3 x 3 matrix
#>   deriv2 : <none>
#>   fixed  : <none>
#> 
#> $treat
#> Parameter vector:
#>   k1 :       3
#>   k2 :  0.367879
#>    A :       2
#> 
#> Attributes:
#>   deriv  : 3 x 4 matrix
#>   deriv2 : <none>
#>   fixed  : <none>
#> 
# }
```
