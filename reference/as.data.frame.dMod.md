# Coerce to a Data Frame

Coerce to a Data Frame

## Usage

``` r
# S3 method for class 'datalist'
as.data.frame(x, ...)

# S3 method for class 'prdlist'
as.data.frame(x, ..., data = NULL, errfn = NULL)
```

## Arguments

- x:

  any R object

- ...:

  not used.

- data:

  object of class `datalist`. Default `NULL`. The condition grid of the
  data adds covariate columns to the output.

- errfn:

  object of class `obsfn`, the error model that predicts the column
  `sigma`. Default `NULL` sets `sigma` to `NaN`.

## Value

a data frame

## Examples

``` r
## Prediction frames of a decay at two rates, without derivatives
times <- 0:4
frame1 <- prdframe(cbind(time = times, A = exp(-times)), parameters = c(k = 1))
frame2 <- prdframe(cbind(time = times, A = exp(-2 * times)), parameters = c(k = 2))
frame1
#>      time          A
#> [1,]    0 1.00000000
#> [2,]    1 0.36787944
#> [3,]    2 0.13533528
#> [4,]    3 0.04978707
#> [5,]    4 0.01831564
#> 
#> The prediction contains 1st-order derivatives: no
#> The prediction contains 2nd-order derivatives: no

## One frame per condition, in long format
pred <- prdlist(slow = frame1, fast = frame2)
as.data.frame(pred)
#>    time name        value condition sigma
#> 1     0    A 1.0000000000      slow   NaN
#> 2     1    A 0.3678794412      slow   NaN
#> 3     2    A 0.1353352832      slow   NaN
#> 4     3    A 0.0497870684      slow   NaN
#> 5     4    A 0.0183156389      slow   NaN
#> 6     0    A 1.0000000000      fast   NaN
#> 7     1    A 0.1353352832      fast   NaN
#> 8     2    A 0.0183156389      fast   NaN
#> 9     3    A 0.0024787522      fast   NaN
#> 10    4    A 0.0003354626      fast   NaN
```
