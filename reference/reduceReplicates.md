# Reduce Replicates to Mean and Standard Error

Replaces the replicates of each condition by their mean and the standard
error of the mean.

## Usage

``` r
reduceReplicates(
  data,
  select = "condition",
  datatrans = NULL,
  keep = NULL,
  weighted = FALSE
)

# S3 method for class 'data.frame'
reduceReplicates(
  data,
  select = "condition",
  datatrans = NULL,
  keep = NULL,
  weighted = FALSE
)
```

## Arguments

- data:

  A data frame with the columns `name` (observable), `time`, `value` and
  `condition`, plus any further columns that define conditions.

- select:

  Names of the columns whose values, together with `name`, `time` and
  `condition`, define a condition. Default `"condition"`.

- datatrans:

  Character, an expression in `x` applied to `value` before the
  reduction, for example `"log(x)"`. Default `NULL`.

- keep:

  Names of columns that are kept although their values differ within a
  condition; the first value is used. Default `NULL`.

- weighted:

  Logical. If `TRUE`, mean and standard error are weighted by
  `1 / sigma^2` from the column `sigma` of `data`. Default `FALSE`.

## Value

A data frame with the columns

- time:

  Time of the measurement.

- value:

  Mean of the replicates.

- sigma:

  Standard error of the mean, `NA` for a single measurement.

- n:

  Number of replicates.

- name:

  Observable.

- condition:

  The values of the `select` columns other than `name` and `time`,
  joined by `"_"`.

followed by the kept columns.

## Details

A column not in `select` is kept when it has the same value for all
replicates of every condition, otherwise it is dropped with a message.

## See also

[`fitErrorModel()`](https://dmodverse.github.io/dMod2/reference/fitErrorModel.md)

## Examples

``` r
data <- data.frame(name = "y", time = rep(c(0, 1, 2), each = 3),
                   value = c(1.0, 1.2, 0.9, 0.6, 0.5, 0.7, 0.3, 0.2, 0.25),
                   condition = "ctrl", replicate = rep(1:3, 3))
reduceReplicates(data)
#> Dropped columns: replicate
#>   time    value      sigma n name condition
#> 1    0 1.033333 0.08819171 3    y      ctrl
#> 4    1 0.600000 0.05773503 3    y      ctrl
#> 7    2 0.250000 0.02886751 3    y      ctrl
```
