# Model Predictions

Evaluates a prediction function for every row of a parameter frame.

## Usage

``` r
# S3 method for class 'prdfn'
predict(object, ..., times, pars, data = NULL)
```

## Arguments

- object:

  object of class `prdfn`.

- ...:

  further arguments going to the prediction function, e.g. `conditions`.

- times:

  numeric vector of time points.

- pars:

  a [parframe](https://dmodverse.github.io/dMod2/reference/parframe.md),
  e.g. `as.parframe(mstrust(...))` or the result of
  [profile()](https://dmodverse.github.io/dMod2/reference/profile.md).

- data:

  object of class `datalist`. Default `NULL`. The condition grid of the
  data adds covariate columns to the output.

## Value

A data frame in long format with one block per row of `pars`. The
non-parameter columns of `pars` are added with a leading dot, e.g.
`.value`. The data, as a data frame, is returned as the attribute
`"data"`.
