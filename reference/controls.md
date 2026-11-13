# List, Get and Set Controls

Reads and changes the settings an object of class `objfn`, `parfn`,
`prdfn` or `obsfn` was built with. `controls(x)` lists the available
controls; with `name`, the value of that control is returned or
replaced. Elements of a list or data frame value can be changed with `$`
or `[]`.

## Usage

``` r
controls(x, ...)

# S3 method for class 'objfn'
controls(x, name = NULL, condition = NULL, ...)

# S3 method for class 'fn'
controls(x, condition = NULL, name = NULL, ...)

controls(x, ...) <- value

# S3 method for class 'objfn'
controls(x, name, condition = NULL, ...) <- value

# S3 method for class 'fn'
controls(x, condition = NULL, name, ...) <- value
```

## Arguments

- x:

  object of class `objfn` or `fn`.

- ...:

  not used.

- name:

  character, the name of the control. Default `NULL` lists the controls.

- condition:

  character, the condition name. Default `NULL`, every condition.

- value:

  the new value of the control.

## Value

`controls(x)` prints the available controls and returns their names
invisibly; for an `fn`, as a list with one entry per function. With
`name`, the value of the control, or `NULL` if no function holds it. The
setter returns `x`.

## Details

A control belongs to the function built with it and is shared by every
sum or composition that contains this function. Changing it on one of
them changes it for all of them.

On an `fn`, the lookup covers the functions `x` is made of, e.g. the
factors of `g * x * p`. `condition` restricts it to the functions
defined for that condition; a function built without a condition counts
for every condition. With `condition = NULL`, the getter returns the
value of the first function that holds the control and the setter
changes every function that holds it. On an `objfn`, including sums,
scaled objectives and objectives composed with a parameter
transformation, the lookup covers the objective terms and the functions
they evaluate: the prediction and error model given to
[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md) and
the parameter transformation of an `objfn * parfn`. The solver
tolerances of an objective can so be changed without the prediction
function it was built from. A control of the objective term itself is
found first. Setting a control that no function holds is an error, as is
an unknown condition.

The second positional argument is `name` for an `objfn` and `condition`
for an `fn`; name both arguments to avoid confusion.

## See also

[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md),
[`P()`](https://dmodverse.github.io/dMod2/reference/P.md),
[`Y()`](https://dmodverse.github.io/dMod2/reference/Y.md)

## Examples

``` r
old <- options(dMod2.outdir = tempdir())

## parfn with condition
p <- P(eqnvec(x = "-a*x"), method = "implicit", condition = "C1")
controls(p)
#> C1:
#> [1] "keepRoot"    "controlsPTC"
controls(p, condition = "C1", name = "keepRoot")
#> [1] TRUE
controls(p, condition = "C1", name = "keepRoot") <- FALSE

## obsfn without condition
g <- Y(g = eqnvec(y = "s*x"), f = NULL, states = "x", parameters = "s")
controls(g)
#> :
#> [1] "attachInput"
controls(g, name = "attachInput")
#> [1] FALSE
controls(g, name = "attachInput") <- TRUE

options(old)
```
