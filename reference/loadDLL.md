# Load Shared Objects for dMod Objects

Loads the shared objects of dMod functions by their model names and the
names of their generated sources, e.g. after restoring a workspace in a
new R session. Searched are the working directory,
`getOption("dMod2.outdir")` and the directories the sources were
generated in. Shared objects already loaded are skipped.

## Usage

``` r
loadDLL(...)
```

## Arguments

- ...:

  objects of class `prdfn`, `obsfn`, `parfn` or `objfn`.

## Value

Character vector of the files loaded by this call, invisibly.

## See also

[`compile()`](https://dmodverse.github.io/dMod2/reference/compile.md),
[`modelname()`](https://dmodverse.github.io/dMod2/reference/modelname.md)

## Examples

``` r
# \donttest{
g <- Y(c(y = "s * x"), parameters = "s", modelname = "loadDLL_ex",
       compile = TRUE, outdir = tempdir())
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC -fopenmp]
loadDLL(g)
# }
```
