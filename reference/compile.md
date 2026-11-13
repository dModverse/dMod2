# Compile Model-Related C/C++ Code

Compiles the generated sources of model objects
([parfn](https://dmodverse.github.io/dMod2/reference/parfn.md),
[obsfn](https://dmodverse.github.io/dMod2/reference/obsfn.md),
[prdfn](https://dmodverse.github.io/dMod2/reference/prdfn.md),
[`odemodel()`](https://dmodverse.github.io/dMod2/reference/odemodel.md))
into shared objects with `R CMD SHLIB` and loads them.

## Usage

``` r
compile(
  ...,
  output = NULL,
  args = NULL,
  cores = detectFreeCores(),
  chunkSize = 100,
  verbose = FALSE
)
```

## Arguments

- ...:

  one or more model objects.

- output:

  character, the name of a combined shared object. Default `NULL`. A
  bare name places it next to the generated sources, a name with a
  directory is taken as given.

- args:

  character, additional compiler and linker flags for every file.
  Default `NULL`.

- cores:

  integer, parallel compilation jobs. Default
  [`detectFreeCores()`](https://dmodverse.github.io/dMod2/reference/detectFreeCores.md).

- chunkSize:

  integer, the maximum number of object files per archiver call when the
  link command would exceed the command-line limit. Default 100.

- verbose:

  logical, print the compiler commands. Default `FALSE`.

## Value

`TRUE`, invisibly.

## Details

Compile and link flags are taken per source from the `"compileInfo"`
attribute of the objects. Objects without it are matched to sources in
the working directory by model name. Unchanged sources are not
recompiled.

Without `output`, each source becomes a shared object named after it,
i.e. after the model name. With `output`, all sources are linked into
one shared object of that name, and the model name of every object
passed as a plain variable is set to `output`, see
[`modelname()`](https://dmodverse.github.io/dMod2/reference/modelname.md).
If a shared object of that name is already loaded, a suffix `_2`, `_3`,
... is appended with a warning. With many sources, set `output`: R
limits the number of loaded shared objects.

Named arguments other than the ones below are an error.

## See also

[`odemodel()`](https://dmodverse.github.io/dMod2/reference/odemodel.md),
[`loadDLL()`](https://dmodverse.github.io/dMod2/reference/loadDLL.md),
[`modelname()`](https://dmodverse.github.io/dMod2/reference/modelname.md)

## Examples

``` r
# \donttest{
g <- Y(c(y = "s * x"), parameters = "s", modelname = "compile_ex_obs",
       outdir = tempdir())
p <- P(c(s = "exp(log_s)"), modelname = "compile_ex_par",
       outdir = tempdir())
compile(g, p, output = "compile_ex", cores = 1)
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC]
modelname(g)
#> [1] "compile_ex"
# }
```
