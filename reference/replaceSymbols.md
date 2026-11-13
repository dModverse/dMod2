# Replace Symbols in Expressions

Substitutes symbols in character expressions. A replacement that is not
a single symbol, a non-negative number, a function call or already
bracketed is put in parentheses, so `x/a` with `a = "b+c"` becomes
`x/(b+c)`. Only variables are replaced, not function names, argument
names or list elements; the rest of the text is left as it is.

## Usage

``` r
replaceSymbols(what, by, x)
```

## Arguments

- what:

  Character vector, the symbols to replace.

- by:

  Character or numeric vector, the replacements, recycled to the length
  of `what`.

- x:

  Character vector of expressions, one per element. Names and attributes
  are kept; `NA` and `"0"` are returned unchanged.

## Value

`x` with the symbols replaced.

## Examples

``` r
replaceSymbols("a", "b+c", c(r1 = "x/a", r2 = "exp(a)*a^2"))
#>                 r1                 r2 
#>          "x/(b+c)" "exp(b+c)*(b+c)^2" 
replaceSymbols(c("a", "b"), c("b", "a"), "a - b")
#> [1] "b - a"
```
