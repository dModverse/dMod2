# Symbols in Expressions

Lists the variables that occur in character expressions, each once, in
order of first appearance. Function names, argument names and `pi`, the
number, are not variables.

## Usage

``` r
getSymbols(char, exclude = NULL)
```

## Arguments

- char:

  Character vector of expressions. Entries `"0"` are skipped.

- exclude:

  Character vector of symbols to drop. Default `NULL`.

## Value

Character vector of the symbols, `NULL` if `char` is `NULL`.

## See also

[`replaceSymbols()`](https://dmodverse.github.io/dMod2/reference/replaceSymbols.md)

## Examples

``` r
getSymbols(c("k1*A/(Km + A)", "exp(-k2*time)"))
#> [1] "k1"   "A"    "Km"   "k2"   "time"
getSymbols("k1*A + k2*B", exclude = "A")
#> [1] "k1" "k2" "B" 
```
