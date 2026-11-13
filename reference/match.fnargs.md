# Match Function Arguments to Choices

Assigns unnamed entries of `arglist` to the entries of `choices` that
are not named in it, in order. A named entry not in `choices`, more
entries than `choices` and a choice left without an entry are errors.

## Usage

``` r
match.fnargs(arglist, choices)
```

## Arguments

- arglist:

  list of arguments, as from `list(...)`.

- choices:

  character, the argument names to match.

## Value

Integer vector of the positions in `arglist` matching `choices`.
