# Overview of the dMod2 Palettes

One row per palette. Pass any of the names as `palette` to
[`dMod2_palette()`](https://dmodverse.github.io/dMod2/reference/dMod2_palette.md),
[`scale_color_dMod2()`](https://dmodverse.github.io/dMod2/reference/scale_color_dMod2.md),
[`scale_color_dMod2_c()`](https://dmodverse.github.io/dMod2/reference/scale_color_dMod2_c.md)
or
[`scale_color_dMod2_div()`](https://dmodverse.github.io/dMod2/reference/scale_color_dMod2_c.md).

## Usage

``` r
dMod2_palettes(type = NULL)
```

## Arguments

- type:

  `NULL` (default) for all palettes, or any of `"qualitative"`,
  `"sequential"` and `"diverging"`.

## Value

A data frame with columns `name`, `type`, `n` (number of seed colors),
`dE` (smallest pairwise CIE2000 color distance under normal and
simulated dichromatic vision) and `colorblind` (`dE >= 10`).

## Examples

``` r
dMod2_palettes()
#>        name        type  n colorblind   dE
#> 16    okabe   diverging  5       TRUE 25.6
#> 17     puor   diverging  5       TRUE 15.7
#> 18     rdbu   diverging  5       TRUE 15.2
#> 19     dMod   diverging  3       TRUE 12.5
#> 20     brbg   diverging  5      FALSE  7.9
#> 5  contrast qualitative  3       TRUE 21.4
#> 3    medium qualitative 12       TRUE 12.1
#> 4      dark qualitative 12       TRUE 11.9
#> 2     muted qualitative 14       TRUE 11.6
#> 1     okabe qualitative 20       TRUE 11.1
#> 8      set1 qualitative  9      FALSE  3.2
#> 7     dark2 qualitative  8      FALSE  2.3
#> 6      dMod qualitative 10      FALSE  2.1
#> 9      set2 qualitative  8      FALSE  1.6
#> 10   paired qualitative 12      FALSE  1.3
#> 12  viridis  sequential  5       TRUE 13.8
#> 13     gray  sequential  2       TRUE 13.2
#> 11    okabe  sequential  6       TRUE 12.0
#> 14     heat  sequential  3       TRUE 10.7
#> 15     dMod  sequential  5      FALSE  5.4
dMod2_palettes("qualitative")
#>        name        type  n colorblind   dE
#> 5  contrast qualitative  3       TRUE 21.4
#> 3    medium qualitative 12       TRUE 12.1
#> 4      dark qualitative 12       TRUE 11.9
#> 2     muted qualitative 14       TRUE 11.6
#> 1     okabe qualitative 20       TRUE 11.1
#> 8      set1 qualitative  9      FALSE  3.2
#> 7     dark2 qualitative  8      FALSE  2.3
#> 6      dMod qualitative 10      FALSE  2.1
#> 9      set2 qualitative  8      FALSE  1.6
#> 10   paired qualitative 12      FALSE  1.3
```
