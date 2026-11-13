# Generate Distinct Colors from a dMod Palette

Returns the first `n` colors of the chosen palette. Beyond the seeds of
the palette the colors come from
[`Polychrome::createPalette()`](https://rdrr.io/pkg/Polychrome/man/createPalette.html)
if Polychrome is installed and from
[`grDevices::hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html)
otherwise. These extra colors are not colorblind-safe; a colorblind-safe
palette warns when it runs out of seeds.

## Usage

``` r
dMod_palette(n, palette = "okabe")
```

## Arguments

- n:

  Integer, number of colors.

- palette:

  Name of a qualitative palette, see
  [`dMod_palettes()`](https://dmodverse.github.io/dMod2/reference/dMod_palettes.md).
  Defaults to `"okabe"`.

## Value

Character vector of length `n` with hex color codes.

## Examples

``` r
dMod_palette(3)
#> [1] "#000000" "#E69F00" "#56B4E9"
dMod_palette(3, palette = "muted")
#> [1] "#332288" "#88CCEE" "#44AA99"
```
