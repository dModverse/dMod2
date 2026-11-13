# Discrete dMod Color Scales

Discrete dMod Color Scales

## Usage

``` r
scale_color_dMod(..., palette = "okabe")

scale_fill_dMod(..., palette = "okabe")
```

## Arguments

- ...:

  Arguments passed to
  [`ggplot2::discrete_scale()`](https://ggplot2.tidyverse.org/reference/discrete_scale.html).

- palette:

  Name of a qualitative palette, see
  [`dMod_palettes()`](https://dmodverse.github.io/dMod2/reference/dMod_palettes.md).
  Defaults to the colorblind-safe `"okabe"`.

## Value

A discrete `ggplot2` scale.

## See also

[`scale_color_dMod_c()`](https://dmodverse.github.io/dMod2/reference/scale_color_dMod_c.md),
[`dMod_palette()`](https://dmodverse.github.io/dMod2/reference/dMod_palette.md)

## Examples

``` r
library(ggplot2)
times <- seq(0, 2*pi, 0.1)
values <- sin(times)
data <- data.frame(
   time = times,
   value = c(values, 1.2*values, 1.4*values, 1.6*values),
   group = rep(c("C1", "C2", "C3", "C4"), each = length(times))
)
ggplot(data, aes(time, value, color = group)) + geom_line() +
   theme_dMod() + scale_color_dMod()
```
