# Continuous dMod2 Color Scales

Sequential (`_c`) and diverging (`_div`) counterparts of the discrete
[`scale_color_dMod2()`](https://dmodverse.github.io/dMod2/reference/scale_color_dMod2.md).
The diverging scales map `mid` to white, also for asymmetric limits.

## Usage

``` r
scale_color_dMod2_c(..., direction = 1, palette = "okabe")

scale_fill_dMod2_c(..., direction = 1, palette = "okabe")

scale_color_dMod2_div(mid = 0, ..., direction = 1, palette = "okabe")

scale_fill_dMod2_div(mid = 0, ..., direction = 1, palette = "okabe")
```

## Arguments

- ...:

  Arguments passed to
  [`ggplot2::scale_color_gradientn()`](https://ggplot2.tidyverse.org/reference/scale_gradient.html)
  or
  [`ggplot2::scale_fill_gradientn()`](https://ggplot2.tidyverse.org/reference/scale_gradient.html).

- direction:

  1 (default) or -1; -1 reverses the ramp.

- palette:

  Name of a sequential palette for `_c`, a diverging one for `_div`, see
  [`dMod2_palettes()`](https://dmodverse.github.io/dMod2/reference/dMod2_palettes.md).
  Defaults to the colorblind-safe `"okabe"`.

- mid:

  Numeric, the value mapped to white. Defaults to 0.

## Value

A continuous `ggplot2` scale.

## See also

[`scale_color_dMod2()`](https://dmodverse.github.io/dMod2/reference/scale_color_dMod2.md)

## Examples

``` r
library(ggplot2)
d <- expand.grid(time = seq(0, 10, 0.1), eps = seq(0, 1, 0.1))
ggplot(d, aes(time, exp(-eps * time), group = eps, color = eps)) +
   geom_line() + theme_dMod2() + scale_color_dMod2_c()
```
