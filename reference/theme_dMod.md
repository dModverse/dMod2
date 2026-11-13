# Standard Plotting Theme of dMod

Standard Plotting Theme of dMod

## Usage

``` r
theme_dMod(base_size = 12, base_family = "", showGrid = FALSE)
```

## Arguments

- base_size:

  Numeric, base font size. Defaults to 12.

- base_family:

  Character, font family. Defaults to `""`.

- showGrid:

  Logical, keep the panel grid. `FALSE` (default) drops it; `TRUE` keeps
  the grid of
  [`ggplot2::theme_bw()`](https://ggplot2.tidyverse.org/reference/ggtheme.html).

## Value

A `ggplot2` theme.

## See also

[`scale_color_dMod()`](https://dmodverse.github.io/dMod2/reference/scale_color_dMod.md)

## Examples

``` r
library(ggplot2)
ggplot(data.frame(x = 1:10, y = (1:10)^2), aes(x, y)) +
  geom_line() + theme_dMod()
```
