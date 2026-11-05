# scRhythmBayes void ggplot theme

This function provides a blank ggplot2 theme used by circular or
annotation-only scRhythmBayes plotting functions.

## Usage

``` r
theme_srb_void(base_size = 12, base_family = "", ...)
```

## Arguments

- base_size:

  A numeric value specifying the base font size. Defaults to `12`.

- base_family:

  A character vector specifying the base font family. Defaults to `""`.

- ...:

  Additional arguments passed to
  [`theme`](https://ggplot2.tidyverse.org/reference/theme.html).

## Value

A ggplot2 theme object.

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)

p <- ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point()

p + theme_srb_void()
} # }
```
