# scRhythmBayes phase palette

This function returns a soft circadian color palette for Zeitgeber time
groups. The default palette is anchored to the six common sampling times
ZT02, ZT06, ZT10, ZT14, ZT18, and ZT22, and interpolated to 24 integer
hours.

## Usage

``` r
palette_srb_phase(zt = NULL, matched = TRUE, names_as = c("key", "label"))
```

## Arguments

- zt:

  A numeric vector of Zeitgeber time values. If `NULL`, colors for
  integer hours 0 to 23 are returned. Defaults to `NULL`.

- matched:

  Logical; if `TRUE`, colors are matched to `zt`. If `FALSE`, the full
  24-hour palette is returned. Defaults to `TRUE`.

- names_as:

  A character vector specifying the name format of returned colors. Must
  be one of `"key"` or `"label"`. `"key"` uses numeric hour names such
  as `"2"`; `"label"` uses names such as `"ZT02"`. Defaults to `"key"`.

## Value

A named character vector of colors.

## Examples

``` r
palette_srb_phase()
#>         0         1         2         3         4         5         6         7 
#> "#D97B93" "#DE7E8A" "#E48282" "#EA867A" "#F08A72" "#F09470" "#F19E6E" "#F2A86C" 
#>         8         9        10        11        12        13        14        15 
#> "#F3B36A" "#E0BA6E" "#CEC272" "#BBC976" "#A9D17A" "#9ECD84" "#94CA8F" "#89C79A" 
#>        16        17        18        19        20        21        22        23 
#> "#7FC4A5" "#7EB1AB" "#7E9EB2" "#7E8BB8" "#7E79BF" "#9479B3" "#AB79A9" "#C27A9E" 
palette_srb_phase(c(2, 6, 10, 14, 18, 22))
#>         2         6        10        14        18        22 
#> "#E48282" "#F19E6E" "#CEC272" "#94CA8F" "#7E9EB2" "#AB79A9" 
```
