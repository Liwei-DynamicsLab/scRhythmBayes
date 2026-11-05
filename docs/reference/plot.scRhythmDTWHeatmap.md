# Plot a scRhythmDTWHeatmap object

This S3 method draws a `scRhythmDTWHeatmap` object returned by
[`plot_dtw_heatmap()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/plot_dtw_heatmap.md).

## Usage

``` r
# S3 method for class 'scRhythmDTWHeatmap'
plot(
  x,
  merge_legends = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  padding = grid::unit(c(6, 4, 4, 4), "mm"),
  ...
)
```

## Arguments

- x:

  A `scRhythmDTWHeatmap` object returned by
  [`plot_dtw_heatmap()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/plot_dtw_heatmap.md).

- merge_legends:

  A logical value indicating whether heatmap and annotation legends
  should be merged. Defaults to `TRUE`.

- heatmap_legend_side:

  A character value specifying the heatmap legend side. Defaults to
  `"right"`.

- annotation_legend_side:

  A character value specifying the annotation legend side. Defaults to
  `"right"`.

- padding:

  A [`unit`](https://rdrr.io/r/grid/unit.html) object specifying padding
  passed to
  [`ComplexHeatmap::draw()`](https://rdrr.io/pkg/ComplexHeatmap/man/draw-dispatch.html).
  Defaults to `grid::unit(c(6, 4, 4, 4), "mm")`.

- ...:

  Additional arguments passed to
  [`ComplexHeatmap::draw()`](https://rdrr.io/pkg/ComplexHeatmap/man/draw-dispatch.html).

## Value

Invisibly returns the drawn ComplexHeatmap object.
