# Plot a scRhythmModuleHeatmap object

This S3 method draws one module-specific heatmap from a
`scRhythmModuleHeatmap` object returned by
[`plot_module_heatmap()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/plot_module_heatmap.md).

## Usage

``` r
# S3 method for class 'scRhythmModuleHeatmap'
plot(
  x,
  module_id,
  title = NULL,
  draw_separators = TRUE,
  separator_lwd = 1.2,
  ...
)
```

## Arguments

- x:

  A `scRhythmModuleHeatmap` object returned by
  [`plot_module_heatmap()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/plot_module_heatmap.md).

- module_id:

  An integer specifying which module to draw.

- title:

  A character value specifying the plot title. If `NULL`, the stored
  module title is used. Defaults to `NULL`.

- draw_separators:

  Logical; if `TRUE`, separators are drawn between phase blocks.
  Defaults to `TRUE`.

- separator_lwd:

  A numeric value specifying separator line width. Defaults to `1.2`.

- ...:

  Additional arguments passed to
  [`ComplexHeatmap::draw()`](https://rdrr.io/pkg/ComplexHeatmap/man/draw-dispatch.html).

## Value

Invisibly returns the drawn ComplexHeatmap object.
