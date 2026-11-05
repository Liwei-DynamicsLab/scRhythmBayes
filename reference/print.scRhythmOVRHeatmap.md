# Plot a scRhythmOVRHeatmap object

This S3 method draws the ComplexHeatmap object stored in a
`scRhythmOVRHeatmap` object returned by
[`plot_ovr_heatmap()`](https://liwei-dynamicslab.github.io/scRhythmBayes/reference/plot_ovr_heatmap.md).

## Usage

``` r
# S3 method for class 'scRhythmOVRHeatmap'
print(x, ...)
```

## Arguments

- x:

  A `scRhythmOVRHeatmap` object.

- ...:

  Additional arguments passed to
  [`ComplexHeatmap::draw()`](https://rdrr.io/pkg/ComplexHeatmap/man/draw-dispatch.html).

## Value

Invisibly returns the drawn ComplexHeatmap object.
