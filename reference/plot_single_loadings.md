# Plot a loading matrix as a heatmap

Creates a ggplot2 heatmap (tile plot) of a numeric loading matrix, with
an optional custom y-axis labeling and an optional colorbar.

## Usage

``` r
plot_single_loadings(
  mat,
  y_label = NULL,
  fill_limits = c(-1.3, 1.5),
  show_colorbar = FALSE
)
```

## Arguments

- mat:

  Numeric matrix. Loading matrix to plot.

- y_label:

  Optional character vector of row labels (length = nrow(mat)). If
  provided, y-axis ticks correspond to row indices and labels are set to
  `y_label`.

- fill_limits:

  Numeric length-2 vector giving the limits for the diverging color
  scale. Use the same limits across plots to make them comparable.

- show_colorbar:

  Logical; whether to display the colorbar legend.

## Value

A ggplot object.
