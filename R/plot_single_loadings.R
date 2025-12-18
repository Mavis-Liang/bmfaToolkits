#' Plot a loading matrix as a heatmap
#'
#' Creates a ggplot2 heatmap (tile plot) of a numeric loading matrix, with an
#' optional custom y-axis labeling and an optional colorbar.
#'
#' @param mat Numeric matrix. Loading matrix to plot.
#' @param y_label Optional character vector of row labels (length = nrow(mat)).
#'   If provided, y-axis ticks correspond to row indices and labels are set to
#'   `y_label`.
#' @param fill_limits Numeric length-2 vector giving the limits for the diverging
#'   color scale. Use the same limits across plots to make them comparable.
#' @param show_colorbar Logical; whether to display the colorbar legend.
#'
#' @return A ggplot object.
#' @export
plot_single_loadings <- function(mat,
                                 y_label = NULL,
                                 fill_limits = c(-1.3, 1.5),
                                 show_colorbar = FALSE) {
  if (is.null(mat)) stop("`mat` is NULL.", call. = FALSE)
  mat <- as.matrix(mat)
  if (!is.numeric(mat)) stop("`mat` must be numeric.", call. = FALSE)
  if (length(fill_limits) != 2 || any(!is.finite(fill_limits))) {
    stop("`fill_limits` must be a numeric vector of length 2.", call. = FALSE)
  }
  
  df <- reshape2::melt(mat)
  
  # Determine the range of the data
  x_range <- range(df$Var2, na.rm = TRUE)
  y_range <- range(df$Var1, na.rm = TRUE)
  
  # X breaks: at most ~6 ticks
  step_x <- max(1, ceiling((x_range[2] - x_range[1]) / 5))
  x_breaks <- seq(x_range[1], x_range[2], by = step_x)
  
  # Y breaks / labels
  if (is.null(y_label)) {
    # Show every 5 rows (reverse axis)
    y_breaks <- seq(y_range[2], y_range[1], by = -5)
    y_labels <- y_breaks
  } else {
    if (length(y_label) != nrow(mat)) {
      stop("`y_label` must have length nrow(mat).", call. = FALSE)
    }
    y_breaks <- seq_along(y_label)
    y_labels <- y_label
  }
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Var2, y = Var1, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_x_continuous(breaks = x_breaks) +
    ggplot2::scale_y_reverse(
      breaks = y_breaks,
      labels = y_labels
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.text.y.left = ggplot2::element_text(
        angle = 0, size = 15,
        margin = ggplot2::margin(r = 8, l = 10)
      ),
      axis.text.x = ggplot2::element_text(
        size = 8,
        margin = ggplot2::margin(t = -15)
      ),
      axis.text.y = ggplot2::element_text(size = 7),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        hjust = 0.5, size = 13,
        margin = ggplot2::margin(b = -10)
      )
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#998ec3",
      mid = "#f7f7f7",
      high = "#f1a340",
      midpoint = 0,
      na.value = "grey50",
      limits = fill_limits,
      name = latex2exp::TeX("Loadings")
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        barwidth = 0.3,
        barheight = 3,
        title.theme = ggplot2::element_text(size = 7),
        label.theme = ggplot2::element_text(size = 6)
      )
    )
  
  if (!isTRUE(show_colorbar)) {
    p <- p + ggplot2::guides(fill = "none")
  }
  
  p
}
