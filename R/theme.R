#' Default node-type -> visNetwork shape map
#'
#' The shapes [plot_network()] gives each node type when no override is
#' supplied. A type not listed here falls back to `theme$shape_default`.
#'
#' @keywords internal
networkplot_default_shapes <- function() {
  c(
    "Kinase" = "dot",
    "Artificial" = "star",
    "Peptide" = "ellipse",
    "Peptide-Kinase" = "box",
    "Sensitivity" = "square",
    "Sensitivity-Kinase" = "hexagon",
    "RNA-Kinase" = "diamond",
    "RNA-Peptide-Kinase" = "star",
    "Protein" = "triangle",
    "Protein-Kinase" = "diamond",
    "Protein-Peptide" = "hexagon",
    "Protein-Peptide-Kinase" = "star",
    "Hidden" = "text"
  )
}

#' Styling configuration for [plot_network()]
#'
#' Every visual choice `plot_network()` makes is a field on the list this
#' returns. Call it with no arguments for the standard look; override a
#' field by name to change just that. `shapes` is *merged* over the default
#' map (an override for `Kinase` leaves every other type's shape alone), so
#' the "config file" is really this one R constructor -- a
#' `read_networkplot_theme()` reader for YAML/JSON could be layered on later
#' without changing `plot_network()`.
#'
#' @param shapes Named character vector, node `type` -> visNetwork shape,
#'   merged over [networkplot_default_shapes()].
#' @param shape_default Shape for a node whose `type` is not in `shapes`.
#' @param logfc_palette Length >= 2 character vector of colours interpolated
#'   into the diverging node-fill scale (low `LogFC` -> ... -> high).
#' @param logfc_limit_fn Function `(numeric) -> numeric(1)`: given the
#'   vector of node `LogFC` values, return the symmetric colour limit
#'   (values are clamped to `c(-lim, lim)` before mapping). Default: the
#'   larger of `|p25 of the negatives|` and `|p75 of the positives|`.
#' @param highlight_degree Nodes with `degree >= highlight_degree` get the
#'   highlight border/width; `NULL` disables highlighting (every node uses
#'   the defaults).
#' @param highlight_border,border_default Node border colour for
#'   high-degree nodes vs. the rest.
#' @param highlight_border_width,border_width_default Node border width for
#'   high-degree nodes vs. the rest.
#' @param edge_colour Colour for every edge.
#' @param edge_width `NULL` for constant-width edges (the default), or a
#'   list `list(column, scale, invert)` mapping a per-edge column to
#'   visNetwork's `width`: `column` is the edge-data column to read
#'   (e.g. `"cost"` or `"weight"`), `scale` a function `(numeric) ->
#'   numeric` producing widths, `invert` (default `FALSE`) negates the
#'   column first so that *small* values render *bold* (use `TRUE` for
#'   `"cost"`, where lower cost = stronger interaction).
#' @param node_font_size,legend_font_size Font sizes.
#' @param layout,seed `visIgraphLayout()` layout name and its random seed.
#'
#' @return A list of class `"networkplot_theme"`.
#' @export
networkplot_theme <- function(shapes = NULL,
                              shape_default = "ellipse",
                              logfc_palette = c("#0000CC", "#BDBDBD", "#C62828"),
                              logfc_limit_fn = NULL,
                              highlight_degree = 5,
                              highlight_border = "#4C9900",
                              border_default = "#454545",
                              highlight_border_width = 2.5,
                              border_width_default = 1,
                              edge_colour = "#000000",
                              edge_width = NULL,
                              node_font_size = 10,
                              legend_font_size = 12,
                              layout = "layout_with_fr",
                              seed = 123) {
  merged_shapes <- networkplot_default_shapes()
  if (!is.null(shapes)) merged_shapes[names(shapes)] <- unname(shapes)

  if (is.null(logfc_limit_fn)) {
    logfc_limit_fn <- function(x) {
      x <- x[is.finite(x)]
      p25_neg <- stats::quantile(x[x < 0], probs = 0.25, na.rm = TRUE)
      p75_pos <- stats::quantile(x[x > 0], probs = 0.75, na.rm = TRUE)
      lim <- max(abs(p25_neg), abs(p75_pos), na.rm = TRUE)
      if (!is.finite(lim) || lim == 0) lim <- max(abs(x), 1, na.rm = TRUE)
      lim
    }
  }

  if (!is.null(edge_width)) {
    edge_width$column <- edge_width$column %||% "cost"
    edge_width$scale <- edge_width$scale %||% function(x) scales::rescale(x, to = c(1, 6))
    edge_width$invert <- isTRUE(edge_width$invert)
  }

  structure(
    list(
      shapes = merged_shapes,
      shape_default = shape_default,
      logfc_palette = logfc_palette,
      logfc_limit_fn = logfc_limit_fn,
      highlight_degree = highlight_degree,
      highlight_border = highlight_border,
      border_default = border_default,
      highlight_border_width = highlight_border_width,
      border_width_default = border_width_default,
      edge_colour = edge_colour,
      edge_width = edge_width,
      node_font_size = node_font_size,
      legend_font_size = legend_font_size,
      layout = layout,
      seed = seed
    ),
    class = "networkplot_theme"
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x
