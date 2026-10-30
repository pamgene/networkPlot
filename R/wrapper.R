#' Enrich a network, then build its interactive HTML widget -- one call
#'
#' Convenience wrapper combining [enrich_network()] and [plot_network()]
#' for the common case: you want both, in sequence, and don't want to wire
#' the enriched `nodes` (the `pathway` column [enrich_network()] joins on)
#' through to `plot_network()` yourself.
#'
#' Equivalent to:
#' ```r
#' enrichment <- enrich_network(nodes, out_dir, label, ...)
#' widget <- plot_network(enrichment$nodes, edges, title, clusters, colour_by, theme)
#' ```
#'
#' @inheritParams enrich_network
#' @inheritParams plot_network
#' @param ... Passed to [enrich_network()] / [enrich_kinases()].
#'
#' @return A list: `enrichment` (the [enrich_network()] result) and `widget`
#'   (the `visNetwork` htmlwidget, built from `enrichment$nodes` so its
#'   `pathway` column drives `plot_network()`'s pathway multi-select and
#'   colouring).
#' @export
enrich_and_plot_network <- function(nodes, edges, title,
                                    out_dir = NULL,
                                    label = NULL,
                                    clusters = NULL,
                                    colour_by = "pathway",
                                    theme = networkplot_theme(),
                                    ...) {
  enrichment <- enrich_network(nodes, out_dir = out_dir, label = label, ...)
  widget <- plot_network(
    enrichment$nodes, edges, title,
    clusters = clusters, colour_by = colour_by, theme = theme
  )
  list(enrichment = enrichment, widget = widget)
}
