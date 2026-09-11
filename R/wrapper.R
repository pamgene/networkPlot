#' Enrich a network, then build its interactive HTML widget -- one call
#'
#' Convenience wrapper combining [enrich_network()] and [plot_network()]
#' for the common case: you want both, in sequence, and don't want to wire
#' the enriched `nodes` (the `pathway` column [enrich_network()] joins on)
#' through to `plot_network()` yourself.
#'
#' Equivalent to:
#' ```r
#' enrichment <- enrich_network(nodes, databases, min_hits, postprocess, out_dir, label, refresh)
#' widget <- plot_network(enrichment$nodes, edges, title, clusters, colour_by, theme)
#' ```
#' Call [enrich_network()] and [plot_network()] separately instead when you
#' need to [reconcile_pathways()] across several comparisons before
#' plotting, or already have an enrichment result from elsewhere.
#'
#' @inheritParams enrich_network
#' @inheritParams plot_network
#'
#' @return A list: `enrichment` (the [enrich_network()] result: `all`,
#'   `short`, `nodes`) and `widget` (the `visNetwork` htmlwidget, built from
#'   `enrichment$nodes` so its `pathway` column drives `plot_network()`'s
#'   pathway multi-select and colouring).
#' @export
enrich_and_plot_network <- function(nodes, edges, title,
                                    databases = "Reactome_Pathways_2024",
                                    min_hits = 2,
                                    postprocess = reactome_postprocess(reactome_refs()),
                                    out_dir = NULL,
                                    label = NULL,
                                    refresh = FALSE,
                                    clusters = NULL,
                                    colour_by = "pathway",
                                    theme = networkplot_theme()) {
  enrichment <- enrich_network(
    nodes, databases = databases, min_hits = min_hits, postprocess = postprocess,
    out_dir = out_dir, label = label, refresh = refresh
  )
  widget <- plot_network(
    enrichment$nodes, edges, title,
    clusters = clusters, colour_by = colour_by, theme = theme
  )
  list(enrichment = enrichment, widget = widget)
}
