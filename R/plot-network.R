#' Build the interactive HTML network widget
#'
#' Ports the interactive `visNetwork` rendering that lived in
#' `Network_generation/R/network_enrichment_and_vis.R::visualize_network_pg()`,
#' with every hardcoded style choice moved onto [networkplot_theme()].
#'
#' Node fill is a diverging colour scale on `LogFC` (clamped to a symmetric
#' limit -- see `theme$logfc_limit_fn`). Node shape is looked up by `type`
#' in `theme$shapes`. Nodes with `degree >= theme$highlight_degree` get a
#' distinct border. Edge boldness is constant unless `theme$edge_width` maps
#' a per-edge column (e.g. `cost`) onto width.
#'
#' @param nodes Data frame with columns `Protein`, `type`, `prize`,
#'   `LogFC`, `degree`; optional `LogFC_all` (hover label) and `pathway`
#'   (comma-joined pathway names -- enables the pathway multi-select when
#'   `colour_by = "pathway"`).
#' @param edges Data frame with columns `from`, `to` (and, for
#'   `theme$edge_width`, whatever column it names -- typically `cost` or
#'   `weight`).
#' @param title Network title, shown as the widget heading.
#' @param clusters Optional data frame with columns `id`, `cluster`
#'   (a `networkGen` result's `wc_df`); merged onto the nodes for reference.
#' @param colour_by `"pathway"` (default) enables the pathway multi-select
#'   dropdown when a `pathway` column is present; anything else leaves it
#'   off.
#' @param theme A [networkplot_theme()] object.
#'
#' @return A `visNetwork` htmlwidget.
#' @export
plot_network <- function(nodes, edges, title, clusters = NULL,
                         colour_by = "pathway", theme = networkplot_theme()) {
  stopifnot(inherits(theme, "networkplot_theme"))
  req <- c("Protein", "type", "prize", "LogFC", "degree")
  missing_cols <- setdiff(req, colnames(nodes))
  if (length(missing_cols)) {
    stop("`nodes` is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  # --- edges ---------------------------------------------------------------
  edges_vis <- edges
  edges_vis$color <- theme$edge_colour
  ew <- theme$edge_width
  if (!is.null(ew)) {
    if (!ew$column %in% colnames(edges_vis)) {
      warning("theme$edge_width$column '", ew$column,
              "' not in `edges`; drawing constant-width edges.", call. = FALSE)
    } else {
      v <- as.numeric(edges_vis[[ew$column]])
      if (ew$invert) v <- -v
      w <- ew$scale(v)
      w[!is.finite(w)] <- 1
      edges_vis$width <- w
    }
  }

  # --- nodes ------------------------------------------------------------------
  nodes_vis <- nodes
  nodes_vis$id <- nodes_vis$Protein
  nodes_vis$group <- nodes_vis$type
  nodes_vis$value <- nodes_vis$prize
  nodes_vis$label <- nodes_vis$id
  nodes_vis <- nodes_vis[!duplicated(nodes_vis$id), , drop = FALSE]

  if (!is.null(clusters)) {
    nodes_vis <- dplyr::left_join(nodes_vis, clusters, by = "id")
  }

  has_pathway <- "pathway" %in% colnames(nodes_vis)
  if (identical(colour_by, "pathway") && has_pathway) {
    nodes_vis$pathway_multi <- ifelse(
      is.na(nodes_vis$pathway) | nodes_vis$pathway == "",
      NA_character_,
      gsub("\\s*,\\s*", ";", gsub('["()]', "", trimws(as.character(nodes_vis$pathway))))
    )
  }

  # shape by type
  nodes_vis$shape <- unname(theme$shapes[as.character(nodes_vis$group)])
  nodes_vis$shape[is.na(nodes_vis$shape)] <- theme$shape_default

  # diverging LogFC fill
  logfc <- as.numeric(nodes_vis$LogFC)
  logfc[is.infinite(logfc)] <- NA
  lim <- theme$logfc_limit_fn(logfc)
  col_limits <- c(-lim, lim)
  col_vals <- grDevices::colorRampPalette(theme$logfc_palette)(201)
  clamped <- pmax(pmin(logfc, col_limits[2]), col_limits[1])
  scaled <- scales::rescale(clamped, to = c(-100, 100), from = col_limits)
  scaled[is.na(scaled)] <- 0
  nodes_vis$color.background <- col_vals[round(scaled) + 101]

  # high-degree border
  if (!is.null(theme$highlight_degree)) {
    hi <- nodes_vis$degree >= theme$highlight_degree
    hi[is.na(hi)] <- FALSE
    nodes_vis$color.border <- ifelse(hi, theme$highlight_border, theme$border_default)
    nodes_vis$borderWidth <- ifelse(hi, theme$highlight_border_width, theme$border_width_default)
  } else {
    nodes_vis$color.border <- theme$border_default
    nodes_vis$borderWidth <- theme$border_width_default
  }

  # --- assemble ------------------------------------------------------------
  net <- visNetwork::visNetwork(nodes = nodes_vis, edges = edges_vis, main = title)
  net <- visNetwork::visPhysics(net)
  net <- visNetwork::visNodes(net, font = list(size = theme$node_font_size))
  net <- visNetwork::visEdges(net, smooth = FALSE)
  net <- visNetwork::visIgraphLayout(net, layout = theme$layout, randomSeed = theme$seed)
  net <- visNetwork::visInteraction(net, multiselect = TRUE)
  net <- visNetwork::visLegend(net, addNodes = build_legend(nodes_vis, theme, lim), useGroups = FALSE)

  if (identical(colour_by, "pathway") && "pathway_multi" %in% colnames(nodes_vis)) {
    net <- visNetwork::visOptions(
      net, highlightNearest = TRUE, nodesIdSelection = TRUE,
      selectedBy = list(variable = "pathway", multiple = TRUE, main = "Select Pathway:")
    )
  } else {
    net <- visNetwork::visOptions(net, highlightNearest = TRUE, nodesIdSelection = TRUE)
  }

  net
}

#' Legend rows: one per present node type, plus the LogFC colour stops,
#' plus an edge-width note when `theme$edge_width` is active.
#' @keywords internal
build_legend <- function(nodes_vis, theme, lim) {
  present <- intersect(names(theme$shapes), unique(as.character(nodes_vis$group)))
  type_legend <- data.frame(
    label = present,
    shape = unname(theme$shapes[present]),
    color = "#BDBDBD",
    borderWidth = 1,
    font.size = theme$legend_font_size,
    stringsAsFactors = FALSE
  )
  ends <- theme$logfc_palette[c(1, length(theme$logfc_palette))]
  mid <- theme$logfc_palette[ceiling(length(theme$logfc_palette) / 2)]
  color_legend <- data.frame(
    label = c(paste0("<= ", round(-lim, 2)), "0 (neutral)", paste0(">= ", round(lim, 2))),
    shape = "dot",
    color = c(ends[1], mid, ends[2]),
    borderWidth = 0,
    font.size = theme$legend_font_size,
    stringsAsFactors = FALSE
  )
  out <- rbind(type_legend, color_legend)
  if (!is.null(theme$edge_width)) {
    note <- if (identical(theme$edge_width$column, "cost")) {
      "thicker edge = stronger interaction"
    } else {
      paste0("edge width = ", theme$edge_width$column)
    }
    out <- rbind(out, data.frame(
      label = note, shape = "dot", color = "#FFFFFF", borderWidth = 0,
      font.size = theme$legend_font_size, stringsAsFactors = FALSE
    ))
  }
  out
}

#' Write a network widget to a self-contained HTML file
#'
#' The `htmlwidgets::saveWidget()` -> `visNetwork::visSave()` fallback that
#' was inlined (three nested `tryCatch`es) in
#' `Network_generation/R/kinograte_PG.R`, in one place.
#'
#' Writes into a `Network_html` subfolder of `out_dir` -- kept separate from
#' `networkGen`'s own output (`nodes_*`/`edges_*`/`wc_df_*`/`missing_nodes_*`
#' CSVs, and [enrich_network()]'s `pathways_*` CSVs), which all land directly
#' in `out_dir` with no subfolder of their own.
#'
#' @param widget A `visNetwork` htmlwidget, e.g. from [plot_network()].
#' @param out_dir The run's main output folder. The file is written to
#'   `out_dir/Network_html/filename`, not directly into `out_dir`.
#' @param filename Output file name, e.g. `"DrugA_vs_DMSO.html"`.
#'
#' @return The full path written (`out_dir/Network_html/filename`), invisibly.
#' @export
save_network_html <- function(widget, out_dir, filename) {
  path <- file.path(out_dir, "Network_html", filename)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  ok <- tryCatch(
    {
      htmlwidgets::saveWidget(widget, file = path, selfcontained = TRUE, background = "white")
      TRUE
    },
    error = function(e) {
      warning("saveWidget failed (", conditionMessage(e), "); retrying with visSave.", call. = FALSE)
      FALSE
    }
  )
  if (!ok) {
    tryCatch(
      visNetwork::visSave(widget, path, selfcontained = FALSE, background = "white"),
      error = function(e) stop("Both saveWidget and visSave failed: ", conditionMessage(e), call. = FALSE)
    )
  }
  invisible(path)
}
