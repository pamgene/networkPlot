#' The `Network_pathways` subfolder of a run's output folder
#'
#' [enrich_network()], [reconcile_pathways_dir()], and
#' [plot_pathway_heatmaps_dir()] all read/write the `pathways_*` CSVs here,
#' kept separate from `networkGen`'s own output
#' (`nodes_*`/`edges_*`/`wc_df_*`/`missing_nodes_*`), which lands directly in
#' `out_dir` with no subfolder of its own -- and from [save_network_html()]'s
#' `Network_html` subfolder.
#'
#' @param out_dir The run's main output folder.
#' @return `file.path(out_dir, "Network_pathways")`.
#' @keywords internal
network_pathways_dir <- function(out_dir) {
  file.path(out_dir, "Network_pathways")
}
