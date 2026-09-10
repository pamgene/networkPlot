#' networkPlot: enrich and visualise PCSF networks built by networkGen
#'
#' Two halves of one workflow:
#'
#' - [enrich_network()] runs Enrichr on a network's nodes and returns a
#'   pathway table plus a per-node `pathway` annotation; the reference
#'   database (Reactome by default) is a pluggable `postprocess` step.
#' - [plot_network()] renders the interactive `visNetwork` HTML, with every
#'   style choice on [networkplot_theme()]; [plot_pathway_heatmaps()] draws
#'   the kinase-by-pathway heatmaps.
#'
#' @keywords internal
#' @importFrom dplyr %>%
#' @importFrom rlang .data
"_PACKAGE"
