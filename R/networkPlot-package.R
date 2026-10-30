#' networkPlot: pathway enrichment and visualisation of kinase networks
#'
#' Two halves of one workflow:
#'
#' - [enrich_network()] runs the pathway enrichment ([enrich_kinases()]:
#'   one-sided Fisher + Benjamini-Hochberg on Reactome v2, against the
#'   kinases UKA can predict) on a network's input kinases and returns the
#'   pathway tables plus a per-node `pathway` annotation; [enrich_uka()]
#'   runs it on UKA data directly, without a network.
#' - [plot_network()] renders the interactive `visNetwork` HTML, with every
#'   style choice on [networkplot_theme()]; [plot_pathway_heatmaps()] draws
#'   the kinase-by-pathway heatmaps.
#'
#' @keywords internal
#' @importFrom dplyr %>%
#' @importFrom rlang .data
"_PACKAGE"
