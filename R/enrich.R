#' Pathway-enrich a PCSF network
#'
#' Runs one Enrichr query (via the `enrichR` package) over the network's
#' nodes and returns the enriched pathways plus a per-node `pathway`
#' annotation. Ports the `per_cluster = FALSE` path of
#' `Network_generation/R/network_enrichment_and_vis.R::do_network_enrichment()`
#' -- the per-cluster / topGO machinery is not carried over.
#'
#' **Which nodes**: every node *except* `type == "Hidden"` (the Steiner
#' connector nodes PCSF adds, which carry no input identity). This widens
#' the original code, which enriched only kinase-typed nodes.
#'
#' **Reference database**: `databases` is passed straight to Enrichr, so any
#' Enrichr library works. Anything Reactome-specific (collapsing redundant
#' terms by shared parent, dropping top-level terms, filtering WikiPathways
#' by ontology tag) lives in `postprocess` -- swap it, or pass `NULL`, for a
#' different database.
#'
#' @param nodes A `networkGen` result's `nodes` data frame (needs
#'   `Protein`, `type`).
#' @param databases Character vector of Enrichr library names.
#' @param min_hits Keep a pathway when its overlap with the gene list is
#'   `>= min_hits` genes.
#' @param postprocess `function(pathway_df) -> pathway_df` applied to the
#'   filtered hits, or `NULL` to skip. Default
#'   `reactome_postprocess(reactome_refs())`.
#' @param out_dir If given, write `pathways_<label>_all.csv` (all candidates)
#'   and `pathways_<label>.csv` (a provisional one-pathway-per-gene-set
#'   pick) here. Requires `label`.
#' @param label Filename stem for `out_dir` output -- use
#'   `"<comparison>_spec<cutoff>"` so [reconcile_pathways()] /
#'   [plot_pathway_heatmaps_dir()] can parse it back.
#' @param refresh If `FALSE` (default) and `out_dir/pathways_<label>.csv`
#'   already exists, read it and skip the Enrichr call. `TRUE` always calls
#'   Enrichr.
#'
#' @return A list: `all` (all filtered + post-processed candidate pathways),
#'   `short` (one pathway per gene set), `nodes` (`nodes` with a `pathway`
#'   column joined on).
#' @export
enrich_network <- function(nodes,
                           databases = "Reactome_Pathways_2024",
                           min_hits = 2,
                           postprocess = reactome_postprocess(reactome_refs()),
                           out_dir = NULL,
                           label = NULL,
                           refresh = FALSE) {
  if (!is.null(out_dir) && is.null(label)) {
    stop("`label` is required when `out_dir` is given.", call. = FALSE)
  }
  short_path <- if (!is.null(out_dir)) file.path(out_dir, paste0("pathways_", label, ".csv")) else NULL
  all_path <- if (!is.null(out_dir)) file.path(out_dir, paste0("pathways_", label, "_all.csv")) else NULL

  if (!refresh && !is.null(short_path) && file.exists(short_path)) {
    short <- readr::read_csv(short_path, show_col_types = FALSE)
    all_df <- if (file.exists(all_path)) readr::read_csv(all_path, show_col_types = FALSE) else short
    return(list(all = all_df, short = short, nodes = join_pathway_column(nodes, short)))
  }

  genes <- unique(nodes$Protein[nodes$type != "Hidden"])
  genes <- genes[!is.na(genes) & nzchar(genes)]

  raw <- enrichr_query(genes, databases)

  filt <- raw %>%
    dplyr::mutate(
      n_hits = as.numeric(stringr::str_extract(.data$Overlap, "^[0-9]+")),
      pathway_size = as.numeric(stringr::str_extract(.data$Overlap, "(?<=/)[0-9]+")),
      overlap = .data$n_hits / .data$pathway_size,
      Adjusted.P.value.mean = signif(as.numeric(.data$Adjusted.P.value), 3),
      Combined.Score.mean = as.numeric(.data$Combined.Score),
      Clusters = "all"
    ) %>%
    dplyr::filter(.data$n_hits >= min_hits) %>%
    dplyr::rename(Pathway = "Term")

  if (nrow(filt) > 0) {
    thr <- stats::quantile(filt$Combined.Score.mean, 0.25, na.rm = TRUE)
    filt <- dplyr::filter(filt, .data$Combined.Score.mean >= thr)
  }

  if (!is.null(postprocess) && nrow(filt) > 0) {
    filt <- postprocess(filt)
  }

  short <- if (nrow(filt) > 0 && "hierarchy_n" %in% colnames(filt)) {
    filt %>%
      dplyr::group_by(.data$Genes) %>%
      dplyr::slice_min(order_by = .data$hierarchy_n, with_ties = FALSE) %>%
      dplyr::ungroup()
  } else if (nrow(filt) > 0) {
    filt %>%
      dplyr::group_by(.data$Genes) %>%
      dplyr::slice_min(order_by = .data$Adjusted.P.value.mean, with_ties = FALSE) %>%
      dplyr::ungroup()
  } else {
    filt
  }

  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    readr::write_csv(filt, all_path)
    readr::write_csv(short, short_path)
  }

  list(all = filt, short = short, nodes = join_pathway_column(nodes, short))
}

#' Join a `pathway` column (comma-joined pathway names per protein) onto `nodes`
#' @keywords internal
join_pathway_column <- function(nodes, short) {
  if (nrow(short) == 0 || !all(c("Genes", "Pathway") %in% colnames(short))) {
    nodes$pathway <- NA_character_
    return(nodes)
  }
  name_col <- if ("Group_Pw" %in% colnames(short)) "Group_Pw" else "Pathway"
  prot2pw <- short %>%
    dplyr::select(pw = dplyr::all_of(name_col), "Genes") %>%
    tidyr::separate_longer_delim("Genes", delim = ";") %>%
    dplyr::distinct() %>%
    dplyr::rename(Protein = "Genes") %>%
    dplyr::group_by(.data$Protein) %>%
    dplyr::summarize(pathway = paste0(.data$pw, collapse = ","), .groups = "drop")
  dplyr::left_join(nodes, prot2pw, by = "Protein")
}

#' Thin wrapper over `enrichR::enrichr()` returning one long data frame
#' @keywords internal
enrichr_query <- function(genes, databases) {
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    stop("`enrichR` is required for enrich_network(). Install it with install.packages('enrichR').", call. = FALSE)
  }
  if (length(genes) == 0) stop("No genes to enrich (every node was type 'Hidden'?).", call. = FALSE)
  res <- enrichR::enrichr(genes, databases)
  out <- purrr::imap(res, function(df, db) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$Database <- db
    df
  })
  out <- Filter(Negate(is.null), out)
  if (length(out) == 0) stop("Enrichr returned no results for these genes.", call. = FALSE)
  dplyr::bind_rows(out)
}
