#' Reconcile per-comparison pathway picks across comparisons
#'
#' [enrich_network()] writes each comparison a *provisional* one-pathway-
#' per-gene-set table, chosen without looking at the other comparisons.
#' Given every comparison's `_all` candidate table, this re-picks -- for
#' each `(comparison, gene set)` -- the pathway that recurs across the most
#' comparisons (tie -> shallower hierarchy -> higher combined score ->
#' name), so the same biology lands under the same label everywhere and the
#' combined heatmaps line up.
#'
#' Ported from
#' `Network_generation/R/network_enrichment_and_vis.R::reconcile_pathway_selection()`,
#' with the folder I/O split into [reconcile_pathways_dir()].
#'
#' @param all_tables Named list of `_all` candidate data frames, one per
#'   comparison (names are the comparison labels). Each needs `Group_Pw`,
#'   `Genes`, `hierarchy_n`, `Combined.Score.mean`, `Pathway`.
#'
#' @return Named list, same names as `all_tables`: each comparison's
#'   reconciled one-pathway-per-gene-set table. Returns `all_tables`
#'   unchanged when fewer than two comparisons are supplied.
#' @export
reconcile_pathways <- function(all_tables) {
  if (length(all_tables) < 2) return(all_tables)

  short_cols <- c("Group_Pw", "Pathway", "Genes", "hierarchy_n", "n_hits",
                  "pathway_size", "overlap", "Adjusted.P.value.mean",
                  "Combined.Score.mean", "Clusters", "hierarchy_ids",
                  "pathway_id", "group")

  all_data <- dplyr::bind_rows(all_tables, .id = "comparison")

  presence <- all_data %>%
    dplyr::distinct(.data$comparison, .data$Group_Pw) %>%
    dplyr::count(.data$Group_Pw, name = "n_comparisons_present")

  reconciled <- all_data %>%
    dplyr::left_join(presence, by = "Group_Pw") %>%
    dplyr::group_by(.data$comparison, .data$Genes) %>%
    dplyr::arrange(
      dplyr::desc(.data$n_comparisons_present),
      .data$hierarchy_n,
      dplyr::desc(.data$Combined.Score.mean),
      .data$Pathway,
      .by_group = TRUE
    ) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  stats::setNames(
    lapply(names(all_tables), function(cmp) {
      out <- reconciled[reconciled$comparison == cmp, , drop = FALSE]
      out[, intersect(short_cols, colnames(out)), drop = FALSE]
    }),
    names(all_tables)
  )
}

#' Reconcile pathway picks for every comparison in a result folder
#'
#' Discovers `pathways_<comparison>_spec<cutoff>_all.csv` files for one
#' `spec_cutoff`, runs [reconcile_pathways()], and **overwrites** each
#' comparison's short `pathways_<comparison>_spec<cutoff>.csv` in place.
#' No-op when fewer than two comparisons are found.
#'
#' @param res_dir Folder holding the `pathways_*` CSVs.
#' @param spec_cutoff The `spec<cutoff>` value whose files to reconcile.
#'
#' @return The reconciled tables (named list), invisibly.
#' @export
reconcile_pathways_dir <- function(res_dir, spec_cutoff) {
  suffix <- paste0("_spec", spec_cutoff, "_all.csv")
  files <- list.files(res_dir, full.names = TRUE)
  bn <- basename(files)
  files <- files[startsWith(bn, "pathways_") & endsWith(bn, suffix)]
  if (length(files) < 2) return(invisible(NULL))

  comparisons <- sub("_all\\.csv$", "", sub("^pathways_", "", basename(files)))
  comparisons <- sub(paste0("_spec", spec_cutoff, "$"), "", comparisons)
  all_tables <- stats::setNames(
    lapply(files, readr::read_csv, show_col_types = FALSE),
    comparisons
  )

  reconciled <- reconcile_pathways(all_tables)
  for (cmp in names(reconciled)) {
    readr::write_csv(
      reconciled[[cmp]],
      file.path(res_dir, paste0("pathways_", cmp, "_spec", spec_cutoff, ".csv"))
    )
  }
  invisible(reconciled)
}
