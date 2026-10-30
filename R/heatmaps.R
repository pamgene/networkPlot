#' Kinase-by-pathway heatmaps
#'
#' Draws the selected pathways (`in_selector`, at least `min_kinases` hits)
#' of each comparison with their hit kinases, filled with the kinase
#' statistic. Two layouts:
#'
#' - `"wrapped"` (default): one small panel per pathway, wrapped into rows;
#'   one figure per comparison and one per pair of comparisons with the
#'   pathways selected in both (see [plot_wrapped_heatmaps()]).
#'
#' With `pair_groups`, pair figures (both layouts) are drawn only for two
#' comparisons of the same group, e.g. the same UKA Supergroup.
#' - `"classic"`: for each comparison one heatmap (pathways x kinases).
#'   Then, across comparisons: a heatmap of the pathways common to *all* of
#'   them; for exactly two comparisons, a "unique to each" heatmap; for three
#'   or more, one pairwise-overlap heatmap per pair (in a
#'   `pairwise_overlaps/` subfolder).
#'
#' @param comparisons Named list, one entry per comparison, each a list
#'   `list(pathways = <pathway data frame>, nodes = <nodes data frame>)`.
#'   `pathways` is [enrich_network()]'s / [enrich_uka()]'s one-pathway-per-
#'   hit-list table (`pw_id`, `pw_name`, `category`, `hits`, `p_value`,
#'   `in_selector`); `nodes` gives the kinase statistic: `Protein` + `LogFC`
#'   (PCSF nodes) or `name` + `kinase_statistic` (SIGNOR nodes CSV,
#'   [enrich_uka()]'s input-kinases CSV).
#' @param out_dir Directory to write the PNGs into.
#' @param layout `"wrapped"` (default) or `"classic"`.
#' @param min_kinases Keep a pathway only if it has at least this many hits.
#' @param max_row_cells Wrapped layout: maximum kinase columns per row of
#'   panels.
#' @param w,h,w_combined Classic layout: PNG sizes in cm.
#' @param pair_groups Optional named character vector, comparison -> group
#'   (e.g. its UKA Supergroup). Pair figures only for comparisons of the
#'   same group; a comparison without a group gets none. `NULL` (default):
#'   every pair.
#'
#' @return `out_dir`, invisibly.
#' @export
plot_pathway_heatmaps <- function(comparisons, out_dir, layout = c("wrapped", "classic"),
                                  min_kinases = 3, max_row_cells = 40,
                                  w = 20, h = 14, w_combined = 25, pair_groups = NULL) {
  layout <- match.arg(layout)
  require_heatmap_pkgs()
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (layout == "wrapped") {
    hits <- dplyr::bind_rows(lapply(names(comparisons), function(nm) {
      selected_hits(comparisons[[nm]]$pathways, comparisons[[nm]]$nodes, min_kinases, nm)
    }))
    if (nrow(hits) == 0) {
      warning("No selected pathways with at least ", min_kinases, " hits; no heatmaps.", call. = FALSE)
      return(invisible(out_dir))
    }
    plot_wrapped_heatmaps(hits, out_dir, max_row_cells, pair_groups = pair_groups)
    return(invisible(out_dir))
  }

  col_fun <- circlize::colorRamp2
  Heatmap <- ComplexHeatmap::Heatmap
  draw <- ComplexHeatmap::draw
  gpar <- grid::gpar
  unit <- grid::unit

  heatmap_data_list <- list()
  for (nm in names(comparisons)) {
    hd <- build_heatmap_data(comparisons[[nm]]$pathways, comparisons[[nm]]$nodes, min_kinases)
    if (is.null(hd)) {
      warning("No usable kinase-pathway data for '", nm, "'; skipping.", call. = FALSE)
      next
    }
    heatmap_data_list[[nm]] <- hd

    m <- hd$matrix
    vals <- as.numeric(m)
    vals <- vals[is.finite(vals)]
    lim <- if (length(vals)) {
      l <- max(abs(stats::quantile(vals[vals < 0], 0.25, na.rm = TRUE)),
               abs(stats::quantile(vals[vals > 0], 0.75, na.rm = TRUE)), na.rm = TRUE)
      if (!is.finite(l) || l == 0) max(abs(vals), na.rm = TRUE) else l
    } else 1
    cf <- col_fun(c(-lim, 0, lim), c("#0000CC", "#e6e2e2", "#C62828"))

    ht <- Heatmap(
      m, name = "LogFC", col = cf,
      cluster_rows = TRUE, cluster_columns = TRUE,
      show_row_names = TRUE, show_column_names = TRUE,
      column_title = paste0("Kinase-Pathway Heatmap:\n", nm),
      column_title_gp = gpar(fontsize = 12, fontface = "bold"),
      row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 8),
      row_names_max_width = unit(8, "cm"),
      rect_gp = gpar(col = "black", lwd = 0.5),
      heatmap_legend_param = list(title = "LogFC", title_gp = gpar(fontsize = 10),
                                  labels_gp = gpar(fontsize = 8))
    )
    fn <- file.path(out_dir, paste0("Kinase_Pathway_heatmap_", gsub("[^A-Za-z0-9_-]", "_", nm), ".png"))
    grDevices::png(fn, width = max(20, ncol(m) * 0.4), height = max(14, nrow(m) * 0.25 + 2),
                   units = "cm", res = 300)
    draw(ht)
    grDevices::dev.off()
  }

  if (length(heatmap_data_list) >= 2) {
    create_combined_heatmap(heatmap_data_list, w_combined, h, out_dir)
    if (length(heatmap_data_list) == 2) {
      create_unique_heatmap_comparison(heatmap_data_list, w_combined, h, out_dir)
    } else {
      create_pairwise_combined_heatmaps(heatmap_data_list, w_combined, h,
                                        file.path(out_dir, "pairwise_overlaps"), pair_groups)
    }
  }
  invisible(out_dir)
}

#' Draw the kinase-pathway heatmaps for one result folder
#'
#' Discovers the selected-pathway tables
#' `Network_pathways/pathways_<comparison>_spec<cutoff>.csv` ([enrich_network()];
#' not `_all` / `_pw_params`) and the matching networkGen nodes CSVs
#' `nodes_<comparison>_spec<cutoff>.csv` in `res_dir`, builds the
#' `comparisons` list, and calls [plot_pathway_heatmaps()]. For UKA mode use
#' [plot_uka_pathway_heatmaps()].
#'
#' @param res_dir The run's main output folder.
#' @param spec_cutoff The `spec<cutoff>` value to render.
#' @param out_dir Where to write the PNGs (default:
#'   `<res_dir>/kinase_pathway_heatmaps/spec_<cutoff>`).
#' @param ... Passed to [plot_pathway_heatmaps()] (`layout`, `min_kinases`,
#'   `max_row_cells`, `w`, `h`, `w_combined`, `pair_groups`).
#'
#' @return `out_dir`, invisibly.
#' @export
plot_pathway_heatmaps_dir <- function(res_dir, spec_cutoff, out_dir = NULL, ...) {
  if (is.null(out_dir)) out_dir <- file.path(res_dir, "kinase_pathway_heatmaps", paste0("spec_", spec_cutoff))
  pw_dir <- network_pathways_dir(res_dir)
  prefix <- "pathways_"
  suffix <- paste0("_spec", spec_cutoff, ".csv")
  files <- list.files(pw_dir, full.names = TRUE)
  bn <- basename(files)
  pw_files <- files[startsWith(bn, prefix) & endsWith(bn, suffix)]
  if (length(pw_files) == 0) stop("No ", prefix, "*", suffix, " files in ", pw_dir, call. = FALSE)

  comparisons <- list()
  for (pf in pw_files) {
    cmp <- sub(paste0("_spec", spec_cutoff, "$"), "", sub(paste0("^", prefix), "", sub("\\.csv$", "", basename(pf))))
    nf <- file.path(res_dir, sub("^pathways_", "nodes_", basename(pf)))
    if (!file.exists(nf)) {
      warning("No nodes file for '", cmp, "' (expected ", basename(nf), "); skipping.", call. = FALSE)
      next
    }
    comparisons[[cmp]] <- list(
      pathways = readr::read_csv(pf, show_col_types = FALSE,
                                 col_types = readr::cols(p_value = "d", n_hits = "i", in_selector = "l")),
      nodes = readr::read_csv(nf, show_col_types = FALSE)
    )
  }
  plot_pathway_heatmaps(comparisons, out_dir, ...)
}

#' Draw the kinase-pathway heatmaps for UKA mode
#'
#' Takes [enrich_uka()]'s result (per comparison its selected pathways and
#' its input kinases with their kinase statistic) and calls
#' [plot_pathway_heatmaps()], with the comparison names as in the UKA data.
#'
#' @param uka_res The list returned by [enrich_uka()].
#' @param out_dir Where to write the PNGs.
#' @param ... Passed to [plot_pathway_heatmaps()] (`layout`, `min_kinases`,
#'   `pair_groups`, ...).
#'
#' @return `out_dir`, invisibly.
#' @export
plot_uka_pathway_heatmaps <- function(uka_res, out_dir, ...) {
  comparisons <- lapply(uka_res, function(r) list(pathways = r$short, nodes = r$nodes))
  names(comparisons) <- sub("_spec[^_]*$", "", sub("^uka_", "", names(uka_res)))
  plot_pathway_heatmaps(comparisons, out_dir, ...)
}

#' Hits of the selected pathways of one comparison, with their kinase statistic
#'
#' @return Long table: `comparison`, `pw_id`, `pw_name`, `category`,
#'   `p_value`, `kinase`, `statistic`.
#' @keywords internal
selected_hits <- function(pathways_df, nodes_df, min_kinases, comparison = NA_character_) {
  pw <- pathways_df[pathways_df$in_selector %in% TRUE & pathways_df$n_hits >= min_kinases, ]
  if (!"category" %in% names(pw)) pw$category <- NA_character_
  stats <- node_statistics(nodes_df)
  pw %>%
    dplyr::transmute(comparison = comparison, .data$pw_id, .data$pw_name, category = as.character(.data$category),
                     .data$p_value, kinase = .data$hits) %>%
    tidyr::separate_longer_delim("kinase", delim = ";") %>%
    dplyr::filter(!is.na(.data$kinase), .data$kinase != "") %>%
    dplyr::mutate(statistic = stats$statistic[match(.data$kinase, stats$name)])
}

#' name + kinase statistic of the nodes: PCSF (`Protein`, `LogFC`) or SIGNOR
#' (`name`, `kinase_statistic`) schema.
#' @keywords internal
node_statistics <- function(nodes_df) {
  if (all(c("Protein", "LogFC") %in% names(nodes_df))) {
    return(tibble::tibble(name = nodes_df$Protein, statistic = as.numeric(nodes_df$LogFC)))
  }
  if (all(c("name", "kinase_statistic") %in% names(nodes_df))) {
    return(tibble::tibble(name = nodes_df$name, statistic = as.numeric(nodes_df$kinase_statistic)))
  }
  stop("Nodes need `Protein` + `LogFC` (PCSF) or `name` + `kinase_statistic` (SIGNOR) for the heatmaps.", call. = FALSE)
}

#' Build one comparison's pathway x kinase matrix of the kinase statistic
#' @keywords internal
build_heatmap_data <- function(pathways_df, nodes_df, min_kinases) {
  expanded <- selected_hits(pathways_df, nodes_df, min_kinases) %>%
    dplyr::transmute(pathway = .data$pw_name, Protein = .data$kinase, LogFC = .data$statistic) %>%
    dplyr::distinct()
  if (nrow(expanded) == 0) return(NULL)

  m <- expanded %>%
    tidyr::pivot_wider(names_from = "Protein", values_from = "LogFC", values_fn = mean, values_fill = NA) %>%
    tibble::column_to_rownames("pathway") %>%
    as.matrix()
  m <- m[!apply(is.na(m), 1, all), !apply(is.na(m), 2, all), drop = FALSE]
  if (nrow(m) == 0 || ncol(m) == 0) return(NULL)
  m[is.na(m)] <- 0
  list(matrix = m, data = expanded)
}

#' @keywords internal
create_combined_heatmap <- function(heatmap_data_list, w_combined, h, save_folder,
                                    filename = "Kinase_Pathway_heatmap_overlapping_comparisons.png") {
  Heatmap <- ComplexHeatmap::Heatmap; draw <- ComplexHeatmap::draw
  gpar <- grid::gpar; unit <- grid::unit; colorRamp2 <- circlize::colorRamp2
  dir.create(save_folder, showWarnings = FALSE, recursive = TRUE)

  all_pathways <- lapply(heatmap_data_list, function(x) rownames(x$matrix))
  common_pathways <- Reduce(intersect, all_pathways)
  if (length(common_pathways) == 0) {
    message("No overlapping pathways across comparisons.")
    return(invisible(NULL))
  }

  vals <- unlist(lapply(heatmap_data_list, function(x) {
    as.numeric(x$matrix[rownames(x$matrix) %in% common_pathways, , drop = FALSE])
  }))
  vals <- vals[is.finite(vals)]
  limit_val <- if (length(vals)) {
    l <- max(abs(stats::quantile(vals[vals < 0], 0.25, na.rm = TRUE)),
             abs(stats::quantile(vals[vals > 0], 0.75, na.rm = TRUE)), na.rm = TRUE)
    if (!is.finite(l) || l == 0) max(abs(vals), na.rm = TRUE) else l
  } else 1
  col_fun <- colorRamp2(c(-limit_val, 0, limit_val), c("#0000CC", "#BDBDBD", "#C62828"))

  individual_heatmaps <- lapply(names(heatmap_data_list), function(comp_name) {
    cm <- heatmap_data_list[[comp_name]]$matrix
    sub <- cm[rownames(cm) %in% common_pathways, , drop = FALSE]
    sub <- sub[common_pathways, , drop = FALSE]
    sub[is.na(sub)] <- 0
    Heatmap(
      sub, name = "LogFC", col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
      column_order = order(colnames(sub)), row_order = order(rownames(sub)),
      show_row_names = TRUE, show_column_names = TRUE, column_title = comp_name,
      column_title_gp = gpar(fontsize = 8, fontface = "bold"),
      row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6),
      row_names_max_width = unit(8, "cm"), rect_gp = gpar(col = "black", lwd = 0.5),
      show_heatmap_legend = FALSE
    )
  })
  individual_heatmaps[[1]]@heatmap_param$show_heatmap_legend <- TRUE
  individual_heatmaps[[1]]@heatmap_param$heatmap_legend_param <- list(
    title = "LogFC", title_gp = gpar(fontsize = 10), labels_gp = gpar(fontsize = 8)
  )
  ht_combined <- Reduce(`+`, individual_heatmaps)

  n_comp <- length(heatmap_data_list)
  wc <- if (!is.null(w_combined)) {
    if (n_comp == 2) 25 else if (n_comp == 3) 35 else 45
  } else w_combined

  fn <- file.path(save_folder, filename)
  grDevices::png(fn, width = wc, height = h, units = "cm", res = 300)
  draw(ht_combined)
  grDevices::dev.off()
  invisible(fn)
}

#' @keywords internal
create_pairwise_combined_heatmaps <- function(heatmap_data_list, w_combined, h, save_folder, pair_groups = NULL) {
  pairs <- comparison_pairs(names(heatmap_data_list), pair_groups)
  if (length(pairs)) dir.create(save_folder, showWarnings = FALSE, recursive = TRUE)
  for (pair in pairs) {
    safe <- gsub("[^A-Za-z0-9_-]", "_", pair)
    create_combined_heatmap(
      heatmap_data_list[pair], w_combined, h, save_folder,
      filename = paste0("Kinase_Pathway_heatmap_overlap_", safe[1], "_vs_", safe[2], ".png")
    )
  }
}

#' @keywords internal
create_unique_heatmap_comparison <- function(heatmap_data_list, w, h, save_folder) {
  Heatmap <- ComplexHeatmap::Heatmap; draw <- ComplexHeatmap::draw
  `%v%` <- ComplexHeatmap::`%v%`
  gpar <- grid::gpar; unit <- grid::unit; colorRamp2 <- circlize::colorRamp2
  if (length(heatmap_data_list) != 2) return(invisible(NULL))

  comp_names <- names(heatmap_data_list)
  mat1 <- heatmap_data_list[[1]]$matrix
  mat2 <- heatmap_data_list[[2]]$matrix
  unique1 <- setdiff(rownames(mat1), rownames(mat2))
  unique2 <- setdiff(rownames(mat2), rownames(mat1))
  if (length(unique1) == 0 && length(unique2) == 0) return(invisible(NULL))

  all_kinases <- sort(union(colnames(mat1), colnames(mat2)))
  extend <- function(mat, pws, cols) {
    if (length(pws) == 0) return(matrix(0, 0, length(cols), dimnames = list(character(0), cols)))
    m <- matrix(0, length(pws), length(cols), dimnames = list(pws, cols))
    ec <- intersect(cols, colnames(mat))
    m[pws, ec] <- mat[pws, ec, drop = FALSE]
    m[is.na(m)] <- 0
    m[sort(rownames(m)), , drop = FALSE]
  }
  ext1 <- extend(mat1, unique1, all_kinases)
  ext2 <- extend(mat2, unique2, all_kinases)

  av <- c(as.numeric(ext1), as.numeric(ext2))
  av <- av[is.finite(av) & av != 0]
  limit_val <- if (length(av)) {
    l <- max(abs(c(stats::quantile(av[av < 0], 0.25, na.rm = TRUE),
                   stats::quantile(av[av > 0], 0.75, na.rm = TRUE))), na.rm = TRUE)
    if (!is.finite(l) || l == 0) max(abs(av), na.rm = TRUE) else l
  } else 1
  col_fun <- colorRamp2(c(-limit_val, 0, limit_val), c("#0000CC", "#BDBDBD", "#C62828"))

  make_ht <- function(mat, comp_name, ht_name, show_legend) {
    Heatmap(
      mat, name = ht_name, col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
      column_order = order(colnames(mat)), show_row_names = TRUE, show_column_names = TRUE,
      column_title = NULL, row_title = comp_name, row_title_side = "left",
      row_title_gp = gpar(fontsize = 8, fontface = "bold"),
      row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6),
      row_names_max_width = unit(8, "cm"), rect_gp = gpar(col = "black", lwd = 0.5),
      show_heatmap_legend = show_legend,
      heatmap_legend_param = list(title = "LogFC", title_gp = gpar(fontsize = 10),
                                  labels_gp = gpar(fontsize = 8))
    )
  }
  ht_combined <- make_ht(ext1, comp_names[1], "LogFC", TRUE) %v%
    make_ht(ext2, comp_names[2], "LogFC_2", FALSE)

  fn <- file.path(save_folder, "Kinase_Pathway_heatmap_unique_comparisons.png")
  grDevices::png(fn, width = w, height = h, units = "cm", res = 300)
  draw(ht_combined)
  grDevices::dev.off()
  invisible(fn)
}

#' @keywords internal
require_heatmap_pkgs <- function() {
  for (p in c("ComplexHeatmap", "circlize")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("`", p, "` is required for the pathway heatmaps. Install it (Bioconductor for ComplexHeatmap).", call. = FALSE)
    }
  }
}

#' The pairs of comparisons to draw overlap figures for: all pairs, or with
#' `pair_groups` (comparison -> group) only pairs within one group.
#' @keywords internal
comparison_pairs <- function(comps, pair_groups = NULL) {
  if (length(comps) < 2) return(list())
  pairs <- utils::combn(comps, 2, simplify = FALSE)
  if (is.null(pair_groups)) return(pairs)
  g <- unname(pair_groups[comps]); names(g) <- comps
  Filter(function(pr) !is.na(g[[pr[1]]]) && identical(g[[pr[1]]], g[[pr[2]]]), pairs)
}
