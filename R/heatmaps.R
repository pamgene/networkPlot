#' Kinase-by-pathway LogFC heatmaps
#'
#' Ported from
#' `Network_generation/R/plotting_functions.R::plot_kinase_pathway_heatmaps()`
#' and its `create_*` helpers, with the folder discovery split into
#' [plot_pathway_heatmaps_dir()]. The dead `create_union_combined_heatmap()`
#' (its call site was already commented out) is dropped.
#'
#' For each comparison it draws one heatmap (pathways x kinases, filled with
#' `LogFC`). Then, across comparisons: a heatmap of the pathways common to
#' *all* of them; for exactly two comparisons, a "unique to each" heatmap;
#' for three or more, one pairwise-overlap heatmap per pair (in a
#' `pairwise_overlaps/` subfolder).
#'
#' @param comparisons Named list, one entry per comparison, each a list
#'   `list(pathways = <pathway data frame>, nodes = <nodes data frame>)`.
#'   `pathways` needs `Group_Pw`, `Genes`; `nodes` needs `Protein`, `type`,
#'   `LogFC`.
#' @param out_dir Directory to write the PNGs into.
#' @param min_kinases Keep a pathway row only if it has at least this many
#'   distinct kinase/artificial nodes.
#' @param w,h,w_combined PNG sizes in cm.
#'
#' @return `out_dir`, invisibly.
#' @export
plot_pathway_heatmaps <- function(comparisons, out_dir, min_kinases = 3,
                                  w = 20, h = 14, w_combined = 25) {
  require_heatmap_pkgs()
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

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
                                        file.path(out_dir, "pairwise_overlaps"))
    }
  }
  invisible(out_dir)
}

#' Draw the kinase-pathway heatmaps for one result folder
#'
#' Discovers `pathways_<comparison>_spec<cutoff>.csv` (excluding `_all.csv`)
#' and the matching `nodes_<comparison>_spec<cutoff>.csv` in `res_dir`,
#' builds the `comparisons` list, and calls [plot_pathway_heatmaps()].
#'
#' @param res_dir Folder holding the `pathways_*` / `nodes_*` CSVs.
#' @param spec_cutoff The `spec<cutoff>` value to render.
#' @param out_dir Where to write the PNGs (default: `<res_dir>/kinase_pathway_heatmaps/spec_<cutoff>`).
#' @param min_kinases,w,h,w_combined Passed to [plot_pathway_heatmaps()].
#'
#' @return `out_dir`, invisibly.
#' @export
plot_pathway_heatmaps_dir <- function(res_dir, spec_cutoff, out_dir = NULL,
                                      min_kinases = 3, w = 20, h = 14, w_combined = 25) {
  if (is.null(out_dir)) {
    out_dir <- file.path(res_dir, "kinase_pathway_heatmaps", paste0("spec_", spec_cutoff))
  }
  suffix <- paste0("_spec", spec_cutoff, ".csv")
  files <- list.files(res_dir, full.names = TRUE)
  bn <- basename(files)
  pw_files <- files[startsWith(bn, "pathways_") & endsWith(bn, suffix) & !endsWith(bn, "_all.csv")]
  if (length(pw_files) == 0) stop("No pathways_*", suffix, " files in ", res_dir, call. = FALSE)

  comparisons <- list()
  for (pf in pw_files) {
    cmp <- sub(paste0("_spec", spec_cutoff, "$"), "", sub("^pathways_", "", sub("\\.csv$", "", basename(pf))))
    nf <- file.path(dirname(pf), sub("^pathways_", "nodes_", basename(pf)))
    if (!file.exists(nf)) {
      warning("No nodes file for '", cmp, "' (expected ", basename(nf), "); skipping.", call. = FALSE)
      next
    }
    comparisons[[cmp]] <- list(
      pathways = readr::read_csv(pf, show_col_types = FALSE),
      nodes = readr::read_csv(nf, show_col_types = FALSE)
    )
  }
  plot_pathway_heatmaps(comparisons, out_dir, min_kinases, w, h, w_combined)
}

#' Build one comparison's pathway x kinase LogFC matrix
#' @keywords internal
build_heatmap_data <- function(pathways_df, nodes_df, min_kinases) {
  expanded <- pathways_df %>%
    dplyr::select("Group_Pw", "Genes") %>%
    tidyr::separate_longer_delim("Genes", delim = ";") %>%
    dplyr::filter(!is.na(.data$Genes), .data$Genes != "", !is.na(.data$Group_Pw), .data$Group_Pw != "") %>%
    dplyr::rename(pathway = "Group_Pw", Protein = "Genes") %>%
    dplyr::left_join(dplyr::select(nodes_df, "Protein", "type", "LogFC"), by = "Protein") %>%
    dplyr::filter(.data$type %in% c("Kinase", "Artificial")) %>%
    dplyr::select("Protein", "pathway", "LogFC") %>%
    dplyr::distinct()
  if (nrow(expanded) == 0) return(NULL)

  expanded <- expanded %>%
    dplyr::group_by(.data$pathway) %>%
    dplyr::filter(dplyr::n_distinct(.data$Protein) >= min_kinases) %>%
    dplyr::ungroup()
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
create_pairwise_combined_heatmaps <- function(heatmap_data_list, w_combined, h, save_folder) {
  dir.create(save_folder, showWarnings = FALSE, recursive = TRUE)
  pairs <- utils::combn(names(heatmap_data_list), 2, simplify = FALSE)
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
