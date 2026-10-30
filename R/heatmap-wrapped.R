#' Wrapped kinase-pathway heatmaps
#'
#' One small panel per pathway -- the pathway's hit kinases as columns, the
#' comparison(s) as rows, filled with the kinase statistic -- placed side by
#' side and wrapped into rows of at most `max_row_cells` kinase columns. A
#' kinase in several pathways appears in each of their panels. Writes:
#'
#' - one figure per comparison, `wrapped_<comparison>.png`: its selected
#'   pathways, one row per panel;
#' - one figure per pair of comparisons,
#'   `pairwise_overlaps/wrapped_overlap_<comparison 1>_vs_<comparison 2>.png`:
#'   the pathways selected in both, two rows per panel; a kinase that is not
#'   a hit in one comparison is grey in that row. Pairs without a shared
#'   pathway get no figure; with `pair_groups`, only pairs within one group.
#'
#' Panels are ordered by pathway category, then by best p. All figures
#' share one colour scale (symmetric, 95th percentile of |kinase statistic|
#' over every drawn cell). Figure width and height are computed from the
#' content (panels, kinase and comparison labels, wrapped pathway titles).
#'
#' @param hits Long table, one row per (comparison, pathway, hit kinase):
#'   `comparison`, `pw_id`, `pw_name`, `category`, `p_value`, `kinase`,
#'   `statistic`.
#' @param out_dir Directory to write the PNGs into.
#' @param max_row_cells Maximum kinase columns per row of panels.
#' @param max_path_chars A pair figure whose full path would be longer
#'   (Windows' limit is 260) is named by the comparisons' numbers instead,
#'   `wrapped_overlap_<i>_vs_<j>.png`, with a message naming the pair.
#' @param pair_groups See [plot_pathway_heatmaps()].
#'
#' @return The written file paths, invisibly.
#' @keywords internal
plot_wrapped_heatmaps <- function(hits, out_dir, max_row_cells = 40, max_path_chars = 250, pair_groups = NULL) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  comps <- unique(hits$comparison)
  vals <- abs(hits$statistic[is.finite(hits$statistic)])
  lim <- if (length(vals)) stats::quantile(vals, 0.95, names = FALSE) else 1
  if (!is.finite(lim) || lim == 0) lim <- 1
  col_fun <- circlize::colorRamp2(c(-lim, 0, lim), c("#0000CC", "#F7F7F7", "#C62828"))
  safe <- function(s) gsub("[^A-Za-z0-9_-]", "_", s)

  files <- character()
  for (cmp in comps) {
    f <- file.path(out_dir, paste0("wrapped_", safe(cmp), ".png"))
    draw_wrapped_figure(hits[hits$comparison == cmp, ], cmp, f, cmp, col_fun, max_row_cells)
    files <- c(files, f)
  }
  pairs <- comparison_pairs(comps, pair_groups)
  if (length(pairs)) {
    pair_dir <- file.path(out_dir, "pairwise_overlaps")
    none <- character()
    for (pr in pairs) {
      shared <- intersect(hits$pw_id[hits$comparison == pr[1]], hits$pw_id[hits$comparison == pr[2]])
      if (!length(shared)) {
        none <- c(none, paste(pr, collapse = " / "))
        next
      }
      dir.create(pair_dir, showWarnings = FALSE, recursive = TRUE)
      f <- file.path(pair_dir, paste0("wrapped_overlap_", safe(pr[1]), "_vs_", safe(pr[2]), ".png"))
      if (nchar(normalizePath(f, mustWork = FALSE)) > max_path_chars) {
        # Windows can't write paths over 260 characters: number the pair instead
        f <- file.path(pair_dir, sprintf("wrapped_overlap_%d_vs_%d.png", match(pr[1], comps), match(pr[2], comps)))
        message("Path too long, overlap figure of '", pr[1], "' and '", pr[2], "' written as ", basename(f))
      }
      draw_wrapped_figure(hits[hits$comparison %in% pr & hits$pw_id %in% shared, ], pr, f,
                          paste0("Pathways selected in both: ", pr[1], "  and  ", pr[2]), col_fun, max_row_cells)
      files <- c(files, f)
    }
    if (length(none)) message("No shared selected pathways, no overlap figure: ", paste(none, collapse = "; "))
  }
  invisible(files)
}

#' Draw one wrapped figure (all lengths in mm)
#' @keywords internal
draw_wrapped_figure <- function(d, cmps, file, title, col_fun, max_row_cells) {
  gpar <- grid::gpar; unit <- grid::unit
  pws <- d %>%
    dplyr::group_by(.data$pw_id, .data$pw_name, .data$category) %>%
    dplyr::summarise(best_p = min(.data$p_value), .groups = "drop") %>%
    dplyr::arrange(.data$category, .data$best_p)
  cell <- 4.5; gap <- 5; row_gap <- 6; title_min <- 24; fs <- 7; lh <- 3.3; margin <- 5; legend_w <- 28

  # measure text on a null device, so the PNG size can be computed before it is opened
  grDevices::pdf(NULL)
  tw <- function(s, face = 1) vapply(s, function(x) grid::convertWidth(
    grid::grobWidth(grid::textGrob(x, gp = gpar(fontsize = fs, fontface = face))), "mm", TRUE), 1, USE.NAMES = FALSE)
  wrap <- function(s, width) {
    words <- strsplit(s, " ")[[1]]; lines <- character(); cur <- ""
    for (w in words) {
      nxt <- if (nzchar(cur)) paste(cur, w) else w
      if (tw(nxt, 2) <= width || !nzchar(cur)) cur <- nxt else { lines <- c(lines, cur); cur <- w }
    }
    c(lines, cur)
  }
  panels <- lapply(seq_len(nrow(pws)), function(i) {
    pd <- d[d$pw_id == pws$pw_id[i], ]
    kin <- sort(unique(pd$kinase))
    bw <- length(kin) * cell; sw <- max(bw, title_min)
    list(cat = dplyr::coalesce(pws$category[i], ""), title = wrap(pws$pw_name[i], sw), kin = kin, d = pd, bw = bw, sw = sw)
  })
  rows <- list(); cur <- integer(); w <- 0   # in order; a panel is never split
  for (i in seq_along(panels)) {
    add <- panels[[i]]$sw + if (length(cur)) gap else 0
    if (length(cur) && w + add > max(max_row_cells * cell, panels[[i]]$sw)) {
      rows[[length(rows) + 1]] <- cur; cur <- integer(); w <- 0; add <- panels[[i]]$sw
    }
    cur <- c(cur, i); w <- w + add
  }
  rows[[length(rows) + 1]] <- cur
  nc <- length(cmps); lab_w <- max(tw(cmps)) + 3; head_h <- 2 * lh + 4
  dims <- lapply(rows, function(r) list(
    width = sum(vapply(panels[r], `[[`, 1, "sw")) + gap * (length(r) - 1),
    title_h = (max(vapply(panels[r], function(p) length(p$title), 1)) + 1) * lh + 1,
    kin_h = max(unlist(lapply(panels[r], function(p) tw(p$kin)))) + 2))
  W <- margin + max(lab_w + max(vapply(dims, `[[`, 1, "width")) + 4 + legend_w, tw(title, 2)) + margin
  body_h <- sum(vapply(dims, function(x) x$title_h + nc * cell + x$kin_h, 1)) + row_gap * (length(rows) - 1)
  H <- margin + head_h + max(body_h, 45) + margin
  grDevices::dev.off()

  grDevices::png(file, width = W, height = H, units = "mm", res = 300)
  on.exit(grDevices::dev.off())
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(xscale = c(0, W), yscale = c(H, 0)))   # mm, y pointing down
  mm <- function(v) unit(v, "native")
  grid::grid.text(title, x = mm(margin), y = mm(margin + lh), just = "left", gp = gpar(fontsize = fs + 2, fontface = 2))
  y <- margin + head_h
  for (ri in seq_along(rows)) {
    x <- margin + lab_w; top <- y + dims[[ri]]$title_h
    grid::grid.text(cmps, x = mm(x - 2), y = mm(top + (seq_len(nc) - 0.5) * cell), just = "right", gp = gpar(fontsize = fs))
    for (i in rows[[ri]]) {
      p <- panels[[i]]; bx <- x + (p$sw - p$bw) / 2; cx <- x + p$sw / 2
      grid::grid.text(p$cat, x = mm(cx), y = mm(top - (length(p$title) + 0.5) * lh - 0.5),
                      gp = gpar(fontsize = fs - 1, col = "grey45", fontface = 3))
      grid::grid.text(p$title, x = mm(cx), y = mm(top - (rev(seq_along(p$title)) - 0.5) * lh - 0.5),
                      gp = gpar(fontsize = fs, fontface = 2))
      for (ci in seq_along(cmps)) {
        v <- p$d$statistic[match(paste(cmps[ci], p$kin), paste(p$d$comparison, p$d$kinase))]
        fill <- ifelse(is.finite(v), col_fun(ifelse(is.finite(v), v, 0)), "#E6E6E6")
        grid::grid.rect(x = mm(bx + (seq_along(p$kin) - 1) * cell), y = mm(top + (ci - 1) * cell),
                        width = mm(cell), height = mm(cell), just = c("left", "bottom"),
                        gp = gpar(fill = fill, col = "white", lwd = 0.6))
      }
      grid::grid.text(p$kin, x = mm(bx + (seq_along(p$kin) - 0.5) * cell), y = mm(top + nc * cell + 1),
                      rot = 90, just = "right", gp = gpar(fontsize = fs))
      x <- x + p$sw + gap
    }
    y <- top + nc * cell + dims[[ri]]$kin_h + row_gap
  }
  grid::upViewport()
  lgd <- ComplexHeatmap::Legend(col_fun = col_fun, title = "Kinase\nstatistic",
                                title_gp = gpar(fontsize = fs, fontface = 2), labels_gp = gpar(fontsize = fs))
  ComplexHeatmap::draw(lgd, x = unit(W - margin - legend_w, "mm"), y = unit(H - margin - head_h, "mm"),
                       just = c("left", "top"))
  if (nc > 1) {
    grid::grid.text("grey: not a hit\nin that comparison", x = unit(W - margin - legend_w, "mm"),
                    y = unit(H - margin - head_h - 35, "mm"), just = c("left", "top"),
                    gp = gpar(fontsize = fs - 1, col = "grey30"))
  }
  invisible(c(width = W, height = H, rows = length(rows)))
}
