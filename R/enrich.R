#' Input kinases and connector nodes of a network
#'
#' Pathway enrichment tests only a network's **input kinases**: the
#' kinases of one UKA comparison that passed the input filter. **Connector
#' nodes** -- the non-input nodes a network method adds to connect them
#' (PCSF hidden nodes, SIGNOR bridges) -- are not drawn from the UKA kinase
#' list, so counting them against a UKAdb background would inflate the hits.
#'
#' Recognised node tables:
#' - a `role` column: `"input"` or `"kinase"` (the SIGNOR nodes CSV) is an
#'   input kinase, anything else (`"connector"`, `"bridge"`,
#'   `"bridge:<type>"`, `"phenotype"`) a connector node;
#' - otherwise a `type` column (networkGen's PCSF nodes): `"Hidden"` is a
#'   connector node, every other type an input.
#'
#' @param nodes A network's node table.
#'
#' @return Character vector, one per row of `nodes`: `"input"` or
#'   `"connector"`.
#' @export
network_roles <- function(nodes) {
  if ("role" %in% names(nodes)) {
    return(ifelse(nodes$role %in% c("input", "kinase"), "input", "connector"))
  }
  if ("type" %in% names(nodes)) {
    return(ifelse(nodes$type %in% "Hidden", "connector", "input"))
  }
  stop("Can't tell input kinases from connector nodes: the node table needs a `role` column ",
       "(\"input\"/\"kinase\" vs anything else, e.g. the SIGNOR nodes CSV) or a `type` column ",
       "(\"Hidden\" = connector, e.g. networkGen's PCSF nodes).", call. = FALSE)
}

#' Pathway enrichment of a list of kinases
#'
#' Over-representation of `kinases` in each pathway of a pathway set,
#' against the kinases that could have been in the list: one-sided Fisher
#' (hypergeometric) per pathway, Benjamini-Hochberg over the tested
#' pathways.
#'
#' - **Kinase universe**: `universe` (UKAdb, the kinases UKA can predict),
#'   restricted to `interactome_ids` when given (the proteins the network
#'   method could use, e.g. STRING for PCSF, SIGNOR-known for SIGNOR).
#' - **Tested pathway**: at least `min_pw_kin` kinase-universe kinases (K).
#' - **Background** (N): kinase-universe kinases in at least one tested
#'   pathway -- a kinase outside every tested pathway can never be a hit.
#' - n = `kinases` in the background, k = `kinases` in the pathway,
#'   p = P(X >= k) = `phyper(k - 1, K, N - K, n, lower.tail = FALSE)`.
#' - **Selector** (`in_selector`): tested, p < `selector_p` and >=
#'   `min_n_hits` hits; when more than `selector_max` pass that, also FDR <
#'   `selector_fdr`.
#'
#' Kinases are matched to the pathway table by UniProt ID. A kinase without
#' a `uniprot` value is looked up by `name` in `universe`; kinases that
#' can't be mapped are listed in `unmapped` and counted in `params`.
#'
#' @param kinases The kinases to test: a character vector of names, or a
#'   data frame with `name` and optionally `uniprot`.
#' @param pathway_table Pathway-to-protein table (`pathway_id`,
#'   `pathway_name`, `uniprot`, optional `category`). Default
#'   [reactome_v2_pathway_table()].
#' @param universe The kinase universe before the interactome restriction,
#'   a data frame with `name` and `uniprot`. Default [ukadb_kinases()].
#' @param interactome_ids Optional names or UniProt IDs of the proteins the
#'   network was built from; the universe is restricted to them. `NULL`
#'   (default, and always in UKA mode): no restriction.
#' @param min_n_hits Minimum hits for a pathway to enter the selector.
#' @param min_pw_kin Minimum kinase-universe kinases for a pathway to be
#'   tested.
#' @param csv_p_cutoff `all` keeps the tested pathways with p below this.
#' @param selector_p,selector_max,selector_fdr The selector rule, see above.
#' @param include_disease,normal_members Filter for pathways with
#'   mutant/defect or infection, applied with this run's kinase universe
#'   (after the interactome restriction); see [filter_disease_pathways()].
#' @param source Free text for `params`: what the kinases are (e.g.
#'   `"network input kinases"`, `"UKA input kinases, no network"`).
#'
#' @return A list:
#'   - `enrichment`: every pathway of `pathway_table` (tested or not), sorted by p;
#'   - `all`: the tested pathways with p < `csv_p_cutoff` -- `pw_id`, `pw_name`,
#'     `category`, `n_hits`, `hits` (";"-joined names), `pw_size_ukadb` (K),
#'     `p_value`, `fdr`, `in_selector`;
#'   - `short`: the selected pathways (`in_selector`), one row per distinct
#'     hit list: the pathway with the smallest p (on equal p the first by
#'     name); the other selected pathways with the same hits in
#'     `also_matches`;
#'   - `params`: `parameter`, `value`, `meaning`;
#'   - `unmapped`: input kinase names without a UniProt ID.
#' @export
enrich_kinases <- function(kinases,
                           pathway_table = reactome_v2_pathway_table(),
                           universe = ukadb_kinases(),
                           interactome_ids = NULL,
                           min_n_hits = 3, min_pw_kin = 3, csv_p_cutoff = 0.1,
                           selector_p = 0.05, selector_max = 20, selector_fdr = 0.3,
                           include_disease = "specific", normal_members = reactome_v2_normal_members(),
                           source = "input kinases") {
  if (is.character(kinases)) kinases <- tibble::tibble(name = kinases)
  if (!"uniprot" %in% names(kinases)) kinases$uniprot <- NA_character_
  kinases$uniprot <- dplyr::coalesce(as.character(kinases$uniprot), universe$uniprot[match(kinases$name, universe$name)])
  unmapped <- unique(kinases$name[is.na(kinases$uniprot)])
  kinases <- dplyr::distinct(kinases[!is.na(kinases$uniprot), c("name", "uniprot")], .data$uniprot, .keep_all = TRUE)

  uni <- universe
  if (!is.null(interactome_ids)) uni <- uni[uni$name %in% interactome_ids | uni$uniprot %in% interactome_ids, ]
  pathway_table <- filter_disease_pathways(pathway_table, uni$uniprot, include_disease, normal_members)
  in_u <- dplyr::distinct(pathway_table[pathway_table$uniprot %in% uni$uniprot, ], .data$pathway_id, .data$uniprot)
  tested_ids <- names(which(table(in_u$pathway_id) >= min_pw_kin))
  bg <- unique(in_u$uniprot[in_u$pathway_id %in% tested_ids])
  hit_bg <- intersect(kinases$uniprot, bg)
  N <- length(bg); n <- length(hit_bg)
  nm <- function(ids) paste(sort(kinases$name[match(ids, kinases$uniprot)]), collapse = ";")
  if (!"category" %in% names(pathway_table)) pathway_table$category <- NA_character_

  res <- pathway_table %>%
    dplyr::group_by(.data$pathway_id, .data$pathway_name, .data$category) %>%
    dplyr::summarise(K = length(intersect(unique(.data$uniprot), bg)),
                     hit_ids = list(intersect(unique(.data$uniprot), hit_bg)), .groups = "drop") %>%
    dplyr::mutate(k = lengths(.data$hit_ids), tested = .data$K >= min_pw_kin,
                  p_value = ifelse(.data$tested, stats::phyper(.data$k - 1, .data$K, N - .data$K, n, lower.tail = FALSE),
                                   NA_real_))
  res$fdr <- NA_real_
  res$fdr[res$tested] <- stats::p.adjust(res$p_value[res$tested], method = "BH")
  enrichment <- res %>%
    dplyr::transmute(pw_id = .data$pathway_id, pw_name = .data$pathway_name, category = .data$category,
                     status = ifelse(.data$tested, "tested", sprintf("not tested (< %d universe kinases)", min_pw_kin)),
                     n_hits = .data$k, hits = vapply(.data$hit_ids, nm, character(1)),
                     pw_size_ukadb = .data$K, p_value = .data$p_value, fdr = .data$fdr) %>%
    dplyr::arrange(.data$p_value, .data$pw_name)
  pass <- enrichment$status == "tested" & enrichment$n_hits >= min_n_hits &
    !is.na(enrichment$p_value) & enrichment$p_value < selector_p
  many <- sum(pass) > selector_max
  enrichment$in_selector <- pass & (!many | enrichment$fdr < selector_fdr)
  selector_rule <- if (many) {
    sprintf("p < %g and >= %d hits: %d pathways (more than %d), so only those with FDR < %g",
            selector_p, min_n_hits, sum(pass), selector_max, selector_fdr)
  } else sprintf("p < %g and >= %d hits", selector_p, min_n_hits)

  all <- enrichment[enrichment$status == "tested" & enrichment$p_value < csv_p_cutoff, ]
  all$status <- NULL
  short <- pick_one_per_hit_list(all[all$in_selector, ])

  params <- tibble::tribble(
    ~parameter, ~value, ~meaning,
    "source", source, "what was tested",
    "n_kin_tested", as.character(n), "n: tested kinases in the background",
    "n_kin_unmapped", as.character(length(unmapped)), "tested kinases without a UniProt ID (left out)",
    "background_size", as.character(N),
      "N: kinase universe (UKAdb, restricted to the interactome when given) in at least one tested pathway",
    "universe", if (is.null(interactome_ids)) "UKAdb kinases (no interactome restriction)" else
      sprintf("UKAdb kinases in the interactome (%d of %d)", nrow(uni), nrow(universe)), "kinase universe",
    "n_pw_tested", as.character(length(tested_ids)), "pathways with at least min_pw_kin universe kinases",
    "min_pw_kin", as.character(min_pw_kin), "minimum universe kinases for a pathway to be tested",
    "min_n_hits", as.character(min_n_hits), "minimum hits for the one-pathway-per-hit-list table and the selector",
    "csv_p_cutoff", format(csv_p_cutoff), "the _all table keeps tested pathways with p below this",
    "selector_rule", selector_rule,
      sprintf("p < %g and >= min_n_hits; when more than %d pass, also FDR < %g", selector_p, selector_max, selector_fdr),
    "test", "one-sided Fisher (hypergeometric)", "p = P(X >= k) = phyper(k - 1, K, N - K, n, lower.tail = FALSE)",
    "correction", "Benjamini-Hochberg", "over the tested pathways",
    "pathway_set", attr(pathway_table, "pathway_set") %||% "user-supplied pathway table", "pathway set",
    "disease_filter", attr(pathway_table, "disease_filter") %||% "not applied", "pathways with mutant/defect or infection",
    "disease_excluded", paste(attr(pathway_table, "disease_excluded"), collapse = "; "),
      "mutant/defect pathways left out: all their universe kinases are in their normal counterpart"
  )
  list(enrichment = enrichment, all = all, short = short, params = params, unmapped = unmapped)
}

#' One pathway per distinct hit list
#'
#' The smallest p wins, on equal p the first by name; the other pathways
#' with the same hits go into `also_matches`. Within one network the same
#' hits mean the same k, n and N, so the smallest p is the smallest
#' (best-covered) pathway -- a pathway-intrinsic rule, so the pick is the
#' same in every comparison.
#' @keywords internal
pick_one_per_hit_list <- function(all) {
  if (nrow(all) == 0) {
    all$also_matches <- character()
    return(all)
  }
  all %>%
    dplyr::group_by(.data$hits) %>%
    dplyr::arrange(.data$p_value, .data$pw_name, .by_group = TRUE) %>%
    dplyr::mutate(also_matches = paste(.data$pw_name[-1], collapse = "; ")) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$p_value, .data$pw_name)
}

#' Pathway-enrich a network
#'
#' Runs [enrich_kinases()] on the network's input kinases ([network_roles()]):
#' connector nodes are never tested. Works on networkGen's PCSF nodes
#' (`Protein`, `type`) and on SIGNOR nodes CSVs (`id`, `name`, `role`).
#'
#' @param nodes The network's node table.
#' @param out_dir If given, write the pathway CSVs to its `Network_pathways`
#'   subfolder ([network_pathways_dir()]): `pathways_<label>_all.csv`,
#'   `pathways_<label>.csv` (the selected pathways, one per hit list) and
#'   `pathways_<label>_pw_params.csv`. Requires `label`.
#' @param label Filename stem -- use `"<comparison>_spec<cutoff>"` so
#'   [plot_pathway_heatmaps_dir()] can parse it back.
#' @param ... Passed to [enrich_kinases()] (`pathway_table`,
#'   `interactome_ids`, `min_n_hits`, ...).
#'
#' @return The [enrich_kinases()] list plus `nodes`: `nodes` with a `pathway`
#'   column, the selector pathways each node belongs to (","-joined), for
#'   [plot_network()]'s pathway multi-select. Connector nodes are annotated
#'   too when their UniProt ID is known (SIGNOR `id`, or a name in UKAdb).
#' @export
enrich_network <- function(nodes, out_dir = NULL, label = NULL, ...) {
  if (!is.null(out_dir) && is.null(label)) stop("`label` is required when `out_dir` is given.", call. = FALSE)
  args <- list(...)
  args$pathway_table <- args$pathway_table %||% reactome_v2_pathway_table()
  args$universe <- args$universe %||% ukadb_kinases()
  args$source <- args$source %||% "network input kinases"

  ids <- node_ids(nodes, args$universe)
  res <- do.call(enrich_kinases, c(list(kinases = ids[network_roles(nodes) == "input", ]), args))
  res$nodes <- join_pathway_column(nodes, ids$uniprot, res$short, args$pathway_table)
  if (!is.null(out_dir)) write_pathway_csvs(res, network_pathways_dir(out_dir), label)
  res
}

#' Pathway-enrich UKA data directly, without a network
#'
#' "UKA mode": per comparison, the input kinases are selected exactly as
#' networkGen does before building a PCSF network
#' ([networkGen::build_network_grid()] with [networkGen::prep_uka()], i.e.
#' [networkGen::uka_top()]: specificity >= `spec_cutoff`, then percentile
#' rank of the kinase statistic >= `perc_cutoff`), and [enrich_kinases()]
#' runs on them. No interactome restriction applies.
#'
#' @param raw_uka Raw UKA data as read by
#'   [networkGen::load_uka_dataset_files()].
#' @param out_dir If given, write per comparison to `<out_dir>/UKA_pathways/`
#'   the three pathway CSVs (`pathways_uka_<comparison>_spec<cutoff>*.csv`).
#'   For heatmaps pass the result to [plot_uka_pathway_heatmaps()].
#' @param comparison_col,spec_cutoff,perc_cutoff,rank_uka_abs Input filter,
#'   as for [networkGen::build_network_grid()]. One value each.
#' @param ... Passed to [enrich_kinases()].
#'
#' @return Named list (`uka_<comparison>_spec<cutoff>`) of [enrich_kinases()]
#'   results, each with `nodes`: the input kinases with their kinase
#'   statistic (`name`, `kinase_statistic`).
#' @export
enrich_uka <- function(raw_uka, out_dir = NULL, comparison_col = "Sgroup_contrast",
                       spec_cutoff, perc_cutoff = 0, rank_uka_abs = TRUE, ...) {
  if (!requireNamespace("networkGen", quietly = TRUE)) {
    stop("`networkGen` is required for enrich_uka() (it holds the UKA input filter).", call. = FALSE)
  }
  cells <- networkGen::build_network_grid(
    raw_uka, clean_fn = function(x) networkGen::prep_uka(x, comparison_col = comparison_col),
    comparison_col = comparison_col, spec_cutoff = spec_cutoff, perc_cutoff = perc_cutoff,
    rank_uka_abs = rank_uka_abs,
    # no network is built: a placeholder keeps the grid from loading networkGen's STRING data
    ppi_network = data.frame(head = character(), tail = character(), cost = numeric())
  )
  args <- list(...)
  args$pathway_table <- args$pathway_table %||% reactome_v2_pathway_table()
  args$source <- args$source %||% "UKA input kinases, no network"
  out <- list()
  for (cell in cells) {
    label <- paste0("uka_", cell$condition, "_spec", cell$spec_cutoff)
    nodes <- tibble::tibble(name = cell$uka_filt$name, kinase_statistic = cell$uka_filt$LogFC)
    res <- do.call(enrich_kinases, c(list(kinases = nodes$name), args))
    res$nodes <- nodes
    if (!is.null(out_dir)) {
      dir <- uka_pathways_dir(out_dir)
      write_pathway_csvs(res, dir, label)
    }
    out[[label]] <- res
  }
  out
}

#' name + uniprot of every node: SIGNOR nodes carry the UniProt ID in `id`,
#' PCSF nodes only a gene symbol (`Protein`), looked up in `universe`.
#' @keywords internal
node_ids <- function(nodes, universe) {
  name <- if ("name" %in% names(nodes)) nodes$name else nodes$Protein
  uniprot <- if ("uniprot" %in% names(nodes)) {
    nodes$uniprot
  } else if (all(c("id", "role") %in% names(nodes))) {
    nodes$id
  } else NA_character_
  uniprot <- dplyr::coalesce(as.character(uniprot), universe$uniprot[match(name, universe$name)])
  tibble::tibble(name = name, uniprot = uniprot)
}

#' Join a `pathway` column (","-joined selector pathway names) onto `nodes`
#'
#' Commas inside a pathway name are dropped: [plot_network()] splits on ",".
#' @keywords internal
join_pathway_column <- function(nodes, uniprot, short, pathway_table) {
  sel <- short[short$in_selector, ]
  nodes$pathway <- vapply(uniprot, function(id) {
    if (is.na(id)) return(NA_character_)
    pw <- sel$pw_name[sel$pw_id %in% pathway_table$pathway_id[pathway_table$uniprot == id]]
    if (length(pw)) paste(gsub(",", "", pw), collapse = ",") else NA_character_
  }, character(1), USE.NAMES = FALSE)
  nodes
}

#' @keywords internal
write_pathway_csvs <- function(res, dir, label) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(res$all, file.path(dir, paste0("pathways_", label, "_all.csv")))
  readr::write_csv(res$short, file.path(dir, paste0("pathways_", label, ".csv")))
  readr::write_csv(res$params, file.path(dir, paste0("pathways_", label, "_pw_params.csv")))
  invisible(dir)
}
