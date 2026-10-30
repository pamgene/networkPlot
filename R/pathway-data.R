#' Reactome v2 pathway table
#'
#' The default pathway set of the pathway enrichment: Reactome with its
#' hierarchy cut at a curated depth per top-level category, so that each
#' pathway is general but not an umbrella category. One row per
#' pathway-protein pair, all human UniProt members. Bundled in
#' `inst/extdata` (see `data-raw/pathway-data.R` for its provenance and the
#' Reactome release). Unfiltered: [enrich_kinases()] applies the filter for
#' pathways with mutant/defect or infection ([filter_disease_pathways()])
#' with its own kinase universe.
#'
#' @return A tibble with `pathway_id`, `pathway_name`, `uniprot`,
#'   `category` (top-level Reactome category), `is_in_disease`,
#'   `disease_group`, `normal_pathway_id` (the normal counterpart(s) of a
#'   mutant/defect pathway, `"; "`-joined; `""` if none). Its
#'   `"pathway_set"` attribute names the set, Reactome release and download
#'   date.
#' @export
reactome_v2_pathway_table <- function() {
  tab <- readr::read_csv(extdata_file("reactome_v2_pathway_table.csv"), show_col_types = FALSE, col_types = "ccc")
  ann <- readr::read_csv(extdata_file("reactome_v2_pathway_annotation.csv"), show_col_types = FALSE, na = "",
                         col_types = readr::cols(.default = "c", is_in_disease = "l"))
  tab <- dplyr::left_join(
    tab,
    dplyr::select(ann, "pathway_id", "category", "is_in_disease", "disease_group", "normal_pathway_id"),
    by = "pathway_id"
  )
  if (anyNA(tab$is_in_disease)) stop("Pathways without disease annotation in reactome_v2_pathway_annotation.csv.", call. = FALSE)
  tab$normal_pathway_id <- dplyr::coalesce(tab$normal_pathway_id, "")
  info <- readLines(extdata_file("reactome_v2_pathway_table_info.txt"), encoding = "UTF-8")
  release <- sub("^Reactome release: *", "", grep("^Reactome release:", info, value = TRUE))
  release <- if (length(release)) sub(" .*", "", release[1]) else "unknown"
  dl <- sub(".* ", "", trimws(grep("UniProt2Reactome_All_Levels.txt", info, value = TRUE)))
  attr(tab, "pathway_set") <- sprintf("Reactome v2 (networkPlot; Reactome release %s, downloaded %s)", release,
                                      if (length(dl)) dl[1] else "unknown")
  tab
}

#' Members of the normal counterparts of Reactome v2's mutant/defect pathways
#'
#' Every human UniProt member (all Reactome levels, also below the Reactome
#' v2 depth cut) of each normal counterpart named in the Reactome v2
#' annotation, for [filter_disease_pathways()].
#'
#' @return A tibble with `normal_pathway_id`, `uniprot`.
#' @export
reactome_v2_normal_members <- function() {
  readr::read_csv(extdata_file("reactome_v2_normal_counterpart_members.csv"), show_col_types = FALSE, col_types = "cc")
}

#' Filter pathways with mutant/defect or infection
#'
#' Reactome flags pathways that only apply with a mutated, fused or
#' defective protein, or a pathogen (`is_in_disease`). With
#' `include_disease = "specific"` (the default) a mutant/defect pathway is
#' dropped when all its kinase-universe kinases (at least one) are also in
#' its normal counterpart: kinase data carry no information on mutations, so
#' it could only duplicate that normal pathway. One with a kinase of its own
#' could point to the mutant situation and is kept; so are pathways without
#' a normal counterpart (infections, some defects). The normal counterpart
#' counts with all its Reactome members (`normal_members`). The rule depends
#' on the kinase universe, so [enrich_kinases()] applies it with the
#' universe of each run; dropped pathways count neither in the background
#' nor in the multiple-testing correction.
#'
#' @param pathway_table A pathway table with `is_in_disease` and
#'   `normal_pathway_id` (e.g. [reactome_v2_pathway_table()]). Without
#'   `is_in_disease` it is returned unfiltered.
#' @param universe_ids UniProt IDs of the kinase universe.
#' @param include_disease `"specific"` (default): the rule above; `"all"`:
#'   keep every pathway; `"none"`: drop every pathway with mutant/defect or
#'   infection.
#' @param normal_members Members of the normal counterparts
#'   (`normal_pathway_id`, `uniprot`). Default [reactome_v2_normal_members()].
#'
#' @return `pathway_table`, filtered. Attributes: `"pathway_set"` (kept),
#'   `"disease_filter"` (what was done) and `"disease_excluded"` (names of
#'   the dropped pathways with at least one universe kinase).
#' @export
filter_disease_pathways <- function(pathway_table, universe_ids, include_disease = "specific",
                                    normal_members = reactome_v2_normal_members()) {
  include_disease <- match.arg(include_disease, c("specific", "all", "none"))
  set <- attr(pathway_table, "pathway_set")
  if (!"is_in_disease" %in% names(pathway_table)) {
    attr(pathway_table, "disease_filter") <- "not applied: the pathway table has no is_in_disease column"
    return(pathway_table)
  }
  if (include_disease == "specific" && !"normal_pathway_id" %in% names(pathway_table)) {
    stop("include_disease = \"specific\" needs a normal_pathway_id column in the pathway table.", call. = FALSE)
  }
  dis <- dplyr::distinct(pathway_table[pathway_table$is_in_disease, , drop = FALSE], .data$pathway_id, .keep_all = TRUE)
  duplicates_normal <- function(i) {
    if (is.na(dis$normal_pathway_id[i]) || !nzchar(dis$normal_pathway_id[i])) return(FALSE)
    normal <- strsplit(dis$normal_pathway_id[i], "; ", fixed = TRUE)[[1]]
    own <- intersect(pathway_table$uniprot[pathway_table$pathway_id == dis$pathway_id[i]], universe_ids)
    length(own) > 0 && all(own %in% normal_members$uniprot[normal_members$normal_pathway_id %in% normal])
  }
  drop <- switch(include_disease, all = character(), none = dis$pathway_id,
                 specific = dis$pathway_id[vapply(seq_len(nrow(dis)), duplicates_normal, logical(1))])
  out <- pathway_table[!pathway_table$pathway_id %in% drop, , drop = FALSE]
  with_kin <- unique(pathway_table$pathway_id[pathway_table$uniprot %in% universe_ids])
  attr(out, "pathway_set") <- set
  attr(out, "disease_excluded") <- sort(unique(dis$pathway_name[dis$pathway_id %in% intersect(drop, with_kin)]))
  attr(out, "disease_filter") <- sprintf(
    "include_disease = %s: %d pathways with mutant/defect or infection (with >= 1 universe kinase) excluded before testing, %d kept",
    include_disease, length(intersect(drop, with_kin)), length(setdiff(intersect(dis$pathway_id, with_kin), drop))
  )
  out
}

#' UKAdb kinases
#'
#' The kinases UKA can predict (UKAdb 0.6.0), with their gene symbol and
#' UniProt ID. They are the kinase universe of the pathway enrichment.
#'
#' @return A tibble with `name` (UKA kinase name, a gene symbol) and
#'   `uniprot`.
#' @export
ukadb_kinases <- function() {
  k <- readr::read_csv(extdata_file("ukadb_0.6.0_kinases.csv"), show_col_types = FALSE)
  dplyr::distinct(tibble::tibble(name = k$Kinase_Name, uniprot = k$Kinase_UniprotID))
}

#' @keywords internal
extdata_file <- function(name) {
  f <- system.file("extdata", name, package = "networkPlot")
  if (!nzchar(f)) stop("Bundled file ", name, " not found in networkPlot's extdata.", call. = FALSE)
  f
}
