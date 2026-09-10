#' Load the Reactome reference tables
#'
#' The pathway list and parent/child relation table
#' [reactome_postprocess()] needs. Defaults to the snapshot bundled with
#' this package (`inst/extdata`); pass paths to use your own, e.g. a fresh
#' download from [download_reactome_refs()].
#'
#' @param pathways_txt Path to a `ReactomePathways.txt`
#'   (`id \t name \t organism`, no header).
#' @param relations_txt Path to a `ReactomePathwaysRelation.txt`
#'   (`parent \t child`, no header).
#'
#' @return A list with `pathways` (`pathway_id`, `pathway_name`, `organism`;
#'   filtered to *Homo sapiens*) and `relations` (`parent`, `child`).
#' @export
reactome_refs <- function(pathways_txt = system.file("extdata", "ReactomePathways.txt", package = "networkPlot"),
                          relations_txt = system.file("extdata", "ReactomePathwaysRelation.txt", package = "networkPlot")) {
  if (!nzchar(pathways_txt) || !file.exists(pathways_txt)) {
    stop("ReactomePathways.txt not found at '", pathways_txt, "'.", call. = FALSE)
  }
  if (!nzchar(relations_txt) || !file.exists(relations_txt)) {
    stop("ReactomePathwaysRelation.txt not found at '", relations_txt, "'.", call. = FALSE)
  }
  pathways <- readr::read_delim(
    pathways_txt, delim = "\t", col_names = c("pathway_id", "pathway_name", "organism"),
    show_col_types = FALSE, progress = FALSE
  )
  pathways <- pathways[pathways$organism == "Homo sapiens", , drop = FALSE]
  relations <- readr::read_delim(
    relations_txt, delim = "\t", col_names = c("parent", "child"),
    show_col_types = FALSE, progress = FALSE
  )
  list(pathways = pathways, relations = relations)
}

#' Fetch the current Reactome reference tables
#'
#' Downloads `ReactomePathways.txt` and `ReactomePathwaysRelation.txt` from
#' reactome.org into `dest`. Run this manually to refresh the bundled
#' snapshot when Reactome cuts a new release (~quarterly) -- it is never
#' called automatically.
#'
#' @param dest Directory to write the two files into.
#'
#' @return The two file paths, invisibly.
#' @export
download_reactome_refs <- function(dest) {
  dir.create(dest, showWarnings = FALSE, recursive = TRUE)
  base <- "https://reactome.org/download/current/"
  files <- c("ReactomePathways.txt", "ReactomePathwaysRelation.txt")
  paths <- file.path(dest, files)
  for (i in seq_along(files)) {
    utils::download.file(paste0(base, files[i]), paths[i], mode = "wb", quiet = TRUE)
  }
  message("Wrote ", paste(paths, collapse = ", "))
  invisible(paths)
}

#' Reactome-specific pathway post-processing
#'
#' Returns the `postprocess` function [enrich_network()] applies to its
#' filtered Enrichr hits: collapse redundant terms by shared Reactome
#' parent, drop top-level (over-general) terms, and -- when WikiPathways is
#' among the queried databases -- filter WikiPathways hits by ontology tag.
#'
#' @param refs A [reactome_refs()] list.
#'
#' @return `function(pathway_df) -> pathway_df`.
#' @export
reactome_postprocess <- function(refs = reactome_refs()) {
  force(refs)
  function(pw_df) {
    if ("Database" %in% colnames(pw_df) && any(pw_df$Database == "WikiPathways_2024_Human")) {
      pw_df <- filter_wps_by_ontology(pw_df)
    }
    add_reactome_hierarchy(pw_df, refs$pathways, refs$relations)
  }
}

#' @keywords internal
id_of <- function(name, rpws) {
  match_idx <- match(gsub("/", " ", tolower(name)), gsub("/", " ", tolower(rpws$pathway_name)))
  rpws$pathway_id[match_idx]
}

#' @keywords internal
get_ancestors <- function(stId, rel, rpws) {
  parents <- rel$parent[rel$child == stId]
  if (length(parents) == 0) return(character(0))
  parent_names <- rpws$pathway_name[match(parents, rpws$pathway_id)]
  unique(c(parent_names, unlist(lapply(parents, get_ancestors, rel = rel, rpws = rpws))))
}

#' Collapse redundant Reactome terms by shared parent; drop top-level terms
#'
#' Ported from
#' `Network_generation/R/network_enrichment_and_vis.R::add_reactome_hierarchy()`.
#' The final `filter(hierarchy_n > 1)` is **deliberate** -- top-level
#' Reactome pathways ("Signal Transduction", "Metabolism", ...) are too
#' generic to be informative in the result -- and is now an explicit
#' `return()` so it is not mistaken for a stray console line.
#'
#' @param pw_df Filtered Enrichr hits (needs `Pathway`, `Genes`, `Clusters`,
#'   `n_hits`, `pathway_size`, `overlap`, `Adjusted.P.value.mean`,
#'   `Combined.Score.mean`).
#' @param rpws,rel [reactome_refs()]'s `pathways` / `relations`.
#'
#' @return `pw_df` with `Group_Pw`, `Pathway` (possibly collapsed to a
#'   parent), `hierarchy_ids`, `hierarchy_n`, `pathway_id`, `group`; rows
#'   whose pathway is top-level (`hierarchy_n <= 1`) removed.
#' @keywords internal
#' @import dplyr
#' @importFrom stringr str_extract
add_reactome_hierarchy <- function(pw_df, rpws, rel) {
  pw_df$pathway_id <- vapply(pw_df$Pathway, function(name) {
    st <- id_of(name, rpws)
    if (is.na(st)) NA_character_ else st
  }, character(1))

  pw_df$hierarchy_ids <- vapply(pw_df$Pathway, function(name) {
    st <- id_of(name, rpws)
    if (is.na(st)) NA_character_ else paste(get_ancestors(st, rel, rpws), collapse = ";")
  }, character(1))

  pw_df <- pw_df %>%
    mutate(
      hierarchy_ids = ifelse(.data$hierarchy_ids == "", NA_character_, .data$hierarchy_ids),
      hierarchy_ids = coalesce(.data$hierarchy_ids, .data$Pathway)
    )
  pw_df$group <- sub("^.*;", "", pw_df$hierarchy_ids)

  top_level_ids <- setdiff(unique(rel$parent), unique(rel$child))
  top_level_names <- rpws$pathway_name[match(top_level_ids, rpws$pathway_id)]
  top_level_names <- top_level_names[!is.na(top_level_names)]
  pw_df <- pw_df %>% filter(!.data$Pathway %in% top_level_names)

  pw_df_n_grouped <- pw_df %>%
    mutate(
      hierarchy_ids = paste0(.data$hierarchy_ids, ";"),
      Immediate_Upstream = str_extract(.data$hierarchy_ids, "^[^;]+")
    ) %>%
    group_by(.data$Clusters, .data$Genes, .data$Immediate_Upstream) %>%
    mutate(N_Redundant = n()) %>%
    ungroup() %>%
    mutate(Final_Pathway = if_else(.data$N_Redundant > 1, .data$Immediate_Upstream, .data$Pathway)) %>%
    relocate("Final_Pathway") %>%
    distinct(.data$Clusters, .data$Genes, .data$Final_Pathway, .keep_all = TRUE) %>%
    select(-"Pathway", -"N_Redundant", -"Immediate_Upstream") %>%
    mutate(Group_Pw = ifelse(!is.na(.data$pathway_id), paste0(.data$group, "_", .data$Final_Pathway), .data$Final_Pathway)) %>%
    relocate("Group_Pw", "Final_Pathway") %>%
    rename(Pathway = "Final_Pathway")

  first_hierarchy_v <- sub(";.*", "", pw_df_n_grouped$hierarchy_ids)
  pw_v <- pw_df_n_grouped$Pathway

  pw_df_n_grouped1 <- pw_df_n_grouped %>%
    mutate(
      first_hierarchy = sub(";.*", "", .data$hierarchy_ids),
      is_upstream = if_else(.data$Pathway %in% first_hierarchy_v, 1, 0),
      is_redundant = if_else(.data$first_hierarchy %in% pw_v, 1, 0)
    ) %>%
    arrange(.data$Genes, desc(.data$is_upstream))

  pw_df_n_grouped2 <- pw_df_n_grouped1 %>%
    filter(.data$is_redundant != 1) %>%
    select(-"first_hierarchy", -"is_upstream", -"is_redundant")
  pw_df_n_grouped2$hierarchy_n <- lengths(strsplit(pw_df_n_grouped2$hierarchy_ids, ";"))

  keep <- c("Group_Pw", "Pathway", "Genes", "hierarchy_n", "n_hits", "pathway_size",
            "overlap", "Adjusted.P.value.mean", "Combined.Score.mean", "Clusters",
            "hierarchy_ids", "pathway_id", "group")
  pw_df_n_grouped2 <- pw_df_n_grouped2[, intersect(keep, colnames(pw_df_n_grouped2)), drop = FALSE]

  # Deliberate: drop top-level / one-deep pathways -- too generic to be
  # useful in the result.
  return(pw_df_n_grouped2[pw_df_n_grouped2$hierarchy_n > 1, , drop = FALSE])
}

#' Filter WikiPathways hits to signalling/metabolic/regulatory ontology tags
#'
#' Ported from
#' `Network_generation/R/network_enrichment_and_vis.R::filter_wps_by_ontology()`.
#' Queries the WikiPathways SPARQL endpoint for each pathway's ontology
#' tags and drops disease/drug/syndrome pathways. Needs `httr` + `jsonlite`;
#' if either is missing, or the query fails, returns `pw_df` unchanged with
#' a warning.
#'
#' @param pw_df A pathway data frame with `Term` and `Database` columns.
#'
#' @return `pw_df`, WikiPathways rows filtered.
#' @keywords internal
#' @import dplyr
#' @importFrom stringr str_detect
#' @importFrom tidyr separate_wider_delim
filter_wps_by_ontology <- function(pw_df) {
  if (!requireNamespace("httr", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE)) {
    warning("httr/jsonlite not installed; skipping WikiPathways ontology filter.", call. = FALSE)
    return(pw_df)
  }

  wp <- pw_df %>%
    filter(.data$Database == "WikiPathways_2024_Human") %>%
    separate_wider_delim("Term", delim = " WP", names = c("Term", "ID")) %>%
    mutate(ID = paste0("WP", .data$ID), Term = tolower(trimws(.data$Term)), Term = sub("-", " ", .data$Term))
  other <- pw_df %>% filter(.data$Database != "WikiPathways_2024_Human") %>% mutate(ID = "")
  if (nrow(wp) == 0) return(pw_df)

  sparql <- paste(
    "PREFIX wp: <http://vocabularies.wikipathways.org/wp#>",
    "PREFIX dcterms: <http://purl.org/dc/terms/>",
    "PREFIX dc: <http://purl.org/dc/elements/1.1/>",
    "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>",
    "SELECT DISTINCT ?wpid ?pathwayLabel ?ontologyLabel WHERE {",
    "?pathway a wp:Pathway ; dcterms:identifier ?wpid ; dc:title ?pathwayLabel ;",
    "wp:ontologyTag ?ontologyTerm ; wp:organismName \"Homo sapiens\" .",
    "?ontologyTerm rdfs:label ?ontologyLabel . }",
    sep = "\n"
  )
  onto <- tryCatch({
    resp <- httr::GET("https://sparql.wikipathways.org/sparql",
                      query = list(query = sparql, format = "json"))
    if (httr::status_code(resp) != 200) return(NULL)
    b <- jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"))$results$bindings
    if (length(b) == 0) return(NULL)
    data.frame(
      wpid = b$wpid$value,
      ontology = tolower(trimws(b$ontologyLabel$value)),
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)

  if (is.null(onto)) {
    warning("WikiPathways SPARQL query failed; returning pathways unfiltered.", call. = FALSE)
    return(pw_df)
  }

  antipattern <- paste(
    "drug", "syndrome", "development", "disorder", "infectious disease", "epilepsy",
    "chromosomal duplication", "parkinson", "alzheimer", "charcot", "pallister",
    "bubonic plague", "impotence", sep = "|"
  )
  onto_sum <- onto %>%
    group_by(.data$wpid) %>%
    summarize(ontologies = paste0(.data$ontology, collapse = "&"), .groups = "drop")
  keep_ids <- onto_sum$wpid[!grepl(antipattern, onto_sum$ontologies)]

  wp_kept <- wp %>% filter(.data$ID %in% trimws(keep_ids))
  bind_rows(wp_kept, other)
}
