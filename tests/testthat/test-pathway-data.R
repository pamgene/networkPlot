test_that("reactome_v2_pathway_table() loads the bundled Reactome v2 set, unfiltered, with disease annotation", {
  pt <- reactome_v2_pathway_table()
  expect_true(all(c("pathway_id", "pathway_name", "uniprot", "category", "is_in_disease", "normal_pathway_id") %in% names(pt)))
  expect_equal(dplyr::n_distinct(pt$pathway_id), 706)
  expect_match(attr(pt, "pathway_set"), "Reactome v2 .*release 97")
  normal_ids <- unique(unlist(strsplit(pt$normal_pathway_id[nzchar(pt$normal_pathway_id)], "; ", fixed = TRUE)))
  expect_gte(mean(normal_ids %in% reactome_v2_normal_members()$normal_pathway_id), 0.95)
})

# N: normal pathway; M_dup / M_own: mutant pathways with N as normal counterpart; INF: infection, no counterpart
disease_fixture <- function() {
  rows <- function(id, members, dis, normal = "") tibble::tibble(pathway_id = id, pathway_name = paste("Pathway", id),
                                                                 uniprot = members, is_in_disease = dis, normal_pathway_id = normal)
  dplyr::bind_rows(rows("N", paste0("U", 1:5), FALSE), rows("M_dup", paste0("U", c(1:3, 90)), TRUE, "N_full"),
                   rows("M_own", paste0("U", c(1:3, 6)), TRUE, "N_full"), rows("INF", paste0("U", 1:3), TRUE))
}
normal_fixture <- function() tibble::tibble(normal_pathway_id = "N_full", uniprot = paste0("U", 1:5))

test_that("filter_disease_pathways() 'specific': the rule depends on the kinase universe", {
  small <- filter_disease_pathways(disease_fixture(), paste0("U", 1:5), "specific", normal_fixture())
  expect_setequal(unique(small$pathway_id), c("N", "INF"))   # both mutants' universe kinases (U1-U3) are all in N
  expect_equal(attr(small, "disease_excluded"), c("Pathway M_dup", "Pathway M_own"))
  big <- filter_disease_pathways(disease_fixture(), paste0("U", c(1:6, 90)), "specific", normal_fixture())
  expect_setequal(unique(big$pathway_id), c("N", "M_dup", "M_own", "INF"))   # U90 / U6: kinases of their own
})

test_that("filter_disease_pathways() never drops a pathway without universe kinases or without a normal counterpart", {
  out <- filter_disease_pathways(disease_fixture(), "U99", "specific", normal_fixture())
  expect_equal(dplyr::n_distinct(out$pathway_id), 4)
  expect_length(attr(out, "disease_excluded"), 0)
})

test_that("filter_disease_pathways() 'all' keeps everything, 'none' drops every flagged pathway, tables without flags pass", {
  expect_equal(dplyr::n_distinct(filter_disease_pathways(disease_fixture(), "U1", "all", normal_fixture())$pathway_id), 4)
  expect_equal(unique(filter_disease_pathways(disease_fixture(), "U1", "none", normal_fixture())$pathway_id), "N")
  pt <- tibble::tibble(pathway_id = "A", pathway_name = "A", uniprot = "U1")
  expect_identical(nrow(filter_disease_pathways(pt, "U1", "none")), 1L)
  expect_match(attr(filter_disease_pathways(pt, "U1"), "disease_filter"), "not applied")
})

test_that("enrich_kinases() applies the disease filter with its own universe and reports it in params", {
  uni <- tibble::tibble(name = paste0("K", c(1:6, 90)), uniprot = paste0("U", c(1:6, 90)))
  r_small <- enrich_kinases(paste0("K", 1:3), pathway_table = disease_fixture(), universe = uni,
                            interactome_ids = paste0("U", 1:5), normal_members = normal_fixture())
  r_big <- enrich_kinases(paste0("K", 1:3), pathway_table = disease_fixture(), universe = uni, normal_members = normal_fixture())
  expect_false(any(c("M_dup", "M_own") %in% r_small$enrichment$pw_id))
  expect_true(all(c("M_dup", "M_own") %in% r_big$enrichment$pw_id))
  expect_equal(r_small$params$value[r_small$params$parameter == "disease_excluded"], "Pathway M_dup; Pathway M_own")
})

test_that("ukadb_kinases() gives name and UniProt of the UKAdb kinases", {
  k <- ukadb_kinases()
  expect_named(k, c("name", "uniprot"))
  expect_true("AKT1" %in% k$name)
  expect_gt(nrow(k), 300)
})
