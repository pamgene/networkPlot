test_that("reactome_refs() loads the bundled snapshot, human-only", {
  refs <- reactome_refs()
  expect_true(all(c("pathway_id", "pathway_name", "organism") %in% colnames(refs$pathways)))
  expect_equal(unique(refs$pathways$organism), "Homo sapiens")
  expect_true(all(c("parent", "child") %in% colnames(refs$relations)))
  expect_gt(nrow(refs$pathways), 1000)
})

test_that("add_reactome_hierarchy() drops top-level pathways and keeps deeper ones", {
  # Tiny hand-made hierarchy: TOP -> MID -> LEAF
  rpws <- data.frame(
    pathway_id = c("R-1", "R-2", "R-3"),
    pathway_name = c("TOP", "MID", "LEAF"),
    organism = "Homo sapiens", stringsAsFactors = FALSE
  )
  rel <- data.frame(parent = c("R-1", "R-2"), child = c("R-2", "R-3"), stringsAsFactors = FALSE)

  pw_df <- data.frame(
    Pathway = c("TOP", "LEAF"),
    Genes = c("K1;K2", "K3;K4"),
    Clusters = "all",
    n_hits = 2, pathway_size = 10, overlap = 0.2,
    Adjusted.P.value.mean = 0.01, Combined.Score.mean = 5,
    stringsAsFactors = FALSE
  )
  out <- add_reactome_hierarchy(pw_df, rpws, rel)
  expect_false("TOP" %in% out$Pathway)   # top-level dropped
  expect_true("LEAF" %in% out$Pathway)   # LEAF has ancestors TOP;MID -> hierarchy_n > 1
  expect_true(all(out$hierarchy_n > 1))
})
