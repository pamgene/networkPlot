nodes_fixture <- function() {
  data.frame(
    Protein = c("AKT1", "MTOR", "STEINER1", "EGFR"),
    type = c("Kinase", "Kinase", "Hidden", "Protein-Kinase"),
    prize = c(0.9, 0.8, 2, 0.7), LogFC = c(2, -1, 0, 1), degree = c(4, 3, 2, 5),
    stringsAsFactors = FALSE
  )
}

test_that("enrich_network() enriches every node except type == 'Hidden'", {
  seen <- NULL
  testthat::local_mocked_bindings(
    enrichr_query = function(genes, databases) {
      seen <<- genes
      data.frame(
        Term = "Some pathway", Overlap = "2/50",
        Adjusted.P.value = 0.001, Combined.Score = 12,
        Genes = "AKT1;MTOR", Database = databases[1], stringsAsFactors = FALSE
      )
    }
  )
  res <- enrich_network(nodes_fixture(), postprocess = NULL, min_hits = 1)
  expect_setequal(seen, c("AKT1", "MTOR", "EGFR"))   # STEINER1 (Hidden) excluded
  expect_true("pathway" %in% colnames(res$nodes))
  expect_equal(res$nodes$pathway[res$nodes$Protein == "AKT1"], "Some pathway")
  expect_true(is.na(res$nodes$pathway[res$nodes$Protein == "STEINER1"]))
})

test_that("enrich_network() keeps pathways with overlap >= min_hits (not strictly >)", {
  testthat::local_mocked_bindings(
    enrichr_query = function(genes, databases) {
      data.frame(
        Term = c("keep", "drop"), Overlap = c("2/40", "1/40"),
        Adjusted.P.value = c(0.01, 0.02), Combined.Score = c(9, 9),
        Genes = c("AKT1;MTOR", "EGFR"), Database = "d", stringsAsFactors = FALSE
      )
    }
  )
  res <- enrich_network(nodes_fixture(), postprocess = NULL, min_hits = 2)
  expect_true("keep" %in% res$all$Pathway)
  expect_false("drop" %in% res$all$Pathway)
})

test_that("enrich_network() resumes from an existing short csv unless refresh = TRUE", {
  d <- file.path(tempdir(), "enrich_resume_test")
  unlink(d, recursive = TRUE); dir.create(d)
  short <- data.frame(Pathway = "cached", Genes = "AKT1;MTOR", Group_Pw = "cached",
                      stringsAsFactors = FALSE)
  readr::write_csv(short, file.path(d, "pathways_condA_spec0.7.csv"))

  called <- FALSE
  testthat::local_mocked_bindings(
    enrichr_query = function(genes, databases) { called <<- TRUE
      data.frame(Term = "fresh", Overlap = "2/9", Adjusted.P.value = 0.01,
                 Combined.Score = 9, Genes = "AKT1;MTOR", Database = "d",
                 stringsAsFactors = FALSE) }
  )
  res <- enrich_network(nodes_fixture(), out_dir = d, label = "condA_spec0.7",
                        postprocess = NULL, refresh = FALSE)
  expect_false(called)
  expect_equal(res$short$Pathway, "cached")

  res2 <- enrich_network(nodes_fixture(), out_dir = d, label = "condA_spec0.7",
                         postprocess = NULL, refresh = TRUE, min_hits = 1)
  expect_true(called)
  expect_true("fresh" %in% res2$all$Pathway)
})
