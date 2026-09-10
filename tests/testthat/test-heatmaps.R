test_that("build_heatmap_data() builds a pathway x kinase matrix and honours min_kinases", {
  pathways <- data.frame(
    Group_Pw = c("PW_big", "PW_big", "PW_big", "PW_small", "PW_small"),
    Genes = c("K1;K2;K3", "K1;K2;K3", "K1;K2;K3", "K4;K5", "K4;K5"),
    stringsAsFactors = FALSE
  )
  # separate_longer_delim expands the ";"-joined Genes; dedup by distinct()
  pathways <- data.frame(
    Group_Pw = c("PW_big", "PW_small"),
    Genes = c("K1;K2;K3", "K4;K5"),
    stringsAsFactors = FALSE
  )
  nodes <- data.frame(
    Protein = c("K1", "K2", "K3", "K4", "K5"),
    type = "Kinase",
    LogFC = c(1, -2, 3, 0.5, -0.5),
    stringsAsFactors = FALSE
  )

  hd <- build_heatmap_data(pathways, nodes, min_kinases = 3)
  expect_type(hd, "list")
  expect_true("PW_big" %in% rownames(hd$matrix))
  expect_false("PW_small" %in% rownames(hd$matrix)) # only 2 kinases < 3
  expect_setequal(colnames(hd$matrix), c("K1", "K2", "K3"))
  expect_equal(hd$matrix["PW_big", "K3"], 3)
})

test_that("build_heatmap_data() returns NULL when nothing survives", {
  pathways <- data.frame(Group_Pw = "PW", Genes = "K1;K2", stringsAsFactors = FALSE)
  nodes <- data.frame(Protein = c("K1", "K2"), type = "Sensitivity", LogFC = c(1, 2))
  expect_null(build_heatmap_data(pathways, nodes, min_kinases = 1)) # no Kinase/Artificial rows
})

test_that("plot_pathway_heatmaps() needs ComplexHeatmap", {
  skip_if(requireNamespace("ComplexHeatmap", quietly = TRUE),
          "ComplexHeatmap installed -- the missing-dep guard can't be exercised")
  expect_error(plot_pathway_heatmaps(list(), tempfile()), "ComplexHeatmap")
})
