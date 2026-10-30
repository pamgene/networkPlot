pathways_fixture <- function() {
  data.frame(
    pw_id = c("PW1", "PW2", "PW3"), pw_name = c("PW big", "PW small", "PW not selected"),
    category = c("Cat A", "Cat B", "Cat A"), n_hits = c(3, 2, 3), hits = c("K1;K2;K3", "K4;K5", "K1;K4;K5"),
    p_value = c(0.001, 0.01, 0.2), in_selector = c(TRUE, TRUE, FALSE), stringsAsFactors = FALSE
  )
}
nodes_fixture_hm <- function() {
  data.frame(Protein = paste0("K", 1:5), type = "Kinase", LogFC = c(1, -2, 3, 0.5, -0.5), stringsAsFactors = FALSE)
}

test_that("build_heatmap_data() builds a pathway x kinase matrix of selected pathways with >= min_kinases hits", {
  hd <- build_heatmap_data(pathways_fixture(), nodes_fixture_hm(), min_kinases = 3)
  expect_equal(rownames(hd$matrix), "PW big")     # PW small: 2 hits; PW not selected: not in selector
  expect_setequal(colnames(hd$matrix), c("K1", "K2", "K3"))
  expect_equal(hd$matrix["PW big", "K3"], 3)
})

test_that("selected_hits() reads the kinase statistic from PCSF and SIGNOR nodes", {
  h <- selected_hits(pathways_fixture(), nodes_fixture_hm(), 2, "cmp")
  expect_equal(nrow(h), 5)
  expect_equal(h$statistic[h$kinase == "K2"], -2)
  signor <- data.frame(name = paste0("K", 1:5), kinase_statistic = 1:5)
  expect_equal(selected_hits(pathways_fixture(), signor, 2)$statistic[1:3], 1:3)
  expect_error(selected_hits(pathways_fixture(), data.frame(x = 1), 2), "kinase_statistic")
})

test_that("plot_pathway_heatmaps() wrapped: one figure per comparison, pair figures only with shared pathways", {
  skip_if_not_installed("ComplexHeatmap")
  pw2 <- pathways_fixture(); pw2$in_selector <- c(TRUE, FALSE, FALSE)
  pw3 <- pathways_fixture(); pw3$in_selector <- c(FALSE, TRUE, FALSE)
  cmps <- list(c1 = list(pathways = pathways_fixture(), nodes = nodes_fixture_hm()),
               c2 = list(pathways = pw2, nodes = nodes_fixture_hm()),
               c3 = list(pathways = pw3, nodes = nodes_fixture_hm()))
  d <- file.path(tempdir(), "wrapped_test"); unlink(d, recursive = TRUE)
  expect_message(plot_pathway_heatmaps(cmps, d, min_kinases = 2), "c2 / c3")
  expect_setequal(list.files(d, "\\.png$"), c("wrapped_c1.png", "wrapped_c2.png", "wrapped_c3.png"))
  expect_setequal(list.files(file.path(d, "pairwise_overlaps")),
                  c("wrapped_overlap_c1_vs_c2.png", "wrapped_overlap_c1_vs_c3.png"))
})

test_that("pair_groups: pair figures only within one group (e.g. UKA Supergroup)", {
  skip_if_not_installed("ComplexHeatmap")
  cmps <- setNames(rep(list(list(pathways = pathways_fixture(), nodes = nodes_fixture_hm())), 4), paste0("c", 1:4))
  d <- file.path(tempdir(), "wrapped_groups_test"); unlink(d, recursive = TRUE)
  plot_pathway_heatmaps(cmps, d, min_kinases = 2, pair_groups = c(c1 = "SgA", c2 = "SgA", c3 = "SgB", c4 = "SgB"))
  expect_setequal(list.files(file.path(d, "pairwise_overlaps")),
                  c("wrapped_overlap_c1_vs_c2.png", "wrapped_overlap_c3_vs_c4.png"))
  expect_length(comparison_pairs(c("c1", "c2"), c(c1 = "SgA")), 0)   # no group -> no pair
  expect_length(comparison_pairs(c("c1", "c2", "c3")), 3)
})

test_that("wrapped figure size follows the content: more pathways -> more rows, never wider than the row limit", {
  skip_if_not_installed("ComplexHeatmap")
  col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  mk <- function(n_pw, n_kin) dplyr::bind_rows(lapply(seq_len(n_pw), function(i) data.frame(
    comparison = "c", pw_id = paste0("P", i), pw_name = paste("A rather long pathway name number", i),
    category = "Cat", p_value = i / 100, kinase = paste0("KINASE", seq_len(n_kin)), statistic = 0.5)))
  f <- tempfile(fileext = ".png")
  one <- draw_wrapped_figure(mk(1, 3), "c", f, "t", col_fun, max_row_cells = 40)
  many <- draw_wrapped_figure(mk(30, 12), "c", f, "t", col_fun, max_row_cells = 40)
  wide <- draw_wrapped_figure(mk(1, 60), "c", f, "t", col_fun, max_row_cells = 40)
  expect_equal(one[["rows"]], 1)
  expect_gt(many[["rows"]], 5)
  expect_gt(many[["height"]], one[["height"]])
  expect_equal(wide[["rows"]], 1)                 # a panel wider than the limit gets its own row
  expect_gt(wide[["width"]], one[["width"]])
})

test_that("plot_pathway_heatmaps() classic layout still writes per-comparison heatmaps", {
  skip_if_not_installed("ComplexHeatmap")
  d <- file.path(tempdir(), "classic_test"); unlink(d, recursive = TRUE)
  plot_pathway_heatmaps(list(c1 = list(pathways = pathways_fixture(), nodes = nodes_fixture_hm())), d,
                        layout = "classic", min_kinases = 3)
  expect_equal(list.files(d, "\\.png$"), "Kinase_Pathway_heatmap_c1.png")
})

test_that("plot_pathway_heatmaps() needs ComplexHeatmap", {
  skip_if(requireNamespace("ComplexHeatmap", quietly = TRUE),
          "ComplexHeatmap installed -- the missing-dep guard can't be exercised")
  expect_error(plot_pathway_heatmaps(list(), tempfile()), "ComplexHeatmap")
})

test_that("plot_pathway_heatmaps_dir() reads network outputs, skipping _all and _pw_params", {
  skip_if_not_installed("ComplexHeatmap")
  d <- file.path(tempdir(), "heatmap_dir_test"); unlink(d, recursive = TRUE)
  dir.create(file.path(d, "Network_pathways"), recursive = TRUE)
  for (f in c("pathways_condA_spec0.7.csv", "pathways_condA_spec0.7_all.csv"))
    readr::write_csv(pathways_fixture(), file.path(d, "Network_pathways", f))
  readr::write_csv(data.frame(parameter = "x", value = "1"), file.path(d, "Network_pathways", "pathways_condA_spec0.7_pw_params.csv"))
  readr::write_csv(nodes_fixture_hm(), file.path(d, "nodes_condA_spec0.7.csv"))   # networkGen output, flat

  out <- plot_pathway_heatmaps_dir(d, spec_cutoff = 0.7, min_kinases = 2)
  expect_equal(list.files(out, "\\.png$"), "wrapped_condA.png")
})

test_that("plot_pathway_heatmaps_dir() reads a comparison without selected pathways (empty CSV)", {
  skip_if_not_installed("ComplexHeatmap")
  d <- file.path(tempdir(), "heatmaps_empty_csv"); unlink(d, recursive = TRUE); dir.create(file.path(d, "Network_pathways"), recursive = TRUE)
  cols <- "pw_id,pw_name,category,n_hits,hits,pw_size_ukadb,p_value,fdr,in_selector,also_matches"
  writeLines(c(cols, "A,Pathway A,Cat,3,K1;K2;K3,10,0.01,0.05,TRUE,"), file.path(d, "Network_pathways", "pathways_c1_spec0.7.csv"))
  writeLines(cols, file.path(d, "Network_pathways", "pathways_c2_spec0.7.csv"))
  for (cmp in c("c1", "c2")) readr::write_csv(data.frame(Protein = c("K1", "K2", "K3"), type = "Kinase", LogFC = 1),
                                              file.path(d, paste0("nodes_", cmp, "_spec0.7.csv")))
  expect_no_error(plot_pathway_heatmaps_dir(d, 0.7, out_dir = file.path(d, "hm"), min_kinases = 3))
})

test_that("plot_uka_pathway_heatmaps() draws enrich_uka()'s result with the UKA comparison names", {
  skip_if_not_installed("ComplexHeatmap")
  uka_res <- list(`uka_A vs B_spec0.7` = list(short = pathways_fixture(),
                                               nodes = data.frame(name = paste0("K", 1:5), kinase_statistic = 1:5)))
  d <- file.path(tempdir(), "uka_heatmap_test"); unlink(d, recursive = TRUE)
  plot_uka_pathway_heatmaps(uka_res, d, min_kinases = 2)
  expect_equal(list.files(d, "[.]png$"), "wrapped_A_vs_B.png")
})
