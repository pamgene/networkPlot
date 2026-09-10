mk <- function(group_pw, genes, hn = 2, cs = 1) {
  data.frame(
    Group_Pw = group_pw, Pathway = group_pw, Genes = genes,
    hierarchy_n = hn, Combined.Score.mean = cs, n_hits = 3, pathway_size = 20,
    overlap = 0.15, Adjusted.P.value.mean = 0.01, Clusters = "all",
    hierarchy_ids = "root;mid", pathway_id = "R-1", group = "root",
    stringsAsFactors = FALSE
  )
}

test_that("reconcile_pathways() is a no-op for fewer than two comparisons", {
  one <- list(A = mk("PW_x", "K1;K2"))
  expect_identical(reconcile_pathways(one), one)
})

test_that("reconcile_pathways() re-picks the pathway present in the most comparisons", {
  # gene set K1;K2: comparison A currently picked PW_shared, B picked PW_onlyB.
  # PW_shared appears in both comparisons' candidates -> both should end on it.
  all_tables <- list(
    A = rbind(mk("PW_shared", "K1;K2", cs = 5), mk("PW_onlyA", "K1;K2", cs = 9)),
    B = rbind(mk("PW_shared", "K1;K2", cs = 1), mk("PW_onlyB", "K1;K2", cs = 9))
  )
  out <- reconcile_pathways(all_tables)
  expect_setequal(names(out), c("A", "B"))
  expect_equal(out$A$Group_Pw, "PW_shared")
  expect_equal(out$B$Group_Pw, "PW_shared")
})

test_that("reconcile_pathways_dir() rewrites each comparison's short csv in place", {
  d <- file.path(tempdir(), "reconcile_dir_test")
  unlink(d, recursive = TRUE)
  dir.create(d)
  readr::write_csv(rbind(mk("PW_shared", "K1;K2", cs = 5), mk("PW_onlyA", "K1;K2", cs = 9)),
                   file.path(d, "pathways_A_spec0.7_all.csv"))
  readr::write_csv(rbind(mk("PW_shared", "K1;K2", cs = 1), mk("PW_onlyB", "K1;K2", cs = 9)),
                   file.path(d, "pathways_B_spec0.7_all.csv"))
  # provisional shorts (wrong picks) that should get overwritten
  readr::write_csv(mk("PW_onlyA", "K1;K2"), file.path(d, "pathways_A_spec0.7.csv"))
  readr::write_csv(mk("PW_onlyB", "K1;K2"), file.path(d, "pathways_B_spec0.7.csv"))

  reconcile_pathways_dir(d, spec_cutoff = 0.7)

  a <- readr::read_csv(file.path(d, "pathways_A_spec0.7.csv"), show_col_types = FALSE)
  b <- readr::read_csv(file.path(d, "pathways_B_spec0.7.csv"), show_col_types = FALSE)
  expect_equal(a$Group_Pw, "PW_shared")
  expect_equal(b$Group_Pw, "PW_shared")
})
