# Synthetic pathway set: universe K1..K40; PW_A = K1..K10, PW_B = K11..K20,
# PW_small = K1, K2 (< 3 universe kinases: not tested), PW_X = K21..K23 + a non-universe protein.
pw_fixture <- function() {
  rows <- function(id, members) tibble::tibble(pathway_id = id, pathway_name = paste("Pathway", id), uniprot = members)
  dplyr::bind_rows(
    rows("A", paste0("U", 1:10)), rows("B", paste0("U", 11:20)),
    rows("small", paste0("U", 1:2)), rows("X", c(paste0("U", 21:23), "NOT_A_KINASE"))
  )
}
universe_fixture <- function() tibble::tibble(name = paste0("K", 1:40), uniprot = paste0("U", 1:40))

test_that("network_roles() splits PCSF and SIGNOR nodes into input kinases and connector nodes", {
  pcsf <- data.frame(Protein = c("AKT1", "STEINER1"), type = c("Kinase", "Hidden"))
  expect_equal(network_roles(pcsf), c("input", "connector"))
  signor <- data.frame(id = c("P1", "P2", "P3", "P4"), name = c("A", "B", "C", "D"),
                       role = c("kinase", "bridge", "bridge:complex", "phenotype"))
  expect_equal(network_roles(signor), c("input", "connector", "connector", "connector"))
  expect_error(network_roles(data.frame(Protein = "AKT1")), "role.*type")
})

test_that("enrich_kinases() runs Fisher against the background of tested pathways", {
  res <- enrich_kinases(c(paste0("K", 1:5), "K11", "K30", "NOPE"), pathway_table = pw_fixture(),
                        universe = universe_fixture(), min_n_hits = 3)
  e <- res$enrichment
  # background N = K1..K23 (in >= 1 tested pathway); K30 is in no tested pathway, NOPE unmapped
  # n = K1..K5 + K11 = 6
  a <- e[e$pw_id == "A", ]
  expect_equal(a$n_hits, 5)
  expect_equal(a$pw_size_ukadb, 10)
  expect_equal(a$p_value, stats::phyper(4, 10, 23 - 10, 6, lower.tail = FALSE))
  expect_equal(a$hits, "K1;K2;K3;K4;K5")
  expect_match(e$status[e$pw_id == "small"], "not tested")
  expect_true(is.na(e$p_value[e$pw_id == "small"]))
  expect_equal(e$pw_size_ukadb[e$pw_id == "X"], 3)   # the non-universe protein doesn't count
  expect_equal(res$unmapped, "NOPE")
  p <- stats::setNames(res$params$value, res$params$parameter)
  expect_equal(p[["background_size"]], "23")
  expect_equal(p[["n_kin_tested"]], "6")
  expect_equal(p[["n_kin_unmapped"]], "1")
  # BH over the 3 tested pathways only
  tested <- e[e$status == "tested", ]
  expect_equal(tested$fdr, stats::p.adjust(tested$p_value, "BH"))
})

test_that("enrich_kinases() restricts the universe to the interactome", {
  res <- enrich_kinases(paste0("K", 1:5), pathway_table = pw_fixture(), universe = universe_fixture(),
                        interactome_ids = c(paste0("K", 1:10), paste0("U", 11:20)))
  p <- stats::setNames(res$params$value, res$params$parameter)
  expect_equal(p[["background_size"]], "20")   # PW_X loses its kinases and is no longer tested
  expect_match(res$enrichment$status[res$enrichment$pw_id == "X"], "not tested")
})

test_that("enrich_kinases() keeps one pathway per hit list: the smallest p, others in also_matches", {
  pt <- dplyr::bind_rows(
    pw_fixture(),
    tibble::tibble(pathway_id = "A_big", pathway_name = "Pathway A_big", uniprot = paste0("U", c(1:10, 24:30)))
  )
  res <- enrich_kinases(paste0("K", 1:6), pathway_table = pt, universe = universe_fixture(), csv_p_cutoff = 1,
                        selector_p = 1.01)
  expect_setequal(res$all$pw_id[res$all$hits == "K1;K2;K3;K4;K5;K6"], c("A", "A_big"))
  s <- res$short[res$short$hits == "K1;K2;K3;K4;K5;K6", ]
  expect_equal(s$pw_id, "A")                     # smaller pathway -> smaller p
  expect_equal(s$also_matches, "Pathway A_big")
})

test_that("enrich_kinases(): `all` keeps every pathway with p < csv_p_cutoff, `short` only the selected ones", {
  res <- enrich_kinases(c(paste0("K", 1:6), "K11"), pathway_table = pw_fixture(), universe = universe_fixture(), csv_p_cutoff = 1)
  expect_true(any(!res$all$in_selector))
  expect_true(nrow(res$short) > 0 && all(res$short$in_selector))
  expect_setequal(res$short$pw_id, res$all$pw_id[res$all$in_selector])
})

test_that("enrich_kinases() switches the selector to FDR when more than selector_max pathways pass", {
  res <- enrich_kinases(paste0("K", 1:20), pathway_table = pw_fixture(), universe = universe_fixture(),
                        selector_p = 1.01, selector_max = 1, selector_fdr = -1)
  expect_false(any(res$enrichment$in_selector))
  expect_match(res$params$value[res$params$parameter == "selector_rule"], "FDR")
})

test_that("enrich_network() tests input kinases only and writes three CSVs", {
  nodes <- data.frame(Protein = c(paste0("K", 1:5), "K11", "K12"), type = c(rep("Kinase", 5), "Kinase", "Hidden"),
                      prize = 1, LogFC = 1, degree = 1)
  d <- file.path(tempdir(), "enrich_network_test"); unlink(d, recursive = TRUE)
  res <- enrich_network(nodes, out_dir = d, label = "condA_spec0.7", pathway_table = pw_fixture(),
                        universe = universe_fixture(), csv_p_cutoff = 1)
  expect_equal(res$enrichment$hits[res$enrichment$pw_id == "B"], "K11")   # K12 is Hidden: not tested
  expect_setequal(list.files(file.path(d, "Network_pathways")),
                  paste0("pathways_condA_spec0.7", c(".csv", "_all.csv", "_pw_params.csv")))
  expect_equal(res$params$value[res$params$parameter == "source"], "network input kinases")
})

test_that("enrich_network() annotates selector pathways on input kinases and on connectors with a UniProt ID", {
  nodes <- data.frame(Protein = c(paste0("K", 1:6), "K7", "STEINER1"),
                      type = c(rep("Kinase", 6), "Hidden", "Hidden"), prize = 1, LogFC = 1, degree = 1)
  res <- enrich_network(nodes, pathway_table = pw_fixture(), universe = universe_fixture())
  expect_true(res$short$in_selector[res$short$pw_id == "A"])
  expect_equal(res$nodes$pathway[res$nodes$Protein == "K1"], "Pathway A")
  expect_equal(res$nodes$pathway[res$nodes$Protein == "K7"], "Pathway A")   # connector, member of A
  expect_true(is.na(res$nodes$pathway[res$nodes$Protein == "STEINER1"]))
})

test_that("enrich_network() reads SIGNOR nodes (id = UniProt, role)", {
  nodes <- data.frame(id = c(paste0("U", 1:6), "Q_BRIDGE"), name = c(paste0("K", 1:6), "BR"),
                      role = c(rep("kinase", 6), "bridge"))
  res <- enrich_network(nodes, pathway_table = pw_fixture(), universe = universe_fixture())
  expect_equal(res$enrichment$n_hits[res$enrichment$pw_id == "A"], 6)
  expect_equal(res$params$value[res$params$parameter == "n_kin_unmapped"], "0")
})

test_that("enrich_network() needs a label with out_dir", {
  expect_error(enrich_network(data.frame(Protein = "K1", type = "Kinase"), out_dir = tempdir()), "label")
})

test_that("enrich_uka() tests exactly networkGen::uka_top()'s input kinases", {
  skip_if_not_installed("networkGen")
  k <- ukadb_kinases()$name[1:30]
  raw <- data.frame(Sgroup_contrast = rep(c("A vs B", "C vs B"), each = 30), `Kinase Name` = rep(k, 2),
                    `Median Kinase Statistic` = c(seq(-1.5, 1.5, length.out = 30), seq(1, -1, length.out = 30)),
                    `Mean Specificity Score` = rep(seq(0, 2.9, by = 0.1), 2), check.names = FALSE)
  d <- file.path(tempdir(), "enrich_uka_test"); unlink(d, recursive = TRUE)
  res <- enrich_uka(raw, out_dir = d, spec_cutoff = 1, perc_cutoff = 0.5)
  expect_named(res, c("uka_A vs B_spec1", "uka_C vs B_spec1"), ignore.order = TRUE)
  prepped <- networkGen::prep_uka(raw)
  expected <- networkGen::uka_top(prepped[prepped$Sgroup_contrast == "A vs B", ], spec_cutoff = 1, perc_cutoff = 0.5)
  expect_setequal(res[["uka_A vs B_spec1"]]$nodes$name, expected$name)
  p <- res[["uka_A vs B_spec1"]]$params
  expect_equal(p$value[p$parameter == "source"], "UKA input kinases, no network")
  expect_true(all(c("pathways_uka_A vs B_spec1.csv", "pathways_uka_A vs B_spec1_all.csv",
                    "pathways_uka_A vs B_spec1_pw_params.csv") %in%
                    list.files(file.path(d, "UKA_pathways"))))
  expect_false(any(startsWith(list.files(file.path(d, "UKA_pathways")), "nodes_")))
})
