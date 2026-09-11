nodes_fixture <- function() {
  data.frame(
    Protein = c("AKT1", "MTOR", "STEINER1"),
    type = c("Kinase", "Kinase", "Hidden"),
    prize = c(0.9, 0.8, 2), LogFC = c(2, -1, 0), degree = c(4, 3, 2),
    stringsAsFactors = FALSE
  )
}
edges_fixture <- function() {
  data.frame(from = c("AKT1", "MTOR"), to = c("MTOR", "STEINER1"),
             cost = c(0.1, 0.4), weight = c(8, 2), stringsAsFactors = FALSE)
}

test_that("enrich_and_plot_network() calls enrich then plot, wiring the enriched nodes through", {
  testthat::local_mocked_bindings(
    enrichr_query = function(genes, databases) {
      data.frame(
        Term = "Some pathway", Overlap = "2/50", Adjusted.P.value = 0.001,
        Combined.Score = 12, Genes = "AKT1;MTOR", Database = databases[1],
        stringsAsFactors = FALSE
      )
    }
  )

  out <- enrich_and_plot_network(nodes_fixture(), edges_fixture(), title = "t",
                                 postprocess = NULL, min_hits = 1)

  expect_named(out, c("enrichment", "widget"))
  expect_true("pathway" %in% colnames(out$enrichment$nodes))
  expect_s3_class(out$widget, "visNetwork")

  # the widget was built from the *enriched* nodes -- its pathway multi-select
  # is wired up, which only happens when a pathway column reaches plot_network()
  expect_true(isTRUE(out$widget$x$byselection$enabled))
})

test_that("enrich_and_plot_network() forwards plotting args (theme, colour_by)", {
  testthat::local_mocked_bindings(
    enrichr_query = function(genes, databases) {
      data.frame(Term = "PW", Overlap = "2/50", Adjusted.P.value = 0.001,
                Combined.Score = 12, Genes = "AKT1;MTOR", Database = databases[1],
                stringsAsFactors = FALSE)
    }
  )
  out <- enrich_and_plot_network(
    nodes_fixture(), edges_fixture(), title = "t", postprocess = NULL, min_hits = 1,
    colour_by = "none", theme = networkplot_theme(shapes = c(Kinase = "diamond"))
  )
  expect_false(isTRUE(out$widget$x$byselection$enabled)) # colour_by = "none" disables it
  expect_equal(unique(out$widget$x$nodes$shape[out$widget$x$nodes$type == "Kinase"]), "diamond")
})
