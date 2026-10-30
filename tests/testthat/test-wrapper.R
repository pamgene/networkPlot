pw_fixture_wrapper <- function() {
  dplyr::bind_rows(tibble::tibble(pathway_id = "A", pathway_name = "Pathway A", uniprot = paste0("U", 1:10)),
                   tibble::tibble(pathway_id = "B", pathway_name = "Pathway B", uniprot = paste0("U", 11:40)))
}
universe_fixture_wrapper <- function() tibble::tibble(name = paste0("K", 1:40), uniprot = paste0("U", 1:40))
nodes_fixture_wrapper <- function() {
  data.frame(
    Protein = c("K1", "K2", "K3", "STEINER1"),
    type = c("Kinase", "Kinase", "Kinase", "Hidden"),
    prize = c(0.9, 0.8, 0.7, 2), LogFC = c(2, -1, 1, 0), degree = c(4, 3, 3, 2),
    stringsAsFactors = FALSE
  )
}
edges_fixture_wrapper <- function() {
  data.frame(from = c("K1", "K2", "K3"), to = c("K2", "STEINER1", "STEINER1"),
             cost = c(0.1, 0.4, 0.2), weight = c(8, 2, 3), stringsAsFactors = FALSE)
}

test_that("enrich_and_plot_network() calls enrich then plot, wiring the enriched nodes through", {
  out <- enrich_and_plot_network(nodes_fixture_wrapper(), edges_fixture_wrapper(), title = "t",
                                 pathway_table = pw_fixture_wrapper(), universe = universe_fixture_wrapper())
  expect_named(out, c("enrichment", "widget"))
  expect_equal(out$enrichment$nodes$pathway[1], "Pathway A")
  expect_s3_class(out$widget, "visNetwork")
  # the widget was built from the *enriched* nodes -- its pathway multi-select
  # is wired up, which only happens when a pathway column reaches plot_network()
  expect_true(isTRUE(out$widget$x$byselection$enabled))
})

test_that("enrich_and_plot_network() forwards plotting args (theme, colour_by)", {
  out <- enrich_and_plot_network(
    nodes_fixture_wrapper(), edges_fixture_wrapper(), title = "t",
    colour_by = "none", theme = networkplot_theme(shapes = c(Kinase = "diamond")),
    pathway_table = pw_fixture_wrapper(), universe = universe_fixture_wrapper()
  )
  expect_false(isTRUE(out$widget$x$byselection$enabled)) # colour_by = "none" disables it
  expect_equal(unique(out$widget$x$nodes$shape[out$widget$x$nodes$type == "Kinase"]), "diamond")
})
