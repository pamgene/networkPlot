make_nodes <- function() {
  data.frame(
    Protein = c("A", "B", "C", "S1"),
    type = c("Kinase", "Kinase", "Hidden", "Sensitivity"),
    prize = c(0.9, 0.8, 2, 0.5),
    LogFC = c(2, -3, 0, 1),
    degree = c(6, 1, 4, 2),
    stringsAsFactors = FALSE
  )
}
make_edges <- function() {
  data.frame(from = c("A", "B", "C"), to = c("B", "C", "S1"),
             cost = c(0.1, 0.4, 0.2), weight = c(8, 2, 5), stringsAsFactors = FALSE)
}

test_that("plot_network() returns a visNetwork widget with type-mapped shapes and degree borders", {
  w <- plot_network(make_nodes(), make_edges(), title = "t", theme = networkplot_theme())
  expect_s3_class(w, "visNetwork")

  n <- w$x$nodes
  expect_equal(n$shape[n$id == "A"], "dot")       # Kinase
  expect_equal(n$shape[n$id == "S1"], "square")   # Sensitivity
  expect_equal(n$shape[n$id == "C"], "text")      # Hidden

  # degree >= 5 gets the highlight border
  expect_equal(n$color.border[n$id == "A"], "#4C9900")
  expect_equal(n$color.border[n$id == "B"], "#454545")
  expect_true(all(!is.na(n$color.background)))
})

test_that("theme shapes override changes only the named type", {
  w <- plot_network(make_nodes(), make_edges(), title = "t",
                    theme = networkplot_theme(shapes = c(Kinase = "diamond")))
  n <- w$x$nodes
  expect_equal(unique(n$shape[n$type == "Kinase"]), "diamond")
  expect_equal(n$shape[n$id == "S1"], "square")
})

test_that("edge_width maps a column onto width, with inversion; missing column warns", {
  w <- plot_network(make_nodes(), make_edges(), title = "t",
                    theme = networkplot_theme(edge_width = list(column = "cost", invert = TRUE)))
  e <- w$x$edges
  # lowest cost (strongest) -> widest
  expect_gt(e$width[e$from == "A"], e$width[e$from == "B"])

  expect_warning(
    plot_network(make_nodes(), make_edges(), title = "t",
                 theme = networkplot_theme(edge_width = list(column = "nope"))),
    "not in .edges."
  )
})

test_that("edge_width = NULL leaves edges without a width column (constant)", {
  w <- plot_network(make_nodes(), make_edges(), title = "t", theme = networkplot_theme())
  expect_false("width" %in% colnames(w$x$edges))
})

test_that("missing required node columns is an error", {
  bad <- make_nodes()[, c("Protein", "type")]
  expect_error(plot_network(bad, make_edges(), title = "t"), "missing required column")
})

test_that("pathway multi-select is wired only when a pathway column is present", {
  n <- make_nodes()
  n$pathway <- c("PW1, PW2", NA, NA, "PW1")
  w <- plot_network(n, make_edges(), title = "t", colour_by = "pathway")
  expect_equal(w$x$byselection$variable, "pathway")
  expect_true(isTRUE(w$x$byselection$enabled))

  w2 <- plot_network(make_nodes(), make_edges(), title = "t", colour_by = "pathway")
  expect_false(isTRUE(w2$x$byselection$enabled))
})
