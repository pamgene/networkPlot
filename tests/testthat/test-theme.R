test_that("networkplot_theme() carries the standard defaults", {
  th <- networkplot_theme()
  expect_s3_class(th, "networkplot_theme")
  expect_equal(th$shapes[["Kinase"]], "dot")
  expect_equal(th$shapes[["Sensitivity"]], "square")
  expect_equal(th$shapes[["Hidden"]], "text")
  expect_equal(th$shape_default, "ellipse")
  expect_equal(th$edge_colour, "#000000")
  expect_null(th$edge_width)
  expect_equal(th$highlight_degree, 5)
})

test_that("shapes override merges over the default map, not replaces it", {
  th <- networkplot_theme(shapes = c(Kinase = "diamond"))
  expect_equal(th$shapes[["Kinase"]], "diamond")   # overridden
  expect_equal(th$shapes[["Sensitivity"]], "square") # untouched
  expect_true("Hidden" %in% names(th$shapes))
})

test_that("edge_width is normalised: defaults column/scale/invert", {
  th <- networkplot_theme(edge_width = list(column = "cost", invert = TRUE))
  expect_equal(th$edge_width$column, "cost")
  expect_true(th$edge_width$invert)
  expect_type(th$edge_width$scale, "closure")

  th2 <- networkplot_theme(edge_width = list())
  expect_equal(th2$edge_width$column, "cost") # default
  expect_false(th2$edge_width$invert)
})

test_that("default logfc_limit_fn is the symmetric p25/p75 rule", {
  th <- networkplot_theme()
  x <- c(-4, -3, -2, -1, 1, 2, 3, 8)
  lim <- th$logfc_limit_fn(x)
  # |p25 of negatives| = 3.25; p75 of positives = 4.25; max = 4.25
  expect_equal(round(lim, 2), 4.25)
})
