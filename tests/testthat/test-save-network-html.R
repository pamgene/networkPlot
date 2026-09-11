test_that("save_network_html() writes into a Network_html subfolder of out_dir", {
  d <- file.path(tempdir(), "save_html_subfolder_test")
  unlink(d, recursive = TRUE)

  n <- data.frame(Protein = c("A", "B"), type = "Kinase", prize = 1,
                  LogFC = c(1, -1), degree = c(2, 1))
  e <- data.frame(from = "A", to = "B", cost = 0.1, weight = 3)
  w <- plot_network(n, e, title = "t")

  out <- save_network_html(w, d, "toy.html")

  expect_equal(out, file.path(d, "Network_html", "toy.html"))
  expect_true(file.exists(out))
  # nothing was written directly into d itself -- only inside Network_html/
  expect_setequal(list.files(d), "Network_html")
})
