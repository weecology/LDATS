context("Check LDA plot functions")

data(rodents)
lda_data <- rodents$document_term_table
ldas <- LDA_set(lda_data, c(2, 4), 2, LDA_controls_list(quiet = TRUE))
lda <- ldas[[1]]
xtime <- rodents$newmoon

test_that("check output from set_LDA_plot_colors", {
  col_default <- set_LDA_plot_colors(x = lda)
  expect_equal(col_default, c("#440154CC", "#BBDF27CC"))
  col_A <- set_LDA_plot_colors(x = lda, option = "A")
  expect_equal(col_A, c("#000004CC", "#FECE91CC"))
  col_grey <- set_LDA_plot_colors(x = lda, cols = "grey")
  expect_equal(col_grey, c("#000000", "#CCCCCC"))
  col_grey_and_A <- set_LDA_plot_colors(x = lda, cols = "grey", option = "A")
  expect_equal(col_grey_and_A, c("#000000", "#CCCCCC"))
  expect_error(set_LDA_plot_colors(x = lda, 1))
})

test_that("check plotting of plot.LDA", {
  plot(x = lda)
  LDA_plot <- recordPlot()
  vdiffr::expect_doppelganger("Base LDA plot", LDA_plot)  
  plot(x = lda, xtime = xtime, xname = "New Moon")
  LDA_plot_xtime <- recordPlot()
  vdiffr::expect_doppelganger("LDA plot with time x", LDA_plot_xtime)  
})

test_that("check plotting of plot.LDA_set", {
  sellda <- select_LDA(ldas)
  plot(x = sellda)
  LDA_set_plot <- recordPlot()
  vdiffr::expect_doppelganger("Base LDA_set selected plot", LDA_set_plot)  
})

test_that("check LDA plot panels", {
  cols <- set_LDA_plot_colors(x = lda, cols = NULL, option = "E")
  LDA_plot_top_panel(x = lda, cols)
  LDA_top_plot <- recordPlot()
  vdiffr::expect_doppelganger("Base LDA plot top panel", LDA_top_plot)  
  LDA_plot_bottom_panel(x = lda, xtime = NULL, xname = NULL, cols)
  LDA_bottom_plot <- recordPlot()
  vdiffr::expect_doppelganger("Base LDA plot bottom panel", LDA_bottom_plot)  
})


