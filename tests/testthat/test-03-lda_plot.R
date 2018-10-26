context("Check LDA plot functions")

data(rodents)
rem <- which(colnames(rodents) %in% c("newmoon", "date", "plots", "traps"))
lda_data <- rodents[ , -rem]
ldas <- LDA_set(lda_data, topics = c(2, 4), nseeds = 2)
lda <- ldas[[1]]
xtime <- rodents$newmoon

test_that("check output from set_LDA_plot_colors", {
  col_default <- set_LDA_plot_colors(x = lda)
  expect_equal(col_default, c("#440154FF", "#FDE725FF"))
  col_A <- set_LDA_plot_colors(x = lda, option = "A")
  expect_equal(col_A, c("#000004FF", "#FCFDBFFF"))
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