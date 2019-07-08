context("Check LDA plot functions")
tenv <- "cran"

data(rodents)
lda_data <- rodents$document_term_table
ldas <- LDA_set(lda_data, c(2, 4), 2, list(quiet = TRUE))
lda <- ldas[[1]]
xtime <- rodents$newmoon

test_that("check output from set_LDA_plot_colors", {
  col_default <- set_LDA_plot_colors(x = lda)
  expect_equal(col_default, c("#0D0887CC", "#FCCE25CC"))
  col_A <- set_LDA_plot_colors(x = lda, option = "A")
  expect_equal(col_A, c("#000004CC", "#FECE91CC"))
  col_grey <- set_LDA_plot_colors(x = lda, cols = "grey")
  expect_equal(col_grey, c("#000000", "#CCCCCC"))
  col_grey_and_A <- set_LDA_plot_colors(x = lda, cols = "grey", option = "A")
  expect_equal(col_grey_and_A, c("#000000", "#CCCCCC"))
  expect_error(set_LDA_plot_colors(x = lda, 1))
  col_extra <- set_LDA_plot_colors(x = lda, cols = 1:3)
  expect_equal(col_extra, 1:2)
})

test_that("check plotting of plot.LDA", {
  if (tenv == "cran"){
    expect_silent(plot(x = lda))
    expect_silent(plot(x = lda, xtime = xtime, xname = "New Moon"))
  } else{
    plot(x = lda)
    LDA_plot <- recordPlot()
    vdiffr::expect_doppelganger("Base LDA plot", LDA_plot)  
    plot(x = lda, xtime = xtime, xname = "New Moon")
    LDA_plot_xtime <- recordPlot()
    vdiffr::expect_doppelganger("LDA plot with time x", LDA_plot_xtime)  
  }
})

test_that("check plotting of plot.LDA_set", {
  sellda <- select_LDA(ldas)
  if (tenv == "cran"){
    expect_silent(plot(x = sellda))
  } else{
    plot(x = sellda)
    LDA_set_plot <- recordPlot()
    vdiffr::expect_doppelganger("Base LDA_set selected plot", LDA_set_plot)  
  }
})

test_that("check LDA plot panels", {
  cols <- set_LDA_plot_colors(x = lda, cols = NULL, option = "E")

  if (tenv == "cran"){
    expect_silent(LDA_plot_top_panel(x = lda, cols))
    expect_silent(LDA_plot_bottom_panel(x = lda, xtime = NULL, 
                                        xname = NULL, cols))
  } else{
    LDA_plot_top_panel(x = lda, cols)
    LDA_top_plot <- recordPlot()
    vdiffr::expect_doppelganger("Base LDA plot top panel", LDA_top_plot)  
    LDA_plot_bottom_panel(x = lda, xtime = NULL, xname = NULL, cols)
    LDA_bottom_plot <- recordPlot()
    vdiffr::expect_doppelganger("Base LDA plot bottom panel", LDA_bottom_plot)  
  }
})


