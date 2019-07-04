context("Check LDA_TS plot functions")
tenv <- "cran"

data(rodents)
lda_data <- rodents$document_term_table
document_term_table <- rodents$document_term_table
document_covariate_table <- rodents$document_covariate_table
   
mod <- LDA_TS(document_term_table, document_covariate_table,
              topics = 2, nseeds = 1, formulas = ~ 1, nchangepoints = 1,
              weights = document_weights(document_term_table), 
              timename = "newmoon",
              control = LDA_TS_controls_list(
                        TS_control = TS_controls_list(nit = 100, seed = 1)))

test_that("check plot for LDA_TS", {
  if (tenv == "cran"){
    expect_silent(plot(mod, 
                       control = LDA_TS_controls_list(
                                       TS_control = TS_controls_list(
                                                     nit = 100, seed = 1)),
                           interactive = FALSE))
  } else{
    plot(mod, interactive = FALSE)
    LDA_TS_set_plot <- recordPlot()
    vdiffr::expect_doppelganger("Base LDA_TS non-interactive plot", 
                                LDA_TS_set_plot) 
  }
})

test_that("check color list creation function", {
  expect_equal(length(set_LDA_TS_plot_cols()), 2)
  expect_equal(names(set_LDA_TS_plot_cols()), c("LDA", "TS"))
  expect_equal(length(set_LDA_TS_plot_cols()[[1]]), 3)
  expect_equal(length(set_LDA_TS_plot_cols()[[2]]), 2)
  expect_equal(names(set_LDA_TS_plot_cols()[[2]]), c("rho", "gamma"))
  expect_equal(length(set_LDA_TS_plot_cols()[[2]][[1]]), 3)
  expect_equal(names(set_LDA_TS_plot_cols()[[2]][[1]]), 
               c("cols", "option", "alpha"))
  expect_equal(names(set_LDA_TS_plot_cols()[[2]][[2]]), 
               c("cols", "option", "alpha"))
})
