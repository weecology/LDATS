context("Check LDA plot functions")

data(rodents)
rem <- which(colnames(rodents) %in% c("newmoon", "date", "plots", "traps"))
lda_data <- rodents[ , -rem]
lda <- LDA_set(lda_data, c(2, 4), nseeds = 2)

test_that("check output from set_LDA_plot_colors", {
  mod <- lda[[1]]
  col_default <- set_LDA_plot_colors(mod)
  expect_equal(col_default, c("#440154FF", "#FDE725FF"))
  col_A <- set_LDA_plot_colors(mod, option = "A")
  expect_equal(col_A, c("#000004FF", "#FCFDBFFF"))
  col_grey <- set_LDA_plot_colors(mod, cols = "grey")
  expect_equal(col_grey, c("#000000", "#CCCCCC"))
  col_grey_and_A <- set_LDA_plot_colors(mod, cols = "grey", option = "A")
  expect_equal(col_grey_and_A, c("#000000", "#CCCCCC"))
  expect_error(set_LDA_plot_colors(mod, 1))
})