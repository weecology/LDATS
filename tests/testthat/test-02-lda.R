context("Check LDA")

data(rodents)
rem <- which(colnames(rodents) %in% c("newmoon", "date", "plots", "traps"))
lda_data <- rodents[ , -rem]
lda <- LDA(lda_data, 2)

test_that("check LDA plot", {
  test_rodent_LDA_plot <- plot(lda)
  expect_doppleganger("rodent LDA plot", test_rodent_LDA_plot)
})