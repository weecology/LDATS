context("Check TS plot functions")
tenv <- "cran"

data(rodents)
lda_data <- rodents$document_term_table
document_term_table <- rodents$document_term_table
document_covariate_table <- rodents$document_covariate_table
topics <- 2
nseeds <- 1
formulas <- ~ 1
nchangepoints <- 1
weights <- document_weights(document_term_table)
control <- LDA_TS_controls_list()
LDAs <- LDA_set(document_term_table, topics, nseeds, control$LDA_control)
LDA_models <- select_LDA(LDAs, control$LDA_control)
control <- TS_controls_list(nit = 1e3, seed = 1)
timename <- "newmoon"
mods <- expand_TS(LDA_models, formulas, nchangepoints)
formula <- mods$formula[[1]]
nchangepoints <- mods$nchangepoints[1]
data <- prep_TS_data(document_covariate_table, LDA_models, mods, 1)
TSmod <- TS(data, formula, nchangepoints, weights, timename, control)


test_that("check rho_hist color generator", {
  rc <- set_TS_summary_plot_cols()$rho
  rho_cols <- set_rho_hist_colors(TSmod$rhos, rc$cols, rc$option, rc$alpha)
  expect_equal(rho_cols, "#44015466")
})

test_that("check pred_gamma color generator", {
  gc <- set_TS_summary_plot_cols()$gamma
  gamma_cols <- set_gamma_colors(TSmod, gc$cols, gc$option, gc$alpha)
  expect_equal(gamma_cols, c("#0D0887CC", "#FCCE25CC"))
})

test_that("check pred_gamma plot", {
  gc <- set_TS_summary_plot_cols()$gamma
  gamma_cols <- set_gamma_colors(TSmod, gc$cols, gc$option, gc$alpha)

  if (tenv == "cran"){
    expect_silent(pred_gamma_TS_plot(TSmod, cols = gamma_cols))

    expect_silent(pred_gamma_TS_plot(TSmod, selection = "mode", 
                  cols = gamma_cols))
  } else{
    TS_gamma_plot <- pred_gamma_TS_plot(TSmod, cols = gamma_cols)
    TS_gamma_plot <- recordPlot()
    vdiffr::expect_doppelganger("Base TS gamma plot", TS_gamma_plot)  
  }
  expect_equal(set_gamma_colors(NULL), NULL)
    expect_error(pred_gamma_TS_plot(TSmod, selection = "ok", 
                  cols = gamma_cols))
})

test_that("check rho_lines", {
  if (tenv == "cran"){
    expect_silent(plot(1, 1, xlim = c(-10, 10), ylim = c(0, 1)))
    expect_silent(rho_lines(1))
  } else{
    plot(1, 1, xlim = c(-10, 10), ylim = c(0, 1))
    rho_lines(1)
    TS_rho_line_plot <- recordPlot()
    vdiffr::expect_doppelganger("rho line plot", TS_rho_line_plot)
  }
  expect_equal(rho_lines(NULL), NULL)
})


test_that("check rho_hist plot", {
  rc <- set_TS_summary_plot_cols()$rho
  rho_cols <- set_rho_hist_colors(TSmod$rhos, rc$cols, rc$option, rc$alpha)
  if (tenv == "cran"){
    expect_silent(rho_hist(TSmod, rho_cols))
  } else{
    TS_rho_plot <- rho_hist(TSmod, rho_cols)
    TS_rho_plot <- recordPlot()
    vdiffr::expect_doppelganger("Base TS rho plot", TS_rho_plot)  
  }
  expect_equal(set_rho_hist_colors(NULL), NULL)
})


test_that("check color list creation function", {
  expect_equal(length(set_TS_summary_plot_cols()), 2)
  expect_equal(names(set_TS_summary_plot_cols()), c("rho", "gamma"))
  expect_equal(length(set_TS_summary_plot_cols()[[1]]), 3)
  expect_equal(length(set_TS_summary_plot_cols()[[2]]), 3)
  expect_equal(names(set_TS_summary_plot_cols()[[2]]), 
               c("cols", "option", "alpha"))
  expect_equal(names(set_TS_summary_plot_cols()[[1]]), 
               c("cols", "option", "alpha"))
})




test_that("check trace_plot", {
  if (tenv == "cran"){
    expect_silent(trace_plot(TSmod$rhos[ , 1]))
  } else{
    TS_trace_plot <- trace_plot(TSmod$rhos[ , 1])
    TS_trace_plot <- recordPlot()
    vdiffr::expect_doppelganger("Base TS trace plot", TS_trace_plot)  
  }
})

test_that("check ecdf_plot", {
  if (tenv == "cran"){
    expect_silent(ecdf_plot(TSmod$rhos[ , 1]))
  } else{
    TS_ecdf_plot <- ecdf_plot(TSmod$rhos[ , 1])
    TS_ecdf_plot <- recordPlot()
    vdiffr::expect_doppelganger("Base TS ecdf plot", TS_ecdf_plot)  
  }
})

test_that("check autocorr_plot", {
  if (tenv == "cran"){
    expect_silent(autocorr_plot(TSmod$rhos[ , 1]))
  } else{
    TS_autocorr_plot <- autocorr_plot(TSmod$rhos[ , 1])
    TS_autocorr_plot <- recordPlot()
    vdiffr::expect_doppelganger("Base TS autocorr plot", TS_autocorr_plot)  
  }
})

test_that("check posterior_plot", {
  if (tenv == "cran"){
    expect_silent(posterior_plot(TSmod$rhos[ , 1]))
  } else{
    TS_posterior_plot <- posterior_plot(TSmod$rhos[ , 1])
    TS_posterior_plot <- recordPlot()
    vdiffr::expect_doppelganger("Base TS posterior plot", TS_posterior_plot)  
  }
})


test_that("check plotting of TS_fit", {
  if (tenv == "cran"){
    expect_silent(plot(TSmod))
    expect_silent(plot(TSmod, plot_type = "diagnostic", interactive = FALSE))
  } else{
    plot(TSmod)
    TS_plot <- recordPlot()
    vdiffr::expect_doppelganger("Base TS plot", TS_plot)  
  }
})

test_that("check TS_diagnostics_plot", {
  TSmod0 <- TS(data, formula, nchangepoints = 0, weights, timename, control)
  if (tenv == "cran"){
    expect_silent(TS_diagnostics_plot(TSmod0, interactive = FALSE))
  } else{
    TS_diagnostics_plot(TSmod0, interactive = FALSE)
    TS_diag_plot <- recordPlot()
    vdiffr::expect_doppelganger("Base TS diagnostic plot", TS_diag_plot)  
  }
})


test_that("check TS_summary_plot", {
  if (tenv == "cran"){
    expect_silent(TS_summary_plot(TSmod, cols = set_TS_summary_plot_cols(),
                        bin_width = 1, xlab = NULL, selection = "median"))
  } else{
    TS_summary_plot(TSmod, cols = set_TS_summary_plot_cols(),
                          bin_width = 1, xlab = NULL, selection = "median")
    TS_summ_plot <- recordPlot()
    vdiffr::expect_doppelganger("Base TS summary plot", TS_summ_plot)  
  }
})


test_that("check rho_diagnostics_plots", {
  if (tenv == "cran"){
    expect_silent(rho_diagnostics_plots(TSmod, interactive = FALSE))
  } else{
    rho_diagnostics_plots(TSmod, interactive = FALSE)
    TS_rdiag_plot <- recordPlot()
    vdiffr::expect_doppelganger("Base TS rho diagnostic plot", TS_rdiag_plot)  
  }
})

test_that("check eta_diagnostics_plots", {
  expect_equal(eta_diagnostics_plots(NULL), NULL)
  if (tenv == "cran"){
    expect_silent(eta_diagnostics_plots(TSmod, interactive = FALSE))
  } else{
    eta_diagnostics_plots(TSmod, interactive = FALSE)
    TS_ediag_plot <- recordPlot()
    vdiffr::expect_doppelganger("Base TS eta diagnostic plot", TS_ediag_plot)  
  }
})
