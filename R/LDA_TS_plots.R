#' @title Plot a LDATS full model
#'
#' @description Generalization of the \code{plot} function to work on fitted
#'   LDA_TS model (class \code{LDA_TS}). 
#' 
#' @param x A \code{LDA_TS} object of a full LDATS model fit by
#'   \code{LDA_TS}.
#' 
#' @param ... Additional arguments to be passed to subfunctions. Not currently
#'   used, just retained for alignment with \code{plot}.
#'
#' @param control Class \code{LDA_TS_controls} list that contains 
#'   \code{LDA_controls}, \code{TS_controls}, and the top-level \code{quiet}.
#'
#' @param bin_width Width of the bins used in the histograms of the summary 
#'   timeseries plot, in units of the x-axis (the time variable used to fit 
#'   the model).
#'
#' @param xlab Label for the x-axis in the summary time series plot.
#'
#' @param selection Indicator of the changepoints to use in the timeseries
#'   summary plot. Currently only defined for "median" and "mode".
#'
#' @param cols \code{list} of elements used to define the colors for the two
#'   panels of the summary plot, as generated simply using 
#'   \code{LDA_TS_summary_cols}. \code{cols} has two elements \code{LDA} and
#'   \code{TS}, each corresponding to its stage in the full model. 
#'   \code{LDA} contains entries \code{cols} and \code{options} (see
#'   \code{set_LDA_plot_colors} for details on usage. \code{TS} contains two
#'   entries, \code{rho} and \code{gamma}, each corresponding to the related 
#'   panel, and each containing default values for entries named \code{cols}, 
#'   \code{option}, and \code{alpha}. See \code{set_gamma_colors} and 
#'   \code{set_rho_hist_colors} for details on usage.
#' 
#' @param interactive \code{logical} indicator, should be set to \code{TRUE}
#'   currently (except for testing).
#'
#' @return Nothing. Plots are generated in the active graphics device.
#' 
#' @examples
#' \dontrun{
#'   data(rodents)
#'   lda_data <- rodents$document_term_table
#'   document_term_table <- rodents$document_term_table
#'   document_covariate_table <- rodents$document_covariate_table
#'   
#'   mod <- LDA_TS(document_term_table, document_covariate_table,
#'                 topics = 2, nseeds = 1, formulas = ~ 1, nchangepoints = 1,
#'                 weights = document_weights(document_term_table), 
#'                 control = LDA_TS_controls_list())
#'   plot(mod)
#' }
#'
#' @export 
#'
plot.LDA_TS <- function(x, ..., control = LDA_TS_controls_list(),
                        cols = LDA_TS_summary_cols(),
                        bin_width = 1, xlab = NULL, selection = "median",
                        interactive = TRUE){
  tname <- control$TS_control$timename
  tvar <- x$"Selected TS model"$data[ , tname]
  if(!is.null(xlab)){
    tname <- xlab
  }
  if (interactive){
    devAskNewPage(TRUE)
  }
  plot.LDA_set(x$"Selected LDA model", xtime = tvar, xname = tname, 
               cols = cols$LDA$cols, option = cols$LDA$option)
  if (interactive){
    devAskNewPage(TRUE)
  }
  plot.TS_fit(x$"Selected TS model", plot_type = "summary", cols = cols$TS,
              bin_width = bin_width, xlab = xlab, selection = selection)
  devAskNewPage(FALSE)
}



#' @title Create the list of colors for the LDATS summary plot
#'
#' @description A default list generator function that produces the options
#'   for the colors controlling the panels of the LDATS summary plots, needed
#'   because the panels should be in different color schemes. See 
#'   \code{set_LDA_plot_colors}, \code{set_gamma_colors}, and 
#'   \code{set_rho_hist_colors} for specific details on usage.
#'
#' @param rho_cols Colors to be used to plot the histograms of changepoints.
#'
#' @param rho_option A character string indicating the colormap option to use
#'   if `rho_cols == NULL`. Four options are available: "magma" (or "A"), 
#'   "inferno" (or "B"), "plasma" (or "C"), "viridis" (or "D", the default
#'   option) and "cividis" (or "E").
#'
#' @param rho_alpha Numeric value [0,1] that indicates the transparency of the 
#'   colors used. Supported only on some devices, see \code{rgb}.
#'
#' @param gamma_cols Colors to be used to plot the LDA topic proportions,
#'   time series of observed topic proportions, and time series of fitted 
#'   topic proportions.
#'
#' @param gamma_option A character string indicating the colormap option to 
#'   use if `gamma_cols == NULL`. Four options are available: "magma" (or 
#'   "A"), "inferno" (or "B", the default option), "plasma" (or "C"), 
#'   "viridis" (or "D") and "cividis" (or "E")
#'
#' @param gamma_alpha Numeric value [0,1] that indicates the transparency of 
#'   the colors used. Supported only on some devices, see \code{rgb}.
#'
#' @return cols \code{list} of elements used to define the colors for the two
#'   panels of the summary plot, as generated simply using 
#'   \code{LDA_TS_summary_cols}. \code{cols} has two elements \code{LDA} and
#'   \code{TS}, each corresponding to its stage in the full model. 
#'   \code{LDA} contains entries \code{cols} and \code{options} (see
#'   \code{set_LDA_plot_colors} for details on usage. \code{TS} contains two
#'   entries, \code{rho} and \code{gamma}, each corresponding to the related 
#'   panel, and each containing default values for entries named \code{cols}, 
#'   \code{option}, and \code{alpha}. See \code{set_gamma_colors} and 
#'   \code{set_rho_hist_colors} for details on usage.
#'
#' @export
#'
LDA_TS_summary_cols <- function(rho_cols = NULL, rho_option = "D", 
                                rho_alpha = 0.4, gamma_cols = NULL, 
                                gamma_option = "C", gamma_alpha = 0.8){
  list(
    LDA = list(cols = gamma_cols, option = gamma_option, alpha = gamma_alpha),
    TS = TS_summary_cols(rho_cols = rho_cols, rho_option = rho_option, 
              rho_alpha = rho_alpha, gamma_cols = gamma_cols, 
              gamma_option = gamma_option, gamma_alpha = gamma_alpha)
  )
}
