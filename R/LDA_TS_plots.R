#' @title Plot the key results from a full LDATS analysis
#'
#' @description 
#'   \code{plot.LDA_TS} generalizes the \code{\link[graphics]{plot}} function
#'     to work on fitted LDA_TS model objects (class \code{LDA_TS}) returned 
#'     by \code{\link{LDA_TS}}). \cr \cr
#'   \code{set_LDA_TS_plot_colors} produces the options for the colors 
#'     controlling the panels of the LDATS summary plots, needed because 
#'     the change point histogram panel should be in a different color scheme
#'     than the LDA and fitted time series model panels, which should be
#'     in a matching color scheme. See \code{\link{set_LDA_plot_colors}},
#'     \code{\link{set_TS_summary_plot_cols}}, \code{\link{set_gamma_colors}},
#'     and \code{\link{set_rho_hist_colors}} for specific details on usage.
#'
#' @param x A \code{LDA_TS} object of a full LDATS model fit by
#'   \code{\link{LDA_TS}}.
#'
#' @param ... Additional arguments to be passed to subfunctions. Not currently
#'   used, just retained for alignment with \code{plot}.
#'
#' @param bin_width Width of the bins used in the histograms of the summary
#'   time series plot, in units of the time variable used to fit the model
#'   (the x-axis).
#'
#' @param xname Label for the x-axis in the summary time series plot. Defaults
#'   to \code{NULL}, which results in usage of the \code{timename} element
#'   of the control list (held in\code{control$TS_control$timename}). To have
#'   no label printed, set \code{xname = ""}.
#'
#' @param border Border for the histogram, default is \code{NA}.
#'
#' @param selection Indicator of the change points to use in the time series
#'   summary plot. Currently only defined for \code{"median"} and
#'   \code{"mode"}.
#'
#' @param rho_cols,gamma_cols Colors to be used in the specific plot. Any 
#'   valid color values (\emph{e.g.}, see
#'   \code{\link[grDevices]{colors}}, \code{\link[grDevices]{rgb}}) can be 
#'   input as with a standard plot. The default (\code{NULL}) triggers use 
#'   of \code{\link[viridis]{viridis}} color options (see
#'   \code{rho_option},\code{gamma_option}).
#'
#' @param rho_option,gamma_option A \code{character} string indicating
#'   the color option from \code{\link[viridis]{viridis}} to use if
#'   "cols == NULL". Four options are available: "magma" (or "A"), 
#'   "inferno" (or "B"), "plasma" (or "C"), "viridis" (or "D", the default 
#'   option) and "cividis" (or "E").
#'
#' @param rho_alpha,gamma_alpha Numeric value [0,1] that indicates the 
#'   transparency of the colors used. Supported only on some devices, see 
#'   \code{\link[grDevices]{rgb}}.
#'
#' @param cols \code{list} of elements used to define the colors for the two
#'   panels of the summary plot, as generated simply using
#'   \code{\link{set_LDA_TS_plot_cols}}. \code{cols} has two elements:
#'   \code{LDA} and \code{TS}, each corresponding the set of plots for
#'   its stage in the full model. \code{LDA} contains entries \code{cols}
#'   and \code{option} (see \code{\link{set_LDA_plot_colors}}). \code{TS}
#'   contains two entries, \code{rho} and \code{gamma}, each corresponding
#'   to the related panel, and each containing default values for entries
#'   named \code{cols}, \code{option}, and \code{alpha} (see
#'   \code{\link{set_TS_summary_plot_cols}}, \code{\link{set_gamma_colors}},
#'   and \code{\link{set_rho_hist_colors}}).
#' 
#' @return 
#'   \code{plot.LDA_TS}: \code{NULL}. \cr \cr
#'   \code{set_LDA_TS_plot_colors}: \code{list} of elements used to define 
#'     the colors for the summary plots, which has two
#'     elements: \code{LDA} and \code{TS}, each corresponding the set of 
#'     plots for its stage in the full model. \code{LDA} contains entries 
#'     \code{cols} and \code{options} (see \code{\link{set_LDA_plot_colors}}). 
#'     \code{TS} contains two entries, \code{rho} and \code{gamma}, each 
#'     corresponding to the related panel, and each containing default values
#'     for entries named \code{cols}, \code{option}, and \code{alpha} (see
#'     \code{\link{set_TS_summary_plot_cols}}, \code{\link{set_gamma_colors}},
#'     and \code{\link{set_rho_hist_colors}}).
#'
#' @name plot.LDA_TS
#'

#' @rdname plot.LDA_TS
#'
#' @export
#'
plot.LDA_TS <- function(x, ..., cols = set_LDA_TS_plot_colors(),
                        bin_width = 1, xname = NULL, border = NA,
                        selection = "median"){
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  tname <- x$"TS models"$"selected_TSs"[[1]]$timename
  tvar <- x$"TS models"$"selected_TSs"[[1]]$data$train$ts_data[ , tname]
  plot.LDA_set(x$"LDA models", xtime = tvar, xname = NULL,
               cols = cols$LDA$cols, option = cols$LDA$option, LDATS = TRUE)
  plot.TS_set(x$"TS models", plot_type = "summary", cols = cols$TS,
              bin_width = bin_width, xname = xname, border = border,
              selection = selection, LDATS = TRUE)
}

#' @rdname plot.LDA_TS
#'
#' @export
#'
set_LDA_TS_plot_colors <- function(rho_cols = NULL, rho_option = "D",
                                   rho_alpha = 0.4, gamma_cols = NULL,
                                   gamma_option = "C", gamma_alpha = 0.8){
  list(
    LDA = list(cols = gamma_cols, option = gamma_option, alpha = gamma_alpha),
    TS = set_TS_summary_plot_cols(rho_cols = rho_cols,
                                  rho_option = rho_option,
                                  rho_alpha = rho_alpha,
                                  gamma_cols = gamma_cols,
                                  gamma_option = gamma_option,
                                  gamma_alpha = gamma_alpha)
  )
}