#' @title Plot an LDATS Time Series model
#'
#' @description 
#'   \code{plot.TS} is a generalization of the \code{\link[graphics]{plot}}
#'     function to work on fitted TS model objects (class \code{TS}) 
#'     returned from \code{\link{TS}}. \cr \cr
#'   \code{plot.TS_set} plots a \code{TS_set} of \code{TS} models, either 
#'     just the \code{selected} models or all. \cr \cr
#'   \code{TS_diagnostics_plot} makes the 4-panel figures (showing trace 
#'     plots, posterior ECDF, posterior density, and iteration 
#'     autocorrelation) for each of the parameters (change point locations 
#'     and regressors) fitted within a compositional time series model (fit
#'     by \code{\link{TS}}). \cr \cr
#'   \code{eta_diagnostics_plots} creates the diagnostic plots
#'     for the regressors (etas) of a time series model. \cr \cr
#'   \code{rho_diagnostics_plots} creates the diagnostic plots
#'     for the change point locations (rho) of a time series model. \cr \cr
#'   \code{trace_plot} produces a trace plot for the parameter of interest 
#'     (rho or eta) as part of \code{\link{TS_diagnostics_plot}}. A 
#'     horizontal line is added to show the median of the posterior. \cr \cr
#'   \code{ecdf_plot} makes a vanilla ECDF (empirical cumulative distribution
#'     function) plot using \code{\link[stats]{ecdf}} for the parameter of 
#'     interest (rho or eta) as part of \code{\link{TS_diagnostics_plot}}.
#'     A horizontal line is added to show the median of the posterior. \cr \cr
#'   \code{autocorr_plot} produces a vanilla ACF plot using 
#'     \code{\link[stats]{acf}} for the parameter of interest (rho or eta)
#'     as part of \code{\link{TS_diagnostics_plot}}.\cr \cr
#'   \code{posterior_plot} makes a vanilla histogram plot using 
#'     \code{\link[graphics]{hist}} for the parameter of interest (rho or eta)
#'     as part of \code{\link{TS_diagnostics_plot}}. A vertical line is added 
#'     to show the median of the posterior. \cr \cr
#'   \code{TS_summary_plot} produces a two-panel figure of [1] the change 
#'     point distributions as histograms over time and [2] the time series of 
#'     the fitted topic proportions over time, based on a selected set of 
#'     change point locations. \cr \cr
#'   \code{pred_gamma_TS_plot} produces a time series of the 
#'     fitted topic proportions over time, based on a selected set of change 
#'     point locations. \cr \cr
#'   \code{rho_hist}: make a plot of the change point distributions as 
#'     histograms over time. \cr \cr
#'   \code{rho_lines} adds vertical lines to the plot of the time series of 
#'     fitted proportions associated with the change points of interest.
#'     \cr \cr
#'   \code{set_gamma_colors} creates the set of colors to be used in
#'     the time series of the fitted gamma (topic proportion) values. \cr \cr
#'   \code{set_rho_hist_colors} creates the set of colors to be used in
#'     the change point histogram. \cr \cr
#'   \code{set_TS_summary_plot_cols} acts as a default \code{list} 
#'     generator function that produces the options for the colors 
#'     controlling the panels of the TS summary plots, so needed
#'     because the panels should be in different color schemes. See 
#'     \code{\link{set_gamma_colors}} and \code{\link{set_rho_hist_colors}} 
#'     for specific details on usage.
#' 
#' @param x In \code{plot.TS}, a \code{TS_fit} object of a multinomial time
#'   series model fit by \code{\link{TS}}. In \code{plot.TS_set}, a 
#'   \code{TS_set} \code{list} of \code{TS} objects.
#' 
#' @param ... Additional arguments to be passed to subfunctions. Not currently
#'   used, just retained for alignment with \code{\link[graphics]{plot}}.
#'
#' @param selected \code{logical} indicator of if only the selected TSs
#'   (the first element in \code{x}) should be plotted or if all the TSs
#'   (the second element in \code{x}) should be plotted.
#'
#' @param spec_rhos \code{numeric} vector indicating the locations along the
#'   x axis where the specific change points being used are located.
#'
#' @param plot_type "diagnostic" or "summary".
#'
#' @param bin_width Width of the bins used in the histograms of the summary 
#'   time series plot, in units of the x-axis (the time variable used to fit 
#'   the model).
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
#' @param cols,rho_cols,gamma_cols 
#'   In \code{plot.TS}, \code{cols} is a \code{list} of elements used to
#'     define the colors for the two panels of the summary plot, as generated 
#'     simply using \code{\link{set_TS_summary_plot_cols}}. 
#'     \code{cols} has two elements \code{rho} and \code{gamma}, each
#'     corresponding to the related panel, and each containing default values
#'     for entries named \code{cols}, \code{option}, and \code{alpha}. \cr
#'   For \code{rho_cols} and \code{gamma_cols} always and for \code{cols} in 
#'     \code{set_rho_hist_colors}, \code{set_gamma_colors}, 
#'     \code{rho_hist}, and \code{pred_gamma_TS_plot}, colors to be used in 
#'     the specific plot. Any valid color values (\emph{e.g.}, see
#'     \code{\link[grDevices]{colors}}, \code{\link[grDevices]{rgb}}) can be 
#'     input as with a standard plot. The default (\code{NULL}) triggers use 
#'     of \code{\link[viridis]{viridis}} color options (see 
#'     \code{option},\code{rho_option},\code{gamma_option}).
#'
#' @param LDATS \code{logical} indicating if the plot is part of a larger 
#'   LDATS plot output.
#'
#' @param interactive \code{logical} input, should be \code{TRUE} unless
#'   testing.
#'
#' @param ylab \code{character} value used to label the y axis.
#'
#' @param draw \code{vector} of parameter values drawn from the posterior 
#'   distribution, indexed to the iteration by the order of the vector.
#'
#' @param xlab \code{character} value used to label the x axis.
#' 
#' @param option,rho_option,gamma_option A \code{character} string indicating
#'   the color option from \code{\link[viridis]{viridis}} to use if
#'   "cols == NULL". Four options are available: "magma" (or "A"), 
#'   "inferno" (or "B"), "plasma" (or "C"), "viridis" (or "D", the default 
#'   option) and "cividis" (or "E").
#'
#' @param alpha,rho_alpha,gamma_alpha Numeric value [0,1] that indicates the 
#'   transparency of the colors used. Supported only on some devices, see 
#'   \code{\link[grDevices]{rgb}}.
#' 
#' @return 
#'   \code{plot.TS},\code{plot.TS_set},\code{TS_diagnostics_plot},
#'     \code{eta_diagnostics_plots},\code{rho_diagnostics_plots},
#'     \code{trace_plot},\code{posterior_plot},\code{autocorr_plot},
#'     \code{ecdf_plot},\code{TS_summary_plot},\code{pred_gamma_TS_plot},
#'     \code{rho_hist},\code{rho_lines}:\code{NULL}. \cr \cr
#'   \code{set_rho_hist_cols},\code{set_gamma_colors}: \code{vector} of 
#'     \code{character} hex codes indicating colors to use.
#'   \code{set_TS_summary_plot_cols}: \code{list} of elements used to define
#'     the colors for the two panels. Contains two elements \code{rho} and 
#'     \code{gamma}, each corresponding to the related panel, and each 
#'     containing default values for entries named \code{cols}, 
#'     \code{option}, and \code{alpha}. 
#'
#' @name plot.TS 
#'


#' @rdname plot.TS
#'
#' @export 
#'
plot.TS <- function(x, ..., plot_type = "summary", interactive = FALSE,
                        cols = set_TS_summary_plot_cols(), bin_width = 1, 
                        xname = NULL, border = NA, selection = "median", 
                        LDATS = FALSE){
  if (plot_type == "diagnostic"){
    TS_diagnostics_plot(x, interactive = interactive)
  } else if (plot_type == "summary"){
    TS_summary_plot(x, cols, bin_width, xname, border, selection, LDATS)
  }
}

#' @rdname plot.TS
#'
#' @export 
#'
plot.TS_set <- function(x, ..., selected = TRUE){
  if(selected){
    x <- x[[1]]
  } else{
    x <- x[[2]]
  }
  on.exit(devAskNewPage(FALSE))
  if (length(x) > 1){
    devAskNewPage(TRUE)
  }
  y <- lapply(x, plot, ...)
  y <- NULL
}

#' @rdname plot.TS
#'
#' @export 
#'
TS_diagnostics_plot <- function(x, interactive = TRUE){
  rho_diagnostics_plots(x, interactive)
  eta_diagnostics_plots(x, interactive)
}

#' @rdname plot.TS
#'
#' @export 
#'
eta_diagnostics_plots <- function(x, interactive){
  on.exit(devAskNewPage(FALSE))
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  etas <- x$etas
  if (is.null(etas)){
    return()
  }
  netas <- ncol(etas) 
  for (i in 1:netas){
    devAskNewPage(interactive)
    par(mfrow = c(2, 2), mar = c(5, 5, 1, 1))
    chead <- colnames(etas)[i]
    spl1 <- strsplit(chead, "_")[[1]]
    seglab <- paste0("Segment ", spl1[1])
    spl2 <- strsplit(spl1[2], ":")[[1]]
    toplab <- paste0(" Topic ", spl2[1])
    coflab <- gsub("\\(Intercept\\)", "Intercept", spl2[2])
    lab <- paste0(seglab, toplab, " ", coflab)
    trace_plot(etas[ , i], lab)
    ecdf_plot(etas[ , i], lab)
    posterior_plot(etas[ , i], lab)
    autocorr_plot(etas[ , i])
  }  
}

#' @rdname plot.TS
#'
#' @export 
#'
rho_diagnostics_plots <- function(x, interactive){
  on.exit(devAskNewPage(FALSE))
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  rhos <- x$focal_rhos
  if (is.null(rhos)){
    return()
  }
  nrhos <- ncol(rhos) 
  for (i in 1:nrhos){
    devAskNewPage(interactive)
    par(mfrow = c(2, 2), mar = c(5, 5, 1, 1))
    lab <- paste0("Change point ", i, " location")
    trace_plot(rhos[ , i], lab)
    ecdf_plot(rhos[ , i], lab)
    posterior_plot(rhos[ , i], lab)
    autocorr_plot(rhos[ , i])
  }
}



#' @rdname plot.TS
#'
#' @export 
#'
trace_plot <- function(draw, ylab = "parameter value"){
  plot(draw, type = "l", lwd = 1, col = 0,
       xlab = "Iteration", ylab = ylab, las = 1, bty = "L")
  ext <- 0.01 * length(draw)
  points(c(-ext, length(draw) + ext), rep(median(draw), 2), 
        type = "l", lwd = 2, lty = 2)
  points(seq_along(draw), draw, type = "l", col = rgb(0.4, 0.4, 0.4, 
         alpha = 0.9))
}


#' @rdname plot.TS
#'
#' @export 
#'
ecdf_plot <- function(draw, xlab = "parameter value"){
  ECDF <- ecdf(draw)
  plot(ECDF, main = "", xlab = xlab, ylab = "%", las = 1, bty = "L")
  abline(a = 0.5, b = 0, lwd = 2, lty = 2)
}



#' @rdname plot.TS
#'
#' @export 
#'
posterior_plot <- function(draw, xlab = "parameter value"){
  hist(draw, las = 1, main = "", xlab = xlab)
  points(rep(median(draw), 2), c(0, 1e5), type = "l", lwd = 2, lty = 2)
}



#' @rdname plot.TS
#'
#' @export 
#'
autocorr_plot <- function(draw){
  acf(draw, las = 1, ylab = "Autocorrelation")
}



#' @rdname plot.TS
#'
#' @export 
#'
set_TS_summary_plot_cols <- function(rho_cols = NULL, rho_option = "D", 
                                     rho_alpha = 0.4, gamma_cols = NULL, 
                                     gamma_option = "C", gamma_alpha = 0.8){
  list(
    rho = list(cols = rho_cols, option = rho_option, alpha = rho_alpha),
    gamma = list(cols = gamma_cols, option = gamma_option, 
                 alpha = gamma_alpha)
  )
}



#' @rdname plot.TS
#'
#' @export 
#'
TS_summary_plot <- function(x, cols = set_TS_summary_plot_cols(), 
                            bin_width = 1, xname = NULL, border = NA, 
                            selection = "median", LDATS = FALSE){
  rho_cols <- set_rho_hist_colors(x = x, cols = cols$rho$cols, 
                                  option = cols$rho$option, 
                                  alpha = cols$rho$alpha)
  rho_hist(x = x, cols = rho_cols, bin_width = bin_width, xname = xname, 
           border = border, together = TRUE, LDATS = LDATS)
  gamma_cols <- set_gamma_colors(x = x, cols = cols$gamma$cols, 
                                 option = cols$gamma$option, 
                                 alpha = cols$gamma$alpha)
  pred_gamma_TS_plot(x = x, selection = selection, cols = gamma_cols, 
                     xname = xname, together = TRUE, LDATS = LDATS)
}

#' @rdname plot.TS
#'
#' @export 
#'
pred_gamma_TS_plot <- function(x, selection = "median", 
                               cols = set_gamma_colors(x),
                               xname = NULL, together = FALSE, LDATS = FALSE){
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  if(LDATS){
    par(fig = c(0, 1, 0, 0.3), new = TRUE)
  } else if(together){
    par(fig = c(0, 1, 0, 0.52), new = TRUE)
  } else{
    par(fig = c(0, 1, 0, 1))
  }    
  rhos <- x$focal_rhos
  nrhos <- ncol(rhos)
  if (!is.null(nrhos)){
    if (selection == "median"){
      spec_rhos <- apply(rhos, 2, median)
    } else if (selection == "mode"){
      spec_rhos <- apply(rhos, 2, modalvalue)
    } else {
      stop("selection input not supported")
    }
  } else{
    spec_rhos <- NULL
  }
  x$control$timename <- NULL # to remove from v0.1.0 model fits

#
# NEEDS TO BE GENERALIZED
#
  seg_mods <- multinom_TS(x$data$train$ts_data, x$formula, spec_rhos,  
                          x$timename, x$weights, x$control)
#
#
#

  nsegs <- length(seg_mods[[1]])
  t1 <- min(x$data$train$ts_data[, x$timename])
  t2 <- max(x$data$train$ts_data[, x$timename])

  if (is.null(xname)){
    xname <- x$timename
  }
  par(mar = c(4, 5, 1, 1))
  plot(1, 1, type = "n", bty = "L", xlab = "", ylab = "", xaxt = "n", 
       yaxt = "n", ylim = c(0, 1), xlim = c(t1 - 1, t2 + 1))
  yax <- round(seq(0, 1, length.out = 5), 3)
  axis(2, at = yax, las = 1)
  axis(1)
  mtext(side = 2, line = 3.5, cex = 1.25, "Proportion")
  mtext(side = 1, line = 2.5, cex = 1.25, xname)
  ntopics <- ncol(as.matrix(x$data$train$ts_data$gamma))
  seg1 <- c(0, spec_rhos[-length(rhos)])
  seg2 <- c(spec_rhos, t2)
  time_obs <- rep(NA, nrow(x$data$train$ts_data))
  pred_vals <- matrix(NA, nrow(x$data$train$ts_data), ntopics)
  sp1 <- 1
  for (i in 1:nsegs){
    mod_i <- seg_mods[[1]][[i]]
    spec_vals <- sp1:(sp1 + nrow(mod_i$fitted.values) - 1)
    pred_vals[spec_vals, ] <- mod_i$fitted.values
    time_obs[spec_vals] <- mod_i$timevals
    sp1 <- sp1 + nrow(mod_i$fitted.values)
  }
  for (i in 1:ntopics){
    points(time_obs, pred_vals[ , i], type = "l", lwd = 3, col = cols[i])
  }
  if(!is.null(spec_rhos)){
    rho_lines(spec_rhos)
  }
}

#' @rdname plot.TS
#'
#' @export 
#'
rho_lines <- function(spec_rhos) {
  if(is.null(spec_rhos)) {
    return()
  }
  for (spec_rho in spec_rhos) {
    points(rep(spec_rho, 2), c(0, 1), type = "l", lwd = 4, 
           col = rgb(0.3, 0.3, 0.3, 0.6))
  }
}


#' @rdname plot.TS
#'
#' @export 
#'
rho_hist <- function(x, cols = set_rho_hist_colors(x$rhos), bin_width = 1, 
                     xname = NULL, border = NA, together = FALSE, 
                     LDATS = FALSE){
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  if(LDATS){
    par(fig = c(0, 1, 0.3, 0.55), mar = c(1.5, 5, 1, 1), new = TRUE)
  } else if(together){
    par(fig = c(0, 1, 0.54, 1), mar = c(1.5, 5, 1, 1))
  } else{
    par(fig = c(0, 1, 0, 1), mar = c(4, 5, 1, 1))
  } 
  rhos <- x$focal_rhos
  nrhos <- ncol(rhos)
  niter <- nrow(rhos)
  timeobs <- x$data$train$ts_data[, x$timename]
  timerange <- range(timeobs)
  timevals <- seq(timerange[1], timerange[2], 1)
  ntimes <- length(timevals) 
  
  binsteps <- seq(timerange[1], timerange[2] + bin_width, bin_width)
  bin1 <- binsteps[1:(length(binsteps)-1)]
  bin2 <- binsteps[2:(length(binsteps))] 
  nbins <- length(bin1)

  maxht <- 1
  if (!is.null(nrhos)){
    hts <- matrix(NA, nbins, nrhos)
    for(i in 1:nbins){
      for(j in 1:nrhos){
        hts[i,j]<- length(which(rhos[,j] >= bin1[i] & rhos[,j] < bin2[i]))
      }
    }
 
    maxht <- max(hts) / niter
    maxht <- ceiling(maxht * 100) / 100
    if (is.null(xname) & !together & !LDATS){
      xname <- x$timename
    }
  }
  plot(1, 1, type = "n", bty = "L", xlab = xname, ylab = "", xaxt = "n", 
       yaxt = "n", ylim = c(0, maxht), xlim = c(bin1[1], bin2[nbins]))
  yax <- round(seq(0, maxht, length.out = 5), 3)
  axis(2, at = yax, las = 1)
  axis(1)
  mtext(side = 2, line = 3.75, cex = 1.25, "Proportion")
  if (!is.null(nrhos)){
    for (i in 1:nbins){
      rho_ord <- order(hts[i,], decreasing = TRUE)
      for(j in 1:nrhos){
        ht_j <- hts[i, rho_ord[j]] / niter
        col_j <- cols[rho_ord[j]]
        rect(bin1[i], 0, bin2[i], ht_j, col = col_j, border = border)      
      }
    }
  }
}



#' @rdname plot.TS
#'
#' @export 
#'
set_rho_hist_colors <- function(x = NULL, cols = NULL, option = "D", 
                                alpha = 1){
  if(is.null(x)){
    return(NULL)
  }

  nrhos <- ncol(x$focal_rhos)
  if (length(cols) == 0){
    cols <- viridis(nrhos, option = option, alpha = alpha, end = 0.9)
  }
  if (length(cols) == 1){
    if (cols == "greys" | cols == "grey" | cols == "grays" | cols == "gray"){
      ggg <- seq(0, 0.8, length.out = nrhos)
      cols <- rep(NA, nrhos)
      for (i in 1:nrhos){
       cols[i] <- rgb(ggg[i], ggg[i], ggg[i], alpha = alpha)
      }
    }
  }
  if (length(cols) > nrhos){
    cols <- cols[1:nrhos]
  }
  if (length(cols) < nrhos){
    nc <- length(cols)
    nt <- nrhos
    msg <- paste0("Fewer colors (", nc, ") provided than change points (", 
                  nt, ")")
    stop(msg)
  }
  cols
}




#' @rdname plot.TS
#'
#' @export 
#'
set_gamma_colors <- function(x, cols = NULL, option = "D", alpha = 1){
  if(is.null(x)){
    return(NULL)
  }

  ntopics <- ncol(as.matrix(x$data$train$ts_data$gamma))
  if (length(cols) == 0){
    cols <- viridis(ntopics, option = option, alpha = alpha, end = 0.9)
  }
  if (length(cols) == 1){
    if (cols == "greys" | cols == "grey" | cols == "grays" | cols == "gray"){
      ggg <- seq(0, 0.8, length.out = ntopics)
      cols <- rep(NA, ntopics)
      for (i in 1:ntopics){
       cols[i] <- rgb(ggg[i], ggg[i], ggg[i], alpha = alpha)
      }
    }
  }
  if (length(cols) > ntopics){
    cols <- cols[1:ntopics]
  }
  if (length(cols) < ntopics){
    nc <- length(cols)
    nt <- ntopics
    msg <- paste0("Fewer colors (", nc, ") provided than number of topics (", 
                  nt, ")")
    stop(msg)
  }
  cols
}

