
plot.TS_fit <- function(x, ..., plot_type = "diagnostic", 
                        cols = TS_summary_cols(),
                        bin_width = 1, xlab = NULL, selection = "median"){
  if (plot_type == "diagnostic"){
    TS_diagnostics_plot(x)
  } else if (plot_type == "summary"){
    TS_summary_plot(x, bin_width, xlab, selection, cols)
  }
}


TS_diagnostics_plot <- function(x){

}

#' @title Create the list of colors for the TS summary plot
#'
#' @description A default list generator function that produces the options
#'   for the colors controlling the panels of the TS summary plots, so needed
#'   because the panels should be in different color schemes. See 
#'   \code{set_gamma_colors} and \code{set_rho_hist_colors} for
#'   specific details on usage.
#'
#' @param rho_cols Colors to be used to plot the histograms of changepoints.
#'
#' @param rho_option A character string indicating the colormap option to use
#'   if `cols == NULL`. Four options are available: "magma" (or "A"), 
#'   "inferno" (or "B"), "plasma" (or "C"), "viridis" (or "D", the default
#'   option) and "cividis" (or "E").
#'
#' @param rho_alpha Numeric value [0,1] that indicates the transparency of the 
#'   colors used. Supported only on some devices, see \code{rgb}.
#'
#' @param gamma_cols Colors to be used to plot the time series of fitted topic 
#'   proportions (gammas).
#'
#' @param gamma_option A character string indicating the colormap option to 
#'   use if `cols == NULL`. Four options are available: "magma" (or "A"),
#'   "inferno" (or "B", the default option), "plasma" (or "C"), "viridis" 
#'   (or "D") and "cividis" (or "E")
#'
#' @param gamma_alpha Numeric value [0,1] that indicates the transparency of 
#'   the colors used. Supported only on some devices, see \code{rgb}.
#'
#' @return \code{list} of elements used to define the colors for the two
#'   panels. Contains two elements \code{rho} and \code{gamma}, each 
#'   corresponding to the related panel, and each containing default values 
#'   for entries named \code{cols}, \code{option}, and \code{alpha}. 
#'
#' @export
#'
TS_summary_cols <- function(rho_cols = NULL, rho_option = "D", 
                            rho_alpha = 0.4, gamma_cols = NULL, 
                            gamma_option = "C", gamma_alpha = 0.8){
  list(
    rho = list(cols = rho_cols, option = rho_option, alpha = rho_alpha),
    gamma = list(cols = gamma_cols, option = gamma_option, 
                 alpha = gamma_alpha)
  )
}


#' @title Create the summary plot for a TS fit to an LDA model
#'
#' @description Produces a two-panel figure of [1] the changepoint 
#'   distributions as histograms over time and [2] the time series of the 
#'   fitted topic proportions (gamma values) over time, based on a selected 
#'   set of changepoint locations.
#'
#' @param x Object of class \code{TS_list}.
#'
#' @param bin_width Width of the bins used in the histograms, in units of the
#'   x-axis (the time variable used to fit the model).
#'
#' @param xlab Label for the x-axis.
#'
#' @param selection Indicator of the changepoints to use. Currently only
#'   defined for \code{"median"}.
#'
#' @param cols \code{list} of elements used to define the colors for the two
#'   panels, as generated simply using \code{TS_summary_cols}. \code{cols}
#'   has two elements \code{rho} and \code{gamma}, each corresponding to the
#'   related panel, and each containing default values for entries named
#'   \code{cols}, \code{option}, and \code{alpha}. See
#'   \code{set_gamma_colors} and \code{set_rho_hist_colors} for
#'   details on usage.
#'
#' @return Nothing. The plot is generated in the active graphics device.
#'
#' @export
#'
TS_summary_plot <- function(x, bin_width, xlab, selection = "median", 
                            cols = TS_summary_cols()){

  par(mfrow = c(2, 1))
  rc <- cols$rho
  rho_cols <- set_rho_hist_colors(x$rhos, rc$cols, rc$option, rc$alpha)
  rho_hist(x, rho_cols, bin_width, xlab = NULL)

  gc <- cols$gamma
  gamma_cols <- set_gamma_colors(x, gc$cols, gc$option, gc$alpha)
  pred_gamma_TS_plot(x, selection, gamma_cols, xlab)

}

#' @title Create the plot of fitted topic proportions over time
#'
#' @description Produces a time series of the fitted topic proportions 
#'   (gamma values) over time, based on a selected set of changepoint 
#'   locations.
#'
#' @param x Object of class \code{TS_list}.
#'
#' @param selection Indicator of the changepoints to use. Currently only
#'   defined for \code{"median"}.
#'
#' @param cols Hex values of the colors to be used to plot the time series.
#'
#' @param xlab Label for the x-axis.
#'
#' @return Nothing. The plot is generated in the active graphics device.
#'
#' @export
#'
pred_gamma_TS_plot <- function(x, selection = "median", cols, xlab){

  rhos <- x$rhos
  nrhos <- ncol(rhos)
  if (selection == "median"){
    spec_rhos <- apply(rhos, 2, median)
  }  
  seg_mods <- multinom_TS(x$data, x$formula, spec_rhos, x$weights, x$control)
  nsegs <- length(seg_mods[[1]])
  t1 <- min(x$data[, x$control$timename])
  t2 <- max(x$data[, x$control$timename])

  if (is.null(xlab)){
    xlab <- x$control$timename
  }
  par(mar = c(4.5, 4, 1, 1))
  plot(1, 1, type = "n", bty = "L", xlab = xlab, ylab = "", xaxt = "n", 
       yaxt = "n", ylim = c(0, 1), xlim = c(t1 - 1, t2 + 1))
  yax <- round(seq(0, 1, length.out = 5), 3)
  axis(2, at = yax, las = 1)
  axis(1)
  ntopics <- ncol(as.matrix(x$data[[x$control$response]]))
  seg1 <- c(0, spec_rhos[-length(rhos)])
  seg2 <- c(spec_rhos, t2)
  time_obs <- rep(NA, nrow(x$data))
  pred_vals <- matrix(NA, nrow(x$data), ntopics)
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
  rho_lines(spec_rhos)

}

#' @title Add changepoint location lines to the time series plot
#'
#' @description Adds vertical lines to the fitted gamma time series plot 
#'   associated with the changepoints of interest.
#'
#' @param spec_rhos \code{numeric} vector indicating the locations along the
#'   x axis where the specific changepoints being used are located.
#'
#' @return Nothing. Lines are added to the active plot object.
#'
#' @export
#'
rho_lines <- function(spec_rhos){
  if(is.null(spec_rhos)){
    return(NULL)
  }
  for (i in 1:length(spec_rhos)){
    points(rep(spec_rhos[i], 2), c(0, 1), type = "l", lwd = 4, 
           col = rgb(0.3, 0.3, 0.3, 0.6))
  }
}

#' @title Create the changepoint histograms plot
#'
#' @description Produces a plot of the changepoint distributions as histograms
#'   over time.
#'
#' @param x Object of class \code{TS_list}.
#'
#' @param cols Hex values of the colors to be used to plot the histograms of 
#'   changepoints.
#'
#' @param bin_width Width of the bins used in the histograms, in units of the
#'   x-axis (the time variable used to fit the model).
#'
#' @param xlab Label for the x-axis.
#'
#' @return Nothing. The plot is generated in the active graphics device.
#'
#' @export
#'
rho_hist <- function(x, cols, bin_width, xlab = NULL){

  rhos <- x$rhos
  nrhos <- ncol(rhos)
  niter <- nrow(rhos)
  timeobs <- x$data[, x$control$timename]
  timerange <- range(timeobs)
  timevals <- seq(timerange[1], timerange[2], 1)
  ntimes <- length(timevals) 
  
  binsteps <- seq(timerange[1], timerange[2] + bin_width, bin_width)
  bin1 <- binsteps[1:(length(binsteps)-1)]
  bin2 <- binsteps[2:(length(binsteps))] 
  nbins <- length(bin1)

  hts <- matrix(NA, nbins, nrhos)
  for(i in 1:nbins){
    for(j in 1:nrhos){
      hts[i,j]<- length(which(rhos[,j] >= bin1[i] & rhos[,j] < bin2[i]))
    }
  }

  maxht <- max(hts) / niter
  maxht <- ceiling(maxht * 100) / 100
  if (is.null(xlab)){
    xlab <- x$control$timename
  }
  par(mar = c(1.5, 4, 1, 1))
  plot(1, 1, type = "n", bty = "L", xlab = xlab, ylab = "", xaxt = "n", 
       yaxt = "n", ylim = c(0, maxht), xlim = c(bin1[1], bin2[nbins]))
  yax <- round(seq(0, maxht, length.out = 5), 3)
  axis(2, at = yax, las = 1)
  axis(1)
  for (i in 1:nbins){
    rho_ord <- order(hts[i,], decreasing = TRUE)
    for(j in 1:nrhos){
      ht_j <- hts[i, rho_ord[j]] / niter
      col_j <- cols[rho_ord[j]]
      rect(bin1[i], 0, bin2[i], ht_j, col = col_j)      
    }
  }

}

#' @title Prepare the colors to be used in the changepoint histogram
#'
#' @description Based on the inputs, create the set of colors to be used in
#'   the changepoint histogram.
#' 
#' @param x Matrix of changepoint locations (element \code{rhos}) from an 
#'   object of class \code{TS_list}.
#'
#' @param cols Colors to be used to plot the histograms of changepoints.
#' 
#' @param option A character string indicating the colormap option to use if 
#'   `cols == NULL`. Four options are available: "magma" (or "A"), "inferno" 
#'   (or "B"), "plasma" (or "C"), "viridis" (or "D", the default option) and 
#'   "cividis" (or "E").
#'
#' @param alpha Numeric value [0,1] that indicates the transparency of the 
#'   colors used. Supported only on some devices, see \code{rgb}.
#'
#' @return Vector of \code{character} hex codes indicating colors to use.
#'
#' @export 
#'
set_rho_hist_colors <- function(x, cols = NULL, option = "D", alpha = 1){
  if(is.null(x)){
    return(NULL)
  }

  nrhos <- ncol(x)
  if (length(cols) == 0){
    cols <- viridis(nrhos, option = option, alpha = alpha)
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
    msg <- paste0("Fewer colors (", nc, ") provided than changepoints (", 
                  nt, ")")
    stop(msg)
  }
  cols
}

#' @title Prepare the colors to be used in the gamma time series
#'
#' @description Based on the inputs, create the set of colors to be used in
#'   the timeseries of the fitted gamma (topic proportion) values.
#' 
#' @param x Object of class \code{TS_list}.
#'
#' @param cols Colors to be used to plot the time series of fitted topic 
#'   proportions (gammas).
#' 
#' @param option A character string indicating the colormap option to use if 
#'   `cols == NULL`. Four options are available: "magma" (or "A"), "inferno" 
#'   (or "B"), "plasma" (or "C"), "viridis" (or "D", the default option) and 
#'   "cividis" (or "E").
#'
#' @param alpha Numeric value [0,1] that indicates the transparency of the 
#'   colors used. Supported only on some devices, see \code{rgb}.
#'
#' @return Vector of \code{character} hex codes indicating colors to use.
#'
#' @export 
#'
set_gamma_colors <- function(x, cols = NULL, option = "D", alpha = 1){
  if(is.null(x)){
    return(NULL)
  }

  ntopics <- ncol(as.matrix(x$data[x$control$response]))
  if (length(cols) == 0){
    cols <- viridis(ntopics, option = option, alpha = alpha)
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

