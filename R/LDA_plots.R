#' @title Plot a set of LDATS LDA models
#'
#' @description Generalization of the \code{\link[graphics]{plot}} function to 
#'   work on a list of LDA topic models (class \code{LDA_set}). 
#' 
#' @param x An \code{LDA_set} object of LDA topic models.
#'
#' @param ... Not used, retained for alignment with base function.
#' 
#' @param selected \code{logical} indicator of if only the selected LDAs
#'   (the first element in \code{x}) should be plotted or if all the LDAs
#'   (the second element in \code{x}) should be plotted.
#'
#' @return \code{NULL}.
#'
#' @export 
#'
plot.LDA_set <- function(x, ..., selected = TRUE){
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

#' @title Plot the results of an LDATS LDA model
#'
#' @description 
#'   \code{plot.LDA} creates an LDATS LDA summary plot, with a top panel 
#'     showing the topic proportions for each word and a bottom panel showing
#'     the topic proportions of each document/over time. \cr \cr
#'   \code{LDA_plot_top_panel} creates an LDATS LDA summary plot 
#'     top panel showing the topic proportions word-by-word. \cr \cr
#'   \code{LDA_plot_bottom_panel} creates an LDATS LDA summary plot
#'     bottom panel showing the topic proportions over time/documents. \cr \cr
#'   \code{set_LDA_plot_colors} creates the set of colors to be used in
#'     the LDA plots based on the variety of argument options.
#' 
#' @param x Object of class \code{LDA}.
#'
#' @param xtime Optional x values used to plot the topic proportions according
#'   to a specific time value (rather than simply the order of observations).
#'
#' @param xname Optional name for the x values used in plotting the topic
#'   proportions (otherwise defaults to "Document"). 
#'
#' @param cols Colors to be used to plot the topics.
#'   Any valid color values (\emph{e.g.}, see \code{\link[grDevices]{colors}},
#'   \code{\link[grDevices]{rgb}}) can be input as with a standard plot. 
#'   The default (\code{cols = NULL}) triggers use of 
#'   \code{\link[viridis]{viridis}} color options (see \code{option}).
#'
#' @param option A \code{character} string indicating the color option
#'   from \code{\link[viridis]{viridis}} to use if `cols == NULL`. Four 
#'   options are available: "magma" (or "A"), "inferno" (or "B"), "plasma" 
#'   (or "C", the default option), "viridis" (or "D") and "cividis" (or "E").
#'
#' @param alpha Numeric value [0,1] that indicates the transparency of the 
#'   colors used. Supported only on some devices, see 
#'   \code{\link[grDevices]{rgb}}.
#'
#' @param together \code{logical} indicating if the subplots are part of a 
#'   larger LDA plot output. 
#'
#' @param LDATS \code{logical} indicating if the LDA plot is part of a larger 
#'   LDATS plot output. 
#'
#' @param ... Not used, retained for alignment with base function.
#' 
#' @return 
#'   \code{plot.LDA},\code{LDA_plot_top_panel},\code{LDA_plot_bottom_panel}: 
#'     \code{NULL}. \cr \cr
#'   \code{set_LDA_plot_colors}: \code{vector} of \code{character} hex codes
#'     indicating colors to use.
#'
#' @name plot.LDA 
#'



#' @rdname plot.LDA
#'
#' @export
#'
plot.LDA <- function(x, ..., xtime = NULL, xname = NULL, cols = NULL, 
                     option = "C", alpha = 0.8, LDATS = FALSE){
  LDA_plot_top_panel(x, cols, option, alpha, TRUE, LDATS)
  LDA_plot_bottom_panel(x, xtime, xname, cols, option, alpha, TRUE, LDATS)
}

#' @rdname plot.LDA 
#' 
#' @export 
#'
LDA_plot_top_panel <- function(x, cols = NULL, option = "C", alpha = 0.8, 
                               together = FALSE, LDATS = FALSE){
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  cols <- set_LDA_plot_colors(x, cols, option, alpha)
  gamma <- x$document_topic_table
  beta <- exp(x$params$beta)
  nobs <- nrow(gamma)
  ntopics <- ncol(gamma)
  nwords <- ncol(beta)
  beta_order <- apply(beta, 2, order)
  beta_sorted <- apply(beta, 2, sort)
  if(length(dim(beta_sorted)) == 0){
    beta_order <- matrix(beta_order, nrow = 1)
    beta_sorted <- matrix(beta_sorted, nrow = 1)
  }

  counter <- 1
  rect_mat <- matrix(NA, nrow = nwords * ntopics, ncol = 4)
  rect_col <- rep(NA, length = nwords * ntopics)
  for (i in 1:nwords){
    x1 <- i - 0.4
    x2 <- i + 0.4
    y1 <- 0
    y2 <- 0
    for (j in 1:ntopics){
      y1 <- y2
      y2 <- y1 + beta_sorted[j, i]
      rect_mat[counter, ] <- c(x1, y1, x2, y2)      
      rect_col[counter] <- cols[beta_order[j, i]]
      counter <- counter + 1
    }
  }

  if (LDATS){
    par(fig = c(0, 0.9, 0.85, 1))
  } else if(together){
    par(fig = c(0, 0.9, 0.7, 1))
  } else{
    par(fig = c(0, 0.9, 0, 1))
  }
  par(mar = c(1, 3.5, 1, 0))
  max_y <- max(rect_mat[,4]) * 1.05
  plot(1, 1, type = "n", bty = "L", xlab = "", ylab = "", las = 1,
       ylim = c(0, max_y), xlim = c(1, nwords), xaxt = "n", cex.axis = 0.75)  
  mtext(side = 2, "Total Proportion", line = 2.5, cex = 0.75)
  for (i in 1:(nwords * ntopics)){
    rect(rect_mat[i, 1], rect_mat[i, 2], rect_mat[i, 3], rect_mat[i, 4],
         col = rect_col[i])
  }
  axis(2, at = seq(0, max_y, 0.1), labels = FALSE, tck = -0.02)
  mtext(side = 1, at = seq(1, nwords, 1), text = x$terms, tck = 0, 
       cex = 0.5, line = 0)

  if (LDATS){
    par(fig = c(0.9, 1, 0.85, 1)) 
  } else if (together){ 
    par(fig = c(0.9, 1, 0.7, 1)) 
  } else{
    par(fig = c(0.9, 1, 0, 1))
  }
  par(mar = c(0, 0, 0, 0), new = TRUE)
  plot(1, 1, type = "n", bty = "n", xlab = "", ylab = "", 
       xaxt = "n", yaxt = "n", ylim = c(0, 1), xlim = c(0,1))

  ypos <- (0.9 / ntopics) * (ntopics:1)
  ttext <- paste("Topic ", 1:ntopics, sep = "")
  for (i in 1:ntopics){
    text(ttext[i], x = 0.25, y = ypos[i], adj = 0, cex = 0.75)
    rect(0.0, ypos[i] - 0.05, 0.15, ypos[i] + 0.05, col = cols[i])
  }
}

#' @rdname plot.LDA
#' 
#' @export 
#'
LDA_plot_bottom_panel <- function(x, xtime = NULL, xname = NULL, cols = NULL,
                                  option = "C", alpha = 0.8, 
                                  together = FALSE, LDATS = FALSE){
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  cols <- set_LDA_plot_colors(x, cols, option, alpha)
  gamma <- x$document_topic_table
  ntopics <- ncol(gamma)

  if (is.null(xtime)){
    xtime <- seq(1, nrow(gamma), 1)
  }
  if (is.null(xname) & !LDATS){
    xname <- "Document"
  }

  if (LDATS){
    par(fig = c(0, 1, 0.53, 0.85), new = TRUE)
  } else if (together){
    par(fig = c(0, 1, 0, 0.7), new = TRUE)
  } else{
    par(fig = c(0, 1, 0, 1))
  }
  par(mar = c(3.25, 5, 1, 1))
  plot(xtime, gamma[ , 1], type = "n", bty = "L", xlab = "", ylab = "", 
       las = 1, ylim = c(0, 1))
  mtext(side = 1, xname, line = 2.2, cex = 1.25)
  mtext(side = 2, "Proportion", line = 3.5, cex = 1.25)
  for (i in 1:ntopics){
    points(xtime, gamma[ , i], col = cols[i], type = "l", lwd = 2)
  }
}

#' @rdname plot.LDA
#'
#' @export
#'
set_LDA_plot_colors <- function(x, cols = NULL, option = "C", alpha = 0.8){

  gamma <- x$document_topic_table
  ntopics <- ncol(gamma)
  if (length(cols) == 0){
    cols <- viridis(ntopics, option = option, alpha = alpha, end = 0.9)
  }
  if (length(cols) == 1){
    if (cols == "greys" | cols == "grey" | cols == "grays" | cols == "gray"){
      ggg <- seq(0, 0.8, length.out = ntopics)
      cols <- rep(NA, ntopics)
      for (i in 1:ntopics){
       cols[i] <- rgb(ggg[i], ggg[i], ggg[i])
      }
    }
  }
  if (length(cols) > ntopics){
    cols <- cols[1:ntopics]
  }
  if (length(cols) < ntopics){
    nc <- length(cols)
    nt <- ntopics
    msg <- paste0("Fewer colors (", nc, ") provided than topics (", nt, ")")
    stop(msg)
  }
  cols
}

