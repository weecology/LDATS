#' @title A simple parallel wrapper to the topicmodels LDA function 
#' 
#' @description Runs each of the number of input topics for the number of 
#'   seeds
#' 
#' @param data matrix of integer data (species counts by period)
#' @param ntopics vector of the topics to evaluate
#' @param nseeds number of seeds (replicate starts) to use for each value of
#'   ntopics
#' @param ncores integer number of cores to use
#' @param ... additional arguments to be passed to the topicmodel LDA function
#' 
#' @return List of LDA models
#' 
#' @examples 
#' \dontrun{
#'   data(rodents)
#'   lda_data <- dplyr::select(rodents, -c(newmoon, date, plots, traps))
#'   r_LDA <- parLDA(data = lda_data, ntopics = 2, nseeds = 2, ncores = 4)
#'                         
#' }
#' @export 
#'
parLDA <- function(data, ntopics = 2, nseeds = 1, ncores = 1, ...) {

  max_cores <- detectCores(logical = TRUE)
  seed_in <- rep(seq(2, nseeds * 2, 2), length(ntopics))
  k_in <- rep(ntopics, each = length(seq(2, nseeds * 2, 2)))
  nruns <- length(seed_in)

  if (ncores > max_cores){
    capped_cores <- floor(0.75 * parallel::detectCores(logical = TRUE))
    msg <- paste(ncores, " is larger than max cores available (", max_cores,
             "); capped at ", capped_cores, " cores.", sep = "")
    warning(msg)
    ncores <- capped_cores
  }
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  i <- 1
  mods <- foreach(i = 1:nruns, .packages = "topicmodels",
                  .errorhandling = "pass") %dopar% {

    topicmodels::LDA(data, k = k_in[i], control = list(seed = seed_in[i]),
                     ...)
  }

  stopCluster(cl)
  names(mods) <- paste("k: ", k_in, ", seed: ", seed_in, sep = "")
  class(mods) <- c("LDA_list", "list")
  return(mods)
}

#' @title Generalization of the AIC function to work on LDA topic models 
#' 
#' @description Using the internal topicmodels calculations of likelihood and 
#'   df
#' 
#' @param object an LDA topic model
#' @param ... additional arguments to be passed to subfunctions
#' @param k penalty per df
#' @param correction AICc correction
#' 
#' @return Named (AIC or AICc) value.
#' 
#' @examples
#' \dontrun{
#'   data(rodents)
#'   lda_data <- select(rodents, -c(newmoon, date, plots, traps))
#'   r_LDA <- LDA(data = lda_data, ntopics = 2, nseeds = 2)
#'   AIC(r_LDA[[1]])
#'   AIC(r_LDA[[2]])
#' }
#'
#' @export 
#'
AIC.LDA <- function(object, ..., k = 2, correction = FALSE){
  val <- logLik(object)
  ll <- as.numeric(val)
  df <- attr(val, "df")
  out <- -2 * ll + k * df
  names(out) <- "AIC"
  if (correction == TRUE){
    nobs <- attr(val, "nobs")
    corr <- ((2 * df^2) + (2 * df)) / (nobs - df - 1)
    out <- -2 * ll + k * df + corr
    names(out) <- "AICc"
  }
  return(out)
}

#' @title Plot a set of LDA models
#'
#' @description Generalization of the plot function to work on a list of LDA 
#'   topic models 
#' 
#' @param x a list of LDA topic model outputs
#' @param ... additional arguments to be passed to subfunctions
#' @return model-by-model plots
#' 
#' @examples 
#' \dontrun{
#'   data(rodents)
#'   lda_data <- select(rodents, -c(newmoon, date, plots, traps))
#'   r_LDA <- parLDA(data = lda_data, ntopics = 2, nseeds = 2)
#'   plot(r_LDA)
#' }
#' @export 
#'
plot.LDA_list <- function(x, ...){
  devAskNewPage(TRUE)
  lapply(x, plot, ...)
  devAskNewPage(FALSE)
}

#' @title Plot an LDA model
#'
#' @description Function used to create an LDA summary plot, called using the 
#'   plot method
#' 
#' @param x an LDA topic model output
#' @param ... additional arguments to be passed to subfunctions
#' @param cols colors to be used to plot the topics
#' @param option A character string indicating the colormap option to use if 
#'   `cols == NULL`. Four options are available: "magma" (or "A"), "inferno" 
#'   (or "B"), "plasma" (or "C"), "viridis" (or "D", the default option) and 
#'   "cividis" (or "E").
#' @return model plots
#' 
#' @examples 
#'   data(rodents)
#'   lda_data <- dplyr::select(rodents, -c(newmoon, date, plots, traps))
#'   lda_models <- parLDA(data = lda_data, ntopics = 4, nseeds = 10)
#'   best_lda <- LDA_select(lda_models)
#'   plot(best_lda, option = "cividis")
#' 
#' @export 
#'
plot.LDA <- function(x, ..., cols = NULL, option = "D"){

  gamma <- x@gamma
  beta <- exp(x@beta)
  nobs <- nrow(gamma)
  ntopics <- ncol(gamma)
  nwords <- ncol(beta)
  beta_order <- apply(beta, 2, order)
  beta_sorted <- apply(beta, 2, sort)

  if (length(cols) == 0){
    cols <- viridis::viridis(ntopics, option = option)
  }
  if (length(cols) == 1){
    if (cols == "greys"){
      ggg <- runif(ntopics, 0, 0.8)
      cols <- rep(NA, ntopics)
      for (i in 1:ntopics){
       cols[i] <- rgb(ggg[i], ggg[i], ggg[i])
      }
    }
  }
  if (length(cols) > ntopics){
    cols <- cols[1:ntopics]
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

  par(fig = c(0, 1, 0, 0.7), mar = c(3.25, 4, 1, 1))
  plot(gamma[ , 1], type = "n", bty = "L", xlab = "", ylab = "", las = 1,
       ylim = c(0, 1))
  mtext(side = 1, "Observation", line = 2.2, cex = 1.25)
  mtext(side = 2, "Proportion", line = 2.8, cex = 1.25)
  for (i in 1:ntopics){
    points(gamma[ , i], col = cols[i], type = "l", lwd = 1)
  }

  par(fig = c(0, 0.85, 0.7, 1), new = TRUE, mar = c(1, 3, 1, 0))
  max_y <- max(rect_mat[,4]) * 1.05
  plot(1, 1, type = "n", bty = "L", xlab = "", ylab = "", las = 1,
       ylim = c(0, max_y), xlim = c(1, nwords), xaxt = "n", cex.axis = 0.75)  
  mtext(side = 2, "Total Proportion", line = 2.125, cex = 0.75)
  for(i in 1:(nwords * ntopics)){
    rect(rect_mat[i, 1], rect_mat[i, 2], rect_mat[i, 3], rect_mat[i, 4],
         col = rect_col[i])
  }
  axis(2, at = seq(0, max_y, 0.1), labels = FALSE, tck = -0.02)
  mtext(side = 1, at = seq(1, nwords, 1), text = x@terms, tck = 0, 
       cex = 0.5, line = 0)

  par(fig = c(0.85, 1, 0.7, 1), new = TRUE, mar = c(0, 0, 0, 0))
  plot(1, 1, type = "n", bty = "n", xlab = "", ylab = "", 
     xaxt = "n", yaxt = "n", ylim = c(0, 1), xlim = c(0,1))

  ypos <- (0.9 / ntopics) * (ntopics:1)
  ttext <- paste("Topic ", 1:ntopics, sep = "")
  for (i in 1:ntopics){
    text(ttext[i], x = 0.1, y = ypos[i], col = cols[i], adj = 0, cex = 0.75)
  }
   
}
