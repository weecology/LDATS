#' @title A simple parallel wrapper to the topicmodels LDA function 
#' 
#' @description Runs each of the number of input topics for the number of 
#'   seeds
#' 
#' @param data Table of integer data (species counts by period)
#' @param ntopics set of topics to evaluate
#' @param nseeds number of seeds to use
#' @param ncores Integer number of cores to use
#' @param ... additional arguments to be passed to the LDA function
#' 
#' @return List of LDA models
#' 
#' @examples 
#'   data(rodents)
#'   lda_data <- dplyr::select(rodents, 
#'                             -c(newmoonnumber, newmoondate, nplots, ntraps))
#'   r_LDA <- LDATS::LDA(data = lda_data, ntopics = 2, nseeds = 2, 
#'                       ncores = 4)
#' @export 
#'
LDA <- function(data, ntopics = 2, nseeds = 1, ncores = 1, ...) {

  max_cores <- parallel::detectCores(logical = TRUE)
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
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  mods <- foreach::foreach(i = 1:nruns, .packages = "topicmodels",
            .errorhandling = "pass") %dopar% {

    k <- k_in[i]
    seed <- seed_in[i]
    mod <- topicmodels::LDA(data, k = k, control = list(seed = seed), ...)
    #class(mod) <- c(class(mod)[1], "LDA")
    mod
  }

  parallel::stopCluster(cl)
  names(mods) <- paste("k: ", k_in, ", seed: ", seed_in, sep = "")
  class(mods) <- c("LDA_list", "list")
  return(mods)
}

#' @title Generalization of the AIC function to work on LDA topic models 
#' 
#' @description Using the internal topicmodels calculations of likelihood and 
#'   df
#' 
#' @param x an LDA topic model
#' @param k penalty per df
#' @param correction whether or not to include the small sample size 
#'   correction (thus generating AICc)
#' 
#' @return Named (AIC or AICc) value.
#' 
#' @examples
#'   data(rodents)
#'   lda_data <- dplyr::select(rodents, 
#'                             -c(newmoonnumber, newmoondate, nplots, ntraps))
#'   r_LDA <- LDATS::LDA(data = lda_data, ntopics = 2, nseeds = 2, 
#'                       ncores = 4)
#'   AIC(r_LDA[[1]])
#'   AIC(r_LDA[[2]])
#' @export 
#'
AIC.LDA <- function(x, k = 2, correction = FALSE){
  val <- topicmodels::logLik(x)
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

#' @title Generalization of the plot function to work on a list of LDA topic 
#'   models 
#' 
#' @param model_list a list of LDA topic model outputs
#' @param cols a vector of color codes to use, one for each of the max topics
#' @return model-by-model plots
#' 
#' @examples 
#'   data(rodents)
#'   lda_data <- dplyr::select(rodents, 
#'                             -c(newmoonnumber, newmoondate, nplots, ntraps))
#'   r_LDA <- LDATS::LDA(data = lda_data, ntopics = 2, nseeds = 2, 
#'                       ncores = 4)
#'   plot(r_LDA)
#' @export 
#'
plot.LDA_list <- function(model_list, cols = NULL){
  devAskNewPage(TRUE)
  lapply(model_list, plot, cols)
  devAskNewPage(FALSE)
}

#' @title Function used to create an LDA summary plot, called using the plot 
#'   method
#' 
#' @param model an LDA topic model output
#' @param cols a vector of color codes, one for each topic
#' @return model plots
#' 
#' @export 
#'
plot.LDA <- function(model, cols = NULL, ...){

  gamma <- model@gamma
  beta <- exp(model@beta)
  nobs <- nrow(gamma)
  ntopics <- ncol(gamma)
  nwords <- ncol(beta)
  beta_order <- apply(beta, 2, order)
  beta_sorted <- apply(beta, 2, sort)

  if (length(cols) == 0){
    cols <- rgb(runif(ntopics), runif(ntopics), runif(ntopics))
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
  mtext(side = 1, at = seq(1, nwords, 1), text = model@terms, tck = 0, 
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
