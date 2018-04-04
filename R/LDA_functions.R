#' A simple parallel wrapper to the topicmodels LDA function 
#' 
#' Runs each of the number of input topics for the number of seeds
#' 
#' @param data Table of integer data (species counts by period)
#' @param ntopics set of topics to evaluate
#' @param nseeds number of seeds to use
#' @param ncores Integer number of cores to use
#' @param ... additional arguments to be passed to the LDA function
#' 
#' @return List of LDA models
#' 
#' @examples NA
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

#' Generalization of the AIC function to work on LDA topic models 
#' 
#' Using the internal topicmodels calculations of likelihood and df
#' 
#' @param x an LDA topic model
#' @param k penalty per df
#' @param correction whether or not to include the small sample size 
#'   correction (thus generating AICc)
#' 
#' @return Named (AIC or AICc) value.
#' 
#' @examples NA
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

