
#' @title Fit a multinomial change point Time Series model
#'
#' @description Fit a set of multinomial regression models to a time series of
#'   of data divided into multiple chunks based on change points. 
#'
#' @param data Class \code{data.frame} object including the predictor and 
#'   response variables.
#'
#' @param formula_RHS Right Hand Side of the continuous time formula as a 
#'   character vector.
#'
#' @param changepoints Numeric vector indicating locations of the change 
#'   points.
#'
#' @param weights Optional class \code{numeric} vector of weights for each 
#'   document. Corresponds to the vector \strong{\eqn{v}} in the math 
#'   description.
#'
#' @param control Class \code{TS_controls} list, holding control parameters
#'   for the Time Series model including the parallel tempering Markov Chain 
#'   Monte Carlo (ptMCMC) controls.
#'
#' @return List of chunk-level model fits and the total log likelihood 
#'   combined across all the chunks.
#'
#' @export 
#'
multinom_TS <- function(data, formula_RHS, changepoints = NULL, 
                        weights = NULL, control = TS_controls_list()){

  if (!check_chunks(data, changepoints)){
    return(list("chunk models" = NA, "logLik" = -Inf))
  }

  chunk_memo <- memoise_fun(multinom_chunk, control$memoise)
  starts <- c(min(data$time) - 1, changepoints)   
  ends <- c(changepoints, max(data$time)) #removed the +1 on the max
  nchunks <- length(changepoints) + 1

  mods <- vector("list", length = nchunks)
  names(mods) <- sprintf("%s %d %s", "chunk", 1:nchunks, "model")
  ll <- rep(0, nchunks)
  for (i in 1:nchunks){
    mods[[i]] <- chunk_memo(data, formula_RHS, starts, ends, weights)
    ll[i] <- logLik(mods)
  }
  list("chunk models" = mods, "logLik" = sum(ll))
}

#' @title Verify the chunks of a multinomial time series model
#'
#' @description Check to verify that a time series can be broken into a set 
#'   of chunks based on input changepoints. 
#'
#' @param data Class \code{data.frame} object including the predictor and 
#'   response variables.
#'
#' @param changepoints Numeric vector indicating locations of the change 
#'   points.
#'
#' @return Logical indicator of the check passing \code{TRUE} or failing
#'   \code{FALSE}.
#'
#' @export 
#'
check_chunks <- function(data, changepoints){
  last_time <- max(data$time)
  time_check <- any(changepoints <= 0) | any(changepoints >= last_time)
  sort_check <- is.unsorted(changepoints, strictly = TRUE)
  check <- time_check | sort_check
  return(check)
}

#' @title Fit a multinomial Time Series model chunk
#'
#' @description Fit a multinomial regression model to a defined chunk of time
#'   \code{(start_time, end_time]} within a time series. The fit is conducted
#'   via \code{\link[nnet]{multinom}}
#'
#' @param data Class \code{data.frame} object including the predictor and 
#'   response variables.
#'
#' @param formula_RHS Right Hand Side of the continuous time formula as a 
#'   character vector.
#'
#' @param start_time Start time for the chunk.
#'
#' @param end_time End time for the chunk.
#'
#' @param weights Optional class \code{numeric} vector of weights for each 
#'   document. Corresponds to the vector \strong{\eqn{v}} in the math 
#'   description.
#' 
#' @return Fitted model for the chunk.
#' 
#' @export 
#'
multinom_chunk <- function(data, formula_RHS, start_time, end_time, 
                           weights = NULL){

  formula <- as.formula(paste("gamma ~", formula))
  chunk <- data$time > start_time & data$time <= end_time
  multinom(formula_RHS, data, weights, subset = chunk, trace = FALSE) 
}