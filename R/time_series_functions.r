#' Fit a multinomial regression model to topic model data for dates within
#'   a defined chunk of time (start_date, end_date] with covariate impacts 
#'   assuming no temporal autocorrelation by default
#'
#' @param formula_RHS Right Hand Side of the continuous time formula as a 
#'   character vector
#' @param data data frame including the predictor and response variables
#' @param start_date start date for the chunk 
#' @param end_date end date for the chunk
#' @param weights weights 
#' @param ... other arguments to be passed to the multinom function
#' @return fitted model for the chunk
#' @examples
#' NA
#' @export 
#'
multinom_chunk <- function(formula_RHS, data, start_date, end_date, weights, 
                           ...) {
  formula <- as.formula(paste("gamma ~", formula_RHS))
  mod <- nnet::multinom(formula, data, weights, 
           subset = date > start_date & date <= end_date, trace = FALSE, ...) 
  return(mod)
}

#' Fit a multinomial regression model with covariate impacts assuming no
#'   temporal autocorrelation 
#'
#' @param formula_RHS Right Hand Side of the continuous time formula as a 
#'   character vector
#' @param data data frame including the predictor and response variables
#' @param changepoints selections of the change points
#' @param weights weights 
#' @param ... other arguments to be passed to subfunctions
#' @return list of chunk-level model fits and the total log likelihood 
#'   combined across all the chunks
#' @examples
#' NA
#' @export 
#'
mutilnom_ts <- function(formula_RHS, data, changepoints = NULL, weights, ...){

  chunk_memo <- memoise::memoise(LDATS::multinom_chunk)

  nchunks <- length(changepoints) + 1
  start_dates <- c(min(data$date) - 1, changepoints)
  end_dates <- c(changepoints, max(data$date) + 1)

  mods <- vector("list", length = nchunks)
  ll <- rep(NA, nchunks)
  for (i in 1:nchunks){
    mods[[i]] <- chunk_memo(formula_RHS, data, start_date = start_dates[i], 
                   end_date = end_dates[i], weights)
    ll[i] <- logLik(mods[[i]])
  }
  total_ll <- sum(ll)
  names(mods) <- sprintf("%s %d %s", "chunk", 1:nchunks, "model")
  output <- list("chunk models" = mods, "logLik" = total_ll)

  return(output)
}

#' Bayesian Time Series analysis of a topic model classification
#'
#'
#' @param model fitted LDA model
#' @param formula formula 
#' @param data 
#'
#'
#'
#'
BTS <- function(model, formula, data){

}







