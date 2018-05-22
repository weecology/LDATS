#' @title Fit multinomial time chunk model
#'
#' @description Fit a multinomial regression model to topic model data for
#'   dates within a defined chunk of time (start_date, end_date] with  
#'   covariate impacts assuming no temporal autocorrelation by default
#'
#' @param data data frame including the predictor and response variables
#' @param formula Right Hand Side of the continuous time formula as a 
#'   character vector
#' @param start_time start time for the chunk 
#' @param end_time end time for the chunk
#' @param weights weights 
#' @param ... other arguments to be passed to the multinom function
#' @return fitted model for the chunk
#' 
#' @export 
#'
multinom_chunk <- function(data, formula, start_time, end_time, weights, 
                           ...) {
  formula <- as.formula(paste("gamma ~", formula))
  mod <- multinom(formula, data, weights, 
           subset = data$time > start_time & data$time <= end_time, 
           trace = FALSE, ...) 
  return(mod)
}

#' @title Fit multinomial time series model
#'
#' @description Fit a multinomial regression model with covariate impacts 
#'   assuming no temporal autocorrelation 
#'
#' @param data data frame including the predictor and response variables
#' @param formula Right Hand Side of the continuous time formula as a 
#'   character vector
#' @param changepoints selections of the change points
#' @param weights weights 
#' @param ... other arguments to be passed to subfunctions
#' @return list of chunk-level model fits and the total log likelihood 
#'   combined across all the chunks
#' @examples
#' NA
#' @export 
#'
multinom_ts <- function(data, formula, changepoints = NULL, weights, ...){

  chunk_memo <- memoise(multinom_chunk)

  nchunks <- length(changepoints) + 1
  start_times <- c(min(data$time) - 1, changepoints)
  end_times <- c(changepoints, max(data$time) + 1)
  last_time <- max(data$time)

  nobs <- length(data$time)
  time_check <- any(changepoints <=0) | any(changepoints >= last_time)
  sort_check <- is.unsorted(changepoints, strictly = TRUE)

  if (time_check | sort_check){
    output <- list("chunk models" = NA, "logLik" = -Inf)
    return(output)
  }

  mods <- vector("list", length = nchunks)
  ll <- rep(NA, nchunks)
  for (i in 1:nchunks){
    mods[[i]] <- chunk_memo(data, formula, start_time = start_times[i], 
                   end_time = end_times[i], weights)
    ll[i] <- logLik(mods[[i]])
  }
  total_ll <- sum(ll)
  names(mods) <- sprintf("%s %d %s", "chunk", 1:nchunks, "model")
  output <- list("chunk models" = mods, "logLik" = total_ll)

  return(output)
}

#' @title Prepare the changepoints for the MTS algorithm
#'
#' @description includes sorting by logLik
#'
#' @param data data frame including the predictor and response variables
#' @param formula formula for the continuous change
#' @param ntemps number of temperatures
#' @param nchangepoints number of change points to include in the model
#' @param weights weights 
#'  
#' @return list of [1] matrix of change points (rows) for each temperature 
#'   (columns) and [2] vector of logLiks
#'
#' @export
#'
prep_changepts<- function(data, formula, ntemps, nchangepoints, weights){

  min_time <- min(data$time)
  max_time <- max(data$time)
  times <- seq(min_time, max_time, 1)
  avail_times <- times[-c(1, length(times))]
  cps <- matrix(NA, nrow = nchangepoints, ncol = ntemps)
  for (i in 1:ntemps){
    cp_times <- sort(sample(avail_times, nchangepoints, replace = FALSE))
    cps[ , i] <- cp_times
  }
  lls <- rep(NA, ntemps)
  for (i in 1:ntemps){
    lls[i] <- multinom_ts(data, formula, cps[ , i], weights)$logLik
  }  
  cps <- cps[ , order(lls, decreasing = TRUE), drop = FALSE]
  lls <- sort(lls, decreasing = TRUE)

  out <- list(cps, lls)
  names(out) <- c("changepts", "lls")
  return(out)
}

#' @title Calculate the proposal distribution for the time series model
#'
#' @description based on the provided magnitude
#'
#' @param nit number of iterations used
#' @param ntemps number of temperatures used
#' @param nchangepoints number of change points to include in the model
#' @param magnitude the inverse of the probability used for the geometic 
#'   distribution to generate the kick magnitudes
#' @return list of two matrices: [1] kicks and [2] which_kicked
#'
#' @export
#'
proposal_dist <- function(nit, ntemps, nchangepoints, magnitude){
  kick_signs <- sample(c(-1, 1), nit * ntemps, replace = TRUE)
  kick_magnitudes <- 1 + rgeom(nit * ntemps, 1 / magnitude)
  kicks <- matrix(kick_signs * kick_magnitudes, nrow = nit)
  if(nchangepoints == 0){
    which_kicked <- matrix(numeric(0), nrow = nit, ncol = ntemps)
  }else{
    which_kicked <- sample.int(nchangepoints, nit * ntemps, replace = TRUE)
    which_kicked <- matrix(which_kicked, nrow = nit)
  }
  out <- list("kicks" = kicks, "which_kicked" = which_kicked)
  return(out)
}

#' @title Multinomial Time Series analysis of a topic model classification
#'
#' @description Full MTS analyses
#'
#' @param data data frame including the predictor and response variables
#' @param formula formula for the continuous change equation
#' @param nchangepoints number of change points to include in the model
#' @param weights weights 
#' @param nit number of iterations used
#' @param magnitude scaling for the kick magnitude used in the proposal dist
#' @param burnin number of iterations to removed from the front-end of the 
#'   simulations prior to summarizing the distributions of change points
#' @param ... additional arguments to be passed to subfunctions
#' @return an object of class MTS
#'
#' @export
#'
MTS <- function(data, formula = ~1, nchangepoints = 1, 
                weights = NULL, nit = 1e3, magnitude = 12, burnin = 1e2, ...){

  character_formula <- as.character(formula)
  formula <- character_formula[length(character_formula)]
  ts_memo <- memoise(multinom_ts)

  if(nchangepoints == 0){
    nit <- 1
  }
  temps <- prep_temps(...)
  betas <- 1 / temps
  ntemps <- length(betas)

  prep_cpts <- prep_changepts(data, formula, ntemps, nchangepoints, 
                 weights)
  changepts <- prep_cpts$changepts
  lls <- prep_cpts$lls

  saved_cps <- array(NA, c(nchangepoints, ntemps, nit))
  saved_lls <- matrix(NA, ntemps, nit)
  saved_ids <- matrix(NA, ntemps, nit)
  saved_accepts <- matrix(NA, ntemps, nit)
  accept_rate <- 0
  temp_ids <- 1:ntemps
  swap_accepted <- matrix(FALSE, nit, ntemps - 1)

  pdist <- proposal_dist(nit, ntemps, nchangepoints, magnitude)
 
  pbform <- "  [:bar] :percent eta: :eta"
  pb <- progress_bar$new(pbform, nit, clear = FALSE, width = 60)

  for (i in 1:nit){
  
    pb$tick()
 
    selection <- cbind(pdist$which_kicked[i, ], 1:ntemps)
    prop_changepts <- changepts
    curr_changepts_s <- changepts[selection]
    prop_changepts_s <- curr_changepts_s + pdist$kicks[i, ]
    if(all(is.na(prop_changepts_s))){
      prop_changepts_s <- NULL
    }
    prop_changepts[selection] <- prop_changepts_s
    prop_lls <- lls

    for (j in 1:ntemps){
      mod <- ts_memo(data, formula, prop_changepts[ , j], weights)
      prop_lls[j] <- mod$logLik
    }

    accepts <- runif(ntemps) < exp((prop_lls - lls) * betas)
    saved_accepts[ , i] <- accepts
    changepts[ , accepts] <- prop_changepts[ , accepts]
    lls[accepts] <- prop_lls[accepts]
    
    revtemps <- seq(ntemps - 1, 1)
    for (j in revtemps){
      cutoff <- exp((betas[j] - betas[j + 1]) * (lls[j + 1] - lls[j]))
      accept_swap <- runif(1) < cutoff
      if (accept_swap == TRUE) {

        swap_accepted[i, j] <- TRUE
        placeholder <- changepts[, j]
        changepts[ , j] <- changepts[, j + 1]
        changepts[ , j + 1] <- placeholder
        
        placeholder <- lls[j]
        lls[j] <- lls[j + 1]
        lls[j + 1] <- placeholder
        
        placeholder <- temp_ids[j]
        temp_ids[j] <- temp_ids[j + 1]
        temp_ids[j + 1] <- placeholder
      }
    }

    saved_cps[ , , i] <- changepts
    saved_ids[ , i] <- temp_ids
    saved_lls[ , i] <- lls
  }

  if (nchangepoints > 0){
    cps <- remove_burnin(saved_cps, burnin)
    ids <- remove_burnin(saved_ids, burnin)
    lls <- remove_burnin(saved_lls, burnin)
    accepts <- remove_burnin(saved_accepts, burnin)
    swap_accept <- t(remove_burnin(t(swap_accepted), burnin))

    cps_t1 <- t(array(cps[ , 1, ], dim = dim(cps)[c(1,3)]))

    cp_summary <- summarize_cps(cps_t1)
    cp_ccmat <- ccmat(cps_t1, lag = 0)

    MCMCsetup <- list(ntemps, temps, magnitude)
    names(MCMCsetup) <- c("ntemps", "temperatures", "magnitude")

    swap_rates <- colMeans(swap_accept)
    accept_rates <- rowMeans(accepts)
    trips <- count_trips(ids)

    MCMCdiagnostics <- list(accept_rates, swap_rates, trips$trip_counts, 
                         trips$trip_rate)
    names(MCMCdiagnostics) <- c("acceptance_rates", "swapping_rates", 
                                "trip_counts", "trip_rates")
  } else{
    lls <- saved_lls
    cps_t1 <- NULL
    cp_summary <- NULL
    cp_ccmat <- NULL
    MCMCsetup <- NULL
    MCMCdiagnostics <- NULL
  }

  out <- list()
  out$call <- match.call()
  out$formula <- formula
  out$nchangepoints <- nchangepoints 
  out$data <- data
  out$weights <- weights
  out$MCMCsetup <- MCMCsetup
  out$lls <- lls[1, ]
  out$lls_full <- lls
  out$cps <- cps_t1
  out$MCMCdiagnostics <- MCMCdiagnostics 
  out$summary <- cp_summary
  out$cor <- cp_ccmat
  class(out) <- c("MTS", "list")
  attr(out, "hidden") <- c("data", "weights", "MCMCsetup", "lls", "lls_full",
                           "cps", "MCMCdiagnostics")

  return(out)
}

#' @export
#'
print.MTS <- function(x, ...){
  hid <- attr(x, "hidden")
  notHid <- !names(x) %in% hid
  print(x[notHid])
}

#' @title Summarize the change point estimations
#'
#' @description Function for summarizing change point estimations
#'
#' @param cps change point estimates
#' @param prob probability used for the interval
#' @return model summary table for the change point locations  
#'
#' @export
#'
summarize_cps <- function(cps, prob = 0.95){

  Mean <- round(apply(cps, 2, mean), 2)
  Median <- round(apply(cps, 2, median), 2)
  SD <- round(apply(cps, 2, sd), 2)
  MCMCerr <- round(SD / sqrt(nrow(cps)), 4)
  HPD <- HPDinterval(as.mcmc(cps), prob = prob)
  Lower <- HPD[ , "lower"]
  Upper <- HPD[ , "upper"]
  AC10 <- t(round(autocorr.diag(as.mcmc(cps), lag = 10), 4))
  ESS <- effectiveSize(cps)
  out <- data.frame(Mean, Median, Lower, Upper, SD, MCMCerr, AC10, ESS)
  colnames(out)[7] <- "AC10"
  rownames(out) <- sprintf("Changepoint_%d", 1:nrow(out))
  return(out)
}

#' @title Measure the cross-correlation among change points
#'
#' @description Measure the cross-correlation among change points
#'
#' @param cps change point estimates
#' @param lag lag to be used in the correlation estimation
#' @return matrix of among-changepoint correlations
#'
#' @export
#'
ccmat <- function(cps, lag = 0){

  CC <- cps %>%
        as.mcmc() %>%
        autocorr(lags = lag) %>%
        round(4)
  out <- matrix(CC[1, , ], dim(CC)[2], dim(CC)[2]) 
  colnames(out) <- sprintf("CP_%d", 1:dim(CC)[2])
  rownames(out) <- sprintf("CP_%d", 1:dim(CC)[2])
  return(out)
}

