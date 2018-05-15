#' @title Prepare the temperatures 
#'
#' @description based on general controls for parallel tempering
#'
#' @param ntemps number of temperatures
#' @param penultimate_temp penultimate temperature (second hottest chain)
#' @param ultimate_temp ultimate temperature (hottest chain)
#' @param k the exponent controlling the temperature sequence: 0 implies 
#'   geometric sequence, 1 implies squaring before exponentiating. Use larger 
#'   values if the cooler chains aren't swapping enough.
#' @param ... additional arguments to be passed to subfunctions
#'  
#' @return temperatures
#'
#' @export
#'
prep_temps <- function(ntemps = 6, penultimate_temp = 2^6, 
                       ultimate_temp = 1E10, k = 0, ...){

  sequence <- seq(0, log2(penultimate_temp), length.out = ntemps - 1)
  log_temps <- sequence^(1 + k) / log2(penultimate_temp)^k
  temps <- 2^(log_temps)
  temps <- c(temps, ultimate_temp) 
  return(temps)
}  

#' @title Remove the first defined number of burn-in samples
#'
#' @description Flexible to the 2- or 3-D structure of objects to be trimmed
#'
#' @param object Object to be trimmed
#' @param burnin number of samples to remove from the beginning
#' @return object with burn-in samples removed
#' 
#' @export 
#'
remove_burnin <- function(object, burnin = 0){

  if (burnin %% 1 !=0 | burnin < 0){
    stop("burnin must be a non-negative integer value.")
  }
  if (burnin > 0){
    to_drop <- -(1:burnin)
    if (length(dim(object)) == 2){
      dobj1 <- dim(object)[1]
      dobj2 <- dim(object)[2]
      dsobjt <- c(dobj1, dobj2 - burnin)
      object_trimmed <- array(object[1:dobj1, to_drop], dim = dsobjt)
    } else if (length(dim(object)) == 3){
      dobj1 <- dim(object)[1]
      dobj2 <- dim(object)[2]
      dobj3 <- dim(object)[3]
      dsobjt <- c(dobj1, dobj2, dobj3 - burnin)
      object_trimmed <- array(object[1:dobj1, 1:dobj2, to_drop], dim = dsobjt)
    }
  } else{
    object_trimmed <- object
  }
  out <- object_trimmed
  return(out)
}

#' @title Count trips
#'
#' Function for counting trips
#'
#' @param ids Matrix of identifiers of the initial temperature for particles
#' @return list of counts and rates of trips by temperature
#' 
#' @export 
#'
count_trips <- function(ids){

  nit <- ncol(ids)
  ntemps <- nrow(ids)
  last_extreme <- NA
  last_extreme_vector <- numeric(nit)
  trips <- numeric(ntemps)
  for(i in 1:ntemps){
    for(j in 1:nit){
      if(ids[1, j] == i){
        last_extreme <- "bottom"
      }
      if(ids[ntemps, j] == i){
        last_extreme <- "top"
      }
      last_extreme_vector[j] <- last_extreme
    }
    first_top <- match("top", last_extreme_vector)
    last_pos <- rle(last_extreme_vector[first_top:nit])$values
    trips[i] <- sum(last_pos == "bottom")
  }
  trip_rates <- trips / nit
  out <- list("trip_counts" = trips, "trip_rates" = trip_rates)
  return(out)
}

