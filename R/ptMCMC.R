
#' @title Summarize the ptMCMC diagnostics
#'
#' @description Summarize the acceptance rates and trip metrics from 
#'   the saved output from the ptMCMC.
#'
#' @param rho_dist List of saved data objects from the ptMCMC estimation of
#'   changepoint locations (unless \code{nchangepoints} is 0, then 
#'   \code{NULL}).
#'
#' @return List of [1] step acceptance rates, [2] swap acceptance rates, [3]
#'   trip counts, and [4] trip rates.
#' 
#' @export 
#'
diagnose_ptMCMC <- function(rho_dist){
  if(is.null(rho_dist)){
    return(NULL)
  }
  trips <- count_trips(rho_dist$ids)
  list(step_acceptance_rate = rowMeans(rho_dist$step_accepts), 
       swap_acceptance_rate = rowMeans(rho_dist$swap_accepts), 
       trip_counts = trips$trip_counts, trip_rates = trips$trip_rates)
}

#' @title Count trips of the ptMCMC particles
#'
#' @description Count the full trips (from one extreme temperature chain to
#'   the other and back again) for each of the ptMCMC particles, as identified
#'   by their id on initialization.
#'
#' @param ids Matrix of identifiers of the particles in each chain for each
#'   iteration of the ptMCMC algorithm.
#'
#' @return List of [1] vector of counts of trips and [2] vector of rates of 
#'   trips by temperature.
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
    if (is.na(first_top)){
      trips[i] <- 0
    } else{
      last_pos <- rle(last_extreme_vector[first_top:nit])$values
      trips[i] <- sum(last_pos == "bottom")
    }
  }
  trip_rates <- trips / nit
  list(trip_counts = trips, trip_rates = trip_rates)
}


#' @title Conduct a set of among-chain swaps for the ptMCMC algorithm
#'
#' @description This function handles the among-chain swapping based on 
#'   temperature and likelihood differential. 
#'
#' @param steps Chain configuration after within-temperature steps.
#'
#' @param inputs Class \code{ptMCMC_inputs} list, containing the static inputs
#'   for use within the ptMCMC algorithm.
#'
#' @param ids The vector of integer chain ids.
#'
#' @return List of updated changepoints, log-likelihoods, and chain ids, as 
#'   well as a vector of acceptance indicators for each swap.
#'
#' @export
#'
swap_chains <- function(steps, inputs, ids){
  temps <- inputs$temps
  itemps <- 1/temps
  ntemps <- length(temps)
  revtemps <- seq(ntemps - 1, 1)
  lls <- steps$lls
  changepts <- steps$changepts
  accept_swap <- rep(FALSE, ntemps - 1)

  for (j in revtemps){
    cutoff <- exp((itemps[j] - itemps[j + 1]) * (lls[j + 1] - lls[j]))
    accept <- runif(1) < cutoff
    if (accept) {

      accept_swap[j] <- TRUE
      placeholder <- changepts[, j]
      changepts[ , j] <- changepts[, j + 1]
      changepts[ , j + 1] <- placeholder
        
      placeholder <- lls[j]
      lls[j] <- lls[j + 1]
      lls[j + 1] <- placeholder
        
      placeholder <- ids[j]
      ids[j] <- ids[j + 1]
      ids[j + 1] <- placeholder
    }
  }
  list(changepts = changepts, lls = lls, ids = ids, accept_swap = accept_swap)
}

#' @title Conduct a within-chain step of the ptMCMC algorithm
#'
#' @description This set of functions steps the chains forward one iteration 
#'   of the within-chain component of the ptMCMC algorithm. \code{step_chains}
#'   is the main function, comprised of a proposal (made by \code{prop_step}),
#'   an evaluation of that proposal (made by \code{eval_step}), and then an 
#'   update of the configuration (made by \code{take_step}). 
#'
#' @param i Integer iteration index.
#'
#' @param cpts Matrix of changepoint locations across chains.
#'
#' @param inputs Class \code{ptMCMC_inputs} list, containing the static inputs
#'   for use within the ptMCMC algorithm.
#'
#' @return \code{step_chains}: the initialized progress bar object.
#'
#' @export
#'
step_chains <- function(i, cpts, inputs){
  prop_step <- propose_step(i, cpts, inputs)
  accept_step <- eval_step(i, cpts, prop_step, inputs)
  take_step(cpts, prop_step, accept_step)
}

#' @rdname step_chains
#'
#' @return \code{propose_step}: List of changepoints and log-likelihood values
#'   for the proposal.
#'
#' @export
#'
propose_step <- function(i, cpts, inputs){

  pdist <- inputs$pdist
  ntemps <- length(inputs$temps)
  selection <- cbind(pdist$which_steps[i, ], 1:ntemps)
  prop_changepts <- cpts$changepts
  curr_changepts_s <- cpts$changepts[selection]
  prop_changepts_s <- curr_changepts_s + pdist$steps[i, ]
  if(all(is.na(prop_changepts_s))){
    prop_changepts_s <- NULL
  }
  prop_changepts[selection] <- prop_changepts_s
  mods <- proposed_step_mods(prop_changepts, inputs)
  lls <- sapply(mods, logLik)
  list(changepts = prop_changepts, lls = lls)
}

#' @rdname step_chains
#'
#' @param prop_step Proposed step output from \code{propose_step}.
#'
#' @return \code{propose_step}: \code{logical} vector indicating if each 
#'   chain's proposal was accepted.
#'
#' @export
#'
eval_step <- function(i, cpts, prop_step, inputs){
  temps <- inputs$temps
  ntemps <- length(temps)
  itemps <- 1 / temps
  runif(ntemps) < exp((prop_step$lls - cpts$lls) * itemps)
}

#' @rdname step_chains
#'
#' @param accept_step \code{logical} indicator of acceptance of each chain's
#'   proposed step.
#'
#' @return \code{take_step}: list of changepoints, log-likelihoods, and
#'   logical indicators of acceptance for each chain.
#'
#' @export
#'
take_step <- function(cpts, prop_step, accept_step){
  changepts <- cpts$changepts
  lls <- cpts$lls
  changepts[ , accept_step] <- prop_step$changepts[ , accept_step]
  lls[accept_step] <- prop_step$lls[accept_step]
  list(changepts = changepts, lls = lls, accept_step = accept_step)
}

#' @title Fit the chunk-level models to a time series, given a set of 
#'   proposed changepoints.
#'
#' @description This function wraps around \code{TS_memo} to provide simpler
#'   interface within the ptMCMC algorithm. 
#'
#' @param prop_changepts Matrix of proposed changepoints across chains.
#'
#' @param inputs Class \code{ptMCMC_inputs} list, containing the static inputs
#'   for use within the ptMCMC algorithm.
#'
#' @return List of models associated with the proposed step.
#'
#' @export
#'
proposed_step_mods <- function(prop_changepts, inputs){ 

  data <- inputs$data
  formula <- inputs$formula
  weights <- inputs$weights
  TS_memo <- inputs$TS_memo
  ntemps <- length(inputs$temps)
  control <- inputs$control
  out <- vector("list", length = ntemps)
  for (i in 1:ntemps){
    out[[i]] <- TS_memo(data, formula, prop_changepts[ , i], weights, control)
  }
  out
}


#' @title Initialize and update the chain ids
#'
#' @description \code{prep_ids} creates and \code{update_ids} updates 
#'   the active vectort of identities (ids) for each of the chains in the 
#'   ptMCMC algorithm.
#'
#' @param control Class \code{TS_controls} list, holding control parameters
#'   for the Time Series model including the parallel tempering Markov Chain 
#'   Monte Carlo (ptMCMC) controls, generated by 
#'   \code{\link{TS_controls_list}}. 
#'
#' @return The vector of chain ids.
#'
#' @export
#'
prep_ids <- function(control){
  if (!is.numeric(control$ntemps) || any(control$ntemps %% 1 != 0)){
    stop("ntemps must be integer-valued")
  }
  1:control$ntemps
}

#' @rdname prep_ids
#'
#' @param ids The existing vector of chain ids.
#'
#' @param swaps Chain configuration after among-temperature swaps.
#'
#' @export
#'
update_ids <- function(ids, swaps){
  swaps$ids
}

#' @title Prepare the inputs for the ptMCMC algorithm
#'
#' @description Package the static inputs (controls and data structures) used
#'   by the ptMCMC algorithm.
#'
#' @param data Class \code{data.frame} object including [1] the time variable
#'   (indicated in \code{control}), [2] the predictor variables (required by
#'   \code{formula}) and [3], the multinomial response variable (indicated
#'   in \code{formula}).  
#'
#' @param formula \code{formula} describing the continuous change. Any 
#'   predictor variable included must also be a column in the
#'   \code{data}.  Any (multinomial) response variable must also be a set of
#'   columns in \code{data}. 
#'
#' @param nchangepoints Integer corresponding to the number of 
#'   change points to include in the model. 0 is a valid input (corresponding
#'   to no change points, so a singular time series model), and the current 
#'   implementation can reasonably include up to 6 change points. The 
#'   number of change points is used to dictate the segementation of the data  
#'   for each continuous model and each LDA model.
#'
#' @param weights Optional class \code{numeric} vector of weights for each 
#'   document. Corresponds to the vector \strong{\eqn{v}} in the math 
#'   description.
#'
#' @param control Class \code{TS_controls} list, holding control parameters
#'   for the Time Series model including the parallel tempering Markov Chain 
#'   Monte Carlo (ptMCMC) controls, generated by 
#'   \code{\link{TS_controls_list}}.
#'
#' @return Class \code{ptMCMC_inputs} list, containing the static inputs for
#'   use within the ptMCMC algorithm.
#'
#' @export
#'
prep_ptMCMC_inputs <- function(data, formula, nchangepoints, weights, 
                               control){
  check_timename(data, control$timename)
  check_formula(data, formula)
  check_weights(weights)
  check_nchangepoints(nchangepoints)
  control$selector <- NULL
  control$measurer <- NULL
  out <- list(control = control, temps = prep_temp_sequence(control), 
              pdist = prep_proposal_dist(nchangepoints, control),
              formula = formula, weights = weights, data = data, 
              TS_memo = memoise_fun(multinom_TS, control$memoise))
  class(out) <- c("ptMCMC_inputs", "list")
  out
}


#' @title Pre-claculate the proposal distribution for the ptMCMC algorithm
#'
#' @description Calculate the proposal distribution in advance of actually
#'   running the ptMCMC algorithm in order to decrease computation time.
#'   The proposal distribution is a joint of three distributions:
#'   [1] a multinomial distribution selecting among the change points within
#'   the chain, [2] a binomial distribution selecting the direction of the 
#'   step of the change point (earlier or later in the time series), and 
#'   [3] a geometric distribution selecting the magnitude of the step.
#'
#' @param nchangepoints Integer corresponding to the number of 
#'   change points to include in the model. 0 is a valid input (corresponding
#'   to no change points, so a singular time series model), and the current 
#'   implementation can reasonably include up to 6 change points. The 
#'   number of change points is used to dictate the segementation of the data  
#'   for each continuous model and each LDA model.
#'
#' @param control Class \code{TS_controls} list, holding control parameters
#'   for the Time Series model including the parallel tempering Markov Chain 
#'   Monte Carlo (ptMCMC) controls, generated by 
#'   \code{\link{TS_controls_list}}. Currently relevant here is 
#'   \code{magnitude} (referenced as \eqn{\kappa} in the math description),
#'   which controls the magnitude of the step size (is the average of the
#'   geometric distribution). 
#'
#' @return List of two matrices: [1] the size of the proposed step for each
#'   iteration of each chain, [2] the identity of the change point location 
#'   to be shifted by the step for each iteration of each chain.
#'
#' @export
#'
prep_proposal_dist <- function(nchangepoints, control = TS_controls_list()){
  check_nchangepoints(nchangepoints)
  check_control(control, "TS_controls")
  ntemps <- control$ntemps
  nit <- control$nit
  if(nchangepoints == 0){
    steps <- matrix(0, nrow = nit, ncol = ntemps)
    which_steps <- matrix(numeric(0), nrow = nit, ncol = ntemps)
  } else{
    magnitude <- control$magnitude 
    step_signs <- sample(c(-1, 1), nit * ntemps, replace = TRUE)
    step_magnitudes <- 1 + rgeom(nit * ntemps, 1 / magnitude)
    steps <- matrix(step_signs * step_magnitudes, nrow = nit)
    which_steps <- sample.int(nchangepoints, nit * ntemps, replace = TRUE)
    which_steps <- matrix(which_steps, nrow = nit)
  }
  list(steps = steps, which_steps = which_steps)
}

#' @title Prepare and update the data structures to save the ptMCMC output
#'
#' @description \code{prep_saves} creates the data structure used to save the 
#'   output from each iteration of the ptMCMC algorithm, which is added via
#'   \code{update_saves}. Once the ptMCMC is complete, the saved data objects
#'   are then processed (burnin iterations are dropped and the remaining
#'   iterations are thinned) via \code{process_saves}.
#'
#' @param nchangepoints Integer corresponding to the number of 
#'   change points to include in the model. 0 is a valid input (corresponding
#'   to no change points, so a singular time series model), and the current 
#'   implementation can reasonably include up to 6 change points. The 
#'   number of change points is used to dictate the segementation of the data  
#'   for each continuous model and each LDA model.
#'
#' @param control Class \code{TS_controls} list, holding control parameters
#'   for the Time Series model including the parallel tempering Markov Chain 
#'   Monte Carlo (ptMCMC) controls, generated by 
#'   \code{\link{TS_controls_list}}.
#'
#' @return List of saved ptMCMC objects: changepoints, log-likelihoods, 
#'   chain ids, step acceptances, and swap acceptances.
#'
#' @export
#'
prep_saves <- function(nchangepoints, control = TS_controls_list()){
  check_nchangepoints(nchangepoints)
  check_control(control, "TS_controls")
  ntemps <- control$ntemps
  nit <- control$nit
  cpts <- array(NA, c(nchangepoints, ntemps, nit))
  lls <- matrix(NA, ntemps, nit)
  ids <- matrix(NA, ntemps, nit)
  step_accepts <- matrix(FALSE, ntemps, nit)
  swap_accepts <- matrix(FALSE, ntemps - 1, nit)
  list(cpts = cpts, lls = lls, ids = ids, step_accepts = step_accepts, 
       swap_accepts = swap_accepts)
}

#' @rdname prep_saves
#'
#' @param i Integer iteration index. 
#'
#' @param saves The existing list of saved data objects.
#'
#' @param steps Chain configuration after within-temperature steps.
#'
#' @param swaps Chain configuration after among-temperature swaps.
#'
#' @export
#'
update_saves <- function(i, saves, steps, swaps){
  saves$cpts[ , , i] <- swaps$changepts
  saves$lls[ , i] <- swaps$lls
  saves$ids[ , i] <- swaps$ids
  saves$step_accepts[ , i] <- steps$accept_step
  saves$swap_accepts[ , i] <- swaps$accept_swap
  saves
}

#' @rdname prep_saves
#'
#' @export
#'
process_saves <- function(saves, control){
  nit <- control$nit
  iters <- 1:nit
  if (control$burnin > 0){
    iters <- iters[-(1:control$burnin)]
  }
  niters <- length(iters)
  thin_interval <- ceiling(1/control$thin_frac)
  iters_thinned <- seq(1, niters, by = thin_interval)
  dims <- c(dim(saves$cpts)[1:2], length(iters_thinned))
  saves$cpts <- array(saves$cpts[ , , iters_thinned], dim = dims)
  saves$lls <- saves$lls[, iters_thinned]
  saves$ids <- saves$ids[, iters_thinned]
  saves$step_accepts <- saves$step_accepts[ , iters_thinned]
  saves$swap_accepts <- saves$swap_accepts[ , iters_thinned]
  saves
}

#' @title Initialize and update the changepoint matrix used in the ptMCMC
#'   algorithm
#' 
#' @description Each of the chains is initialized by \code{prep_cpts} using a 
#'   draw from the available times (i.e. assuming a uniform prior), the best 
#'   fit (by likelihood) draw is put in the focal chain with each subsequently 
#'   worse fit placed into the subsequently hotter chain. \code{update_cpts}
#'   updates the changepoints after every iteration in the ptMCMC algorithm.
#'
#' @param data Class \code{data.frame} object including [1] the time variable
#'   (indicated in \code{control}), [2] the predictor variables (required by
#'   \code{formula}) and [3], the multinomial response variable (indicated
#'   in \code{formula}). 
#'
#' @param formula \code{formula} describing the continuous change. Any 
#'   predictor variable included must also be a column in the
#'   \code{data}.  Any (multinomial) response variable must also be a set of
#'   columns in \code{data}. 
#'
#' @param nchangepoints Integer corresponding to the number of 
#'   change points to include in the model. 0 is a valid input (corresponding
#'   to no change points, so a singular time series model), and the current 
#'   implementation can reasonably include up to 6 change points. The 
#'   number of change points is used to dictate the segementation of the data  
#'   for each continuous model and each LDA model.
#'
#' @param weights Optional class \code{numeric} vector of weights for each 
#'   document. Corresponds to the vector \strong{\eqn{v}} in the math 
#'   description.
#'
#' @param control Class \code{TS_controls} list, holding control parameters
#'   for the Time Series model including the parallel tempering Markov Chain 
#'   Monte Carlo (ptMCMC) controls, generated by 
#'   \code{\link{TS_controls_list}}.
#'
#' @return List of [1] matrix of change points (rows) for each temperature 
#'   (columns) and [2] vector of log-likelihood values for each of the chains.
#'
#' @export
#'
prep_cpts <- function(data, formula, nchangepoints, weights, 
                            control = TS_controls_list()){

  check_formula(data, formula)
  check_nchangepoints(nchangepoints)
  check_weights(weights)
  check_control(control, "TS_controls")
  temps <- prep_temp_sequence(control)
  ntemps <- length(temps)
  min_time <- min(data[ , control$timename])
  max_time <- max(data[ , control$timename])
  times <- seq(min_time, max_time, 1)
  avail_times <- times[-c(1, length(times))]
  cps <- matrix(NA, nrow = nchangepoints, ncol = ntemps)
  for (i in 1:ntemps){
    cp_times <- sort(sample(avail_times, nchangepoints, replace = FALSE))
    cps[ , i] <- cp_times
  }
  lls <- rep(NA, ntemps)
  for (i in 1:ntemps){
    lls[i] <- multinom_TS(data, formula, cps[ , i], weights, control)$logLik
  }  
  cps <- cps[ , order(lls, decreasing = TRUE), drop = FALSE]
  lls <- sort(lls, decreasing = TRUE)

  out <- list(cps, lls)
  names(out) <- c("changepts", "lls")
  out
}

#' @rdname prep_cpts
#'
#' @param cpts The existing matrix of changepoints.
#'
#' @param swaps Chain configuration after among-temperature swaps.
#'
#' @export
#'
update_cpts <- function(cpts, swaps){
  list(changepts = swaps$changepts, lls = swaps$lls)
}

#' @title Prepare the ptMCMC temperature sequence
#'
#' @description Create the series of temperatures used in the ptMCMC algorithm
#'
#' @param control Class \code{TS_controls} list, holding control parameters
#'   for the Time Series model including the parallel tempering Markov Chain 
#'   Monte Carlo (ptMCMC) controls, generated by 
#'   \code{\link{TS_controls_list}}.
#'
#' @return Vector of temperatures.
#'
#' @export
#'
prep_temp_sequence <- function(control = TS_controls_list()){
  check_control(control, "TS_controls")
  ntemps <- control$ntemps
  penultimate_temp <- control$penultimate_temp
  ultimate_temp <- control$ultimate_temp
  q <- control$q
  sequence <- seq(0, log2(penultimate_temp), length.out = ntemps - 1)
  log_temps <- sequence^(1 + q) / log2(penultimate_temp)^q
  c(2^(log_temps), ultimate_temp) 
}
