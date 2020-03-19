

prep_cpts <- function(TS, control = list()){

  data <- TS$data$train$ts_data
  temps <- prep_temp_sequence(control)
  ntemps <- length(temps)
  min_time <- min(data[ , timename])
  max_time <- max(data[ , timename])
  times <- seq(min_time, max_time, 1)
  avail_times <- times[-c(1, length(times))]
  cps <- matrix(NA, nrow = TS$nchangepoints, ncol = ntemps)
  for (i in 1:ntemps){
    cp_times <- sort(sample(avail_times, TS$nchangepoints, replace = FALSE))
    cps[ , i] <- cp_times
  }
  lls <- rep(NA, ntemps)
  for (i in 1:ntemps){
    fun <- eval(parse(text = paste0(TS$response, "_TS")))
    fun <- memoise_fun(fun, control$memoise)
    args <- list(data = data, formula = TS$formula, changepoints = cps[ , i], 
                 timename = TS$timename, weights = TS$weights, 
                 control = control)
    modfit <- soft_call(fun, args, TRUE)
    lls[i] <- modfit$logLik
  }  
  cps <- cps[ , order(lls, decreasing = TRUE), drop = FALSE]
  lls <- sort(lls, decreasing = TRUE)

  out <- list(cps, lls)
  names(out) <- c("changepts", "lls")
  out
}


prep_temp_sequence <- function(control = list()){
  ntemps <- control$ntemps
  penultimate_temp <- control$penultimate_temp
  ultimate_temp <- control$ultimate_temp
  q <- control$q
  sequence <- seq(0, log2(penultimate_temp), length.out = ntemps - 1)
  log_temps <- sequence^(1 + q) / log2(penultimate_temp)^q
  c(2^(log_temps), ultimate_temp) 
}


prep_saves <- function(TS, control = list()){
  nchangepoints <- TS$nchangepoints
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



prep_ptMCMC_inputs <- function(TS, control = list()){
    
  out <- list(control = control, 
              temps = prep_temp_sequence(control), 
              pdist = prep_proposal_dist(TS$nchangepoints, control),
              formula = TS$formula, 
              weights = TS$weights, 
              data = TS$data$train$ts_data, 
              response = TS$response,
              timename = TS$timename)
  class(out) <- c("ptMCMC_inputs", "list")
  out
}



prep_proposal_dist <- function(nchangepoints, control = list()){
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


#' @title Calculate ptMCMC summary diagnostics
#'
#' @description Summarize the step and swap acceptance rates as well as trip
#'   metrics from the saved output of a ptMCMC estimation. 
#'
#' @details Within-chain step acceptance rates are averaged for each of the
#'   chains from the raw step acceptance histories 
#'   (\code{ptMCMCout$step_accepts}) and between-chain swap acceptance rates
#'   are similarly averaged for each of the neighboring pairs of chains from
#'   the raw swap acceptance histories (\code{ptMCMCout$swap_accepts}).
#'   Trips are defined as movement from one extreme chain to the other and
#'   back again (Katzgraber \emph{et al.} 2006). Trips are counted and turned 
#'   to per-iteration rates using \code{\link{count_trips}}.
#'   \cr \cr 
#'   This function was first designed to work within \code{\link{TS}} and 
#'   process the output of \code{\link{est_changepoints}}, but has been 
#'   generalized and would work with any output from a ptMCMC as long as 
#'   \code{ptMCMCout} is formatted properly.
#'
#' @param ptMCMCout Named \code{list} of saved data objects from a ptMCMC 
#'   estimation including elements named \code{step_accepts} (matrix of 
#'   \code{logical} outcomes of each step; rows: chains, columns: iterations),
#'   \code{swap_accepts} (matrix of \code{logical} outcomes of each swap;
#'   rows: chain pairs, columns: iterations), and \code{ids} (matrix of 
#'   particle identifiers; rows: chains, columns: iterations). 
#'   \code{ptMCMCout = NULL} indicates no use of ptMCMC and so the function
#'   returns \code{NULL}.
#'
#' @return \code{list} of [1] within-chain average step acceptance rates 
#'   (\code{$step_acceptance_rate}), [2] average between-chain swap acceptance
#'   rates (\code{$swap_acceptance_rate}), [3] within particle trip counts 
#'   (\code{$trip_counts}), and [4] within-particle average trip rates 
#'   (\code{$trip_rates}).
#' 
#' @references 
#'   Katzgraber, H. G., S. Trebst, D. A. Huse. And M. Troyer. 2006. 
#'   Feedback-optimized parallel tempering Monte Carlo. \emph{Journal of 
#'   Statistical Mechanics: Theory and Experiment} \strong{3}:P03018
#'   \href{https://bit.ly/2LICGXh}{link}.
#'
#' @examples 
#' \donttest{
#'   data(rodents)
#'   document_term_table <- rodents$document_term_table
#'   document_covariate_table <- rodents$document_covariate_table
#'   LDA_models <- LDA_set(document_term_table, topics = 2)[[1]]
#'   data <- document_covariate_table
#'   data$gamma <- LDA_models@gamma
#'   weights <- document_weights(document_term_table)
#'   data <- data[order(data[,"newmoon"]), ]
#'   rho_dist <- est_changepoints(data, gamma ~ 1, 1, "newmoon", 
#'                                weights, TS_control())
#'   diagnose_ptMCMC(rho_dist)
#' }
#' 
#' @export 
#'
diagnose_ptMCMC <- function(ptMCMCout){
  if(is.null(ptMCMCout)){
    return(NULL)
  }
  trips <- count_trips(ptMCMCout$ids)
  list(step_acceptance_rate = rowMeans(ptMCMCout$step_accepts), 
       swap_acceptance_rate = rowMeans(ptMCMCout$swap_accepts), 
       trip_counts = trips$trip_counts, trip_rates = trips$trip_rates)
}

#' @title Count trips of the ptMCMC particles
#'
#' @description Count the full trips (from one extreme temperature chain to
#'   the other and back again; Katzgraber \emph{et al.} 2006) for each of the
#'   ptMCMC particles, as identified by their id on initialization.
#'   \cr \cr
#'   This function was designed to work within \code{\link{TS}} and process
#'   the output of \code{\link{est_changepoints}} as a component of 
#'   \code{\link{diagnose_ptMCMC}}, but has been generalized
#'   and would work with any output from a ptMCMC as long as \code{ids}
#'   is formatted properly.
#'
#' @param ids \code{matrix} of identifiers of the particles in each chain for 
#'   each iteration of the ptMCMC algorithm (rows: chains, 
#'   columns: iterations).
#'
#' @return \code{list} of [1] \code{vector} of within particle trip counts 
#'   (\code{$trip_counts}), and [2] \code{vector} of within-particle average 
#'   trip rates (\code{$trip_rates}).
#' 
#' @references 
#'   Katzgraber, H. G., S. Trebst, D. A. Huse. And M. Troyer. 2006. 
#'   Feedback-optimized parallel tempering Monte Carlo. \emph{Journal of 
#'   Statistical Mechanics: Theory and Experiment} \strong{3}:P03018
#'   \href{https://bit.ly/2LICGXh}{link}.
#'
#' @examples 
#' \donttest{
#'   data(rodents)
#'   document_term_table <- rodents$document_term_table
#'   document_covariate_table <- rodents$document_covariate_table
#'   LDA_models <- LDA_set(document_term_table, topics = 2)[[1]]
#'   data <- document_covariate_table
#'   data$gamma <- LDA_models@gamma
#'   weights <- document_weights(document_term_table)
#'   data <- data[order(data[,"newmoon"]), ]
#'   rho_dist <- est_changepoints(data, gamma ~ 1, 1, "newmoon", weights,
#'                                TS_control())
#'   count_trips(rho_dist$ids)
#' }
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
#'   temperatures and likelihood differentials.  
#'   \cr \cr
#'   This function was designed to work within \code{\link{TS}} and 
#'   specifically \code{\link{est_changepoints}}. It is still hardcoded to do
#'   so, but has the capacity to be generalized to work with any estimation
#'   via ptMCMC with additional coding work.
#'
#' @details The ptMCMC algorithm couples the chains (which are 
#'   taking their own walks on the distribution surface) through "swaps", 
#'   where neighboring chains exchange configurations (Geyer 1991, Falcioni 
#'   and Deem 1999) following the Metropolis criterion (Metropolis 
#'   \emph{et al.} 1953). This allows them to share information and search the
#'   surface in combination (Earl and Deem 2005). 
#'
#' @param chainsin Chain configuration to be evaluated for swapping.
#'
#' @param inputs Class \code{ptMCMC_inputs} list, containing the static inputs
#'   for use within the ptMCMC algorithm.
#'
#' @param ids The vector of integer chain ids.
#'
#' @return \code{list} of updated change points, log-likelihoods, and chain
#'   ids, as well as a vector of acceptance indicators for each swap.
#'
#' @references
#'   Earl, D. J. and M. W. Deem. 2005. Parallel tempering: theory, 
#'   applications, and new perspectives. \emph{Physical Chemistry Chemical 
#'   Physics} \strong{7}: 3910-3916.
#'   \href{https://rsc.li/2XkxPCm}{link}.
#'
#'   Falcioni, M. and M. W. Deem. 1999. A biased Monte Carlo scheme for 
#'   zeolite structure solution.  \emph{Journal of Chemical Physics}
#'   \strong{110}: 1754-1766.
#'   \href{https://aip.scitation.org/doi/10.1063/1.477812}{link}.
#' 
#'   Geyer, C. J. 1991. Markov Chain Monte Carlo maximum likelihood. \emph{In
#'   Computing Science and Statistics: Proceedings of the 23rd Symposium on 
#'   the Interface}. pp 156-163. American Statistical Association, New York,
#'   USA. \href{https://www.stat.umn.edu/geyer/f05/8931/c.pdf}{link}.
#'  
#'   Metropolis, N., A. W. Rosenbluth, M. N. Rosenbluth, A. H. Teller, and E.
#'   Teller. 1953. Equations of state calculations by fast computing machines.
#'   \emph{Journal of Chemical Physics} \strong{21}: 1087-1092.
#'   \href{https://bayes.wustl.edu/Manual/EquationOfState.pdf}{link}.
#'
#' @examples 
#' \donttest{
#'   data(rodents)
#'   document_term_table <- rodents$document_term_table
#'   document_covariate_table <- rodents$document_covariate_table
#'   LDA_models <- LDA_set(document_term_table, topics = 2)[[1]]
#'   data <- document_covariate_table
#'   data$gamma <- LDA_models@gamma
#'   weights <- document_weights(document_term_table)
#'   data <- data[order(data[,"newmoon"]), ]
#'   saves <- prep_saves(1, TS_control())
#'   inputs <- prep_ptMCMC_inputs(data, gamma ~ 1, 1, "newmoon", weights, 
#'                                TS_control())
#'   cpts <- prep_cpts(data, gamma ~ 1, 1, "newmoon", weights, TS_control())
#'   ids <- prep_ids(TS_control())
#'   for(i in 1:TS_control()$nit){
#'     steps <- step_chains(i, cpts, inputs)
#'     swaps <- swap_chains(steps, inputs, ids)
#'     saves <- update_saves(i, saves, steps, swaps)
#'     cpts <- update_cpts(cpts, swaps)
#'     ids <- update_ids(ids, swaps)
#'   }
#' }
#'
#' @export
#'
swap_chains <- function(chainsin, inputs, ids){
  temps <- inputs$temps
  itemps <- 1/temps
  ntemps <- length(temps)
  revtemps <- seq(ntemps - 1, 1)
  lls <- chainsin$lls
  changepts <- chainsin$changepts
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
#'   \cr \cr
#'   This set of functions was designed to work within \code{\link{TS}} and 
#'   specifically \code{\link{est_changepoints}}. They are still hardcoded to
#'   do so, but have the capacity to be generalized to work with any 
#'   estimation via ptMCMC with additional coding work.
#'
#' @details For each iteration of the ptMCMC algorithm, all of the chains
#'   have the potential to take a step. The possible step is proposed under
#'   a proposal distribution (here for change points we use a symmetric
#'   geometric distribution), the proposed step is then evaluated and either
#'   accepted or not (following the Metropolis-Hastings rule; Metropolis,
#'   \emph{et al.} 1953, Hasting 1960, Gupta \emph{et al.} 2018), and then
#'   accordingly taken or not (the configurations are updated). 
#'
#' @param i \code{integer} iteration index.
#'
#' @param cpts \code{matrix} of change point locations across chains.
#'
#' @param inputs Class \code{ptMCMC_inputs} \code{list}, containing the 
#'   static inputs for use within the ptMCMC algorithm.
#'
#' @param prop_step Proposed step output from \code{propose_step}.
#'
#' @param accept_step \code{logical} indicator of acceptance of each chain's
#'   proposed step.
#'
#' @return 
#'   \code{step_chains}: \code{list} of change points, log-likelihoods, 
#'   and logical indicators of acceptance for each chain. \cr \cr
#'   \code{propose_step}: \code{list} of change points and 
#'   log-likelihood values for the proposal. \cr \cr
#'   \code{eval_step}: \code{logical} vector indicating if each 
#'   chain's proposal was accepted. \cr \cr
#'   \code{take_step}: \code{list} of change points, log-likelihoods, 
#'   and logical indicators of acceptance for each chain.
#'
#' @references
#'   Gupta, S., L. Hainsworth, J. S. Hogg, R. E. C. Lee, and J. R. Faeder. 
#'   2018. Evaluation of parallel tempering to accelerate Bayesian parameter
#'   estimation in systems biology. 
#'   \href{https://arxiv.org/abs/1801.09831}{link}.
#'
#'   Hastings, W. K. 1970. Monte Carlo sampling methods using Markov Chains
#'   and their applications. \emph{Biometrika} \strong{57}:97-109.
#'   \href{https://doi.org/10.2307/2334940}{link}.
#'  
#'   Metropolis, N., A. W. Rosenbluth, M. N. Rosenbluth, A. H. Teller, and E.
#'   Teller. 1953. Equations of state calculations by fast computing machines.
#'   \emph{Journal of Chemical Physics} \strong{21}: 1087-1092.
#'   \href{https://bayes.wustl.edu/Manual/EquationOfState.pdf}{link}.
#'
#' @examples 
#' \donttest{
#'   data(rodents)
#'   document_term_table <- rodents$document_term_table
#'   document_covariate_table <- rodents$document_covariate_table
#'   LDA_models <- LDA_set(document_term_table, topics = 2)[[1]]
#'   data <- document_covariate_table
#'   data$gamma <- LDA_models@gamma
#'   weights <- document_weights(document_term_table)
#'   data <- data[order(data[,"newmoon"]), ]
#'   saves <- prep_saves(1, TS_control())
#'   inputs <- prep_ptMCMC_inputs(data, gamma ~ 1, 1, "newmoon", weights, 
#'                                TS_control())
#'   cpts <- prep_cpts(data, gamma ~ 1, 1, "newmoon", weights, TS_control())
#'   ids <- prep_ids(TS_control())
#'   for(i in 1:TS_control()$nit){
#'     steps <- step_chains(i, cpts, inputs)
#'     swaps <- swap_chains(steps, inputs, ids)
#'     saves <- update_saves(i, saves, steps, swaps)
#'     cpts <- update_cpts(cpts, swaps)
#'     ids <- update_ids(ids, swaps)
#'   }
#'   # within step_chains()
#'   cpts <- prep_cpts(data, gamma ~ 1, 1, "newmoon", weights, TS_control())
#'   i <- 1
#'   prop_step <- propose_step(i, cpts, inputs)
#'   accept_step <- eval_step(i, cpts, prop_step, inputs)
#'   take_step(cpts, prop_step, accept_step)
#' }
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
  lls <- vapply(mods, logLik, 0)
  list(changepts = prop_changepts, lls = lls)
}

#' @rdname step_chains
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
#'   proposed change points within the ptMCMC algorithm
#'
#' @description This function wraps around \code{TS_memo} 
#'   (optionally memoised \code{\link{multinom_TS}}) to provide a
#'   simpler interface within the ptMCMC algorithm and is implemented within
#'   \code{\link{propose_step}}. 
#'
#' @param prop_changepts \code{matrix} of proposed change points across 
#'   chains.
#'
#' @param inputs Class \code{ptMCMC_inputs} list, containing the static inputs
#'   for use within the ptMCMC algorithm.
#'
#' @return List of models associated with the proposed step, with an element
#'   for each chain.
#'
#' @examples 
#' \donttest{
#'   data(rodents)
#'   document_term_table <- rodents$document_term_table
#'   document_covariate_table <- rodents$document_covariate_table
#'   LDA_models <- LDA_set(document_term_table, topics = 2)[[1]]
#'   data <- document_covariate_table
#'   data$gamma <- LDA_models@gamma
#'   weights <- document_weights(document_term_table)
#'   data <- data[order(data[,"newmoon"]), ]
#'   saves <- prep_saves(1, TS_control())
#'   inputs <- prep_ptMCMC_inputs(data, gamma ~ 1, 1, "newmoon", weights,
#'                                TS_control())
#'   cpts <- prep_cpts(data, gamma ~ 1, 1, "newmoon", weights, TS_control())
#'   i <- 1
#'   pdist <- inputs$pdist
#'   ntemps <- length(inputs$temps)
#'   selection <- cbind(pdist$which_steps[i, ], 1:ntemps)
#'   prop_changepts <- cpts$changepts
#'   curr_changepts_s <- cpts$changepts[selection]
#'   prop_changepts_s <- curr_changepts_s + pdist$steps[i, ]
#'   if(all(is.na(prop_changepts_s))){
#'     prop_changepts_s <- NULL
#'   }
#'   prop_changepts[selection] <- prop_changepts_s
#'   mods <- proposed_step_mods(prop_changepts, inputs)
#' }
#'
#' @export
#'
proposed_step_mods <- function(prop_changepts, inputs){ 

  data <- inputs$data
  formula <- inputs$formula
  weights <- inputs$weights
  TS_function <- inputs$TS_function
  ntemps <- length(inputs$temps)
  control <- inputs$control
  timename <- inputs$timename
  out <- vector("list", length = ntemps)
  for (i in 1:ntemps){
    out[[i]] <- TS_function(data, formula, prop_changepts[ , i], timename, 
                            weights, control)
  }
  out
}


#' @title Initialize and update the chain ids throughout the ptMCMC algorithm
#'
#' @description \code{prep_ids} creates and \code{update_ids} updates 
#'   the active vector of identities (ids) for each of the chains in the 
#'   ptMCMC algorithm. These ids are used to track trips of the particles
#'   among chains.
#'   \cr \cr
#'   These functions were designed to work within \code{\link{TS}} and 
#'   specifically \code{\link{est_changepoints}}, but have been generalized
#'   and would work within any general ptMCMC as long as \code{control},
#'   \code{ids}, and \code{swaps} are formatted properly.
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   Time Series model including the parallel tempering Markov Chain 
#'   Monte Carlo (ptMCMC) controls. Values not input assume defaults set by 
#'   \code{\link{TS_control}}.
#'
#' @param ids The existing vector of chain ids.
#'
#' @param swaps Chain configuration after among-temperature swaps.
#'
#' @return The vector of chain ids.
#'
#' @examples
#'   prep_ids()
#' \donttest{
#'   data(rodents)
#'   document_term_table <- rodents$document_term_table
#'   document_covariate_table <- rodents$document_covariate_table
#'   LDA_models <- LDA_set(document_term_table, topics = 2)[[1]]
#'   data <- document_covariate_table
#'   data$gamma <- LDA_models@gamma
#'   weights <- document_weights(document_term_table)
#'   data <- data[order(data[,"newmoon"]), ]
#'   saves <- prep_saves(1, TS_control())
#'   inputs <- prep_ptMCMC_inputs(data, gamma ~ 1, 1, "newmoon", weights,
#'                                TS_control())
#'   cpts <- prep_cpts(data, gamma ~ 1, 1, "newmoon", weights, TS_control())
#'   ids <- prep_ids(TS_control())
#'   for(i in 1:TS_control()$nit){
#'     steps <- step_chains(i, cpts, inputs)
#'     swaps <- swap_chains(steps, inputs, ids)
#'     saves <- update_saves(i, saves, steps, swaps)
#'     cpts <- update_cpts(cpts, swaps)
#'     ids <- update_ids(ids, swaps)
#'   }
#' }
#'
#' @export
#'
prep_ids <- function(control = list()){
  control <- do.call("TS_control", control)
  if (!is.numeric(control$ntemps) || any(control$ntemps %% 1 != 0)){
    stop("ntemps must be integer-valued")
  }
  1:control$ntemps
}

#' @rdname prep_ids
#'
#' @export
#'
update_ids <- function(ids, swaps){
  swaps$ids
}

#' @title Prepare the inputs for the ptMCMC algorithm estimation of 
#'   change points
#'
#' @description Package the static inputs (controls and data structures) used
#'   by the ptMCMC algorithm in the context of estimating change points. 
#'   \cr \cr
#'   This function was designed to work within \code{\link{TS}} and 
#'   specifically \code{\link{est_changepoints}}. It is still hardcoded to do
#'   so, but has the capacity to be generalized to work with any estimation
#'   via ptMCMC with additional coding work.
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
#'   number of change points is used to dictate the segmentation of the data  
#'   for each continuous model and each LDA model.
#'
#' @param timename \code{character} element indicating the time variable
#'   used in the time series. Defaults to \code{"time"}. The variable must be
#'   integer-conformable or a \code{Date}. If the variable named
#'   is a \code{Date}, the input is converted to an integer, resulting in the
#'   timestep being 1 day, which is often not desired behavior.
#'
#' @param weights Optional class \code{numeric} vector of weights for each 
#'   document. Defaults to \code{NULL}, translating to an equal weight for
#'   each document. When using \code{multinom_TS} in a standard LDATS 
#'   analysis, it is advisable to weight the documents by their total size,
#'   as the result of \code{\link[topicmodels]{LDA}} is a matrix of 
#'   proportions, which does not account for size differences among documents.
#'   For most models, a scaling of the weights (so that the average is 1) is
#'   most appropriate, and this is accomplished using 
#'   \code{\link{document_weights}}.
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   Time Series model including the parallel tempering Markov Chain 
#'   Monte Carlo (ptMCMC) controls. Values not input assume defaults set by 
#'   \code{\link{TS_control}}.
#'
#' @return Class \code{ptMCMC_inputs} \code{list}, containing the static 
#'   inputs for use within the ptMCMC algorithm for estimating change points. 
#'
#' @examples 
#' \donttest{
#'   data(rodents)
#'   document_term_table <- rodents$document_term_table
#'   document_covariate_table <- rodents$document_covariate_table
#'   LDA_models <- LDA_set(document_term_table, topics = 2)[[1]]
#'   data <- document_covariate_table
#'   data$gamma <- LDA_models@gamma
#'   weights <- document_weights(document_term_table)
#'   data <- data[order(data[,"newmoon"]), ]
#'   saves <- prep_saves(1, TS_control())
#'   inputs <- prep_ptMCMC_inputs(data, gamma ~ 1, 1, "newmoon", weights, 
#'                                TS_control())
#' }
#' @export
#'
prep_ptMCMC_inputsx <- function(data, formula, nchangepoints, timename, 
                               weights = NULL, control = list()){
  check_timename(data, timename)
  check_formula(data, formula)
  check_weights(weights)
  check_nchangepoints(nchangepoints)
  check_control(control)
  control <- do.call("TS_control", control)
  control$selector <- NULL
  control$measurer <- NULL
  out <- list(control = control, temps = prep_temp_sequence(control), 
              pdist = prep_proposal_dist(nchangepoints, control),
              formula = formula, weights = weights, data = data, 
              TS_function = TS_fun(control),
              timename = timename)
  class(out) <- c("ptMCMC_inputs", "list")
  out
}


#' @title Pre-calculate the change point proposal distribution for the ptMCMC 
#'   algorithm
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
#'   number of change points is used to dictate the segmentation of the data  
#'   for each continuous model and each LDA model.
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   Time Series model including the parallel tempering Markov Chain 
#'   Monte Carlo (ptMCMC) controls. Values not input assume defaults set by 
#'   \code{\link{TS_control}}. Currently relevant here is 
#'   \code{magnitude}, which controls the magnitude of the step size (is the 
#'   average of the geometric distribution). 
#'
#' @return \code{list} of two \code{matrix} elements: [1] the size of the 
#'   proposed step for each iteration of each chain and [2] the identity of 
#'   the change point location to be shifted by the step for each iteration of
#'   each chain.
#'
#' @examples
#'   prep_proposal_dist(nchangepoints = 2)
#'
#' @export
#'
prep_proposal_distx <- function(nchangepoints, control = list()){
  check_nchangepoints(nchangepoints)
  check_control(control)
  control <- do.call("TS_control", control)
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
#'   are then processed (burn-in iterations are dropped and the remaining
#'   iterations are thinned) via \code{process_saves}.
#'   \cr \cr
#'   This set of functions was designed to work within \code{\link{TS}} and 
#'   specifically \code{\link{est_changepoints}}. They are still hardcoded to
#'   do so, but have the capacity to be generalized to work with any
#'   estimation via ptMCMC with additional coding work.
#'
#' @param nchangepoints \code{integer} corresponding to the number of 
#'   change points to include in the model. 0 is a valid input (corresponding
#'   to no change points, so a singular time series model), and the current 
#'   implementation can reasonably include up to 6 change points. The 
#'   number of change points is used to dictate the segmentation of the data  
#'   for each continuous model and each LDA model.
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   Time Series model including the parallel tempering Markov Chain 
#'   Monte Carlo (ptMCMC) controls. Values not input assume defaults set by 
#'   \code{\link{TS_control}}.
#'
#' @param i \code{integer} iteration index. 
#'
#' @param saves The existing list of saved data objects.
#'
#' @param steps Chain configuration after within-temperature steps.
#'
#' @param swaps Chain configuration after among-temperature swaps.
#'
#' @return \code{list} of ptMCMC objects: change points (\code{$cpts}), 
#'   log-likelihoods (\code{$lls}), chain ids (\code{$ids}), step acceptances
#'   (\code{$step_accepts}), and swap acceptances (\code{$swap_accepts}).
#'
#' @examples 
#' \donttest{
#'   data(rodents)
#'   document_term_table <- rodents$document_term_table
#'   document_covariate_table <- rodents$document_covariate_table
#'   LDA_models <- LDA_set(document_term_table, topics = 2)[[1]]
#'   data <- document_covariate_table
#'   data$gamma <- LDA_models@gamma
#'   weights <- document_weights(document_term_table)
#'   data <- data[order(data[,"newmoon"]), ]
#'   saves <- prep_saves(1, TS_control())
#'   inputs <- prep_ptMCMC_inputs(data, gamma ~ 1, 1, "newmoon", weights, 
#'                                TS_control())
#'   cpts <- prep_cpts(data, gamma ~ 1, 1, "newmoon", weights, TS_control())
#'   ids <- prep_ids(TS_control())
#'   for(i in 1:TS_control()$nit){
#'     steps <- step_chains(i, cpts, inputs)
#'     swaps <- swap_chains(steps, inputs, ids)
#'     saves <- update_saves(i, saves, steps, swaps)
#'     cpts <- update_cpts(cpts, swaps)
#'     ids <- update_ids(ids, swaps)
#'   }
#'   process_saves(saves, TS_control())
#' }
#'
#' @export
#'
prep_savesx <- function(nchangepoints, control = list()){
  check_nchangepoints(nchangepoints)
  check_control(control)
  control <- do.call("TS_control", control)
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
process_saves <- function(saves, control = list()){
  control <- do.call("TS_control", control)
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

#' @title Initialize and update the change point matrix used in the ptMCMC
#'   algorithm
#' 
#' @description Each of the chains is initialized by \code{prep_cpts} using a 
#'   draw from the available times (i.e. assuming a uniform prior), the best 
#'   fit (by likelihood) draw is put in the focal chain with each subsequently 
#'   worse fit placed into the subsequently hotter chain. \code{update_cpts}
#'   updates the change points after every iteration in the ptMCMC algorithm.
#'
#' @param data \code{data.frame} including [1] the time variable (indicated 
#'   in \code{timename}), [2] the predictor variables (required by
#'   \code{formula}) and [3], the multinomial response variable (indicated in
#'   \code{formula}) as verified by \code{\link{check_timename}} and 
#'   \code{\link{check_formula}}. Note that the response variables should be
#'   formatted as a \code{data.frame} object named as indicated by the 
#'   \code{response} entry in the \code{control} list, such as \code{gamma} 
#'   for a standard TS analysis on LDA output.  
#'
#' @param formula \code{formula} defining the regression relationship between
#'   the change points, see \code{\link[stats]{formula}}. Any 
#'   predictor variable included must also be a column in 
#'   \code{data} and any (multinomial) response variable must be a set of
#'   columns in \code{data}, as verified by \code{\link{check_formula}}.
#'
#' @param nchangepoints \code{integer} corresponding to the number of 
#'   change points to include in the model. 0 is a valid input (corresponding
#'   to no change points, so a singular time series model), and the current 
#'   implementation can reasonably include up to 6 change points. The 
#'   number of change points is used to dictate the segmentation of the data  
#'   for each continuous model and each LDA model.
#'
#' @param timename \code{character} element indicating the time variable
#'   used in the time series. Defaults to \code{"time"}. The variable must be
#'   integer-conformable or a \code{Date}. If the variable named
#'   is a \code{Date}, the input is converted to an integer, resulting in the
#'   timestep being 1 day, which is often not desired behavior.
#'
#' @param weights Optional class \code{numeric} vector of weights for each 
#'   document. Defaults to \code{NULL}, translating to an equal weight for
#'   each document. When using \code{multinom_TS} in a standard LDATS 
#'   analysis, it is advisable to weight the documents by their total size,
#'   as the result of \code{\link[topicmodels]{LDA}} is a matrix of 
#'   proportions, which does not account for size differences among documents.
#'   For most models, a scaling of the weights (so that the average is 1) is
#'   most appropriate, and this is accomplished using 
#'   \code{\link{document_weights}}.
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   Time Series model including the parallel tempering Markov Chain 
#'   Monte Carlo (ptMCMC) controls. Values not input assume defaults set by 
#'   \code{\link{TS_control}}.
#'
#' @param cpts The existing matrix of change points.
#'
#' @param swaps Chain configuration after among-temperature swaps.
#'
#' @return \code{list} of [1] \code{matrix} of change points (rows) for 
#'   each temperature (columns) and [2] \code{vector} of log-likelihood 
#'   values for each of the chains.
#'
#' @examples 
#' \donttest{
#'   data(rodents)
#'   document_term_table <- rodents$document_term_table
#'   document_covariate_table <- rodents$document_covariate_table
#'   LDA_models <- LDA_set(document_term_table, topics = 2)[[1]]
#'   data <- document_covariate_table
#'   data$gamma <- LDA_models@gamma
#'   weights <- document_weights(document_term_table)
#'   data <- data[order(data[,"newmoon"]), ]
#'   saves <- prep_saves(1, TS_control())
#'   inputs <- prep_ptMCMC_inputs(data, gamma ~ 1, 1, "newmoon", weights,
#'                                TS_control())
#'   cpts <- prep_cpts(data, gamma ~ 1, 1, "newmoon", weights, TS_control())
#'   ids <- prep_ids(TS_control())
#'   for(i in 1:TS_control()$nit){
#'     steps <- step_chains(i, cpts, inputs)
#'     swaps <- swap_chains(steps, inputs, ids)
#'     saves <- update_saves(i, saves, steps, swaps)
#'     cpts <- update_cpts(cpts, swaps)
#'     ids <- update_ids(ids, swaps)
#'   }
#' }
#'
#' @export
#'
prep_cptsx <- function(data, formula, nchangepoints, timename, weights, 
                            control = list()){

  check_formula(data, formula)
  check_nchangepoints(nchangepoints)
  check_weights(weights)
  check_timename(data, timename)
  check_control(control)
  control <- do.call("TS_control", control)
  temps <- prep_temp_sequence(control)
  ntemps <- length(temps)
  min_time <- min(data[ , timename])
  max_time <- max(data[ , timename])
  times <- seq(min_time, max_time, 1)
  avail_times <- times[-c(1, length(times))]
  cps <- matrix(NA, nrow = nchangepoints, ncol = ntemps)
  for (i in 1:ntemps){
    cp_times <- sort(sample(avail_times, nchangepoints, replace = FALSE))
    cps[ , i] <- cp_times
  }
  lls <- rep(NA, ntemps)
  for (i in 1:ntemps){
    modfit <- multinom_TS(data, formula, cps[ ,i], timename, weights, control)
    lls[i] <- modfit$logLik
  }  
  cps <- cps[ , order(lls, decreasing = TRUE), drop = FALSE]
  lls <- sort(lls, decreasing = TRUE)

  out <- list(cps, lls)
  names(out) <- c("changepts", "lls")
  out
}

#' @rdname prep_cpts
#'
#' @export
#'
update_cpts <- function(cpts, swaps){
  list(changepts = swaps$changepts, lls = swaps$lls)
}

#' @title Prepare the ptMCMC temperature sequence
#'
#' @description Create the series of temperatures used in the ptMCMC 
#'   algorithm.
#'   \cr \cr
#'   This function was designed to work within \code{\link{TS}} and 
#'   \code{\link{est_changepoints}} specifically, but has been generalized
#'   and would work with any ptMCMC model as long as \code{control}
#'   includes the relevant control parameters (and provided that the 
#'   \code{\link{check_control}} function and its use here are generalized).
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   Time Series model including the parallel tempering Markov Chain 
#'   Monte Carlo (ptMCMC) controls. Values not input assume defaults set by 
#'   \code{\link{TS_control}}.
#'
#' @return \code{vector} of temperatures.
#'
#' @examples
#'   prep_temp_sequence()
#'
#' @export
#'
prep_temp_sequencex <- function(control = list()){
  check_control(control)
  control <- do.call("TS_control", control)
  ntemps <- control$ntemps
  penultimate_temp <- control$penultimate_temp
  ultimate_temp <- control$ultimate_temp
  q <- control$q
  sequence <- seq(0, log2(penultimate_temp), length.out = ntemps - 1)
  log_temps <- sequence^(1 + q) / log2(penultimate_temp)^q
  c(2^(log_temps), ultimate_temp) 
}
