#' @title Estimate changepoints using the LDATS classic ptMCMC method
#'
#' @description Uses the LDATS classic ptMCMC method to fit a changepoint
#'  model. 
#'
#' @param TS \code{list} time series model object. 
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   time series model via the LDATS classic ptMCMC method. Values not input 
#'   assume defaults set by \code{\link{ldats_classic_control}}.
#'
#' @return \code{list} of changepoint locations, log likelihoods, and model
#'  diagnostics.
#'
#' @export
#'
ldats_classic <- function(TS, control = list()){
  control <- do.call("ldats_classic_control", control)
  saves <- prep_saves(TS = TS, control = control)

  inputs <- prep_ptMCMC_inputs(TS = TS, control = control)
  cpts <- prep_cpts(TS = TS, control = control)
  ids <- prep_ids(TS = TS, control)
  pbar <- prep_pbar(control = control, bar_type = "rho")

  for(i in 1:control$nit){
    update_pbar(pbar = pbar, control = control)
    steps <- step_chains(TS = TS, i = i, cpts = cpts, inputs = inputs)
    swaps <- swap_chains(chainsin = steps, inputs = inputs, ids = ids)
    saves <- update_saves(i = i, saves = saves, steps = steps, swaps = swaps)
    cpts <- update_cpts(cpts = cpts, swaps = swaps)
    ids <- update_ids(ids = ids, swaps = swaps)
  }

  process_saves(saves = saves, control = control)

}


#' @title Create the controls list for the classic LDATS ptMCMC sampler
#'
#' @description This function provides a simple creation and definition of a
#'   list used to control time series model fitting via the
#'   \code{\link{ldats_classic}} sampler. 
#'
#' @param ntemps \code{integer} number of temperatures (chains) to use in the 
#'   ptMCMC algorithm.
#'
#' @param penultimate_temp Penultimate temperature in the ptMCMC sequence.
#'
#' @param ultimate_temp Ultimate temperature in the ptMCMC sequence.
#'
#' @param q Exponent controlling the ptMCMC temperature sequence from the 
#'   focal chain (reference with temperature = 1) to the penultimate chain. 0
#'   (default) implies a geometric sequence. 1 implies squaring before 
#'   exponentiating.
#'
#' @param nit \code{integer} number of iterations (steps) used in the ptMCMC
#'   algorithm.
#'
#' @param magnitude Average magnitude (defining a geometric distribution)
#'   for the proposed step size in the ptMCMC algorithm.
#'
#' @param burnin \code{integer} number of iterations to remove from the 
#'   beginning of the ptMCMC algorithm.
#'
#' @param thin_frac Fraction of iterations to retain, must be \eqn{(0, 1]},
#'   and the default value of 1 represents no thinning.
#'
#'
#' @param quiet \code{logical} indicator of whether the model should run 
#'   quietly (if \code{FALSE}, a progress bar and notifications are printed).
#'
#' @param memoise \code{logical} indicator of whether the response 
#'   function should be memoised (via \code{\link[memoise]{memoise}}). 
#'
#' @return \code{list}, with named elements corresponding to the arguments.
#'
#' @export
#'
ldats_classic_control <- function(ntemps = 6, penultimate_temp = 2^6, 
                                  ultimate_temp = 1e10, q = 0, 
                                  nit = 1e4, magnitude = 12, 
                                  burnin = 0, thin_frac = 1, 
                                  memoise = TRUE, quiet = FALSE){
  list(ntemps = ntemps, penultimate_temp = penultimate_temp, 
       ultimate_temp = ultimate_temp, q = q, nit = nit, 
       magnitude = magnitude, burnin = burnin, thin_frac = thin_frac, 
       memoise = memoise, quiet = quiet)
}


#' @title Count trips of the ptMCMC particles in a classic_ldats estimation
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


#' @title Initialize and update the change point matrix used in the LDATS
#'   classic ptMCMC algorithm
#' 
#' @description Each of the chains is initialized by \code{prep_cpts} using a 
#'   draw from the available times (i.e. assuming a uniform prior), the best 
#'   fit (by likelihood) draw is put in the focal chain with each subsequently 
#'   worse fit placed into the subsequently hotter chain. \code{update_cpts}
#'   updates the change points after every iteration in the ptMCMC algorithm.
#'
#' @param TS \code{list} time series model object. 
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   time series model via the LDATS classic ptMCMC method. Values not input 
#'   assume defaults set by \code{\link{ldats_classic_control}}.
#'
#' @param cpts The existing matrix of change points.
#'
#' @param swaps Chain configuration after among-temperature swaps.
#'
#' @return \code{list} of [1] \code{matrix} of change points (rows) for 
#'   each temperature (columns) and [2] \code{vector} of log-likelihood 
#'   values for each of the chains.
#'
#'
#' @export
#'
prep_cpts <- function(TS, control = list()){
  control <- do.call("ldats_classic_control", control)
  data <- TS$data$train$ts_data
  temps <- prep_temp_sequence(TS = TS, control = control)
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
    fun <- TS$control$response
    fun <- memoise_fun(fun, TS$control$memoise)
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


#' @rdname prep_cpts
#'
#' @export
#'
update_cpts <- function(cpts, swaps){
  list(changepts = swaps$changepts, lls = swaps$lls)
}

#' @title Prepare the ptMCMC temperature sequence for the LDATS classic 
#'  method
#'
#' @description Create the series of temperatures used in the ptMCMC 
#'   algorithm.
#'   \cr \cr
#'   This function was designed to work within \code{\link{TS}} and 
#'   \code{\link{est_changepoints}} specifically, but has been generalized
#'   and would work with any ptMCMC model as long as \code{control}
#'   includes the relevant control parameters.
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   time series model via the LDATS classic ptMCMC method. Values not input 
#'   assume defaults set by \code{\link{ldats_classic_control}}.
#'
#' @param TS \code{list} time series model object. 
#'
#' @return \code{vector} of temperatures.
#'
#' @examples
#'   prep_temp_sequence()
#'
#' @export
#'
prep_temp_sequence <- function(TS, control = list()){
  control <- do.call("ldats_classic_control", control)
  ntemps <- control$ntemps
  penultimate_temp <- control$penultimate_temp
  ultimate_temp <- control$ultimate_temp
  q <- control$q
  sequence <- seq(0, log2(penultimate_temp), length.out = ntemps - 1)
  log_temps <- sequence^(1 + q) / log2(penultimate_temp)^q
  c(2^(log_temps), ultimate_temp) 
}


#' @title Prepare and update the data structures to save the LDATS classic 
#'   ptMCMC output
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
#' @param TS \code{list} time series model object. 
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   time series model via the LDATS classic ptMCMC method. Values not input 
#'   assume defaults set by \code{\link{ldats_classic_control}}.
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
#'
#' @export
#'
prep_saves <- function(TS, control = list()){
  control <- do.call("ldats_classic_control", control)
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
  control <- do.call("ldats_classic_control", control)
  nit <- control$nit
  iters <- 1:nit
  if (control$burnin > 0){
    iters <- iters[-(1:control$burnin)]
  }
  niters <- length(iters)
  thin_interval <- ceiling(1/control$thin_frac)
  iters_thinned <- seq(1, niters, by = thin_interval)
  dims <- c(dim(saves$cpts)[1:2], length(iters_thinned))

  trips <- count_trips(saves$ids)
  diagnostics <- list(step_acceptance_rate = rowMeans(saves$step_accepts), 
                      swap_acceptance_rate = rowMeans(saves$swap_accepts), 
                      trip_counts = trips$trip_counts, 
                      trip_rates = trips$trip_rates)


  saves$cpts <- array(saves$cpts[ , , iters_thinned], dim = dims)
  saves$lls <- saves$lls[, iters_thinned]
  saves$ids <- saves$ids[, iters_thinned]
  saves$step_accepts <- saves$step_accepts[ , iters_thinned]
  saves$swap_accepts <- saves$swap_accepts[ , iters_thinned]

  saves$diagnostics <- diagnostics


  saves
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
#' @param TS \code{list} time series model object. 
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   time series model via the LDATS classic ptMCMC method. Values not input 
#'   assume defaults set by \code{\link{ldats_classic_control}}.
#'
#' @return \code{list} containing the static 
#'   inputs for use within the ptMCMC algorithm for estimating change points. 
#'
#' @export
#'
prep_ptMCMC_inputs <- function(TS, control = list()){
  control <- do.call("ldats_classic_control", control)
  fun <- TS$control$response
  fun <- memoise_fun(fun, control$memoise)
  list(control = control, 
              temps = prep_temp_sequence(TS = TS, control = control), 
              pdist = prep_proposal_dist(TS = TS, control = control),
              formula = TS$formula, 
              weights = TS$weights, 
              data = TS$data$train$ts_data, 
              response = TS$response,
              timename = TS$timename,
              fun = fun)
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
#' @param TS \code{list} time series model object. 
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   time series model via the LDATS classic ptMCMC method. Values not input 
#'   assume defaults set by \code{\link{ldats_classic_control}}.
#'
#' @return \code{list} of two \code{matrix} elements: [1] the size of the 
#'   proposed step for each iteration of each chain and [2] the identity of 
#'   the change point location to be shifted by the step for each iteration of
#'   each chain.
#'
#' @export
#'
prep_proposal_dist <- function(TS, control = list()){
  control <- do.call("ldats_classic_control", control)
  nchangepoints <- TS$nchangepoints
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


#' @title Initialize and update the chain ids throughout the LDATS classic
#'   ptMCMC algorithm
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
#' @param TS \code{list} time series model object. 
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   time series model via the LDATS classic ptMCMC method. Values not input 
#'   assume defaults set by \code{\link{ldats_classic_control}}.
#'
#' @param ids The existing vector of chain ids.
#'
#' @param swaps Chain configuration after among-temperature swaps.
#'
#' @return The vector of chain ids.
#'
#' @export
#'
prep_ids <- function(TS, control = list()){
  control <- do.call("ldats_classic_control", control)
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
#' @export
#'
step_chains <- function(TS, i, cpts, inputs){
  prop_step <- propose_step(TS = TS, i = i, cpts = cpts, inputs = inputs)
  accept_step <- eval_step(i = i, cpts = cpts, prop_step = prop_step, 
                           inputs = inputs)
  take_step(cpts = cpts, prop_step = prop_step, accept_step = accept_step)
}

#' @rdname step_chains
#'
#' @export
#'
propose_step <- function(TS, i, cpts, inputs){

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
  mods <- proposed_step_mods(TS = TS, prop_changepts = prop_changepts, 
                             inputs = inputs)
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
#' @export
#'
proposed_step_mods <- function(TS, prop_changepts, inputs){ 

  data <- inputs$data
  formula <- inputs$formula
  weights <- inputs$weights
  TS_function <- inputs$TS_function
  ntemps <- length(inputs$temps)
  control <- inputs$control
  timename <- inputs$timename
  out <- vector("list", length = ntemps)
  for (i in 1:ntemps){

    fun <- inputs$fun
    args <- list(data = data, formula = TS$formula, 
                 changepoints = prop_changepts[ , i], 
                 timename = TS$timename, weights = TS$weights, 
                 control = control)
    out[[i]] <- soft_call(fun, args, TRUE)
  }
  out
}



