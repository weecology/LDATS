
ldats_classic <- function(TS, control = list()){

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


process_saves <- function(saves, control = list()){
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


update_saves <- function(i, saves, steps, swaps){
  saves$cpts[ , , i] <- swaps$changepts
  saves$lls[ , i] <- swaps$lls
  saves$ids[ , i] <- swaps$ids
  saves$step_accepts[ , i] <- steps$accept_step
  saves$swap_accepts[ , i] <- swaps$accept_swap
  saves
}
update_cpts <- function(cpts, swaps){
  list(changepts = swaps$changepts, lls = swaps$lls)
}

update_ids <- function(ids, swaps){
  swaps$ids
}

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



prep_cpts <- function(TS, control = list()){

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


prep_temp_sequence <- function(TS, control = list()){
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
              temps = prep_temp_sequence(TS = TS, control = control), 
              pdist = prep_proposal_dist(TS = TS, control = control),
              formula = TS$formula, 
              weights = TS$weights, 
              data = TS$data$train$ts_data, 
              response = TS$response,
              timename = TS$timename)
  class(out) <- c("ptMCMC_inputs", "list")
  out
}



prep_proposal_dist <- function(TS, control = list()){
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

prep_ids <- function(TS, control = list()){
  if (!is.numeric(control$ntemps) || any(control$ntemps %% 1 != 0)){
    stop("ntemps must be integer-valued")
  }
  1:control$ntemps
}


step_chains <- function(TS, i, cpts, inputs){
  prop_step <- propose_step(TS = TS, i = i, cpts = cpts, inputs = inputs)
  accept_step <- eval_step(i = i, cpts = cpts, prop_step = prop_step, 
                           inputs = inputs)
  take_step(cpts = cpts, prop_step = prop_step, accept_step = accept_step)
}


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

    fun <- eval(parse(text = paste0(TS$response, "_TS")))
    fun <- memoise_fun(fun, control$memoise)
    args <- list(data = data, formula = TS$formula, 
                 changepoints = prop_changepts[ , i], 
                 timename = TS$timename, weights = TS$weights, 
                 control = control)
    out[[i]] <- soft_call(fun, args, TRUE)
  }
  out
}



eval_step <- function(i, cpts, prop_step, inputs){
  temps <- inputs$temps
  ntemps <- length(temps)
  itemps <- 1 / temps
  runif(ntemps) < exp((prop_step$lls - cpts$lls) * itemps)
}


take_step <- function(cpts, prop_step, accept_step){
  changepts <- cpts$changepts
  lls <- cpts$lls
  changepts[ , accept_step] <- prop_step$changepts[ , accept_step]
  lls[accept_step] <- prop_step$lls[accept_step]
  list(changepts = changepts, lls = lls, accept_step = accept_step)
}



