multinom_chunk <- function(formula_RHS, data, start_time, end_time, weights, 
                           ...) {
  formula <- as.formula(paste("gamma ~", formula_RHS))
  mod <- nnet::multinom(formula, data, weights, 
           subset = data$time > start_time & data$time <= end_time, 
           trace = FALSE, ...) 
  return(mod)
}


formula_RHS<-"1"

chunk_memo <- memoise::memoise(multinom_chunk)

system.time(for(W in 1:1e3){
chunk_memo(formula_RHS, data, start_time = start_times[i], 
  end_time = end_times[i], weights)
})

system.time(for(W in 1:1e3){
multinom_chunk(formula_RHS, data, start_time = start_times[i], 
  end_time = end_times[i], weights)
})




multinom_ts <- function(formula_RHS, data, changepoints = NULL, weights, 
                        ...){

  chunk_memo <- memoise::memoise(multinom_chunk)

  max_time <- max(data$time)
  unsorted <- is.unsorted(changepoints, strictly = TRUE)
  if (any(changepoints <= 0) | any(changepoints >= max_time) | unsorted){
    return(-Inf)
  }

  nchunks <- length(changepoints) + 1
  start_times <- c(min(data$time) - 1, changepoints)
  end_times <- c(changepoints, max(data$time) + 1)

  mods <- vector("list", length = nchunks)
  ll <- rep(NA, nchunks)
  for (i in 1:nchunks){
    mods[[i]] <- chunk_memo(formula_RHS, data, start_time = start_times[i], 
                   end_time = end_times[i], weights)
    ll[i] <- logLik(mods[[i]])
  }
  total_ll <- sum(ll)
  names(mods) <- sprintf("%s %d %s", "chunk", 1:nchunks, "model")
  output <- list("chunk models" = mods, "logLik" = total_ll)

  return(output)
}



ts_memo <- memoise::memoise(multinom_ts)


system.time(for(W in 1:1e3){
changepoints <- sample(4, 2:300, replace = F)
x1<-multinom_ts(formula_RHS, data, changepoints, weights)
})

system.time(for(W in 1:1e3){
changepoints <- sample(4, 2:300, replace = F)
x2<-ts_memo(formula_RHS, data, changepoints, weights)
})



#
# the massive difference in time is lost at this step in the hierachy
#


MTS <- function(formula, data, nchangepoints, weights, nit = 1e4, 
                ...){
  ts_memo <- memoise::memoise(multinom_ts)

  temps <- prep_temps(...)
  betas <- 1 / temps
  ntemps <- length(betas)

  prep_cpts <- prep_changepts(formula, data, ntemps, nchangepoints, 
                 weights)
  changepts <- prep_cpts$changepts
  lls <- prep_cpts$lls

  saved_cps <- array(NA, c(nchangepoints, ntemps, nit))
  saved_lls <- matrix(NA, ntemps, nit)
  saved_ids <- matrix(NA, ntemps, nit)
  accept_rate <- 0
  temp_ids <- 1:ntemps
  swap_accepted <- matrix(FALSE, nit, ntemps - 1)

  pdist <- proposal_dist(nit, ntemps, nchangepoints, magnitude = 12)
 
  pbform <- "  [:bar] :percent eta: :eta"
  pb <- progress::progress_bar$new(pbform, nit, clear = FALSE, width = 60)

  for (i in 1:nit){
  
    pb$tick()
 
    selection <- cbind(pdist$which_kicked[i, ], 1:ntemps)
    prop_changepts <- changepts
    curr_changepts_s <- changepts[selection]
    prop_changepts_s <- curr_changepts_s + pdist$kicks[i, ]
    prop_changepts[selection] <- prop_changepts_s
    prop_lls <- lls

    for (j in 1:ntemps){
      mod <- ts_memo(formula, data, prop_changepts[ , j], weights)
      prop_lls[j] <- tryCatch(mod$logLik, error = function(x) {-Inf})
    }

    accepts <- runif(ntemps) < exp((prop_lls - lls) * betas)
    accept_rate <- accept_rate + accepts / nit
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
    saved_lls[ , i] <- temp_ids
    saved_ids[ , i] <- lls
  }

  swap_rates <- colMeans(swap_accepted)
  out <- list(temps, saved_cps, saved_lls, saved_ids, swap_rates, accept_rate)
  names(out) <- c("temps", "changepts", "lls", "ids", "swap_rates", 
                  "accept_rate")
  return(out)
}


xxMTS <- function(formula, data, nchangepoints, weights, nit = 1e4, 
                ...){
  ts_memo <- multinom_ts

  temps <- prep_temps(...)
  betas <- 1 / temps
  ntemps <- length(betas)

  prep_cpts <- prep_changepts(formula, data, ntemps, nchangepoints, 
                 weights)
  changepts <- prep_cpts$changepts
  lls <- prep_cpts$lls

  saved_cps <- array(NA, c(nchangepoints, ntemps, nit))
  saved_lls <- matrix(NA, ntemps, nit)
  saved_ids <- matrix(NA, ntemps, nit)
  accept_rate <- 0
  temp_ids <- 1:ntemps
  swap_accepted <- matrix(FALSE, nit, ntemps - 1)

  pdist <- proposal_dist(nit, ntemps, nchangepoints, magnitude = 12)
 
  pbform <- "  [:bar] :percent eta: :eta"
  pb <- progress::progress_bar$new(pbform, nit, clear = FALSE, width = 60)

  for (i in 1:nit){
  
    pb$tick()
 
    selection <- cbind(pdist$which_kicked[i, ], 1:ntemps)
    prop_changepts <- changepts
    curr_changepts_s <- changepts[selection]
    prop_changepts_s <- curr_changepts_s + pdist$kicks[i, ]
    prop_changepts[selection] <- prop_changepts_s
    prop_lls <- lls

    for (j in 1:ntemps){
      mod <- ts_memo(formula, data, prop_changepts[ , j], weights)
      prop_lls[j] <- tryCatch(mod$logLik, error = function(x) {-Inf})
    }

    accepts <- runif(ntemps) < exp((prop_lls - lls) * betas)
    accept_rate <- accept_rate + accepts / nit
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
    saved_lls[ , i] <- temp_ids
    saved_ids[ , i] <- lls
  }

  swap_rates <- colMeans(swap_accepted)
  out <- list(temps, saved_cps, saved_lls, saved_ids, swap_rates, accept_rate)
  names(out) <- c("temps", "changepts", "lls", "ids", "swap_rates", 
                  "accept_rate")
  return(out)
}


system.time(
x1<-MTS(formula, data, nchangepoints, weights, nit = 1e3)
)

system.time(
x2<-xxMTS(formula, data, nchangepoints, weights, nit = 1e3)
)






