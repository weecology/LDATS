
MTS <- function(data, formula = ~1, nchangepoints = 1, weights = NULL, 
                nit = 1e3, magnitude = 12, burnin = 1e2, ...){



  if(nchangepoints == 0){ # consider handling this cleaner
    nit <- 1
  }   	

  temps <- prep_temps() #...
  betas <- 1 / temps
  ntemps <- length(betas)

  prep_cpts <- prep_changepts(data, formula, ntemps, nchangepoints, weights)
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

