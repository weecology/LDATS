library(magrittr)
data(rodents)
lda_data <- dplyr::select(rodents, -c(newmoon, date, plots, traps))
ts_data <- data.frame("time" = data.frame(rodents)[ , "newmoon"])

document_term_matrix = lda_data
document_covariate_matrix = ts_data
nchangepoints = 1

formula = c~time

nit = 100
weights <- LDATS::doc_weights(document_term_matrix)
ldas <- LDATS::LDA(data = document_term_matrix, ntopics = 2:5, nseeds = 2, 
                   ncores = 4)
selected <- LDATS::LDA_select(ldas) 
mtss <- selected %>%
          LDATS::MTS_prep(document_covariate_matrix) 
data = mtss[[1]]
magnitude = 12


nit<-1e5
############################


prep_temps <- function (ntemps = 6, ultimate_temp = 2^6, k = 0, ...){
    sequence <- seq(0, log2(ultimate_temp), length = ntemps)
    log_temps <- sequence^(1 + k)/log2(ultimate_temp)^k
    temps <- 2^(log_temps)
    return(temps)
}

MTS <- function(data, formula = ~1, nchangepoints = 1, 
                weights = NULL, nit = 1e3, magnitude = 12, ...){

  character_formula <- as.character(formula)
  formula <- character_formula[length(character_formula)]
  ts_memo <- memoise::memoise(LDATS::multinom_ts)

  if(nchangepoints == 0){
    nit <- 1
  }
  temps <- prep_temps()#...)
  betas <- 1 / temps
  ntemps <- length(betas)

  MCMCsetup <- list(temps, magnitude)
  names(MCMCsetup) <- c("temperatures", "magnitude")

  prep_cpts <- LDATS::prep_changepts(data, formula, ntemps, nchangepoints, 
                 weights)
  changepts <- prep_cpts$changepts
  lls <- prep_cpts$lls

  saved_cps <- array(NA, c(nchangepoints, ntemps, nit))
  saved_lls <- matrix(NA, ntemps, nit)
  saved_ids <- matrix(NA, ntemps, nit)
  accept_rate <- 0
  temp_ids <- 1:ntemps
  swap_accepted <- matrix(FALSE, nit, ntemps - 1)

  pdist <- LDATS::proposal_dist(nit, ntemps, nchangepoints, magnitude)
 
  pbform <- "  [:bar] :percent eta: :eta"
  pb <- progress::progress_bar$new(pbform, nit, clear = FALSE, width = 60)

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
    saved_ids[ , i] <- temp_ids
    saved_lls[ , i] <- lls
  }

  swap_rates <- colMeans(swap_accepted)

  
  last_extreme <- NA
  last_extreme_vector <- numeric(nit)
  trips <- numeric(ntemps)
  for(i in 1:ntemps){
    for(j in 1:nit){
      if(saved_ids[1, j] == i){
        last_extreme <- "bottom"
      }
      if(saved_ids[ntemps, j] == i){
        last_extreme <- "top"
      }
      last_extreme_vector[j] <- last_extreme
    }
    first_top <- match("top", last_extreme_vector)
    last_pos <- rle(last_extreme_vector[first_top:nit])$values
    trips[i] <- sum(last_pos == "bottom")
  }
  trip_rate <- trips / nit

  MCMCdiagnostics <- list(accept_rates, swap_rates, trips, trip_rate)
  names(MCMCdiagnostics) <- c("acceptance_rates", "swapping_rates", 
                              "trips", "trip_rate")

  out <- list()
  out$call <- match.call()
  out$data <- data
  out$weights <- weights
  out$MCMCsetup <- MCMCsetup
  out$lls <- saved_lls[1, ]
  out$lls_full <- saved_lls
  out$MCMCdiagnostics <- MCMCdiagnostics 

  class(out) <- c("MTS", "list")
  attr(out, "hidden") <- c("data", "weights", "MCMCsetup", "lls", "lls_full",
                           "MCMCdiagnostics")



  out <- list(saved_cps, saved_lls, saved_ids, swap_rates, accept_rate)
  names(out) <- c("changepts", "lls", "ids", "swap_rates", 
                  "accept_rate")
  return(out)
}