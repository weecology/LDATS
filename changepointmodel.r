if (packageVersion("memoise") <= "1.0.0") {
  devtools::install_github("hadley/memoise")
}

library(memoise)     # For avoiding redundante computations
library(lubridate)   # For dates
library(progress)    # For progress bar
library(topicmodels) # For LDA
library(ggplot2)
library(viridis)
library(nnet)        # For multinomial model (part of changepoint analysis)
library(RColorBrewer)

my_colors = brewer.pal(5, name = "Set2")

# Begin changepoint model -------------------------------------------------


# Fit a model to dates in (start,end] and return the log likelihood.
# Memoization is a trick that lets us save the output for a given chunk and
# avoid finding the answer more than once.
fit_chunk_non_memoized = function(ldamodel, x, start, end, make_plot = FALSE, ...) {
  # Weights average to 1, & are proportional to total rodents caught that month
  m = multinom(
    ldamodel@gamma ~ year_continuous + sin_year + cos_year, 
    data = x,
    maxit = 1E5,
    weights = rowMeans(dat) / mean(as.matrix(dat)),
    subset = x$year_continuous > start & x$year_continuous <= end,
    trace = FALSE
  )
  
  if (make_plot) {
    plotfun = ifelse(start == -Inf, matplot, matlines)
    plotfun(
      x$year_continuous[x$year_continuous > start & x$year_continuous <= end], 
      fitted(m),
      ylim = c(0, 1), 
      xlim = range(x$year_continuous),
      type = "l",
      lty = 1,
      col = my_colors,
      ...
    )
    abline(v = start)
  }
  
  logLik(m)
}


# Get the log-likelihood associated with a set of breakpoints
get_ll_non_memoized = function(ldamodel, x, changepoints, make_plot = FALSE, ...){
  if (make_plot) {
    fit_chunk = fit_chunk_non_memoized
  }
  
  if (any(changepoints <= 0) | any(changepoints >= length(x$year_continuous)) | 
      is.unsorted(changepoints, strictly = TRUE)) {
    return(-Inf)
  }
  
  
  changedates = c(-Inf, x$year_continuous[changepoints], Inf)
  sum(
    sapply(
      seq_len(length(changedates) - 1),
      function(i){
        fit_chunk(ldamodel, x, changedates[i], changedates[i + 1], make_plot = make_plot, ...)
      }
    )
  )
}


# Saving the caches as hidden folders to prevent silly Mac computers
# (and RStudio) from wasting resources trying to index them
fit_chunk = memoise(fit_chunk_non_memoized, 
                    cache = cache_filesystem(".cache_chunk"))
get_ll = memoise(get_ll_non_memoized, 
                 cache = cache_filesystem(".cache_ll"))





# k is the exponent controlling the temperature sequence: 0 implies
# geometric sequence, 1 implies squaring before exponentiating.
# Use larger values if the cooler chains aren't swapping enough.
changepoint_model = function(ldamodel,
                             x,
                             n_changepoints,
                             N_temps = 6,
                             maxit = 1E4,
                             penultimate_temp = 2^6,
                             k = 0){
  
  # Temperature sequence
  sequence = seq(0, log2(penultimate_temp), length.out = N_temps - 1)
  log_temps = sequence^(1 + k) / log2(penultimate_temp)^k
  temps = 2^(log_temps)
  temps = c(temps, 1E10) # Highest temperature
  betas = 1/temps # "inverse temperature"
  
  # Initialize randomly, with the best starting values in the coldest chain
  changepoints = replicate(N_temps, sort(sample.int(length(x$year_continuous), n_changepoints)))
  lls = sapply(1:N_temps, function(j){get_ll(ldamodel, x, changepoints[ , j])})
  changepoints = changepoints[ , order(lls, decreasing = TRUE)]
  lls = sort(lls, decreasing = TRUE)
  
  saved = array(NA, c(n_changepoints, N_temps, maxit))
  saved_lls = matrix(NA, N_temps, maxit)
  saved_ids = saved_lls
  
  accept_rate = 0
  ids = 1:N_temps
  swap_accepted = matrix(FALSE, maxit, N_temps - 1)
  
  # Pre-calculate proposal distributions
  kick_signs = sample(c(-1, 1), maxit * N_temps, replace = TRUE)
  kick_magnitudes = 1 + rgeom(maxit * N_temps, 1/12)
  kicks = matrix(kick_signs * kick_magnitudes, nrow = maxit)
  which_kicked = matrix(
    sample.int(n_changepoints, maxit * N_temps, replace = TRUE),
    nrow = maxit
  )
  
  pb = progress_bar$new(format = "  [:bar] :percent eta: :eta",
                        total = maxit, clear = FALSE, width = 60)
  for (i in 1:maxit) {
    pb$tick()
    # Make proposals for each temperature
    proposed_changepoints = changepoints
    proposed_changepoints[cbind(which_kicked[i, ], 1:N_temps)] = 
      changepoints[cbind(which_kicked[i, ], 1:N_temps)] + kicks[i, ]
    proposed_lls = sapply(1:N_temps, function(j){get_ll(ldamodel, x, proposed_changepoints[ , j])})
    
    # Accept some proposals via Metropolis rule; update the changepoints
    # and the associated log-likelihoods
    accepts = runif(N_temps) < exp((proposed_lls - lls) * betas)
    accept_rate = accept_rate + accepts / maxit
    changepoints[ , accepts] = proposed_changepoints[ , accepts]
    lls[accepts] = proposed_lls[accepts]
    
    for (j in seq(N_temps - 1, 1)) {
      # Propose a swap between temperature j and temperature j+1
      accept_swap = runif(1) < exp((betas[j] - betas[j + 1]) * (lls[j + 1] - lls[j]))
      if (accept_swap) {
        swap_accepted[i, j] = TRUE
        
        # Swap changepoint vectors between MCMC replicas
        placeholder = changepoints[, j]
        changepoints[ , j] = changepoints[, j + 1]
        changepoints[ , j + 1] = placeholder
        
        # Swap the associated log-likelihood values
        placeholder = lls[j]
        lls[j] = lls[j + 1]
        lls[j + 1] = placeholder
        
        placeholder = ids[j]
        ids[j] = ids[j + 1]
        ids[j + 1] = placeholder
      }
    }
    saved[,,i] = changepoints
    saved_ids[,i] = ids
    saved_lls[,i] = lls
  }
  
  list(
    temps = temps,
    saved = saved,
    saved_ids = saved_ids,
    saved_lls = saved_lls,
    swap_rates = colMeans(swap_accepted),
    accept_rate = accept_rate
  )
}



# Functions for viewing/diagnosing the changepoints -----------------------

# Number of times the particle went from hottest chain to the coldest one,
# indicating good mixing.  Should probably be in the mid-hundreds or low
# thousands if we want to be really confident about the results.
count_trips = function(results){
  N_temps = length(results$accept_rate)
  maxit = ncol(results$saved_lls)
  sapply(
    1:N_temps,
    function(k){
      last_extreme = NA
      last_extreme_vector = numeric(maxit)
      for (i in 1:maxit) {
        if (results$saved_ids[1, i] == k) {
          last_extreme = "bottom"
        }
        if (results$saved_ids[N_temps, i] == k) {
          last_extreme = "top"
        }
        last_extreme_vector[i] = last_extreme
      }
      
      first_top = match("top", last_extreme_vector)
      
      sum(rle(last_extreme_vector[first_top:maxit])$values == "bottom")
    }
  )
}


# Histogram showing percentage of MCMC samples that contained
# a changepoint in a given year.
annual_hist = function(results, year_continuous){
  if (missing(year_continuous)) {
    year_continuous = 1:max(results$saved) / 12
  }
  hist(year_continuous[results$saved[,1,]], 
       breaks = seq(0, 3000), xlim = range(year_continuous), 
       axes = FALSE, yaxs = "i", ylim = c(0, 1.04 * length(results$saved[1,1,])))
  axis(2, seq(0, 1, 0.25) * length(results$saved[1,1,]), 
       seq(0, 1, .25))
  axis(1)
}

# # =========================================================================
# # Run the model
# 
# nstart = 20 # For the final analysis, maybe do 1000
# ldamodel2 = LDA(dat,2,control=list(estimate.alpha=F,alpha=1, nstart = nstart),method="VEM")
# 
#results3_3 = changepoint_model(ldamodel3, x, 3)
#annual_hist(results,year_continuous)
#get_ll_non_memoized(ldamodel3,x,2,)