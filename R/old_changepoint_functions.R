#' Fit a model to dates in (start,end] and return the log likelihood.
#' 
#' @param ldamodel ldamodel
#' @param x x
#' @param start start
#' @param end end 
#' @param make_plot make_plot
#' @param weights weights
#' @examples
#' NA
#' @export 

fit_chunk_non_memoized = function(ldamodel, x, start, end, make_plot = FALSE, 
                                  weights, ...) {
  # Weights average to 1, & are proportional to total rodents caught that month
  m = nnet::multinom(
    ldamodel@gamma ~ sin_year + cos_year, 
    data = x,
    maxit = 1E5,
    weights = weights,
    subset = x$year_continuous > start & x$year_continuous <= end,
    trace = FALSE
  )
  
  
  logLik(m)
}


#' Fit a model to dates in (start,end] and return data frame for plotting
#' 
#' @param ldamodel output of LDA()
#' @param x same x used in changepoint_model()
#' @param start value: start of the section to be fit
#' @param end value: end of the section to be fit
#' @param weights same weights used in changepoint_model()
#' @export 

fit_section = function(ldamodel, x, start, end, weights, ...) {
  m = nnet::multinom(
    ldamodel@gamma ~ sin_year + cos_year, 
    data = x,
    maxit = 1E5,
    weights = weights,
    subset = x$year_continuous > start & x$year_continuous <= end,
    trace = FALSE
  )
  
  section_df = as.data.frame(fitted(m))
  section_df$date = format(date_decimal(x$year_continuous[x$year_continuous > start & x$year_continuous <= end]), '%Y-%m-%d') %>% as.Date()
  section = melt(section_df,id.var='date')
  return(section)
}


#' Get the log-likelihood associated with a set of breakpoints
#' @param ldamodel ldamodel
#' @param x x
#' @param changepoints changepoints
#' @param make_plot make_plot
#' @param weights weights
#' @export 

get_ll_non_memoized = function(ldamodel, x, changepoints, make_plot = FALSE, 
                               weights, ...){

  fit_chunk = memoise::memoise(fit_chunk_non_memoized, 
                      cache = memoise::cache_filesystem(".cache_chunk"))
  
  changedates = c(-Inf, x$year_continuous[changepoints], Inf)
  sum(
    sapply(
      seq_len(length(changedates) - 1),
      function(i){
        fit_chunk(ldamodel, x, changedates[i], changedates[i + 1], 
                  make_plot = make_plot, weights = weights, ...)
      }
    )
  )
}

#' Model to locate given number of change points in an LDA model output
#' 
#' @param ldamodel output object from LDA model using VEM
#' @param x x
#' @param n_changepoints integer number of changepoints to look for
#' @param maxit maximum iterations
#' @param penultimate_temp penultimate_temp
#' @param k the exponent controlling the temperature sequence: 0 implies geometric sequence, 
#'          1 implies squaring before exponentiating. Use larger values if the cooler chains aren't swapping enough.
#'          
#' @return 
#' 
#' @author Dave Harris
#' @export 

changepoint_model = function(ldamodel,
                             x,
                             n_changepoints,
                             N_temps = 6,
                             maxit = 1E5,
                             penultimate_temp = 2^6,
                             k = 0,
                             weights){
  
  file.remove(dir(".cache_chunk/", full.names = TRUE))
  file.remove(dir(".cache_ll/", full.names = TRUE))
  
  get_ll = memoise::memoise(get_ll_non_memoized, 
                   cache = memoise::cache_filesystem(".cache_ll"))
  
  
  # Temperature sequence
  sequence = seq(0, log2(penultimate_temp), length.out = N_temps - 1)
  log_temps = sequence^(1 + k) / log2(penultimate_temp)^k
  temps = 2^(log_temps)
  temps = c(temps, 1E10) # Highest temperature
  betas = 1/temps # "inverse temperature"
  # Initialize randomly, with the best starting values in the coldest chain
  changepoints = matrix(
    replicate(N_temps, sort(sample.int(length(x$year_continuous) - 1, 
                                       n_changepoints, replace = FALSE))),
    ncol = N_temps
  )
  lls = sapply(1:N_temps, 
               function(j){get_ll(ldamodel, x, changepoints[ , j], 
                                  weights = weights)})
  changepoints = changepoints[ , order(lls, decreasing = TRUE), drop = FALSE]
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
  
  pb = progress::progress_bar$new(format = "  [:bar] :percent eta: :eta",
                        total = maxit, clear = FALSE, width = 60)
  for (i in 1:maxit) {
    pb$tick()
    # Make proposals for each temperature
    proposed_changepoints = changepoints
    proposed_changepoints[cbind(which_kicked[i, ], 1:N_temps)] = 
      changepoints[cbind(which_kicked[i, ], 1:N_temps)] + kicks[i, ]
    proposed_lls = sapply(1:N_temps, function(j){get_ll(ldamodel, x, proposed_changepoints[ , j],
                                                        weights = weights)})
    
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



#' Run a suite of changepoint models for the same data set
#' 
#' @param data data set to use
#' @param time time of each sample (row in data)
#' @param ntopics number of topics to use in the baseline LDA
#' @param SEED seed to use in the baseline LDA
#' @param weights weights to use through time for the data points
#' @param maxit maximum iterations
#' @param maxcp maximum number of change points to use in the models
#'          
#' @return list of change point model results
#' 
#' @author Juniper Simonis
#' @export 

cp_models <- function(data = NULL, dates = NULL, ntopics = NULL, 
                      SEED = NULL, weights = NULL, maxit = NULL, 
                      maxcps = NULL){

  # run baseline lda model 

    bl_lda <- topicmodels::LDA(data, ntopics, 
                               control = list(seed = SEED), method = 'VEM')

  # Change point model

    # set up time for model

      cds <- as.Date(as.character(dates))
      year_continuous <- 1970 + as.integer(julian(cds)) / 365.25
      x <- data.frame(year_continuous = year_continuous,
                               sin_year = sin(year_continuous * 2 * pi),
                               cos_year = cos(year_continuous * 2 * pi))


    # run models with 0:maxcps changepoints

      output <- vector("list", 1 + maxcps)

      # 0 change points

        m <- nnet::multinom(bl_lda@gamma ~ sin_year + cos_year, 
                            data = x, maxit = maxit, weights = weights,
                            trace = FALSE) 
        output[[1]] <- m
        names(output)[1] <- "0 changepoints"

      # 1:maxcps changepoints

        for(i in 1:maxcps){
  
          print(paste("Running model with ", i, " changepoint(s).", sep = ""))
          cpm <- changepoint_model(bl_lda, x, i, maxit = maxit, 
                                   weights = weights)
          output[[i + 1]] <- cpm
          names(output)[i + 1] <- paste(i, " changepoints", sep = "")
        }


  return(output)
}

# Functions for viewing/diagnosing the changepoints -----------------------

#' Number of times the particle went from hottest chain to the coldest one,
#' indicating good mixing.  Should probably be in the mid-hundreds or low
#' thousands if we want to be really confident about the results.
#' @param results results
#' @export 

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


#' Histogram showing percentage of MCMC samples that contained
#' a changepoint in a given year.
#' @param results results
#' @param year_continuous year_continuous
#' @export 
annual_hist = function(results, year_continuous){
  if (missing(year_continuous)) {
    year_continuous = 1:max(results$saved) / 12
  }
  hist(year_continuous[results$saved[,1,]], 
       breaks = seq(0, 3000), xlim = range(year_continuous), 
       axes = FALSE, yaxs = "i", ylim = c(0, 1.04 * length(results$saved[1,1,])),main='',xlab='')
  axis(2, seq(0, 1, 0.25) * length(results$saved[1,1,]), 
       seq(0, 1, .25))
  axis(1)
}



#' Find changepoint locations
#' 
#' @param results results object output from changepoint_model
#' 
#' @param return vector of changepoint locations
#' 
#' @author Erica Christensen
#' 
#' @export 
find_changepoint_location = function(results) {
  cpts = c()
  for (n in seq(dim(results$saved)[1])) {
    cp = results$saved[n,1,]
    x = unique(cp)[which.max(tabulate(match(cp, unique(cp))))]
    cpts = append(cpts,x)
  }
  return(cpts)
}


