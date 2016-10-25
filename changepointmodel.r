if (packageVersion("memoise") <= "1.0.0") {
  devtools::install_github("hadley/memoise")
}

library(memoise)     # For avoiding redundante computations
library(lubridate)   # For dates
library(progress)
library(topicmodels) # For LDA
library(ggplot2)
library(viridis)
library(nnet)        # For multinomial model (part of changepoint analysis)
library(RColorBrewer)


my_colors = brewer.pal(5, name = "Set2")

# load csv of abundance of each sp per plot
dat = read.csv('Rodent_table_dat.csv',na.strings = '',as.is=T)

# how many species
#nsp = length(names(dat))

# dates of trapping periods
perdat = read.csv('Period_dates_single.csv')
perdat$date = as.Date(perdat$date,format='%m/%d/%Y')

#date_dat = dat
#date_dat$year = perdat$yr[1:length(date_dat[,1])]

#yearly_dat = aggregate(date_dat[,1:21],by=list(year=date_dat$year),FUN=mean)
#yearly_dat = round(yearly_dat[-1],0)

#=======================
# LDA

k = 4
ldamodel = LDA(dat,k,control=list(seed=30,estimate.alpha=F,alpha=1),method="VEM")
set.seed(12345)

n_changepoints = 4 

d = ymd(as.Date(perdat$date[1:length(dat[,1])]))
yday = yday(d)
year_continuous = 1970 + as.integer(julian(d)) / 365.25
x = data.frame(
  sin_year = sin(year_continuous * 2 * pi),
  cos_year = cos(year_continuous * 2 * pi)
)


# Begin changepoint model -------------------------------------------------


changepoint_model = function(ldamodel,n_changepoints) {
  
  # Simulated annealing uses a "cooling schedule" or "annealing schedule" that 
  # drops the temperature from very high (i.e. very random, very easy to kick the 
  # changepoints into bad configurations) to very low (where the model freezes up 
  # an optimum). When the termperature is 1, we're basically doing MCMC on the
  # original model. Schedules of the form a/(x+b) have good theoretical properties
  # for convergence, but they fall too quickly at the start of the schedule, so I
  # use a more convenient scheudle to get us down to a temperature of 1 and then 
  # do a/(a + x) from there to the end.
  {
    a = 1E4
    switch = 5E5
    start_temp = 10
    cooling = function(x, shift, scale){
      out = a / (a + x - switch)
      ifelse(x < switch, (start_temp + 1) - x/switch * start_temp, out)
    }
    par(mfrow = c(1, 2))
    curve(cooling(x, shift, scale), to = 1E6, n = 1E4, log = "y", type = "p")
    curve(cooling(x, shift, scale), to = 1E6, n = 1E4, log = "", yaxs = "i", ylim = c(0, 20), xaxs = "i")
    par(mfrow = c(1, 1))
  }
  
  # Fit a model to dates in (start,end] and return the log likelihood.
  # Memoization is a trick that lets us save the output for a given chunk and
  # avoid finding the answer more than once.
  if (!exists("fit_chunk")) {
    fit_chunk_non_memoized = function(start, end, make_plot = FALSE, ...) {
      # Weights average to 1, & are proportional to total rodents caught that month
      m = multinom(
        ldamodel@gamma ~ sin_year + cos_year, 
        data = x,
        maxit = 1E5,
        weights = rowMeans(dat) / mean(as.matrix(dat)),
        subset = year_continuous > start & year_continuous <= end,
        trace = FALSE
      )
      
      if (make_plot) {
        plotfun = ifelse(start == -Inf, matplot, matlines)
        plotfun(
          year_continuous[year_continuous > start & year_continuous <= end], 
          fitted(m),
          ylim = c(0, 1), 
          xlim = range(year_continuous),
          type = "l",
          lty = 1,
          col = my_colors,
          ...
        )
        abline(v = start)
      }
      
      logLik(m)
    }
    fit_chunk = memoise(fit_chunk_non_memoized, cache = cache_filesystem("cache_chunk"))
  }
  
  # Get the log-likelihood associated with a set of breakpoints
  if (!exists("get_ll")) {
    get_ll_non_memoized = function(changepoints, make_plot = FALSE, ...){
      if(make_plot){
        fit_chunk = fit_chunk_non_memoized
      }
      
      if (any(changepoints <= 0) | any(changepoints >= length(year_continuous)) | 
          is.unsorted(changepoints, strictly = TRUE)) {
        return(-Inf)
      }
      
      
      changedates = c(-Inf, year_continuous[changepoints], Inf)
      sum(
        sapply(
          seq_len(length(changedates) - 1),
          function(i){
            fit_chunk(changedates[i], changedates[i + 1], make_plot = make_plot, ...)
          }
        )
      )
    }
    get_ll = memoise(get_ll_non_memoized, cache = cache_filesystem("cache_ll"))
    
  }
  
  # Set initial changepoints randomly
  changepoints = sort(sample.int(length(year_continuous), n_changepoints))
  ll = get_ll(changepoints)
  
  maxit = 1E6
  cp = matrix(NA, n_changepoints, maxit / 1E3)
  
  kick_signs = sample(c(-1, 1), maxit, replace = TRUE)
  kick_magnitudes = 1 + rbinom(maxit, size = 11, prob = 1/4)
  
  which_kicked = sample.int(n_changepoints, maxit, replace = TRUE)
  runifs = runif(maxit)
  
  max_observed_ll = -Inf
  
  temperatures = cooling(1:maxit, shift, scale)
  
  
  #plot(temperatures, type = "l", log = "y")
  #plot(temperatures, type = "l", log = "")
  
  
  for (i in 1:maxit) {
    T = temperatures[i]
    
    # Choose a changepoint at random and either kick it one step to the left or
    # to the right.
    proposed_changepoints = changepoints
    proposed_changepoints[which_kicked[i]] = 
      proposed_changepoints[which_kicked[i]] + kick_signs[i] * kick_magnitudes[i]
    
    # Metropolis decision rule for accepting proposals
    proposed_ll = get_ll(proposed_changepoints)
    energy_difference = (proposed_ll - ll) / T
    if (exp(energy_difference) > runifs[i]) {
      changepoints = proposed_changepoints
      ll = proposed_ll
    }
    
    # Print status
    if (i %% 1000 == 0) {
      print(paste("i:", i, "  T:", T, "  ll:", ll))
      cp[ ,i / 1000] = changepoints
    }
    
    # Save the best value observed so far
    if (ll > max_observed_ll) {
      max_observed_ll = ll
      best_changepoints = changepoints
    }
  }
  
  
  
  # Plot changepoint model & LDA output over time
  par(mfrow = c(1, 2))
  get_ll_non_memoized(best_changepoints, make_plot = TRUE, lwd = 2)
  matplot(year_continuous, ldamodel@gamma, type = "l", lty = 1, col = my_colors, lwd = 2)
  abline(v = year_continuous[best_changepoints], lwd = 2)
  par(mfrow = c(1, 1))
  
  
  matplot(
    seq(1E3, maxit, 1E3),
    matrix(year_continuous[t(cp)], ncol = n_changepoints), 
    col = brewer.pal(5, "Dark2"), 
    type = "l", 
    lty = 1,
    xlab = "iterations",
    ylim = range(year_continuous),
    yaxs = "i",
    ylab = "year of breakpoint",
    lwd = 1/2
  )
  
  #beepr::beep()
}
