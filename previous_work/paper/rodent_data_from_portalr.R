library(tidyverse)
library(portalr)

get_exclosure_rodents = function(cont_or_exclosures = 'control', time_or_plots = 'time') {
  if (cont_or_exclosures == 'exclosure') {
  
  if (tolower(time_or_plots) == 'plots') {
    length = 'all'
    startperiod = 118
    standardeffort = 8
  }
  if(tolower(time_or_plots) == 'time') {
    length = 'longterm'
    startperiod = 1
    standardeffort = 4
  }
  }
  if(cont_or_exclosures == 'control') {
    length = 'all'
    startperiod = 1
    standardeffort = 8
  }
  
  
  dat <- abundance(path = "repo", clean = FALSE, 
                   level = 'Plot', type = "Rodents", length = length, 
                   unknowns = FALSE, fill_incomplete = F, shape = 'crosstab',
                   time = 'period', effort = TRUE, min_plots = 0)
  
  
  if (cont_or_exclosures == 'exclosure') {
    dat2 <- dat %>%
    filter(treatment == cont_or_exclosures, period %in% startperiod:436,
           ntraps >= 1) %>%
    mutate(effort = 1) %>%
    group_by(period) %>%
    summarise_at(c(colnames(dat)[5:25], 'effort'), sum)
  } else {
    dat2 <- dat %>%
      filter(plot %in% c(2,4,8,11,12,14,17,22), 
             period %in% startperiod:436,
             ntraps >= 1) %>%
      mutate(effort = 1) %>%
      group_by(period) %>%
      summarise_at(c(colnames(dat)[5:25], 'effort'), sum)
  }
  
  datsums = vector(length =nrow(dat2))
  
  for(i in 1:nrow(dat2)) {
    thiseffort = dat2[i, 'effort']
    for (j in 2:22) {
      dat2[i,j] = round((dat2[i, j] / thiseffort) * standardeffort)
    }
    datsums[i] = sum(dat2[i,2:22])
  }
  
  dat2 = dat2[ which(datsums >= 1), 1:22]

  return(dat2)
}