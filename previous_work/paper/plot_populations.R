library(dplyr)
library(tidyr)
library(ggplot2)

dat = read.csv('Rodent_table_dat.csv')
# dates
moondat = read.csv(text=RCurl::getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/moon_dates.csv"),stringsAsFactors = F)
moondat$date = as.Date(moondat$censusdate)

period_dates = filter(moondat,period %in% rownames(dat)) %>% select(period,date)
dates = period_dates$date

dat$date = dates

longdat = gather(dat, species, abundance, BA:SO, factor_key=TRUE)


ggplot(longdat, aes(x=date,y=log(abundance),colour=species)) +
  geom_line(size=1.5)
