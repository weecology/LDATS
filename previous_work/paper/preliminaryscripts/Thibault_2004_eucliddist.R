
# Trying to re-create Thibault et al 2004


# read in latest Rodent data, make a table of species/time
rdat = read.csv('C:/Users/EC/Desktop/git/PortalData/Rodents/Portal_rodent.csv',stringsAsFactors = F,na.strings = '')
granivores = c('BA','DM','DO','DS','PB','PH','PI','PP','PF','PE','PM','RM','RO','RF')

# read in trapping data and find censuses where all control plots were trapped
trapdat = read.csv('C:/Users/EC/Desktop/git/PortalData/Rodents/Portal_rodent_trapping.csv')
allplots = aggregate(trapdat$Sampled,by=list(period=trapdat$Period),FUN=sum) %>% filter(x>21)

# create table of species abundances
absabund = rdat %>% filter(period>7, period< 295, period %in% allplots$period, species %in% granivores, plot %in% c(1,2,4,8,9,11,12,14,17,22)) %>% select(period,species) %>% table()

# load trapping date data
pdat = read.csv('C:/Users/EC/Desktop/git/PortalData/Rodents/moon_dates.csv',stringsAsFactors = F)
pdat$CensusDate = as.Date(pdat$CensusDate)

# set up for yearly aggregation
absdat_dates = merge(absabund,pdat,by.x='period',by.y='Period')
absdat_dates$year = format(absdat_dates$CensusDate,'%Y')

# make column of total abundance per period
period_totalabund = aggregate(absdat_dates$Freq,by=list(period=absdat_dates$period),FUN=sum)
absdat_dates= merge(absdat_dates,period_totalabund)
absdat_dates$relabund = absdat_dates$Freq/absdat_dates$x

# relative abundance by year
reldat_yr = aggregate(absdat_dates$relabund,by=list(year = absdat_dates$year,species=absdat_dates$species),FUN=mean)

reltable = data.frame(reshape(reldat_yr,timevar='species',idvar='year',direction='wide'))

ED_yr = as.matrix(vegdist(reltable[,-1],method='euclidean'))


#yearly_table = data.frame(reshape(absdat_yr,timevar='species',idvar='year',direction='wide'))
#yearly_avg_total = data.frame(year = as.numeric(yearly_table$year),total = rowSums(yearly_table[,-1]))

#relative_yr = yearly_table[,-1]/yearly_avg_total$total
# adjacent pairs
ED_yr = as.matrix(vegdist(relative_yr,method='euclidean'))

bylag = c()
for (r in seq(nrow(ED_yr))) {
  for (c in r:nrow(ED_yr)) {
    bylag = rbind(bylag,c(c-r,ED_yr[r,c]))
  }
}
df = as.data.frame(bylag)

#====================
# all the lag times (Thibault 2004) (only used the 14 granivores)


ggplot(data.frame(df),aes(x=V1,y=V2,group=V1)) +
  geom_boxplot() +
  scale_y_continuous(name='Euclidean Distance') +
  scale_x_continuous(name='Lag Years')
