#Script for making the euclidean distance boxplot (Figure 1)
library(vegan)
library(dplyr)
library(ggplot2)
library(multipanelfigure)

# read in latest Rodent data, make a table of species/time
rdat = read.csv('C:/Users/EC/Desktop/git/PortalData/Rodents/Portal_rodent.csv',stringsAsFactors = F,na.strings = '')
granivores = c('BA','DM','DO','DS','PB','PH','PI','PP','PF','PE','PM','RM','RO','RF')

# read in trapping data and find censuses where all control plots were trapped
trapdat = read.csv('C:/Users/EC/Desktop/git/PortalData/Rodents/Portal_rodent_trapping.csv')
allplots = aggregate(trapdat$sampled,by=list(period=trapdat$period),FUN=sum) %>% filter(x>21)

# create table of species abundances
absabund = rdat %>% filter(period>0, period< 437, period %in% allplots$period, species %in% granivores, plot %in% c(2,4,8,11,12,14,17,22)) %>% select(period,species) #%>% table()

# load trapping date data
pdat = read.csv('C:/Users/EC/Desktop/git/PortalData/Rodents/moon_dates.csv',stringsAsFactors = F)
pdat$censusdate = as.Date(pdat$censusdate)

# set up for yearly aggregation
absdat_dates = merge(absabund,pdat,by.x='period',by.y='period')
absdat_dates$year = format(absdat_dates$censusdate,'%Y')

absdat_yr = aggregate(absdat_dates$Freq,by=list(year = absdat_dates$year,species=absdat_dates$species),FUN=mean)

yearly_table = data.frame(reshape(absdat_yr,timevar='species',idvar='year',direction='wide'))
yearly_avg_total = data.frame(year = as.numeric(yearly_table$year),total = rowSums(yearly_table[,-1]))

relative_yr = yearly_table[,-1]/yearly_avg_total$total
# adjacent pairs
ED_yr = as.matrix(vegdist(relative_yr,method='euclidean'))


#====================
# all the lag times (Thibault 2004) (only used the 14 granivores)
bylag = c()
for (n in 1:37) {
  bylag = rbind(bylag,cbind(diag(ED_yr[,-(1:n)]),rep(n)))
}

df = data.frame(bylag)
boxfigure = ggplot(data.frame(bylag),aes(x=X2,y=X1,group=X2)) +
  geom_boxplot() +
  scale_y_continuous(name='Euclidean Distance') +
  scale_x_continuous(name='Lag Years')
boxfigure
#boxplot(bylag[,1]~bylag[,2],xlab='lag years',ylab='Euclidean dist', main='by lag time')


# =======================================
# multi panel figure
(figure1 <- multi_panel_figure(
  width = c(150,10),
  height = c(70,70),
  column_spacing = 0,
  panel_label_type = "lower-alpha"))
figure1 %<>% fill_panel(
  boxfigure,
  row = 1, column = 1)
figure1 %<>% fill_panel(
  'C:/Users/EC/Dropbox/Photos/Plot1beforeafter.tif',
  row = 2, column = 1:2)

figure1
