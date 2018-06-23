library(fitdistrplus)
library(stats)
library(dplyr)
library(RCurl)
library(ggplot2)
# scripts for making supporting figures in the extreme events/ LDA project



# ==============================================
# Abundance figure in manuscript (Figure 2)

# total abundance
dat = read.csv('Rodent_table_dat.csv',na.strings = '',as.is=T)

# find censuses that were complete
trappingdat = read.csv(text=getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent_trapping.csv"),
                       na.strings=c(""), as.is=T, stringsAsFactors = FALSE)
trapdat = aggregate(trappingdat$sampled,by=list(period=trappingdat$period),FUN=sum)
fullcensus = trapdat[trapdat$x>20,]
perioddates = read.csv(text=getURL('https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/moon_dates.csv'),as.is=T,na.strings = '')
perioddates$censusdate = as.Date(perioddates$censusdate)
fullcensus = merge(fullcensus,perioddates)

# data frame of abundance by period - restricted to only complete censuses
abund_dat = data.frame(period = 1:436, n = rowSums(dat))
abund_dat = merge(abund_dat,fullcensus[,c('period','censusdate')])
abund_dat$density = abund_dat$n/2 # density is in rodents/hectare: each plot is 1/4 hectare, there are 8 control plots, so total area covered by control plots is 2 hectares

# finding data in the lowest .15 fraction of the data. 
fitdistrplus::descdist(abund_dat$n, discrete=TRUE)
fit.negbin = fitdistrplus::fitdist(abund_dat$n, "nbinom")
plot(fit.negbin)
dist_size = fit.negbin$estimate[[1]]
dist_mu = fit.negbin$estimate[[2]]
crit_value = qnbinom(.15, size=dist_size, mu=dist_mu)
abund_dat$extreme = ifelse(abund_dat$n < crit_value, 1,0)



# in ggplot
chpts = data.frame(x1=c(as.Date('1983-12-01'),as.Date('1988-10-01'),as.Date('1998-09-01'),as.Date('2009-06-01')),
                   x2=c(as.Date('1984-07-01'),as.Date('1996-01-01'),as.Date('1999-12-01'),as.Date('2010-09-01')),
                   y1=c(-20,-20,-20,-20),y2=c(130,130,130,130))
ggplot(abund_dat,aes(x=censusdate,y=density)) +
  coord_cartesian(ylim=c(-10,110)) +
  geom_rect(data=chpts, inherit.aes=F,aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), alpha=0.5) +
  geom_line(size=1) +  
  geom_point(size=1.6, aes(fill = as.factor(abund_dat$extreme)),colour='black',pch=21) +
  geom_hline(yintercept = mean(abund_dat$density),linetype=2) +
  labs(x='',y='rodent density\n(rodents per hectare)') +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=12)) +
  scale_fill_manual(values = c('#000000',"#56B4E9")) +
  annotate("text", x=as.Date('1994-01-01'), y=-3.5, label= "Drought",fontface=2,size=3) +
  annotate("text", x=as.Date('1983-08-01'), y=-3.5, label= "Storm", fontface=2,size=3) +
  annotate("text", x=as.Date('1999-08-01'), y=-3.5, label= "Storm", fontface=2,size=3) +
  annotate("text", x=as.Date('2009-06-01'), y=-3.5, label= "Drought",fontface=2,size=3) +
  geom_segment(aes(x=as.Date('1993-09-01'), xend=as.Date('1994-10-01'), y=-10, yend=-10),size=1) +
  geom_segment(aes(x=as.Date('1993-09-01'), xend=as.Date('1993-09-01'), y=-12, yend=-8), size=1) +
  geom_segment(aes(x=as.Date('1994-10-01'), xend=as.Date('1994-10-01'), y=-12, yend=-8), size=1) +
  geom_segment(aes(x=as.Date('2009-01-01'), xend=as.Date('2009-12-31'), y=-10, yend=-10),size=1) +
  geom_segment(aes(x=as.Date('2009-01-01'), xend=as.Date('2009-01-01'), y=-12, yend=-8), size=1) +
  geom_segment(aes(x=as.Date('2009-12-31'), xend=as.Date('2009-12-31'), y=-12, yend=-8), size=1) +
  geom_point(size=2, inherit.aes=F,aes(x=as.Date('1983-10-01'),y=-10),pch=15) +
  geom_point(size=2, inherit.aes=F,aes(x=as.Date('1999-08-15'),y=-10),pch=15) +
  theme(legend.position = 'none') 

ggsave(filename='Figure2.tiff',width=6,height=2.8,units='in',dpi=600,compression='lzw')
