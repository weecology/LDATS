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

# exclude censuses that were incomplete
trappingdat = read.csv(text=getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent_trapping.csv"),
                       na.strings=c(""), as.is=T, stringsAsFactors = FALSE)
trapdat = aggregate(trappingdat$sampled,by=list(period=trappingdat$period),FUN=sum)
fullcensus = trapdat[trapdat$x>20,]
perioddates = read.csv(text=getURL('https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/moon_dates.csv'),as.is=T,na.strings = '')
perioddates$censusdate = as.Date(perioddates$censusdate)
fullcensus = merge(fullcensus,perioddates)

# data frame of abundance by period
abund_dat = data.frame(period = 1:436, n = rowSums(dat))
abund_dat = merge(abund_dat,fullcensus[,c('period','censusdate')])
abund_dat$density = abund_dat$n/2

# finding data in the lowest .15 fraction of the data. 
descdist(abund_dat$n, discrete=TRUE)
fit.negbin = fitdist(abund_dat$n, "nbinom")
plot(fit.negbin)
dist_size = fit.negbin$estimate[[1]]
dist_mu = fit.negbin$estimate[[2]]
crit_value = qnbinom(.15, size=dist_size, mu=dist_mu)
abund_dat$extreme = ifelse(abund_dat$n < crit_value, 1,0)

# plot
plot(abund_dat$censusdate,abund_dat$density,xlab='',ylab='Total Abundundance',pch=19,ylim=c(0,150),cex.axis=1.5,cex.lab=1.5)
# gray boxes showing changepoint 95% intervals
rect(xleft = as.Date('1983-12-01'),xright = as.Date('1984-07-01'),ytop = 250,ybottom=-10,col='gray',border=NA)
rect(xleft = as.Date('1988-10-01'),xright = as.Date('1996-01-01'),ytop = 250,ybottom=-10,col='gray',border=NA)
rect(xleft = as.Date('1998-09-01'),xright = as.Date('1999-12-01'),ytop = 250,ybottom=-10,col='gray',border=NA)
rect(xleft = as.Date('2009-06-01'),xright = as.Date('2010-09-01'),ytop = 250,ybottom=-10,col='gray',border=NA)
lines(abund_dat$censusdate,abund_dat$density)
points(abund_dat$censusdate,abund_dat$density,pch=16,col=as.factor(abund_dat$extreme))
abline(h=mean(abund_dat$density))
box(which='plot')

# line segments showing changepoint 95% intervals
rect(xleft = as.Date('1999-07-01'),xright = as.Date('1999-10-01'),ytop = 250,ybottom=0,col='gray',border=NA)
rect(xleft = as.Date('1983-08-01'),xright = as.Date('1983-11-01'),ytop = 250,ybottom=0,col='gray',border=NA)
rect(xleft = as.Date('1993-09-01'),xright = as.Date('1994-10-01'),ytop = 250,ybottom=0,col='gray',border=NA)
rect(xleft = as.Date('2009-01-01'),xright = as.Date('2009-12-31'),ytop = 250,ybottom=0,col='gray',border=NA)
segments(as.Date('1983-12-01'),200 , as.Date('1984-07-01'), 200, col='red', lwd=3, xpd = FALSE)
segments(as.Date('1988-10-01'),200 , as.Date('1996-01-01'), 200, col='red', lwd=3, xpd = FALSE)
segments(as.Date('1998-09-01'),200 , as.Date('1999-12-01'), 200, col='red', lwd=3, xpd = FALSE)
segments(as.Date('2009-06-01'),200 , as.Date('2010-09-01'), 200, col='red', lwd=3, xpd = FALSE)

# in ggplot
chpts = data.frame(x1=c(as.Date('1983-12-01'),as.Date('1988-10-01'),as.Date('1998-09-01'),as.Date('2009-06-01')),
                   x2=c(as.Date('1984-07-01'),as.Date('1996-01-01'),as.Date('1999-12-01'),as.Date('2010-09-01')),
                   y1=c(-20,-20,-20,-20),y2=c(130,130,130,130))
ggplot(abund_dat,aes(x=censusdate,y=density)) +
  coord_cartesian(ylim=c(-10,110)) +
  geom_rect(data=chpts, inherit.aes=F,aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), alpha=0.5) +
  geom_line(size=1.5) +  
  geom_point(size=3, aes(fill = as.factor(abund_dat$extreme)),colour='black',pch=21) +
  geom_hline(yintercept = mean(abund_dat$density),linetype=2) +
  labs(x='',y='rodent density (rodents per hectare)') +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14)) +
  scale_fill_manual(values = c('#000000',"#56B4E9")) +
  annotate("text", x=as.Date('1994-01-01'), y=-4, label= "Drought",fontface=2) +
  annotate("text", x=as.Date('1983-08-01'), y=-4, label= "Storm", fontface=2) +
  annotate("text", x=as.Date('1999-08-01'), y=-4, label= "Storm", fontface=2) +
  annotate("text", x=as.Date('2009-06-01'), y=-4, label= "Drought",fontface=2) +
  geom_segment(aes(x=as.Date('1993-09-01'), xend=as.Date('1994-10-01'), y=-10, yend=-10),size=1.5) +
  geom_segment(aes(x=as.Date('1993-09-01'), xend=as.Date('1993-09-01'), y=-12, yend=-8), size=1.5) +
  geom_segment(aes(x=as.Date('1994-10-01'), xend=as.Date('1994-10-01'), y=-12, yend=-8), size=1.5) +
  geom_segment(aes(x=as.Date('2009-01-01'), xend=as.Date('2009-12-31'), y=-10, yend=-10),size=1.5) +
  geom_segment(aes(x=as.Date('2009-01-01'), xend=as.Date('2009-01-01'), y=-12, yend=-8), size=1.5) +
  geom_segment(aes(x=as.Date('2009-12-31'), xend=as.Date('2009-12-31'), y=-12, yend=-8), size=1.5) +
  geom_point(size=3, inherit.aes=F,aes(x=as.Date('1983-10-01'),y=-10),pch=15) +
  geom_point(size=3, inherit.aes=F,aes(x=as.Date('1999-08-15'),y=-10),pch=15) +
  theme(legend.position = 'none') 



geom_segment(aes(x=as.Date('1983-10-01'), xend=as.Date('1983-10-01'), y=60, yend=50), 
             arrow = arrow(length = unit(0.35, "cm"))) +
  geom_segment(aes(x=as.Date('1999-08-15'), xend=as.Date('1999-08-15'), y=60, yend=50), 
               arrow = arrow(length = unit(0.35, "cm"))) +
  
  
  par(xpd=F)


# ======================================================================
# make this plot piece-wise for presentation figure
plot(abund_dat$censusdate[1:96],log(abund_dat$n[1:96]),xlab='',ylab='Log Total Abund.',pch=19,ylim=c(2,6),xlim=range(abund_dat$censusdate),cex.axis=1.5,cex.lab=1.5)
rect(xleft = as.Date('1983-08-01'),xright = as.Date('1983-11-01'),ytop = 250,ybottom=0,col='gray',border=NA)
lines(abund_dat$censusdate[1:96],log(abund_dat$n[1:96]))
points(abund_dat$censusdate[1:96],log(abund_dat$n[1:96]),pch=16)
box(which='plot')

# line segments showing changepoint 95% intervals
segments(as.Date('1983-12-01'),5.8 , as.Date('1984-07-01'), 5.8, col='red', lwd=3, xpd = FALSE)
segments(as.Date('1988-10-01'),5.8 , as.Date('1996-01-01'), 5.8, col='red', lwd=3, xpd = FALSE)
segments(as.Date('1998-09-01'),5.8 , as.Date('1999-12-01'), 5.8, col='red', lwd=3, xpd = FALSE)
segments(as.Date('2009-06-01'),5.8 , as.Date('2010-09-01'), 5.8, col='red', lwd=3, xpd = FALSE)
