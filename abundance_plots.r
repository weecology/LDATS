# scripts for making supporting figures in the extreme events/ LDA project

# Figure 5 in manuscript

# ==============================================
# total abundance

dat = read.csv('Rodent_table_dat.csv',na.strings = '',as.is=T)
rdat = read.csv('../PortalData/Rodents/Portal_rodent.csv',as.is=T)
trappingdat = read.csv('../PortalData/Rodents/Portal_rodent_trapping.csv',as.is=T,na.strings = '')
trapdat = aggregate(trappingdat$sampled,by=list(period=trappingdat$period),FUN=sum)
fullcensus = trapdat[trapdat$x>20,]
perioddates = read.csv('../PortalData/Rodents/moon_dates.csv',as.is=T,na.strings = '')
perioddates$censusdate = as.Date(perioddates$censusdate)
fullcensus = merge(fullcensus,perioddates)

# data frame
abund_dat = data.frame(period = 1:436, n = rowSums(dat))
abund_dat = merge(abund_dat,fullcensus[,c('period','censusdate')])

# plot
plot(abund_dat$censusdate,log(abund_dat$n),xlab='',ylab='Log Total Abund.',pch=19,ylim=c(2,6),cex.axis=1.5,cex.lab=1.5)
rect(xleft = as.Date('1999-07-01'),xright = as.Date('1999-10-01'),ytop = 250,ybottom=0,col='gray',border=NA)
rect(xleft = as.Date('1983-08-01'),xright = as.Date('1983-11-01'),ytop = 250,ybottom=0,col='gray',border=NA)
rect(xleft = as.Date('1993-09-01'),xright = as.Date('1994-10-01'),ytop = 250,ybottom=0,col='gray',border=NA)
rect(xleft = as.Date('2009-03-01'),xright = as.Date('2010-01-01'),ytop = 250,ybottom=0,col='gray',border=NA)
lines(abund_dat$censusdate,log(abund_dat$n))
points(abund_dat$censusdate,log(abund_dat$n),pch=16)
abline(h=log(mean(abund_dat$n)))
box(which='plot')

# line segments showing changepoint 95% intervals
segments(as.Date('1983-12-01'),5.8 , as.Date('1984-07-01'), 5.8, col='red', lwd=3, xpd = FALSE)
segments(as.Date('1988-10-01'),5.8 , as.Date('1996-01-01'), 5.8, col='red', lwd=3, xpd = FALSE)
segments(as.Date('1998-09-01'),5.8 , as.Date('1999-12-01'), 5.8, col='red', lwd=3, xpd = FALSE)
segments(as.Date('2009-06-01'),5.8 , as.Date('2010-09-01'), 5.8, col='red', lwd=3, xpd = FALSE)


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
