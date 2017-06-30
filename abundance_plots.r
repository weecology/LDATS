# scripts for making supporting figures in the extreme events/ LDA project



# ==============================================
# total abundance

dat = read.csv('Rodent_table_dat.csv',na.strings = '',as.is=T)
rdat = read.csv('../PortalData/Rodents/Portal_rodent.csv',as.is=T)
trappingdat = read.csv('../PortalData/Rodents/Portal_rodent_trapping.csv',as.is=T,na.strings = '')
trapdat = aggregate(trappingdat$Sampled,by=list(Period=trappingdat$Period),FUN=sum)
fullcensus = trapdat[trapdat$x>20,]
perioddates = read.csv('../PortalData/Rodents/moon_dates.csv',as.is=T,na.strings = '')
perioddates$CensusDate = as.Date(perioddates$CensusDate)
fullcensus = merge(fullcensus,perioddates)

# data frame
abund_dat = data.frame(Period = 1:436, n = rowSums(dat))
abund_dat = merge(abund_dat,fullcensus[,c('Period','CensusDate')])

# plot
plot(abund_dat$CensusDate,log(abund_dat$n),xlab='',ylab='Log Total Abundance',pch=19,ylim=c(2,6))
rect(xleft = as.Date('1999-07-01'),xright = as.Date('1999-10-01'),ytop = 250,ybottom=0,col='gray',border=NA)
rect(xleft = as.Date('1983-08-01'),xright = as.Date('1983-11-01'),ytop = 250,ybottom=0,col='gray',border=NA)
rect(xleft = as.Date('1993-09-01'),xright = as.Date('1994-10-01'),ytop = 250,ybottom=0,col='gray',border=NA)
rect(xleft = as.Date('2009-03-01'),xright = as.Date('2010-01-01'),ytop = 250,ybottom=0,col='gray',border=NA)
lines(abund_dat$CensusDate,log(abund_dat$n))
points(abund_dat$CensusDate,log(abund_dat$n),pch=16)
abline(h=log(mean(abund_dat$n)))
box(which='plot')

# line segments showing changepoint 95% intervals
segments(as.Date('1983-12-01'),5.8 , as.Date('1984-07-01'), 5.8, col='red', lwd=3, xpd = FALSE)
segments(as.Date('1988-10-01'),5.8 , as.Date('1996-01-01'), 5.8, col='red', lwd=3, xpd = FALSE)
segments(as.Date('1998-09-01'),5.8 , as.Date('1999-12-01'), 5.8, col='red', lwd=3, xpd = FALSE)
segments(as.Date('2009-06-01'),5.8 , as.Date('2010-09-01'), 5.8, col='red', lwd=3, xpd = FALSE)


par(xpd=F)

