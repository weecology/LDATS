
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



# ==================================================================================
# histogram to show low abundances
abund_dat$ratio = abund_dat$n/8
histo = hist(abund_dat$ratio,breaks = seq(0,30,.5),xlab='#animals per plot',ylab='freqency',main='')
abline(v=1.125,col='red',lwd=2) # Sep 1999
abline(v=1.5,col='blue',lwd=2) # Jan 2010
abline(v=2.25,col='forestgreen',lwd=2) # Sep 1993
abline(v=3.25,col='orange',lwd=2) # 1984 after octave


# stacked bar plot version, with colors for events
abund_dat$color = rep(0)
abund_dat$color[257] = 1      # Sep 1999
abund_dat$color[377:385] = 2  # Sep 2009 - May 2010 (< 4 animals/plot)
abund_dat$color[184:199] = 3  # Sep 1993 - Sep 1994
abund_dat$color[77:82] = 4    # Mar 1984 - Aug 1984

# expanded definition of droughts
#abund_dat$color = rep(0)
#abund_dat$color[257] = 1      # Sep 1999
#abund_dat$color[377:391] = 2  # Sep 2009 - Mar 2011 (< 4 animals/plot)
#abund_dat$color[178:209] = 3  # Dec 1992 - Jun 1995
#abund_dat$color[77:82] = 4    # Mar 1984 - Aug 1984

abund_dat$bin = rep(0)
for (n in seq(length(abund_dat$Period))) {
  for (b in seq(length(histo$breaks))) {
    if (abund_dat$ratio[n] > histo$breaks[b] & abund_dat$ratio[n] <= histo$breaks[b+1]) {
      abund_dat$bin[n] = b
    }
  }
}

qplot(abund_dat$bin,data=abund_dat, geom='bar',fill=factor(abund_dat$color)) + 
  geom_bar() + 
  scale_fill_manual(breaks = c("0", "1", "2","3","4"), 
                    values=c("gray", "red", "blue","forestgreen","orange"),
                    labels=c("","1999 flood","2009 drought","1993 drought","1983 hurricane"),
                    name = "") +
  scale_x_continuous(breaks=seq(0,60,10),
                     labels=seq(0,30,5),
                     limits = c(0,60),
                     name = '')

# =======================================================================================
# histogram of ndvi
ndvidat = read.csv('C:/Users/EC/Documents/NDVI/Monthly_Landsat_NDVI.csv')
ndvi = ndvidat[!is.na(ndvidat$NDVI),]
ndvi$n = seq(length(ndvi$Date))

#histogram of all monthly ndvi averages
histo = hist(ndvi$NDVI,breaks=seq(.05,.5,.01),xlab='NDVI')
ndvi$color = rep(0)
ndvi$color[271:281] = 1  # Sep 2009 - May 2010 
ndvi$color[86:98] = 2  # Sep 1993 - Sep 1994

ndvi$bin = rep(0)
for (n in seq(length(ndvi$Date))) {
  for (b in seq(length(histo$breaks))) {
    if (ndvi$NDVI[n] > histo$breaks[b] & ndvi$NDVI[n] <= histo$breaks[b+1]) {
      ndvi$bin[n] = b
    }
  }
}

qplot(ndvi$bin,data=ndvi, geom='bar',fill=factor(ndvi$color)) + 
  geom_bar() + 
  scale_fill_manual(breaks = c("0", "1", "2"), 
                    values=c("gray", "blue","forestgreen"),
                    labels=c("","2009 drought","1993 drought"),
                    name = "") 

#histogram of just summer peaks
summerndvi_peak = data.frame(year = c(),ndvi = c())
summerndvi = data.frame(year = c(),ndvi=c())
for (n in seq(1986,2015,1)) {
  summ_months = c(paste(n,'-07',sep=''),paste(n,'-08',sep=''))
  summ_dat = ndvi[ndvi$Date %in% summ_months,]
  summ_peak = max(summ_dat$NDVI,na.rm=T)
  summerndvi_peak = rbind(summerndvi_peak,c(n,summ_peak))
  summerndvi = rbind(summerndvi,summ_dat)
}
names(summerndvi_peak) = c('year','NDVI')
hist(summerndvi_peak$NDVI)
histo = hist(summerndvi$NDVI,breaks = seq(.05,.5,.02))
summerndvi$bin = rep(0)
for (n in seq(length(summerndvi$Date))) {
  for (b in seq(length(histo$breaks))) {
    if (summerndvi$NDVI[n] > histo$breaks[b] & summerndvi$NDVI[n] <= histo$breaks[b+1]) {
      summerndvi$bin[n] = b
    }
  }
}

qplot(summerndvi$bin,data=summerndvi, geom='bar',fill=factor(summerndvi$color)) + 
  geom_bar() + 
  scale_fill_manual(breaks = c("0", "1", "2"), 
                    values=c("gray", "blue","forestgreen"),
                    labels=c("","2009 drought","1993 drought"),
                    name = "") 

# ==========================================================================================================
# look at 2009-2014
plot(dates[369:433],abund_dat$ratio[369:433],main='rodents',xlab='',ylab='animals',type = 'l',lwd=2)
plot(weather_ndvi$date[224:294],weather_ndvi$ppt[224:294],xlab='',ylab='precip',type='l',lwd=2,main='monthly precip')
plot(weather_ndvi$date[224:294],weather_ndvi$NDVI[224:294],xlab='',ylab='NDVI',type='l',lwd=2,main='monthly NDVI')

az = read.csv('C:/Users/EC/Desktop/Drought Index/CDODiv8736526894107_AZ.csv')
nm = read.csv('C:/Users/EC/Desktop/Drought Index/CDODiv9809606894112_NM.csv')

plot(weather_ndvi$date[224:294],az$PDSI[385:455],main='PDSI',xlab='',ylab='PDSI',type = 'l',lwd=2)

histo=hist(az$ZNDX,breaks=seq(-5,9,.5))
az$bin = rep(0)
for (n in seq(length(az$YearMonth))) {
  for (b in seq(length(histo$breaks))) {
    if (az$ZNDX[n] > histo$breaks[b] & az$ZNDX[n] <= histo$breaks[b+1]) {
      az$bin[n] = b
    }
  }
}
az$color = rep(0)
az$color[273:281] = 1  # Sep 2009 - May 2010 
az$color[86:98] = 2  # Sep 1993 - Sep 1994

qplot(az$bin,data=az, geom='bar',fill=factor(az$color)) + 
  geom_bar() + 
  scale_fill_manual(breaks = c("0", "1", "2"), 
                    values=c("gray", "blue","forestgreen"),
                    labels=c("","2009 drought","1993 drought"),
                    name = "") 

