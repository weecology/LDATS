
# Plot NDVI over time to show drought
library(RCurl)
library(dplyr)
library(ggplot2)

ndvi = read.csv('Monthly_Landsat_NDVI.csv',
                   na.strings=c(""), stringsAsFactors = FALSE)

ndvi$Date = as.Date(ndvi$Date,format='%Y-%M')
ndvi$NDVI = as.numeric(ndvi$NDVI)
ndvi$year = format(ndvi$date,'%Y') %>% as.numeric()

plot(ndvi$date,ndvi$NDVI)
lines(ndvi$date,ndvi$NDVI)


# ===============
# yearly avg - from 1984 there are >12 images per year (should go back and remove excessively cloudy images too)

ndviyr = aggregate(ndvi$NDVI,by=list(year = ndvi$year),FUN=mean,na.rm=T) %>% filter(year>1983)


# ggplot version - with long term mean
ggplot(ndviyr,aes(x=as.integer(year),y=x)) +
  geom_line(size=1.5) +  
  geom_point(size=3) +
  geom_hline(yintercept = mean(ndviyr$x),linetype=2) +
  labs(x='',y='NDVI (yearly mean)') +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14))
  
