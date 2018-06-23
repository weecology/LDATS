
# Plot NDVI over time to show drought
#library(RCurl)
library(dplyr)
library(ggplot2)

ndvi = read.csv('Monthly_Landsat_NDVI.csv',
                   na.strings=c(""), stringsAsFactors = FALSE)

ndvi$NDVI = as.numeric(ndvi$NDVI)
ndvi$year = as.integer(substr(ndvi$Date,1,4))



# ===============
# yearly avg - from 1984 there are >12 images per year (should go back and remove excessively cloudy images too)

ndviyr = aggregate(ndvi$NDVI,by=list(year = ndvi$year),FUN=mean,na.rm=T) %>% filter(year>1983)


# ggplot version - with long term mean
ggplot(ndviyr,aes(x=as.integer(year),y=x)) +
  geom_line(size=1) +  
  geom_point(size=2) +
  geom_hline(yintercept = mean(ndviyr$x),linetype=2) +
  labs(x='',y='NDVI (yearly mean)') +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))
  
ggsave(filename='FigureB-6.tiff',width=6,height=2.8,units='in',dpi=600,compression='lzw')
