# LDA analysis (with figures) for the extreme events/ LDA project

# LDA =================================
library(topicmodels)
library(ggplot2)

# load csv of abundance of each sp per plot
dat = read.csv('Rodent_table_dat.csv',na.strings = '',as.is=T)

# dates of trapping periods
perdat= read.csv('period_dates_single.csv')

perdat$date = as.Date(perdat$date,format='%m/%d/%Y')

# LDA
k = 5
ldamodel = LDA(dat,k,control=list(seed=100,estimate.alpha=F,alpha=1),method="VEM")

# composition of component communities, to identify most important species
structure(round(exp(ldamodel@beta), 3), dimnames = list(NULL, ldamodel@terms))

# plots like in Valle et al 2014
#get parameter estimates
z=posterior(ldamodel)
commun.plot=z$topics
commun.spp=z$terms

dates = as.Date(perdat$date[1:length(dat[,1])])


# ===========================================
# figures

#plot relative abundance of component communities for each sampling unit
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(dates[1],dates[436]),ylim=c(0,1),xlab='',ylab='Relative abundance',xaxt='n')
for (i in 1:k){
  lines(dates,commun.plot[,i],col=i)
}
axis(1,at=as.Date(c('1980-01-01','1985-01-01','1990-01-01','1995-01-01','2000-01-01','2005-01-01','2010-01-01','2015-01-01')),
     labels=c('1980','1985','1990','1995','2000','2005','2010','2015'))

# ggplot version
ldaplot = data.frame(date = dates,comm1=commun.plot[,1],comm2=commun.plot[,2])
ggplot(ldaplot, aes(date)) + 
  geom_line(aes(y = comm1, colour = "Community A"),size=1) + 
  geom_line(aes(y = comm2, colour = "Community B"),size=1) +
  scale_x_date(name='') +
  scale_y_continuous(name='Percent Similarity') +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18),
        legend.text=element_text(size=18),
        legend.title=element_text(size=18))

