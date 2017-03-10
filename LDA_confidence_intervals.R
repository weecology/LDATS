# This script does an LDA analysis on rodent data, using the gibbs samp method so we can get credible intervals
# 3/2017

library(ggplot2)


setwd('C:/Users/EC/Desktop/git/Extreme-events-LDA')
source('gibbs functions.R')

# ================================================================
# load data and run LDA model
dat.agg=data.matrix(read.csv('Rodent_table_dat.csv',as.is=T))

ngibbs=1000 # has to be greater than 200
ncommun=3   # number of topics
results=gibbs.samp(dat.agg=dat.agg,ngibbs=ngibbs,ncommun=ncommun,a.betas=1,a.theta=1)

#plot results
plot(1:ngibbs,results$logL,type='l',ylim=range(results$log[50:ngibbs]))

coda::effectiveSize(results$logL) 

nplots=nrow(dat.agg)
nspp=ncol(dat.agg)

beta1=matrix(apply(results$beta,2,mean),ncommun,nspp)
par(mfrow=c(ncommun,1))
for (i in 1:ncommun) {plot(1:nspp,beta1[i,],type='l',xaxt='none',xlab='',ylim=c(0,.8),col=i,lwd=2) 
  axis(1,at=1:length(colnames(dat.agg)),labels = colnames(dat.agg))}

theta1=matrix(apply(results$theta,2,mean),ncommun,nplots)
par(mfrow=c(1,1))
plot(NA,NA,ylim=c(0,1),xlab='',xlim=c(0,436),xaxt='none',ylab='Relative Composition')
for (i in 1:ncommun) {lines(1:nplots,theta1[i,],type='l',xlab='',col=i,ylim=c(0,1))
  axis(1,at=c(30,86,142,207,261,320,382,434),
       labels=c('1980','1985','1990','1995','2000','2005','2010','2015'))}

# =================================================
# credible intervals

betasd = matrix(apply(results$beta,2,sd),ncommun,nspp)
thetasd = matrix(apply(results$theta,2,sd),ncommun,nplots)

# ===================================================================
# ggplot version

thetadf = data.frame(grp1 = theta1[1,],grp2 = theta1[2,],grp3=theta1[3,],sd1=thetasd[1,],sd2=thetasd[2,],sd3=thetasd[3,])

ggplot(thetadf) +
  geom_ribbon(data = thetadf, mapping = aes_string(x = seq(436), ymin = thetadf$grp1-1.96*thetadf$sd1, ymax = thetadf$grp1+1.96*thetadf$sd1), fill = "grey") +
  geom_line(aes(y = thetadf$grp1,x=seq(436))) +
  geom_ribbon(data = thetadf, mapping = aes_string(x = seq(436), ymin = thetadf$grp2-1.96*thetadf$sd2, ymax = thetadf$grp2+1.96*thetadf$sd2), fill = "pink") +
  geom_line(aes(y = thetadf$grp2,x=seq(436)),color='red') +
  geom_ribbon(data = thetadf, mapping = aes_string(x = seq(436), ymin = thetadf$grp3-1.96*thetadf$sd3, ymax = thetadf$grp3+1.96*thetadf$sd3), fill = "lightgreen") +
  geom_line(aes(y = thetadf$grp3,x=seq(436)),color='forestgreen')
