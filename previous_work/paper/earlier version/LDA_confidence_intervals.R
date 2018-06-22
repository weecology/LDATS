# This script does an LDA analysis on rodent data, using the gibbs samp method so we can get credible intervals
# 3/2017

library(ggplot2)


setwd('C:/Users/EC/Desktop/git/Extreme-events-LDA')
source('gibbs_functions.R')

# ================================================================
# load data and run LDA model
dat.agg=data.matrix(read.csv('Rodent_table_dat.csv',as.is=T))

ngibbs=1000 # has to be greater than 200
ncommun=2  # number of topics
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

# make data frame of theta estimate and sd

thetadf = data.frame(grp1 = theta1[1,],grp2 = theta1[2,],sd1=thetasd[1,],sd2=thetasd[2,],grp3=theta1[3,],sd3=thetasd[3,])#,grp4=theta1[4,],sd4=thetasd[4,])

ggplot(betadf) +
  geom_ribbon(data = thetadf, mapping = aes_string(x = seq(436), ymin = thetadf$grp1-1.96*thetadf$sd1, ymax = thetadf$grp1+1.96*thetadf$sd1,alpha=.4), fill = "grey") +
  geom_line(aes(y = thetadf$grp1,x=seq(21))) +
  geom_ribbon(data = thetadf, mapping = aes_string(x = seq(436), ymin = thetadf$grp2-1.96*thetadf$sd2, ymax = thetadf$grp2+1.96*thetadf$sd2,alpha=.4), fill = "pink") +
  geom_line(aes(y = thetadf$grp2,x=seq(21)),color='red') 
  geom_ribbon(data = thetadf, mapping = aes_string(x = seq(436), ymin = thetadf$grp3-1.96*thetadf$sd3, ymax = thetadf$grp3+1.96*thetadf$sd3,alpha=.4), fill = "lightgreen") +
  geom_line(aes(y = thetadf$grp3,x=seq(436)),color='forestgreen') 
  geom_ribbon(data = thetadf, mapping = aes_string(x = seq(436), ymin = thetadf$grp4-1.96*thetadf$sd4, ymax = thetadf$grp4+1.96*thetadf$sd4,alpha=.4), fill = "plum") +
  geom_line(aes(y = thetadf$grp4,x=seq(436)),color='purple4') 
  
  
  # plot credible intervals around sp comp
  betadf = data.frame(grp1 = beta1[1,],grp2 = beta1[2,],sd1=betasd[1,],sd2=betasd[2,])
  
  ggplot(betadf) +
    geom_ribbon(data = betadf, mapping = aes_string(x = seq(21), ymin = betadf$grp1-1.96*betadf$sd1, ymax = betadf$grp1+1.96*betadf$sd1,alpha=.4), fill = "grey") +
    geom_line(aes(y = betadf$grp1,x=seq(21))) +
    geom_ribbon(data = betadf, mapping = aes_string(x = seq(21), ymin = betadf$grp2-1.96*betadf$sd2, ymax = betadf$grp2+1.96*betadf$sd2,alpha=.4), fill = "pink") +
    geom_line(aes(y = betadf$grp2,x=seq(21)),color='red') 
