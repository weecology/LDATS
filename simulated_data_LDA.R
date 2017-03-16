# Use simulated data to investigate the results of the LDA model with known dynamics:
#    1. fast transition between 2 topics with steady-state before and after
#    2. slow transition between 2 topics with little to no steady-state before/after
#    3. fixed proportion of 2 topics through time (no transition)


# Started EMC 3/2017

# ====================================================================================
# Create simulated data:
#  beta = matrix of species composition of the topics
#  gamma = matrix of prevalence of topics through time
#  assume even species composition and no overlap of species between topics

nspecies = 20
topics = 2
tsteps = 400 #I think of these as monthly time steps
N = 200      # total number of "animals" -- LDA model depends on integer counts as data

beta = matrix(rep(0,topics*nspecies),nrow=topics,ncol=nspecies)
beta[1,] = c(rep(.1,nspecies/2),rep(0,nspecies/2))
beta[2,] = c(rep(0,nspecies/2),rep(.1,nspecies/2))

# gamma for a constant topic prevalence through time: topic1 at 80% and topic2 at 20%
gamma_constant = matrix(rep(0,tsteps*topics),nrow=tsteps,ncol=topics)
gamma_constant[,1] = rep(.8,tsteps)
gamma_constant[,2] = rep(.2,tsteps)

# gamma for a fast transition from topic1 to topic2 (one year/12 time steps)
gamma_fast = matrix(rep(0,tsteps*topics),nrow=tsteps,ncol=topics)
# proportions are constant for first 200 time steps
gamma_fast[1:200,1] = rep(.9)
gamma_fast[1:200,2] = rep(.1)
# fast transition from tstep 201-212
gamma_fast[201:212,1] = seq(12)*(-.8/12)+.9
gamma_fast[201:212,2] = seq(12)*(.8/12)+.1
# proportions are constant for rest of time series
gamma_fast[213:400,1] = rep(.1)
gamma_fast[213:400,2] = rep(.9)

# gamma for a slow transition from topic1 to topic2 that takes the entire time series
gamma_slow = matrix(rep(0,tsteps*topics),nrow=tsteps,ncol=topics)
gamma_slow[,1] = seq(400)*(-.8/400)+.9
gamma_slow[,2] = seq(400)*(.8/400)+.1


# ==================================================================================
# plot beta and gammas
par(mfrow=c(topics,1))
for (i in 1:topics) {plot(1:nspecies,beta[i,],type='l',xlab='',ylim=c(0,.2),col=i,lwd=2) }

plot(gamma_constant[,1],type='l',xlab='time',ylab='topic',ylim=c(0,1),main='constant')
lines(gamma_constant[,2],col='red')
plot(gamma_fast[,1],type='l',xlab='time',ylab='topic',ylim=c(0,1),main='fast')
lines(gamma_fast[,2],col='red')
plot(gamma_slow[,1],type='l',xlab='time',ylab='topic',ylim=c(0,1),main='slow')
lines(gamma_slow[,2],col='red')


# create data set from beta and gamma; data must be in integer form
dataset1 = round(as.data.frame(gamma_fast %*% beta) *N,digits=0)
dataset2 = round(as.data.frame(gamma_constant %*% beta) *N,digits=0)
dataset3 = round(as.data.frame(gamma_slow %*% beta) *N,digits=0)


# =================================================================================
# run LDA model
nstart = 20 # For the final analysis, maybe do 1000
ldamodel = LDA(dataset3,2,control=list(estimate.alpha=F,alpha=.5, nstart = nstart),method="VEM")
plot_component_communities(ldamodel,2,seq(400))

