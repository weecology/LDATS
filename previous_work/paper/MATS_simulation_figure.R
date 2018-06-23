
library(multipanelfigure)
library(ggplot2)




###############################################################################
# Functions


#' create simulated data for demonstrating LDATS
#' 2 topics, nonuniform distribution of species in two community-types
#'     
#' @param tsteps = number of [monthly] time steps
#' 
#' @return 
#'    beta = matrix of species composition of the groups
#'    gamma = matrix of topic composition over time
#'            3 simulations of gamma: uniform, slow transition, and fast transitions
create_sim_data_ldats = function(tsteps=400) {
  
  topics = 2
  nspecies = 12
  
  # beta: species composition of topics
  # I calculated this distribution by taking the average of each Portal sampling sp distribution (periods 1:436)
  distribution = c(27,13,7, 5, 3, 2, 1, 1, 1, 0, 0, 0)
  # simple permutation of the first distribution
  distribution2 = c(3,1, 0, 1, 0, 13,2, 0, 1,27, 5, 7)
  
  beta = matrix(rep(0,topics*nspecies),nrow=topics,ncol=nspecies)
  beta[1,] = distribution/sum(distribution)
  beta[2,] = distribution2/sum(distribution2)
  
  # gamma for a constant topic prevalence through time: topic1 at 90% and topic2 at 10%
  gamma_constant = matrix(rep(0,tsteps*topics),nrow=tsteps,ncol=topics)
  gamma_constant[,1] = rep(.9,tsteps)
  gamma_constant[,2] = rep(.1,tsteps)
  
  # gamma for a slow gradual transition from topic1 to topic2 
  slow_series = matrix(rep(0,tsteps*topics),nrow=tsteps,ncol=topics)
  # proportions are constant before change
  slow_series[1:50,1] = rep(1)
  slow_series[1:50,2] = rep(0)
  # duration of change
  slow_series[(50+1):350,1] = seq(300)*(-1/300)+1
  slow_series[(50+1):350,2] = seq(300)*(1/300)+0
  # proportions are constant for rest of time series
  slow_series[(350+1):400,1] = rep(0)
  slow_series[(350+1):400,2] = rep(1)
  
  
  # gamma for two fast transitions
  fast_series = matrix(rep(0,tsteps*topics),nrow=tsteps,ncol=topics)
  # proportions are constant before change
  fast_series[1:100,1] = rep(1)
  fast_series[1:100,2] = rep(0)
  # change 1: lasts 12 timesteps
  fast_series[(100+1):112,1] = seq(12)*(-.6/12)+1
  fast_series[(100+1):112,2] = seq(12)*(.6/12)+0
  # stable at .6 and .4 for a while
  fast_series[113:250,1] = rep(.4)
  fast_series[113:250,2] = rep(.6)
  # change 2: lasts 12 timesteps
  fast_series[(250+1):262,1] = seq(12)*(.3/12)+.4
  fast_series[(250+1):262,2] = seq(12)*(-.3/12)+.6
  # constant for the rest
  fast_series[263:400,1] = rep(.7)
  fast_series[263:400,2] = rep(.3)
  
  return(list(beta,gamma_constant,slow_series,fast_series))
}


#' @param ldamodel lda model output
#' @param sim_dates xaxis values
#' @param x_labels T/F whether you want to plot x ticks
#' 
#' 
gg_plot_gamma = function(ldamodel,sim_dates,x_labels=F) {
  z = posterior(ldamodel)
  xticks=sim_dates
  ldaplot = data.frame(date=c(),relabund = c(), community = c())
  for (t in 1:2) {
    ldaplot = rbind(ldaplot,data.frame(date=xticks,relabund=z$topics[,t],community = as.factor(rep(t,length(z$topics[,1])))))
  }
  if (x_labels==F) {x_text = c('','','','');xname=''} else {x_text=c('1980','1990','2000','2010');xname='Time'}
  g = ggplot(ldaplot, aes(x=date,y=relabund,colour=community)) + 
    geom_line(size=1.5) +
    scale_y_continuous(limits=c(0,1),name='') +
    scale_x_date(name=xname,breaks=c(as.Date('1980-01-01'),as.Date('1990-01-01'),as.Date('2000-01-01'),as.Date('2010-01-01')),labels=x_text) +
    theme(axis.text.x=element_text(size=8),
          panel.border=element_rect(colour='black',fill=NA),
          legend.position='none') +
    scale_colour_manual(breaks=as.character(seq(2)),
                        values=cbPalette[c(1,3)],
                        guide=FALSE)
  return(g)
}


##################################################################################

cbPalette <- c( "#e19c02","#999999", "#56B4E9", "#0072B2", "#D55E00", "#F0E442", "#009E73", "#CC79A7")


N = 200   # total number of individuals
set.seed(1)

n <- 400
N <- cumsum(sample(c(-10:10), n, TRUE)) + 150

output = create_sim_data_ldats()

# distribution of species in the sample communities follows a species abundance distribution derived from Portal data (average of sampling periods 1:436)
# the distribution of species in the second sample community is a simple permutation of the first
# the two communities both contain 9 out of 12 species
beta = as.matrix(as.data.frame(output[1]))
colnames(beta) <- list('S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11','S12')

# plot three types of simulated dynamics for the 2 sample communities
gamma_constant = as.matrix(as.data.frame(output[2]))
gamma_slow = as.matrix(as.data.frame(output[3]))
gamma_fast = as.matrix(as.data.frame(output[4]))

sim_dates = seq.Date(from=as.Date('1977-01-01'),by=30,length.out = 400) 

slow =  data.frame(date = rep(sim_dates,dim(gamma_slow)[2]),
                   relabund = as.vector(gamma_slow),
                   community = as.factor(c(rep(1,dim(gamma_slow)[1]),rep(2,dim(gamma_slow)[1]))))
fast =  data.frame(date = rep(sim_dates,dim(gamma_fast)[2]),
                   relabund = as.vector(gamma_fast),
                   community = as.factor(c(rep(1,dim(gamma_fast)[1]),rep(2,dim(gamma_fast)[1]))))
const = data.frame(date = rep(sim_dates,dim(gamma_constant)[2]),
                   relabund = as.vector(gamma_constant),
                   community = as.factor(c(rep(1,dim(gamma_constant)[1]),rep(2,dim(gamma_constant)[1]))))


# create data sets from beta and gamma; data must be in integer form (simulating species counts)
dataset1 = round(as.data.frame(gamma_slow %*% beta) *N,digits=0)
dataset2 = round(as.data.frame(gamma_fast %*% beta) *N,digits=0)
dataset3 = round(as.data.frame(gamma_constant %*% beta) *N,digits=0)

dataset1$time = sim_dates
dataset2$time = sim_dates
dataset3$time = sim_dates

# plot species pop over time for 3 datasets
data1 <- tidyr::gather(dataset1, species, n, S1:S12, factor_key=TRUE)
data2 <- tidyr::gather(dataset2, species, n, S1:S12, factor_key=TRUE)
data3 <- tidyr::gather(dataset3, species, n, S1:S12, factor_key=TRUE)

datebreaks=c(as.Date('1980-01-01'),as.Date('1990-01-01'),as.Date('2000-01-01'),as.Date('2010-01-01'))

pop1 = ggplot(data1,aes(x=time,y=n,colour=species)) +
  geom_line(size=1.5) +
  ylab('') +
  scale_x_date(name='',breaks=datebreaks,labels=c('','','','')) +
  scale_y_continuous(limits=c(0,140)) +
  theme(legend.position="none",
        axis.text.x=element_text(size=8),
        panel.border=element_rect(color='black',fill=NA)) 
pop2 = ggplot(data2,aes(x=time,y=n,colour=species)) +
  geom_line(size=1.5) +
  scale_x_date(name='',breaks=datebreaks,labels=c('','','','')) +
  scale_y_continuous(name='',limits=c(0,140)) +
  theme(legend.position="none",
        axis.text.x=element_text(size=8),
        panel.border=element_rect(color='black',fill=NA))
pop3 = ggplot(data3,aes(x=time,y=n,colour=species)) +
  geom_line(size=1.5) +
  ylab('') +
  scale_x_date(name='Time',breaks=datebreaks,labels=c('1980','1990','2000','2010')) +
  scale_y_continuous(limits=c(0,140)) +
  theme(legend.position="none",
        axis.text.x=element_text(size=8),
        panel.border=element_rect(color='black',fill=NA)) 
gridExtra::grid.arrange(pop1,pop2,pop3,nrow=1)


#############################################################################################
# Run LDA
SEED  = 1

ldamodel1 = LDA(dataset1[,-13],k=2, control = list(seed = SEED,estimate.alpha=F,alpha=.1),method='VEM')
ldamodel2 = LDA(dataset2[,-13],k=2, control = list(seed = SEED,estimate.alpha=F,alpha=10),method='VEM')
ldamodel3 = LDA(dataset3[,-13],k=2, control = list(seed = SEED,estimate.alpha=F,alpha=.1),method='VEM')


#  plot gammas
g1 = gg_plot_gamma(ldamodel1,sim_dates)
g2 = gg_plot_gamma(ldamodel2,sim_dates)
g3 = gg_plot_gamma(ldamodel3,sim_dates,T)
gridExtra::grid.arrange(g1,g2,g3,nrow=1)


# species composition bar plots
spcomp1 = data.frame(topic=c(rep('t1',12),rep('t2',12)),
                     species=c(colnames(beta1),colnames(beta1)),
                     percent=c(beta1[1,],beta1[2,]))
spcomp1$species <- factor(spcomp1$species,levels = colnames(beta1))
gg_spcomp1 = ggplot(spcomp1) +
  geom_bar(stat='identity',aes(x=species,y=percent,fill=topic),
           position = position_dodge(width = .33)) +
  scale_x_discrete(name='',breaks=levels(spcomp1$species),labels=rep('',12)) +
  scale_y_continuous(name='') +
  scale_fill_manual(values=cbPalette[c(1,3)]) +
  theme(legend.position = 'none',
        panel.border=element_rect(color='black',fill=NA),
        axis.text.x=element_text(size=8))


spcomp2 = data.frame(topic=c(rep('t1',12),rep('t2',12)),
                     species=c(colnames(beta2),colnames(beta2)),
                     percent=c(beta2[1,],beta2[2,]))
spcomp2$species <- factor(spcomp2$species,levels = colnames(beta2))
gg_spcomp2 = ggplot(spcomp2) +
  geom_bar(stat='identity',aes(x=species,y=percent,fill=topic),
           position = position_dodge(width = .33)) +
  scale_x_discrete(name='',breaks=levels(spcomp2$species),labels=rep('',12)) +
  scale_y_continuous(name='') +
  scale_fill_manual(values=cbPalette[c(1,3)]) +
  theme(legend.position = 'none',
        panel.border=element_rect(color='black',fill=NA),
        axis.text.x=element_text(size=8))


spcomp3 = data.frame(topic=c(rep('t1',12),rep('t2',12)),
                     species=c(colnames(beta3),colnames(beta3)),
                     percent=c(beta3[1,],beta3[2,]))
spcomp3$species <- factor(spcomp3$species,levels = colnames(beta3))
gg_spcomp3 = ggplot(spcomp3) +
  geom_bar(stat='identity',aes(x=species,y=percent,fill=topic),
           position = position_dodge(width = .33)) +
  scale_x_discrete(name='Species') +
  scale_y_continuous(name='') +
  scale_fill_manual(values=cbPalette[c(1,3)]) +
  theme(legend.position = 'none',
        panel.border=element_rect(color='black',fill=NA),
        axis.text.x = element_text(size=8))
gridExtra::grid.arrange(gg_spcomp1,gg_spcomp2,gg_spcomp3,nrow=1)

#############################################################################################
# Big figure

(figure1 <- multi_panel_figure(
  width = c(80,80,80),
  height = c(80,80,80),
  panel_label_type = "none",
  column_spacing = 0,
  row_spacing = 0))
figure1 %<>% fill_panel(
  pop1,
  row = 1, column = 1)
figure1 %<>% fill_panel(
  pop2,
  row = 2, column = 1)
figure1 %<>% fill_panel(
  pop3,
  row = 3, column = 1)
figure1 %<>% fill_panel(
  gg_spcomp1,
  row = 1, column = 2)
figure1 %<>% fill_panel(
  gg_spcomp2,
  row = 2, column = 2)
figure1 %<>% fill_panel(
  gg_spcomp3,
  row = 3, column = 2)
figure1 %<>% fill_panel(
  g1,
  row = 1, column = 3)
figure1 %<>% fill_panel(
  g2,
  row = 2, column = 3)
figure1 %<>% fill_panel(
  g3,
  row = 3, column = 3)


figure1

